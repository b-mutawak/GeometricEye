import numpy as np
import json
from scipy import optimize
from math import ceil, cos, sin, sqrt, pi, atan, asin, atan2, radians
from numpy import linalg
from copy import deepcopy
import maya.cmds as cmds
import math_utils
import curve_utils
import surface_utils
from copy import copy
import maya.OpenMaya as OpenMaya

# Define rectus insertion points
rectus_insertions = {'LR':12.45, 'MR':10.85, 'SR':12.65, 'IR':11.85}
oblique_insertions = {'SO':[20,2], 'IO':[19.5,2]}
muscle_total_lengths = {'LR':43.5, 'MR':40.7, 'SR':41.6, 'IR':41.7, 'SO':35.5, 'IO':38.5}
tendon_lengths = {'SO':15.5, 'MR':3.0, "IR":4.7, "LR": 7.2, "SR":4.3, 'IO':-1}
muscle_lengths = {'LR':36.3, 'MR':37.7, 'SR':37.3, 'IR':37.0, 'SO':20.0, 'IO':-1}

muscle_color =  'red'
tendon_color =  'yellow'

#cmds.curveRGBColor( '*muscle_bspline', muscle_color[0], muscle_color[1], muscle_color[2] )
#cmds.curveRGBColor('*tendon_bspline', tendon_color[0], tendon_color[1], tendon_color[2])

def mat_svd(matrix):
    U, sdiag, VH = linalg.svd(matrix)
    S = np.zeros(matrix.shape)
    np.fill_diagonal(S, sdiag)
    V = VH.T.conj()
    return U, S, V

def projsplx(y):

    m = len(y)
    bget = False
    
    s = y[::-1].sort()
    tmpsum = 0
    
    for ii in range(m):
        tmpsum = tmpsum + s[ii]
        tmax = (tmpsum - 1)/ii
        if tmax >= s[ii+1]:
            bget = True
            break
        
    if not bget:
        tmax = (tmpsum + s[m-1] - 1) / m
    
    x = max(y-tmax,0)
    
    return x


def Project_on_B(q0):

    Q0= np.vstack(([q0[0],q0[3]/sqrt(2),q0[4]/sqrt(2)], 
                  [q0[3]/sqrt(2),q0[1],q0[5]/sqrt(2)], 
                  [q0[4]/sqrt(2),q0[5]/sqrt(2),q0[2]]))
    U,S0,V = mat_svd(Q0)
    e1=np.vstack((1, 1, 1, 0, 0, 0, 0, 0, 0, -1))
    e1=e1/np.linalg.norm(e1)
    q=q0.reshape(len(q0), 1)-((np.dot(q0,e1)-1) * e1)
    return q



def ellipsoid_fit_matlab(x, nit):
    """
        Python translation of Ellipsoid3D_Fitting_OBLIQUE.m
        Taken from: https://github.com/pierre-weiss/FitEllipsoid
    """
    # Grab number of points
    num_points = x.shape[1]
    
    # First block, D array
    D = np.zeros((10, num_points))
    D[:3, :] = np.power(x[:3, :], 2)         # Up to line 15
    D[3, :] = sqrt(2) * x[0, :] * x[1, :]
    D[4, :] = sqrt(2) * x[0, :] * x[2, :]
    D[5, :] = sqrt(2) * x[1, :] * x[2, :]
    D[6:9, :] = x[:3, :]
    D[9, :] = 1                              # Up to line 22
    
    # k?
    k = np.matmul(D, D.T)                    # line 24
    [_, s1, v1] = mat_svd(k)
    a = v1[:, 9]
    a = a / np.sum(a[:3])
    if a[0] < 0:
        a = -a
    
    # A? Line 32
    A = np.vstack(( [a[0], a[3] / sqrt(2), a[4] / sqrt(2)],
                    [a[3] / sqrt(2), a[1], a[5] / sqrt(2)],
                    [a[4] / sqrt(2), a[5] / sqrt(2), a[2]] ))
    if np.min(np.linalg.eig(A)[0]) < 0:
        print "Beware, solution is not an ellipse"
    
    
    # SM, line 37
    sm = np.vstack((a[6],
                    a[7],
                    a[8]))
    c = -0.5 * np.linalg.solve(A, sm)     # will need to make sure this is the same.
    r = np.sqrt(1e-16 - a[9] + np.dot(np.matmul(A, c).T, c))
    
    # SVD, line 41
    uA, sA, vA = mat_svd(A)
    Ademi = np.matmul(np.matmul(uA, np.sqrt(np.maximum(sA, 1e-16))) , vA.T)  # check this one...
    xx = np.matmul((1/r) * Ademi, x - np.tile(c, (1, num_points)))
    
    # Second D array, line 45
    D = np.zeros((10, num_points))
    D[:3, :] = np.power(xx[:3, :], 2)        
    D[3, :] = sqrt(2) * xx[0, :] * xx[1, :]
    D[4, :] = sqrt(2) * xx[0, :] * xx[2, :]
    D[5, :] = sqrt(2) * xx[1, :] * xx[2, :]
    D[6:9, :] = xx[:3, :]
    D[9, :] = 1                              
    
    k = np.matmul(D, D.T)                 # line 57
    
    # The objective is now to solve min <q,Kq>, Tr(Q)=1, Q>=0 using a Douglas-Rachford algorithm
    # Begin with a sphere. Line 60
    c1 = np.mean(xx[0, :])
    c2 = np.mean(xx[1, :])
    c3 = np.mean(xx[2, :])
    r2 = np.var(xx[0, :]) + np.var(xx[1, :]) + np.var(xx[2, :])
    u= np.vstack(((1.0/3), (1.0/3),(1.0/3),(0), (0), (0), (-2*c1*1.0/3), (-2*c2*1.0/3), (-2*c3*1.0/3), ((1.0*c1**2+c2**2+c3**2-r2)/3)));
    
    # Now for the iterative douglas-rachford algorithm
    gamma = 10
    m = gamma * k + np.identity(len(k))
    proxf1 = lambda q: linalg.solve(m, q)
    proxf2 = lambda q: Project_on_B(q)
    u1, s1, v1 = mat_svd(k)
    p = np.real(v1[:, 9])
    p = p / np.sum(p[0:3])
    cf = np.zeros((nit + 1, 1))
    p = np.reshape(p, [10, 1])
    # Line 79
    for i in range(nit):
        q = Project_on_B(np.reshape(p, [10]))
        p = p + proxf1(2 * q - p) - q
        cf[i] = np.matmul(np.matmul(0.5*q.T, k), q)
    
    q = proxf2(np.reshape(q, [10]))
    q2 = q
    cf[nit] = np.matmul(np.matmul(0.5*q.T, k), q)
    
    qs = np.zeros((10, 1))
    
    # Line 92
    b = np.vstack(([q2[0], q2[3] / sqrt(2), q2[4] / sqrt(2)],
                   [q2[3] / sqrt(2), q2[1], q2[5]/sqrt(2)],
                   [q2[4]/sqrt(2), q2[5] / sqrt(2), q2[2]])).reshape([3, 3])
    
    sm = np.vstack((q2[6], q2[7], q2[8]))
    cprime = -0.5 * np.linalg.solve(b, sm)
    rprime = np.sqrt(-1 * q2[9] + np.dot(np.matmul(b, cprime).T, cprime))
    qs = np.zeros((10, 1))
    m = np.matmul(np.matmul(Ademi, b), Ademi)
    qs[:6] = (1 / r**2) * np.vstack((m[0, 0], m[1, 1], m[2, 2], sqrt(2) * m[1, 0], sqrt(2) * m[2, 0], sqrt(2) * m[2, 1]))
    qs[6:9] = np.matmul(np.matmul(-(2/r)*Ademi, b), np.matmul((1/r)*Ademi, c) + cprime)
    qs[9] = np.dot(np.matmul(b, np.matmul((1/r) * Ademi, c) + cprime).T,np.matmul((1/r) * Ademi, c) + cprime) - rprime**2
    qs = qs / np.sum(qs[:3])
    qs = np.real(qs)
    
    return getRealEllipsoidParameters(qs)
    
    
def getRealEllipsoidParameters(q):
    # Returns radii and center from q matrix
    
    matSR = np.zeros((3,3))
    matSR[0,0] = q[0]
    matSR[0,1] = q[3]
    matSR[0,2] = q[4]
    matSR[1,0] = q[3]
    matSR[1,1] = q[1]
    matSR[1,2] = q[5]
    matSR[2,0] = q[4]
    matSR[2,1] = q[5]
    matSR[2,2] = q[2]
    b = np.zeros((3, 1))
    b[0, 0] = q[6] / 2.0
    b[1, 0] = q[7] / 2.0
    b[2, 0] = q[8] / 2.0
    b = b * -1
    centerMat = np.linalg.solve(matSR, b)
    scalarMat = np.matmul(matSR, centerMat).T
    scalarMat = np.matmul(scalarMat, centerMat)
    r2 = scalarMat[0,0] - q[9]
    tmp = np.copy(matSR)
    u, s, v = mat_svd(tmp)
    xSemiLength = sqrt(r2 / s[0,0])
    ySemiLength = sqrt(r2 / s[1, 1])
    zSemiLength = sqrt(r2 / s[2, 2])
    
    return (xSemiLength, ySemiLength, zSemiLength), centerMat.T



def sphereFit(points):
    """
        Least-squares fit of 3D points to a sphere.
        Returns the radius and center of sphere.
        
        Credit to:
            https://jekel.me/2015/Least-Squares-Sphere-Fit/
    
    """
    spX = np.array([])
    spY = np.array([])
    spZ = np.array([])
    
    # Grab the points
    for z_index, values in points.items():
        values = np.array(values)
        spX = np.append(spX, values[:, 0])
        spY = np.append(spY, values[:, 1])
        spZ = np.append(spZ, values[:, 2])

    A = np.zeros((len(spX),4))
    A[:,0] = spX*2
    A[:,1] = spY*2
    A[:,2] = spZ*2
    A[:,3] = 1

    #   Assemble the f matrix
    f = np.zeros((len(spX),1))
    f[:,0] = (spX*spX) + (spY*spY) + (spZ*spZ)
    C, residules, rank, singval = np.linalg.lstsq(A,f)

    #   solve for the radius
    t = (C[0]*C[0])+(C[1]*C[1])+(C[2]*C[2])+C[3]
    radius = sqrt(t)

    return radius, C[0], C[1], C[2]



def circleFit(x_in, y_in):
    """ 
        Fits a circle to the set of X and Y points given.
        
        Returns center and radius. Adapted from:
            https://scipy-cookbook.readthedocs.io/items/Least_Squares_Circle.html
    """
    
    def calc_R(x, y, xc, yc):
        """ calculate the distance of each 2D points from the center (xc, yc) """
        x = np.asarray(x)
        y = np.asarray(y)
        xc = np.asarray(xc)
        yc = np.asarray(yc)
        return np.sqrt((x-xc)**2 + (y-yc)**2)
    
    def f_2(c, x, y):
        """ calculate the algebraic distance between the data points and the mean circle centered at c=(xc, yc) """
        
        Ri = calc_R(x, y, *c)
        return np.asarray(Ri) - np.asarray(Ri).mean()
    
    # Use least square fit to find circle
    x_in = np.asarray(x_in)
    y_in = np.asarray(y_in)
    center_estimate = np.mean(x_in), np.mean(y_in)
    center, ier = optimize.leastsq(f_2, center_estimate, args=(x_in, y_in))

    xc, yc = center
    Ri_2       = calc_R(x_in, y_in,*center)
    R        = Ri_2.mean()
    return xc, yc, R 


def getRotationMatrix(axis, angle):
    """
        Returns 3D rotation matrix along axis of angle
    """
    
    # Convert to radians
    angle = np.radians(angle)
    
    if axis == "x":
        rotMatrix = np.array([[1,       0,              0],
                              [0, cos(angle), -sin(angle)],
                              [0, sin(angle), cos(angle)]])
    elif axis == "y":
        rotMatrix = np.array([[cos(angle),  0, sin(angle)],
                              [0,       1,              0],
                              [-sin(angle), 0, cos(angle)]])
    elif axis == "z":
        rotMatrix = np.array([[cos(angle), -sin(angle), 0],
                              [sin(angle), cos(angle),  0],
                              [0,           0,          1]])
    
    return rotMatrix


def is_number(x):
    
    isnum = False
    try:
        float(x)
        isnum = True
    except:
        pass
    
    return isnum

def rotatePoints(axis, degrees, points):
    
    # Get rotation
    rotation = getRotationMatrix(axis, degrees)
    points = np.array(points)
    
    # Define axis of rotation to be the center
    center_axis = np.mean(points, axis=0)
        
    # Apply rotation                
    points = rotation.dot((points - center_axis).T).T + center_axis
    
    return points.tolist()

def correctRotation(data):
    """
        Corrects MR image rotation of muscles.
        Simply rotates OS and OD point clouds around the superior-inferior axis +- 23 degrees according to literature.
    
    """

    # Create rotation matrices
    os_rotation = getRotationMatrix("y", 23)
    od_rotation = getRotationMatrix("y", -23)
    
    # Copy data
    data_copy = deepcopy(data)
    
    # To map back
    data_lengths = {}
    
    # Data structures to hold os and od points
    data_points = {'OS':np.zeros((0,3)), 'OD':np.zeros((0,3))}
    data_keys = {'OS':[], 'OD':[]}
    
    # Loop through each plane/eye
    for point_set_name, point_sets in data_copy.items():
        
        # grab key for eye
        eye = 'OS' if 'OS' in point_set_name else 'OD'

        # add to keys. This structure will be something like [OD Sagittal Prim, LR, MR, ... OD Coronal Prim, LR, SR...]
        # when iterating back through to modify things, we simply iterate through this key list, grabbing the first key always.
        data_keys[eye].append(point_set_name)
        data_lengths[point_set_name] = {}
                
        
        # Iterate through each muscle in the set
        for muscle_name, points in point_sets.items():
            
            # Grab the points
            sorted_points = data[point_set_name][muscle_name]
            
            # add to keys
            data_keys[eye].append(muscle_name)
            
            data_lengths[point_set_name][muscle_name] = {}
            
            # Loop through each z-index points
            for z_index, points in sorted_points.items():
            
                # Grab the points
                points = np.array(points)
                
                # add to keys
                data_keys[eye].append(z_index)
                
                # Update lengths
                data_lengths[point_set_name][muscle_name][z_index] = points.shape[0]
                
                # Add to eye structure
                data_points[eye] = np.append(data_points[eye], points, axis=0)
                
    
    # Now apply rotation to each eye
    for eye, points in data_points.items():
        
        
        # Check which eye this is
        rotation = os_rotation if 'OS' in eye else od_rotation 
        
        
        # Define axis of rotation to be the center
        center_axis = np.mean(points, axis=0)
        
        # Apply rotation                
        points = rotation.dot((points - center_axis).T).T + center_axis
        
        data_points[eye] = points
    
    
    # Now map things back
    for eye, keys in data_keys.items():
        
        # This is the eye keys
        major_key = None
        
        # This is the muscle key
        muscle_key = None
                
        # Iterate as long as we have a key available
        while not len(keys) == 0:
            
            # Grab first key
            key = keys.pop(0)
            
            
            # check if this is a major key
            if 'OS' in key or 'OD' in key:
                major_key = key
                continue
            
            # Check if this is a muscle key
            if not is_number(key):
                muscle_key = key
                continue
                 
            # Grab length of this feature, reassign
            data_length = data_lengths[major_key][muscle_key][key]
            data_copy[major_key][muscle_key][key] = np.copy(data_points[eye][:data_length, :])
            
            # Remove the data from copy
            data_points[eye] = np.delete(data_points[eye], np.s_[:data_length], axis=0)
            
            
            
    # return
    return data_copy


def translateSortedPoints(sorted_points, center):
    """
        Translates the given sorted points to center
    """
    
    # Loop through each z-index points
    for z_index, points in sorted_points.items():
    
        # Grab the points
        points = np.array(points)
        
        # Translate
        points = points - np.array(center)
        
        # Reassign
        sorted_points[z_index] = points
    
    return sorted_points


def globeCenter(muscles):
    """
        Translates all objects to globe center
    """
    
    # Grabe OD, OS globes
    OD_globe = muscles['Coronal']['OD']['Eyeball']
    OS_globe = muscles['Coronal']['OS']['Eyeball']
    globes = [OD_globe, OS_globe]

    # Iterate through each globe
    for globe in globes:
        
        # Grab center, plane, eye
        center = np.array(globe.center)
        eye = globe.eye
        
        # Apply transform to all other points
        for plane in muscles.keys():
            
            # Grab the eye
            eye_muscles = muscles[plane][eye]

            # Iterate through each muscle
            for muscle, muscle_object in eye_muscles.items():
                
                # Translate object if applicable
                if not muscle_object.object is None:
                    cmds.select(muscle_object.object)
                    cmds.xform(cp=True)
                    cmds.move(center[0] * -1, center[1]*-1, center[2]*-1, r=True, ws=True)
                    cmds.select(cl=True)
                    muscle_object.setCenter(queryCoords(muscle_object.object, ws=1))
                
                # Translate locators if applicable
                cmds.select(muscle_object.locators)
                cmds.xform(cp=True)
                cmds.move(center[0] * -1, center[1]*-1, center[2]*-1, r=True, ws=True)
                cmds.select(cl=True)
                
                # Adjust raw points too
                sorted_points = muscle_object.rawPoints
                for z_index, points in sorted_points.items():
                
                    # Grab the points
                    sorted_points[z_index] = np.array(points) - center
                    
                    
                    
                
    return muscles

def setRGBColor(ctrl, color = (1,1,1)):
    
    rgb = ("R","G","B")
    
    cmds.setAttr(ctrl + ".overrideEnabled",1)
    cmds.setAttr(ctrl + ".overrideRGBColors",1)
    
    for channel, color in zip(rgb, color):
        
        cmds.setAttr(ctrl + ".overrideColor%s" %channel, color)
    
    return
            
def getCoords(radius, phi, theta):
    x = radius * cos(phi) * sin(theta)
    y = radius * sin(phi) * sin(theta)
    z = radius * cos(theta)
    return [x,y,z]
        
        
def generateInsertions(muscles):
    """
        Generates muscle insertions points
    
    """
    
    # Grabe OD, OS globes
    OD_globe = muscles['Coronal']['OD']['Eyeball']
    OS_globe = muscles['Coronal']['OS']['Eyeball']
    globes = [OD_globe, OS_globe]
    
    # For each eye
    for globe in globes:
        
        # L/R flip for OD / OS
        flip = 1 if 'OD' in globe.eye else -1
        
        # Create group
        eye_group = cmds.group(em=True, name=globe.eye)
        
        # Grab radius
        radius = globe.radius
        
        # MR Insertion
        arclength = rectus_insertions['MR']
        phi = 0
        theta = arclength / radius
        pos = getCoords(radius, phi, theta)
        pos[0] *= flip
        mr_insert = cmds.spaceLocator(p=(pos[0], pos[1], pos[2]))
        cmds.xform(cp=True)
        cmds.parent(mr_insert[0], eye_group)
        setRGBColor(mr_insert[0], color=(1, 0, 0))
        
        # LR Insertion
        arclength = rectus_insertions['LR']
        phi = 0
        theta = arclength / radius * -1
        pos = getCoords(radius, phi, theta)
        pos[0] *= flip
        lr_insert = cmds.spaceLocator(p=(pos[0], pos[1], pos[2]))
        cmds.xform(cp=True)
        cmds.parent(lr_insert[0], eye_group)
        setRGBColor(lr_insert[0], color=(1, 0, 0))

        
        # SR Insertion
        arclength = rectus_insertions['SR']
        phi = 1.5708 # 90 degrees in radians
        theta = arclength / radius * -1
        pos = getCoords(radius, phi, theta)
        sr_insert = cmds.spaceLocator(p=(pos[0], pos[1], pos[2]))
        cmds.xform(cp=True)
        cmds.parent(sr_insert[0], eye_group)
        setRGBColor(sr_insert[0],color= (1, 0, 0))
        
        # IR Insertion
        arclength = rectus_insertions['IR']
        phi = 1.5708 # 90 degrees in radians
        theta = arclength / radius
        pos = getCoords(radius, phi, theta)
        ir_insert = cmds.spaceLocator(p=(pos[0], pos[1], pos[2]))
        cmds.xform(cp=True)
        cmds.parent(ir_insert[0], eye_group)
        setRGBColor(ir_insert[0], color=(1, 0, 0))
        
        # SO
        vert_length, horizontal_length = oblique_insertions['SO']
        if 'OD' in globe.eye:
            phi = 1.5708  - (horizontal_length / radius)# 90 degrees in radians
        else:
            phi = 1.5708  + (horizontal_length / radius)# 90 degrees in radians
            
        theta = vert_length / radius * -1
        pos = getCoords(radius, phi, theta)        
        so_insert = cmds.spaceLocator(p=(pos[0], pos[1], pos[2]))
        cmds.xform(cp=True)
        cmds.parent(so_insert[0], eye_group)
        setRGBColor(so_insert[0],color= (1, 0, 0))
        
        # IO
        vert_length, horizontal_length = oblique_insertions['IO']
        phi = horizontal_length / radius * -1
        theta = vert_length / radius * -1
        pos = getCoords(radius, phi, theta)
        pos[0] *= flip
        io_insert = cmds.spaceLocator(p=(pos[0], pos[1], pos[2]))
        cmds.xform(cp=True)
        cmds.parent(io_insert[0], eye_group)
        setRGBColor(io_insert[0], color=(1, 0, 0))
        
        globe.setMuscleInsertions({'MR':mr_insert[0], 'LR': lr_insert[0], 'SR': sr_insert[0], 'IR': ir_insert[0], 'SO':so_insert[0], 'IO':io_insert[0]})
    
    return muscles

def extendMuscles(muscles):
    """
        Extends muscle bsplines based on lengths.    
    """
    muscle_names = ['MR', 'LR', 'SR', 'IR', 'SO', 'IO']
        
    # Iterate through each eye
    for eye_muscles in muscles['Coronal'].values():
        
        for muscle_name in muscle_names:
            
            # Grab muscle
            muscle = eye_muscles[muscle_name]

            # Grab length
            muscle_len = cmds.arclen(muscle.bspline)
            
            # Compute amount of extension required
            req_extension = muscle_total_lengths[muscle_name] - muscle_len
            
            # Extend inplace
            cmds.extendCurve(muscle.bspline, rpo=True, s=1, em=0, d=req_extension)
            

    
    return muscles

    

def generateBSplines(muscle_objs, save_points=False):
    """
        Generates BSplines for all muscles.
        Interpolates from muscle insertion to bspline if required.
    """

    # Temporarily, we only care about simple stuff
    planes = ["Coronal"]
    muscle_names = ['MR', 'LR', 'SR', 'IR', 'SO', 'IO']
    
    # Go through each plane
    for plane in planes:
        
        # Go through each eye, grab muscles
        for eye, muscles in muscle_objs[plane].items():
            
            # Iterate through priority muscles
            for muscle in muscle_names:
                
                # Grab control points
                muscle_obj = muscles[muscle]
                muscle_group = muscle_obj.objectGroup
                cps = muscle_obj.controlPoints
                
                # Create group for control points
                ctrl_group = cmds.group(n='{}_{}_control_points'.format(muscle_obj.eye, muscle_obj.muscle_name), parent=muscle_group)
                
 
                
                # Create curve
                curve = cmds.curve(p=cps, d=3)
                bspline = cmds.fitBspline(n='{}_{}_bspline'.format(muscle_obj.eye, muscle_obj.muscle_name))[0]
                cmds.parent(curve, muscle_group)
                cmds.parent(bspline, muscle_group)
                cmds.hide(curve)
                
                # Plot control points
                for pt in cps:
                    loc = cmds.spaceLocator(p=(pt[0], pt[1], pt[2]))
                    cmds.parent(loc[0], ctrl_group)
                    
                # Save curve
                muscle_objs[plane][eye][muscle].bspline = bspline
                muscle_objs[plane][eye][muscle].curve = curve
                
                # Figure out why this bspline is not showing up
                if save_points:
                    
                    # Save control vertices for debugging
                    with open(r"C:\Users\admin\Dropbox\Maya Model Generator\Project\Maya_New\Maya2020\bin\plug-ins\Model_Generator\export/{}_{}_{}_cv.txt".format(plane, muscle_obj.eye, muscle), 'w') as out_file:
                        cvs = cmds.getAttr("{}.cv[*]".format(bspline))
                        for cv in cvs:
                            out_file.write("{}\t{}\t{}\n".format(cv[0], cv[1], cv[2]))
                    
                    # Save data points for debugging
                    with open(r"C:\Users\admin\Dropbox\Maya Model Generator\Project\Maya_New\Maya2020\bin\plug-ins\Model_Generator\export/{}_{}_{}_dp.txt".format(plane, muscle_obj.eye, muscle), 'w') as out_file:
                        for cp in cps:
                            out_file.write("{}\t{}\t{}\n".format(cp[0], cp[1], cp[2]))
    
    return muscle_objs


def closest_node(node, nodes, m):
    """
        Finds the closest point to node within nodes.
        Returns the closest point and the euclidean distance
    """
    deltas = nodes[:, :m] - node[:m]
    dist_2 = np.einsum('ij,ij->i', deltas, deltas)
    min_index = np.argmin(dist_2)
    return [nodes[min_index], dist_2[min_index], min_index]

def queryCoords(obj, ws=0):
    """
        Returns world space coordinates of the given object
    """
    return cmds.xform(obj,q=1,ws=ws, rp=1)



# Split this function. First: have a before function to generate the two sets
# of SO control points. Then call this function. 
def generateControlPointsSO(muscle, muscles, num_skip, sort_lambda, threshold=False, threshold_val=-1.0):
    """
        Generates control points for SO muscle
    
    """

    # Grab slice means
    center_means = []
    rawPoints = muscle.rawPoints
    
    # Group
    muscle_grp = cmds.group(em=True, name='SO_{}'.format(muscle.eye))
    
    # Grab globe properties
    globe = muscles['Coronal'][muscle.eye]['Eyeball']
    globe_radius = globe.radius
    globe_center = globe.center
    
    # Populate control points
    for i, (z_index, points) in enumerate(rawPoints.items()):               
        avged_point = tuple(np.average(points, axis=0))
        center_means.append(avged_point)
    
    # Sort
    center_means = sorted(center_means, key=sort_lambda)
    
        # Plot
    for i,p in enumerate(center_means):
        loc = cmds.spaceLocator(p=p, n='center_means_SO_{}_{}'.format(muscle.eye, i))[0]
        cmds.parent(loc, muscle_grp)
        
    # Fit a line to control points
    center_means = np.array(center_means)
    centers = center_means
    radii_points = rawPoints
    cps_mean = np.mean(center_means, axis=0)
    _, _, vv = np.linalg.svd(center_means - cps_mean)
    
    # Range to produce paramtric points on
    points_range = np.arange(-1 * globe_radius * 2, globe_radius*2, 0.25)
        
    # Find the points closest to globe
    print "finding points for: ", muscle.eye
    pt1, pt2, parametric_points = math_utils.findClosestPoints(cps_mean,vv[0], points_range, globe_radius, globe_center)
    for pt in parametric_points:
        loc = cmds.spaceLocator(p=pt)
        cmds.parent(loc[0], muscle_grp)
        
    
    # Convert both points to lat/lon
    geo_pt1, geo_pt2 = math_utils.pointCart2Geo(globe_radius, pt1), math_utils.pointCart2Geo(globe_radius, pt2)
    geo_pt2 = list(geo_pt2)
    
    # Compute bearing
    bearing = math_utils.bearing(geo_pt1, geo_pt2)

    # Grab length of sagittal segment
    center_means = center_means.tolist()
    dist = math_utils.pointDist(center_means[0], center_means[-1])
    
    # Subtract from length of tendon
    dist = tendon_lengths['SO'] - dist
    
    # Compute destination point
    SO_insert = math_utils.pointGeo2Cart(globe_radius, math_utils.destinationPoint(geo_pt1, bearing, dist, globe_radius))
    cps = center_means
    SO_insert = list(SO_insert)
#    SO_insert[1] = max(SO_insert[1], 0)
    
    # Remove some of those center points
    if muscle.plane == 'Sagittal':
        
        num_to_skip = 1
    else:
        num_to_skip = 3
    ctr = num_to_skip
    pts_copy = []
    for i in range(len(cps)):
        if ctr == num_to_skip:
            ctr = 0
            pts_copy.append(cps[i])
        else:
            ctr += 1
    cps = pts_copy
    cps.append(tuple(pt1))
    
    # Normalize
    pt1 = np.array(pt1) / globe_radius
    SO_insert = np.array(SO_insert) / globe_radius
    
    # Interpolate points on sphere
    pts = math_utils.interpolate(pt1, SO_insert, steps=20) * globe_radius
    for i in range(len(pts)):
        cps.append(pts[i])
        
    # Plot CPS
    for i, cp in enumerate(cps):
        loc = cmds.spaceLocator(p=cp)
        cmds.parent(loc, muscle.objectGroup)
        
    # Grab the coronal verison of this muscle for logging
    logging_muscle = muscles['Coronal'][muscle.eye][muscle.muscle_name]
    
    # Append these points if we already generated others
    if not len(logging_muscle.controlPoints) == 0:
        cps = np.append(cps, logging_muscle.controlPoints, axis=0)
        logging_muscle.rawPoints.update(muscle.rawPoints)

    if not len(logging_muscle.radii_points.keys()) == 0:
        for key in logging_muscle.radii_points.keys():
            radii_points[key] = logging_muscle.radii_points[key]

    if not len(logging_muscle.centers) == 0:
        centers = np.append(centers, logging_muscle.centers, axis=0)   
    
    # Set control points
    logging_muscle.setControlPoints(cps)
    logging_muscle.centers = centers
    logging_muscle.radii_points = radii_points
    return


def generateMeans(muscle, sort_lambda):
    """
    Computes the mean of each slice.
    """
    
    # Grab slice means
    center_means = []
    rawPoints = muscle.rawPoints
    
    # Populate control points
    for i, (z_index, points) in enumerate(rawPoints.items()):               
        avged_point = tuple(np.average(points, axis=0))
        center_means.append(avged_point)
    
    # Sort
    center_means = sorted(center_means, key=sort_lambda)
    return center_means
    

def getPositionsFromGroup(group):
    """
    Returns the location of all objects in a group
    """
    
    cmds.select(cl=True)
    cmds.select(group, hi=True)
    names = cmds.ls(sl=True)
    names.pop(0)
    
    locs = []
    for name in names:
        locs.append(cmds.objectCenter(name))
    
    return locs

def deletePointsFromGroup(group):
    """
    Deletes all points in a given group
    """
    cmds.select(cl=True)
    cmds.select(group, hi=True)
    names = cmds.ls(sl=True)
    cmds.select(names[0], toggle=True)
    names.pop(0)
    
    cmds.delete()
    return

def generateControlPointsSOSplit(muscle, muscles, num_skip, muscle_grps, sort_lambda, threshold=False, threshold_val=-1.0):
    """
    Generates control points for a split SO
    """
    
    # Group
    muscle_grp = cmds.group(em=True, name='SO_{}'.format(muscle.eye))
    
    # Grab globe properties
    globe = muscles['Coronal'][muscle.eye]['Eyeball']
    globe_radius = globe.radius
    globe_center = globe.center
    
    # Iterate through each muscle grp
    for grp in muscle_grps:
        
        # Grab all locator positions
        center_means = getPositionsFromGroup(grp)
        
        # Sort 
        center_means = sorted(center_means, key=sort_lambda)
            
        # Fit a line to control points
        center_means = np.array(center_means)
        cps_mean = np.mean(center_means, axis=0)
        _, _, vv = np.linalg.svd(center_means - cps_mean)
        
        # Range to produce paramtric points on
        points_range = np.arange(-1 * globe_radius, globe_radius, 0.25)
            
        # Find the points closest to globe
        pt1, pt2, parametric_points = math_utils.findClosestPoints(cps_mean,vv[0], points_range, globe_radius, globe_center)
        for pt in parametric_points:
            loc = cmds.spaceLocator(p=pt)
            cmds.parent(loc[0], muscle_grp)
            
        
        # Convert both points to lat/lon
        geo_pt1, geo_pt2 = math_utils.pointCart2Geo(globe_radius, pt1), math_utils.pointCart2Geo(globe_radius, pt2)
        geo_pt2 = list(geo_pt2)
        
        # plot
        loc = cmds.spaceLocator(p=pt1, n='geo_pt1')
        cmds.parent(loc[0], muscle_grp)
        loc = cmds.spaceLocator(p=pt2, n='geo_pt2')
        cmds.parent(loc[0], muscle_grp)
        
        # Compute bearing
        bearing = math_utils.bearing(geo_pt1, geo_pt2)
    
        # Grab length of sagittal segment
        center_means = center_means.tolist()
        dist = math_utils.pointDist(center_means[0], center_means[-1])
        
        # Subtract from length of tendon
        dist = tendon_lengths['SO'] - dist
        
        # Compute destination point
        SO_insert = math_utils.pointGeo2Cart(globe_radius, math_utils.destinationPoint(geo_pt1, bearing, dist, globe_radius))
        cps = center_means
        
        # Remove some of those center points
        num_to_skip = 1
        ctr = num_to_skip
        pts_copy = []
        for i in range(len(cps)):
            if ctr == num_to_skip:
                ctr = 0
                pts_copy.append(cps[i])
            else:
                ctr += 1
        cps = pts_copy
        cps.append(tuple(pt1))
        
        # Normalize
        pt1 = np.array(pt1) / globe_radius
        SO_insert = np.array(SO_insert) / globe_radius
        
        # Interpolate points on sphere
        pts = math_utils.interpolate(pt1, SO_insert, steps=20) * globe_radius
        for i in range(len(pts)):
            cps.append(pts[i])
        
        
        # Remove previous points
        deletePointsFromGroup(grp)
        
        # Plot CPS
        for i, cp in enumerate(cps):
            loc = cmds.spaceLocator(p=cp)
            cmds.parent(loc, grp)
            
        # Grab the coronal verison of this muscle for logging
        logging_muscle = muscles['Coronal'][muscle.eye][muscle.muscle_name]
        
        # Append these points if we already generated others
        if not len(logging_muscle.controlPoints) == 0:
            cps = np.append(cps, logging_muscle.controlPoints, axis=0)
                    
        # Set control points
        logging_muscle.setControlPoints(cps)
    
    return

    
def generateSliceControlPoints(muscles, desired_muscles, desired_plane, num_skip, SOP=False, threshold=False, threshold_val=-1.0):
    """
        Generates control points for all slices of a given muscle in the
        desired plane.
        
        Can control how many slices are included by num_skip.
    """
    # How to sort points
    sort_lambda = {'OS':{'Sagittal':lambda cp: cp[0], 'Coronal':lambda cp: cp[2]},
                   'OD':{'Sagittal':lambda cp: -cp[0], 'Coronal':lambda cp: cp[2]}}
    
    # Go through each eye muscle set, only grab the desired_muscle
    for eye in muscles[desired_plane].values():
        
        # Go through all desired muscles
        for muscle_name in desired_muscles:
            
            muscle = eye[muscle_name]

            
            # Special condition for SO
            if muscle_name == 'SO' and desired_plane == 'Sagittal' and SOP:
                
                # Grab center means
                center_means = generateMeans(eye[muscle_name], sort_lambda[muscle.eye][muscle.plane])
                
                # Make a group
                so1_grp = cmds.group(n='{}_{}_SO1'.format(muscle.plane, muscle.eye), p=muscle.objectGroup, em=True)
#                print so1_grp, muscle.objectGroup
                # Make locators, add to group
                for mean in center_means:
                    loc = cmds.spaceLocator(p=mean)[0]
                    cmds.parent(loc, so1_grp)
                
                # Grab pivot point
                pivot_loc = center_means[0]
                
                # Duplicate
                so2_grp = cmds.duplicate(so1_grp)[0]
                so2_grp = cmds.rename(so2_grp, so2_grp.replace('1', '2'))

                # Rotate each one                
                cmds.select(cl=True)
                cmds.select(so1_grp)
                cmds.rotate('-10deg',  y=True, p=pivot_loc, os=True)
                cmds.select(cl=True)
                cmds.select(so2_grp)
                cmds.rotate('10deg',  y=True, p=pivot_loc, os=True)
                
                # Generate control points
                generateControlPointsSOSplit(eye[muscle_name], muscles, num_skip, [so1_grp, so2_grp], sort_lambda[muscle.eye][muscle.plane], threshold=False, threshold_val=-1.0)
                
            elif muscle_name == 'SO' and desired_plane == 'Sagittal':
                generateControlPointsSO(eye[muscle_name], muscles, num_skip, sort_lambda[muscle.eye][muscle.plane], threshold=False, threshold_val=-1.0)
                continue
            
            # Grab the muscle, points
            rawPoints = muscle.rawPoints
            cps = []
            
            # Special case for SO coronal
            if muscle_name == 'SO' and desired_plane == 'Coronal':
                skip_amount = 3
                
                # Populate control points
                slice_ctr = skip_amount
                radii_points = {}
                for i, (z_index, points) in enumerate(rawPoints.items()):               
                    if slice_ctr == skip_amount or i == len(rawPoints.items()) - 1:
                        radii_points[z_index] = points
                        avged_point = tuple(np.average(points, axis=0))
                        cps.append(avged_point)
                        slice_ctr = 0
                    else:
                        slice_ctr += 1
            
            else:
            # Populate control points
                slice_ctr = num_skip
                radii_points = {}
                for i, (z_index, points) in enumerate(rawPoints.items()):               
                    if slice_ctr == num_skip or i == len(rawPoints.items()) - 1:
                        radii_points[z_index] = points
                        avged_point = tuple(np.average(points, axis=0))
                        cps.append(avged_point)
                        slice_ctr = 0
                    else:
                        slice_ctr += 1
            # Sort
            cps = sorted(cps, key=sort_lambda[muscle.eye][muscle.plane])
            cps = np.array(cps)
            centers = cps
#                        
            # Threshold undesired slices if needed
            if threshold and not muscle_name == 'SO':
                to_delete = []
                for i, cp in enumerate(cps):
                    if cp[2] > threshold_val:
                        to_delete.append(i)
                cps = np.delete(cps, to_delete, axis=0)    
            

            print "Radii: ", len(radii_points), " cps: ", len(cps)
            # Grab the coronal verison of this muscle for logging
            logging_muscle = muscles['Coronal'][muscle.eye][muscle.muscle_name]
            print '{}_{}_{}_centroids'.format(logging_muscle.plane, logging_muscle.eye, muscle_name)
            # Plot CPS
            centroid_group = cmds.group(n='{}_{}_{}_centroids'.format(logging_muscle.plane, logging_muscle.eye, muscle_name), p=logging_muscle.objectGroup, em=False)
            for i, cp in enumerate(cps):
                loc = cmds.spaceLocator(p=cp)
                cmds.parent(loc, centroid_group)
                

            # Append these points if we already generated others
            if not len(logging_muscle.controlPoints) == 0:
                cps = np.append(cps, logging_muscle.controlPoints, axis=0)
#            logging_muscle.rawPoints.update(muscle.rawPoints)
            
            if not len(logging_muscle.radii_points.keys()) == 0:
                for key in logging_muscle.radii_points.keys():
                    radii_points[key] = logging_muscle.radii_points[key]
        
            if not len(logging_muscle.centers) == 0:
                centers = np.append(centers, logging_muscle.centers, axis=0)   
            
            # Set control points
            logging_muscle.setControlPoints(cps)
            logging_muscle.centers = centers
            logging_muscle.radii_points = radii_points
            
    return muscles

def interpolateControlPoints(muscles, desired_muscles, desired_plane, globes,
                             sort_lambda, delete_closest_point=False):
    """
        Interpolates points between a muscle's control points and the insertion
        on the globe.
    
    """
    

    
    # Go through each eye
    for eye_name, eye_muscles in muscles['Coronal'].items():
        
        # Grab globe, insertions
        globe = globes[eye_name]
        insertions = globe.muscleInsertions
        globe_center = globe.center
        radius = globe.radius
        
        # Go through each muscle
        for muscle_name in desired_muscles:
            
            # Grab cps
            muscle = eye_muscles[muscle_name]
            cps = muscle.controlPoints
            
            if len(cps) == 0:
                continue
            
            
            # Find the control point closest to muscle insertion
            insert = np.array(queryCoords(insertions[muscle_name]))
            
            # Find closest point, remove from consideration
            closest_point, dist, idx = closest_node(insert, cps, cps.shape[1])
            if delete_closest_point:
                cps = np.delete(cps, idx, axis=0)
            
            # Find points of insertion on the globe from closest muscle cp
            intersection1, intersection2 = findIntersectionPoints(insert, closest_point, globe_center, radius)
            cps = np.append(cps, np.array([np.array(intersection1)]), axis=0)
            intersection = cmds.spaceLocator(p=(intersection1[0], intersection1[1], intersection1[2]), n="{}_intersect_point".format(muscle.muscle_name))
            cmds.parent(intersection, muscle.objectGroup)
            
            # Interpolate points between insertion and intersection
            intersection1 = np.array(intersection1) / radius # Unit distance
            insertion = insert / radius
            pts = math_utils.interpolate(intersection1, insertion) * radius
            
            # Delete some of those interpolated points
            num_to_skip = 3
            ctr = 3
            pts_copy = []
            for i in range(len(pts)):
                if ctr == num_to_skip:
                    ctr = 0
                    pts_copy.append(pts[i])
                else:
                    ctr += 1
            
            # Convert to numpy, plot and group
            pts = pts_copy
            pts = np.array(pts)
            grp = cmds.group(parent=muscle.objectGroup, name="{}_{}_interp_pts".format(globe.eye, muscle.muscle_name))
            cps = np.append(cps, pts, axis=0)
            for pt in pts:
                loc = cmds.spaceLocator(p=(pt[0], pt[1], pt[2]))
                cmds.parent(loc, grp)
                
            # Set control points again
            muscle.setControlPoints(cps)    
    
    return muscles

def generateControlPoints(muscles, SOP=False):
    """
        Generates muscle bspline control points
    """

    # Grab globes
    OD_globe = muscles['Coronal']['OD']['Eyeball']
    OS_globe = muscles['Coronal']['OS']['Eyeball']
    globes = {'OD':OD_globe, 'OS':OS_globe}
    
    # Plane/muscle pairs
    plane_pairs = {'Coronal':['MR', 'LR', 'SR', 'IR', 'SO'],
                   'Sagittal':['SO', 'IO']}

    # Slice skip amount
    num_skip = 1
    slice_threshold = -2.0
    
        
    # Iterate once to generate control points
    for plane, muscle_names in plane_pairs.items():

        # Switch based on plane
        if plane == 'Coronal':
            threshold = True
            threshold_val = slice_threshold
        elif plane == 'Sagittal':
            threshold = False
            threshold_val = -1

        
        # Generate slice control points
        muscles = generateSliceControlPoints(muscles, muscle_names, plane, 
                                             num_skip, threshold=threshold,
                                             threshold_val=threshold_val,
                                             SOP=SOP)


    # Remove sagittal from interpolation
    plane_pairs['Sagittal'].remove('SO')
    plane_pairs['Coronal'].remove('SO')
    
    # Iterate again to interpolate. This must be separate due to one muscle 
    # being present multiple times.
    for plane, muscle_names in plane_pairs.items():

        # Switch based on plane
        if plane == 'Coronal':
            sort_lambda = lambda cp: cp[2]
        else:
            sort_lambda = lambda cp: cp[0]

        
        muscles = interpolateControlPoints(muscles, muscle_names, plane, globes,
                                           sort_lambda, delete_closest_point=True)
        
    return muscles

def filterPoints(muscles):
    """
        Filters out points not needed.
        
        Muscles filtered:
            - SO
            
    """
    
    # Threshold to filter out SO points from sagittal plane. This is a multiplier w.r.t globe radius,
    # so any points posterior to the globe by this multipler will be removed.
    SO_filter_threshold_multiplier = -0.5
    
    # Go through each plane
    for plane, eyes in muscles.items():
        
        for eye in eyes.values():
            
            # Go through each eye
            for muscle in eye.values():
                if muscle.muscle_name == 'SO' and plane == 'Sagittal':
                    # Filter out SO bad points
                    globe_radius = muscles['Coronal'][muscle.eye]['Eyeball'].radius
                    limit = SO_filter_threshold_multiplier * globe_radius
                    for z_index, points in muscle.rawPoints.items():    
                        points = np.array(points)
                        points = points[points[:, 2] > limit]
#                        muscle.rawPoints[z_index] = points
                        x_points_to_include = np.logical_and(points[:, 0] < globe_radius, points[:, 0] > globe_radius*-1)
                        if not np.any(x_points_to_include):
                            del muscle.rawPoints[z_index]
                            continue
                        else:
                            points = points[x_points_to_include]
                            muscle.rawPoints[z_index] = points


    return muscles

def findIntersectionPointsVector(direction_vec, cps_mean, globe_center, globe_radius):
    """
        Computes the points of insertion (if any) of the line segment
        by the direction vector and the sphere at center globe_center and radius.
        
        Math from: http://www.ambrsoft.com/TrigoCalc/Sphere/SpherLineIntersection_.htm
    
    """
    # Quadratic line equation: at^2 + bt + c = 0
    # a = xt * xt + yt * yt + zt * zt
    # b = 2 * (xt * (x1 - xc) + yt * (y1 - yc) + zt * (z1 - zc))
    # c = (xc - x1) * (xc - x1) + (yc - y1) * (yc - y1) + (zc - z1) * (zc - z1) - r * r
    
    # Get line equation variables
    x1, y1, z1 = direction_vec
    xt, yt, zt = cps_mean
    a = (xt**2) + (yt**2) + (zt**2)
    b = 2 * (xt * (x1 - globe_center[0]) + yt * (y1 - globe_center[1]) + zt * (z1 - globe_center[2]))
    c = ((globe_center[0] - x1)**2) + ((globe_center[1] - y1)**2) + ((globe_center[2] - z1)**2) - globe_radius**2
    
    distance = (b**2) - (4*a*c)
    pt1 = None
    pt2 = None
    if distance > 0:
        # There is intersection
        t1 = (-b + sqrt(distance)) / (2*a)
        t2 = (-b - sqrt(distance)) / (2*a)
        pt1 = (x1 + t1 * xt, y1 + t1 * yt, z1 + t1 * zt)
        pt2 = (x1 + t2 * xt, y1 + t2 * yt, z1 + t2 * zt)
    elif distance == 0:
        # Line is tangent
        t1 = -b / (2 * a)
        pt1 = (x1 + t1 * xt, y1 + t1 * yt, z1 + t1 * zt)
    else:
        # No intersection
        pt1 = None
        pt2 = None
    
    return pt1, pt2



def findIntersectionPoints(pt1, pt2, globe_center, globe_radius):
    """
        Computes the points of insertion (if any) of the line segment
        by pt1 -> pt2 and the sphere at center globe_center and radius.
        
        Math from: http://www.ambrsoft.com/TrigoCalc/Sphere/SpherLineIntersection_.htm
    """
    
    # Quadratic line equation: at^2 + bt + c = 0
    # a = xt * xt + yt * yt + zt * zt
    # b = 2 * (xt * (x1 - xc) + yt * (y1 - yc) + zt * (z1 - zc))
    # c = (xc - x1) * (xc - x1) + (yc - y1) * (yc - y1) + (zc - z1) * (zc - z1) - r * r
    
    # Get line equation variables
    x1 = pt1[0]
    y1 = pt1[1]
    z1 = pt1[2]
    xt = pt2[0] - x1
    yt = pt2[1] - y1
    zt = pt2[2] - z1
    a = (xt**2) + (yt**2) + (zt**2)
    b = 2 * (xt * (x1 - globe_center[0]) + yt * (y1 - globe_center[1]) + zt * (z1 - globe_center[2]))
    c = ((globe_center[0] - x1)**2) + ((globe_center[1] - y1)**2) + ((globe_center[2] - z1)**2) - globe_radius**2
    
    distance = (b**2) - (4*a*c)
    pt1 = None
    pt2 = None
    if distance > 0:
        # There is intersection
        t1 = (-b + sqrt(distance)) / (2*a)
        t2 = (-b - sqrt(distance)) / (2*a)
        pt1 = (x1 + t1 * xt, y1 + t1 * yt, z1 + t1 * zt)
        pt2 = (x1 + t2 * xt, y1 + t2 * yt, z1 + t2 * zt)
    elif distance == 0:
        # Line is tangent
        t1 = -b / (2 * a)
        pt1 = (x1 + t1 * xt, y1 + t1 * yt, z1 + t1 * zt)
    else:
        # No intersection
        pt1 = None
        pt2 = None
    
    return pt1, pt2


def createOriginPlane():
    """
        Defines the origin plane for all muscles
    
    """
    
    plane = cmds.nurbsPlane(w=50, n='originPlane')[0]
    cmds.select(plane)
    cmds.rotate('90deg', y=True)
    cmds.move(-34, z=True)
    cmds.select(cl=True)
    return plane

    
        
def splitMuscleTendon(muscles):
    """
        Splits muscle bsplines into tendon, muscle
    """

    # Temporarily, we only care about simple stuff
    planes = ["Coronal"]
    muscle_names = ['MR', 'LR', 'SR', 'IR', 'SO', 'IO']
    
    # Go through each plane
    for plane in planes:
        
        # Go through each eye, grab muscles
        for eye, eye_muscles in muscles[plane].items():
            
            # Iterate through priority muscles
            for muscle_name in muscle_names:
                
                # Grab object
                muscle_obj = eye_muscles[muscle_name]
                
                # Grab muscle bspline
                bspline = muscle_obj.bspline
                
                # Length from muscle origin to split at
                split_length = muscle_total_lengths[muscle_name] - tendon_lengths[muscle_name]
                if muscle_name == 'IO':
                    split_length = -1
                
                # Get the split parameter based on length
                split_param = curve_utils.getParamFromLength(bspline, split_length)
                
                # Split into muscle, tendon
                muscle_bspline, tendon_bspline = cmds.detachCurve(bspline, p=split_param, rpo=False)
                
                # Rename
                muscle_bspline = cmds.rename(muscle_bspline, '{}_{}_muscle_bspline'.format(muscle_obj.eye, muscle_name))
                tendon_bspline = cmds.rename(tendon_bspline, '{}_{}_tendon_bspline'.format(muscle_obj.eye, muscle_name))
                
                # Group
                cmds.parent(muscle_bspline, muscle_obj.objectGroup)
                cmds.parent(tendon_bspline, muscle_obj.objectGroup)
                
                # Recolor?
                curve_utils.colorize(muscle_bspline, muscle_color)
                curve_utils.colorize(tendon_bspline, tendon_color)
                
                # Set 
                muscle_obj.setMuscleBspline(muscle_bspline)
                muscle_obj.setTendonBspline(tendon_bspline)
                
                # hide
                cmds.hide(bspline)
                
    
    return muscles

def generatePulleyLocations(muscles):
    """
        Generates pulley locations for simple muscles.
    """
    
    muscle_names = ['MR', 'LR', 'SR', 'IR', 'SO']
    
    # Iterate through each eye
    for eye_muscles in muscles['Coronal'].values():
        
        for muscle_name in muscle_names:
            
            # Grab muscle
            muscle = eye_muscles[muscle_name] 
            
            # Query the muscle/tendon split point
            pulley_loc = cmds.getAttr('{}.cv[0]'.format(muscle.tendonBspline))[0]
            
            # Create locator
            locator = cmds.spaceLocator(p=pulley_loc, n='{}_{}_pulley'.format(muscle.eye, muscle_name))[0]
            setRGBColor(locator, color=(0.4, 0.2, 1))
            cmds.parent(locator, muscle.objectGroup)
            
            # Set
            muscle.setPulleyLocation(locator)
            
    return

            

def extendToOrigin(muscles):
    
    """
        Extends muscle origins to origin Plane
    """
    
    muscle_names = ['MR', 'LR', 'SR', 'IR', 'SO']
    plane_z = -34
    
    # Iterate through each eye
    for eye, eye_muscles in muscles['Coronal'].items():
                
        for muscle_name in muscle_names:
            
            # Grab muscle
            muscle = eye_muscles[muscle_name] 
                            
            # Query for origin point on muscle
            origin_point = cmds.getAttr('{}.cv[0]'.format(muscle.muscleBspline))[0]

            print origin_point, muscle_name
            # Flat distance value
#            dist = 15
            dist = math_utils.pointDist(origin_point, (0, 0, plane_z))
            
            
            # Extend
            if origin_point[2] > plane_z:
                # Extend based on distance
                
                
                cmds.extendCurve(muscle.muscleBspline, rpo=True, s=1, em=0, d=dist, et=2)
                print "here with: ", muscle_name
                loc1 = cmds.spaceLocator(p=origin_point, n='Origin')[0]
                cmds.parent(loc1, muscle.objectGroup)
                
                # Now grab the new origin
                new_origin = cmds.getAttr('{}.cv[0]'.format(muscle.muscleBspline))[0]
            
                loc = cmds.spaceLocator(p=origin_point, n='{}_{}_{}_origin_point1'.format(muscle.eye, muscle.plane, muscle_name))[0]
                cmds.parent(loc, muscle.objectGroup)
            else:
                new_origin = cmds.getAttr('{}.cv[6]'.format(muscle.muscleBspline))[0]
            
            # Create vector from old -> new
            vec = (np.array(new_origin) - np.array(origin_point)).tolist()
            
            # Grab intersection point
            intersect = surface_utils.intersect('originPlane', origin_point, vec)
            
            # Plot point
            loc = cmds.spaceLocator(p=intersect, n='{}_{}_origin_point'.format(muscle.eye, muscle_name))[0]
            cmds.parent(loc, muscle.objectGroup)
            
            new_origin = cmds.getAttr('{}.cv[0]'.format(muscle.muscleBspline))[0]

            # Grab this new distance
            dist = math_utils.pointDist(intersect, new_origin)
                        
            # Get the split parameter based on length
            split_param = curve_utils.getParamFromLength(muscle.muscleBspline, dist)
            
            # Split into muscle, tendon
            extra_bspline, muscle_bspline = cmds.detachCurve(muscle.muscleBspline, p=split_param, rpo=False)
            
            # Delete old
            cmds.delete(muscle.muscleBspline)
            
            # Rename
            muscle_bspline = cmds.rename(muscle_bspline, '{}_{}_muscle_bspline'.format(muscle.eye, muscle_name))
            
            # Group
            cmds.parent(muscle_bspline, muscle.objectGroup)
            
            # Recolor?
            curve_utils.colorize(muscle_bspline, muscle_color)
            
            # Set 
            muscle.setMuscleBspline(muscle_bspline)
            
    return muscles


def generateRadius(muscles):
    """
    Computes slice-radius for each muscle of interest.
    """
    muscle_names = ['MR', 'LR', 'SR', 'IR', 'SO', 'IO']

    # Iterate through each eye
    for eye, eye_muscles in muscles['Coronal'].items():
        
        for muscle_name in muscle_names:
            
            # Grab muscle
            muscle = eye_muscles[muscle_name]
            
            # Grab raw points
            pts = muscle.rawPoints
            radii = {}
            max_radius = 0
            min_radius = 10000
            
            # Iterate through each slice, compute average radius
            for z_index, points in pts.items():    
                
                # Grab points, computer center and distance from center (radius)
                points = np.array(points)
                center = np.mean(points)
                radius = np.mean(np.power(points - center, 2))
                
                radii[z_index] = [center, radius]
                if radius < min_radius:
                    min_radius = radius
                
            muscle.radius = min_radius
            muscle.slice_radii = radii
    
    return muscles

                
def generateBiomechanicalProperties(muscles):
    
    muscle_names = ['MR', 'LR', 'SR', 'IR', 'SO', 'IO']

    # Iterate through each eye
    for eye, eye_muscles in muscles['Coronal'].items():
        
        for muscle_name in muscle_names:
            
            # Grab muscle
            muscle = eye_muscles[muscle_name]
            
            #Take area, devide by pi and take the sqrt
            # Grab radius, compute other properties
            
            muscleR = muscle.radius
            muscle.PCSA = muscleR ** 2 * pi
            muscle.musclePassive = 0.8036 / muscle.PCSA
            muscle.muscleActive = 0.9878 / muscle.PCSA
            muscle.tendonPassive = 80.36 / muscle.PCSA
            muscle.muscleDamping = 2e-3
            muscle.tendonDamping = 2e-3
            muscle.initialStrain = -0.01
            muscle.activeFL = "data/forces/active.xml"
            muscle.passiveFL = "data/forces/passive.xml"
            muscle.tendonFL = "data/forces/passive.xml"
            muscle.add = True
            
            
    
    return muscles


def addBiomechanicalProperties(muscleData, muscleR):
    muscleData['PCSA'] = muscleR ** 2 * pi
    muscleData['musclePassive'] = 0.8036 / muscleData['PCSA']
    muscleData['muscleActive'] = 0.9878 / muscleData['PCSA']
    muscleData['tendonPassive'] = 80.36 / muscleData['PCSA']
    muscleData['muscleDamping'] = 2e-3
    muscleData['tendonDamping'] = 2e-3
    muscleData['initialStrain'] = -0.01
    muscleData['activeFL'] = "data/forces/active.xml"
    muscleData['passiveFL'] = "data/forces/passive.xml"
    muscleData['tendonFL'] = "data/forces/passive.xml"
    muscleData['add'] = True
    return muscleData
                
def generateOrbitalBspline(muscles):
    """
    Generates the orbital bsline for each simple muscle. Simply a linear
    interpolation from pulley locations to origin locations
    """

    muscle_names = ['MR', 'LR', 'SR', 'IR']

    # Iterate through each eye
    for eye, eye_muscles in muscles['Coronal'].items():
        
        for muscle_name in muscle_names:
            
            # Grab muscle
            muscle = eye_muscles[muscle_name]
            
            # Grab origin and pulley locations
            origin_loc = cmds.objectCenter('{}_{}_origin_point'.format(muscle.eye, muscle_name))
            pulley_loc = cmds.objectCenter('{}_{}_pulley'.format(muscle.eye, muscle_name))
            
            # Create vector from pulley -> origin
            vec =  np.array(pulley_loc) - np.array(origin_loc)
            
            # Generate 6 points along this line
            globe_layer_pts = []
            num_points = 6
            prev_pt = np.array(origin_loc)
            orbital_group = cmds.group(n='{}_{}_{}_orbital_points'.format(muscle.plane, muscle.eye, muscle_name), parent=muscle.objectGroup, em=True)
            loc = cmds.spaceLocator(p=prev_pt, n='{}_{}_orbital_point_0'.format(muscle.eye, muscle_name))[0]
            globe_layer_pts.append(origin_loc)
            cmds.parent(loc, orbital_group)
            diff = vec / num_points
            for i in range(num_points):
                new_pt = prev_pt + diff
                globe_layer_pts.append(new_pt)
                prev_pt = new_pt
                loc = cmds.spaceLocator(p=prev_pt, n='{}_{}_orbital_point_{}'.format(muscle.eye, muscle_name, i+1))[0]
                cmds.parent(loc, orbital_group)
            
            cmds.hide(orbital_group)
            # Create bspline
            curve = cmds.curve(p=globe_layer_pts, d=3)
            bspline = cmds.fitBspline(n='{}_{}_orbital_bspline'.format(muscle.eye, muscle.muscle_name))[0]
            cmds.parent(curve, muscle.objectGroup)
            cmds.parent(bspline, muscle.objectGroup)
            cmds.hide(curve)
            cmds.hide(bspline)
            muscle.globalLayerPoints = globe_layer_pts
            
    return muscles

def splitSOMuscle(muscles):
    
    # Iterate through each eye
    for eye, eye_muscles in muscles['Coronal'].items():
        
        # Grab the pivot point, this is the point where muscle+tendon connect.
        # Grab the location of this.
        muscle = eye_muscles['SO']
        pulley_name = muscle.eye + '_SO_pulley'
        bspline_name = muscle.eye + '_SO_tendon_bspline'
        new_name = muscle.eye + '_SO_tendon_bspline2'
        
        pivot_loc = cmds.objectCenter(pulley_name)
        
        # Make a copy of the tendon bspline.
        cmds.select(bspline_name)
        duplicated_obj = cmds.duplicate()[0]
        cmds.select(cl=True)
        
        # Select the tendon bespline to be copied.
        cmds.select(duplicated_obj)
        
        # Rotate
        cmds.rotate('30deg',  y=True, p=pivot_loc, os=True)
        
        cmds.select(cl=True)
    
    return muscles


def polygon_area(poly):
    if isinstance(poly, list):
        poly = np.array(poly)
    edges = poly[1:] - poly[0:1]
    cross_product = np.cross(edges[:-1],edges[1:], axis=1)
    area = np.linalg.norm(cross_product, axis=1)/2
    return sum(area)

def getSliceRadii(points):
    
    radii = []
    for z_index, points in points.items():    
        points = np.array(points)
        radii.append(sqrt(polygon_area(points) / pi))
    return np.array(radii)

def getCentroids(points):
    
    centroids = []
    for z_index, points in points.items():    
        centroids.append(np.mean(points, axis=0).tolist())
    return np.array(centroids)

def getSlices(bspline, muscle):
    """
        brute force find muscle slices closest to bspline control vertices
    """
    
    # Grab control points
    control_points = cmds.getAttr("{}.controlPoints[*]".format(bspline))
    
    # Grab first/last 
    start_pt, end_pt = control_points[0], control_points[-1]
    
    # For each point, find the closest slice based on centroid
    slices = muscle.rawPoints.values()
    z_indices = muscle.rawPoints.keys()
    centroids = []
    for slice_points in slices:
        centroids.append(tuple(np.average(slice_points, axis=0)))
    centroids = np.array(centroids)
    
    # Get closest indices
    closest_start_slice_idx = np.argmin(np.linalg.norm(np.array(start_pt) - centroids,
                                                       axis=1))
    closest_end_slice_idx = np.argmin(np.linalg.norm(np.array(end_pt) - centroids,
                                                       axis=1))
    
    if closest_start_slice_idx >= closest_end_slice_idx:
        temp = closest_end_slice_idx
        closest_end_slice_idx = closest_start_slice_idx
        closest_start_slice_idx = temp
    
    out_slices = {}
    for i in range(closest_start_slice_idx, closest_end_slice_idx + 1):
        out_slices[str(z_indices[i])] = slices[i]
    
    return out_slices

def getCurvePoints(curve):
    points = []
    for i in range(len(cmds.getAttr(curve + ".cv[*]"))):
        points.append(cmds.pointPosition("{}.cv[{}]".format(curve, str(i))))
    return points

        
def export(muscles, out_file_path="./exported_data.json"):
    """
    Generates a json to house the whole scene.
    """
    
    data = {}
    muscle_names = ['SO', 'IO', 'IR', 'SR', 'LR', 'MR', 'Eyeball']
    orbital_names =  ['MR', 'LR', 'SR', 'IR']
    scale = 1000
    rotation_axis = 'x'
    rotation_degrees = 0
    
    # Select all
    cmds.select(all=True)
    cmds.rotate('90deg', x=True, ws=True)
    cmds.rotate('180deg', y=True, ws=True)
    cmds.select(cl=True)
    
    # Iterate through each eye
    for eye, eye_muscles in muscles['Coronal'].items():
        data[eye] = {'muscles':{}}
        
        # Iterate through each muscle
        for muscle_name in muscle_names:
            
            # Grab the muscle
            muscle = eye_muscles[muscle_name]
            
            # Check if eyeball
            if muscle_name == 'Eyeball':
                
                # Rebuild spans
                cmds.rebuildSurface('{}_{}_Eyeball'.format('Coronal', eye), 
                                    spansU=15, spansV=10)
                # Grab eyeball control points
                eyeball_cvs, spansU, spansV = getSurfacePoints('{}_{}_Eyeball'.format('Coronal', eye))
                eyeball_cvs = np.array(eyeball_cvs) / scale
                
                # Duplicate first and final rings
                first_ring = eyeball_cvs[:10]
                final_ring = eyeball_cvs[-10:]
                eyeball_cvs = np.vstack((first_ring,
                                         first_ring,
                                         eyeball_cvs,
                                         final_ring,
                                         final_ring))
                #eyeball_cvs = rotatePoints(rotation_axis, rotation_degrees, eyeball_cvs)
                # Grab radius
                radius = muscle.radius / scale
                
                eye_data = {'eyeball_cvs':eyeball_cvs.tolist(), 'radius':radius,
                            'spansU':spansU + 4, 'spansV':spansV}
                data[eye]['Eyeball'] = eye_data
                continue
            
            
            elif muscle_name == "SO":
                
                # Split tendon into two pieces.
                tendon_name = '{}_{}_tendon_bspline'.format(muscle.eye, muscle_name)
                tendon_length = cmds.arclen(tendon_name)
                print "Tendon length: ", tendon_length
                # How much we want SO2 vs SO3. 
                split_proportion = 0.66
                SO2_length = split_proportion * tendon_length
                
                # Get param on the bspline
                split_param = curve_utils.getParamFromLength(tendon_name, SO2_length)
                so2_tendon, so3_tendon = cmds.detachCurve(tendon_name, p=split_param, rpo=False)
                so2_tendon = cmds.rename(so2_tendon, "{}_SO2_tendon_bspline".format(muscle.eye))
                so3_tendon = cmds.rename(so3_tendon, "{}_SO3_tendon_bspline".format(muscle.eye))

                cmds.parent(so2_tendon, muscle.objectGroup)
                cmds.parent(so3_tendon, muscle.objectGroup)
                cmds.rebuildCurve(so3_tendon,kr=4, spans=15, kep=True)
                
                # Get slices
                muscle_slices = getSlices('{}_{}_muscle_bspline'.format(muscle.eye, muscle_name), muscle)
                so2_slices = getSlices(so2_tendon, muscle)
                so3_slices = getSlices(so3_tendon, muscle)
                
                # Get centroids
                muscle_centroids = getCentroids(muscle_slices) / scale
                muscle_centroids = muscle_centroids.tolist()
                so2_centroids = getCentroids(so2_slices) / scale
                so2_centroids = so2_centroids.tolist()
                so3_centroids = getCentroids(so3_slices) / scale
                so3_centroids = so3_centroids.tolist()
                
                # Get radii
                muscle_radii = getSliceRadii(muscle_slices) / scale
                muscle_radii = muscle_radii.tolist()
                so2_radii = getSliceRadii(so2_slices) / scale
                so2_radii = so2_radii.tolist()
                so3_radii = getSliceRadii(so3_slices) / scale
                so3_radii = so3_radii.tolist()
                
                # Get control vertices
                #muscle_cvs = cmds.getAttr('{}_{}_muscle_bspline.cv[*]'.format(muscle.eye, muscle_name))
                #so2_cvs = cmds.getAttr("{}.cv[*]".format(so2_tendon))
                #so3_cvs = cmds.getAttr("{}.cv[*]".format(so3_tendon))
                muscle_cvs = getCurvePoints('{}_{}_muscle_bspline'.format(muscle.eye, muscle_name))
                muscle_cvs = np.array(muscle_cvs)
#                muscle_cvs = np.vstack((muscle_cvs[0], muscle_cvs, muscle_cvs[-1])).tolist()
                muscle_cvs = np.vstack((muscle_cvs[0] + np.random.normal(0, 1e-13),
                                        muscle_cvs[0] + np.random.normal(0, 1e-13),
                                        muscle_cvs[0] + np.random.normal(0, 1e-13),
                                        muscle_cvs[0] + np.random.normal(0, 1e-13), 
                                        muscle_cvs,
                                        muscle_cvs[-1] + np.random.normal(0, 1e-13),
                                        muscle_cvs[-1] + np.random.normal(0, 1e-13),
                                        muscle_cvs[-1] + np.random.normal(0, 1e-13),
                                        muscle_cvs[-1] + np.random.normal(0, 1e-13))).tolist()
                so2_cvs = getCurvePoints(so2_tendon)
                so2_cvs = np.array(so2_cvs)
#                so2_cvs = np.vstack((so2_cvs[0], so2_cvs, so2_cvs[-1])).tolist()
                so2_cvs = np.vstack((so2_cvs[0] + np.random.normal(0, 1e-13),
                                        so2_cvs[0] + np.random.normal(0, 1e-13),
                                        so2_cvs[0] + np.random.normal(0, 1e-13),
                                        so2_cvs[0] + np.random.normal(0, 1e-13), 
                                        so2_cvs,
                                        so2_cvs[-1] + np.random.normal(0, 1e-13),
                                        so2_cvs[-1] + np.random.normal(0, 1e-13),
                                        so2_cvs[-1] + np.random.normal(0, 1e-13),
                                        so2_cvs[-1] + np.random.normal(0, 1e-13))).tolist()
                so3_cvs = getCurvePoints(so3_tendon)
                so3_cvs = np.array(so3_cvs)
#                so3_cvs = np.vstack((so3_cvs[0], so3_cvs, so3_cvs[-1])).tolist()
                so3_cvs = np.vstack((so3_cvs[0] + np.random.normal(0, 1e-13),
                        so3_cvs[0] + np.random.normal(0, 1e-13),
                        so3_cvs[0] + np.random.normal(0, 1e-13),
                        so3_cvs[0] + np.random.normal(0, 1e-13), 
                        so3_cvs,
                        so3_cvs[-1] + np.random.normal(0, 1e-13),
                        so3_cvs[-1] + np.random.normal(0, 1e-13),
                        so3_cvs[-1] + np.random.normal(0, 1e-13),
                        so3_cvs[-1] + np.random.normal(0, 1e-13))).tolist()
                
                # Rotate
                #muscle_cvs = rotatePoints(rotation_axis, rotation_degrees, muscle_cvs)
                #so2_cvs = rotatePoints(rotation_axis, rotation_degrees, so2_cvs)
                #so3_cvs = rotatePoints(rotation_axis, rotation_degrees, so3_cvs)
                muscle_cvs = (np.array(muscle_cvs) / scale).tolist()
                so2_cvs = (np.array(so2_cvs) / scale).tolist()
                so3_cvs = (np.array(so3_cvs) / scale).tolist()
                
                # Plot for temp
                muscle_data = {'muscle_centroids':muscle_centroids, 'muscle_bspline_cvs':muscle_cvs, 
                               'muscle_radii':muscle_radii, 'so2_centroids':so2_centroids,
                               'so2_radii':so2_radii, 'so2_bspline_cvs':so2_cvs,
                               'so3_centroids':so3_centroids, 'so3_radii':so3_radii,
                               'so3_bspline_cvs':so3_cvs, 'muscle_init_strain':0.5,
                               'so2_init_strain':0.5, 'so3_init_strain':0.5,
                               'name':muscle_name}
                
                # Add properties
                radius = np.mean(muscle_radii)
                muscle_data = addBiomechanicalProperties(muscle_data, radius)
                data[eye]['muscles'][muscle_name] = muscle_data
                continue

            # Grab muscle control vertices
            cmds.rebuildCurve('{}_{}_tendon_bspline'.format(muscle.eye, muscle_name),kr=4, spans=15, kep=True)
            #tendon_cvs = cmds.getAttr('{}_{}_tendon_bspline.cv[*]'.format(muscle.eye, muscle_name))
            tendon_cvs = getCurvePoints('{}_{}_tendon_bspline'.format(muscle.eye, muscle_name))
            tendon_cvs = np.array(tendon_cvs)
#            tendon_cvs = np.vstack((tendon_cvs[0], tendon_cvs, tendon_cvs[-1])).tolist()
            tendon_cvs = np.vstack((tendon_cvs[0] + np.random.normal(0, 1e-13),
                    tendon_cvs[0] + np.random.normal(0, 1e-13),
                    tendon_cvs[0] + np.random.normal(0, 1e-13),
                    tendon_cvs[0] + np.random.normal(0, 1e-13), 
                    tendon_cvs,
                    tendon_cvs[-1] + np.random.normal(0, 1e-13),
                    tendon_cvs[-1] + np.random.normal(0, 1e-13),
                    tendon_cvs[-1] + np.random.normal(0, 1e-13),
                    tendon_cvs[-1] + np.random.normal(0, 1e-13))).tolist()
            #tendon_cvs = rotatePoints(rotation_axis, rotation_degrees, tendon_cvs)
            tendon_cvs = (np.array(tendon_cvs) / scale).tolist()
            
            
            # Grab tendon control vertices
            if not muscle_name == 'IO':
                cmds.rebuildCurve('{}_{}_muscle_bspline'.format(muscle.eye, muscle_name),kr=4, spans=6, kep=True)
                muscle_cvs = getCurvePoints('{}_{}_muscle_bspline'.format(muscle.eye, muscle_name))
                muscle_cvs = np.array(muscle_cvs)
#                muscle_cvs = np.vstack((muscle_cvs[0],muscle_cvs, muscle_cvs[-1])).tolist()
                muscle_cvs = np.vstack((muscle_cvs[0] + np.random.normal(0, 1e-13),
                        muscle_cvs[0] + np.random.normal(0, 1e-13),
                        muscle_cvs[0] + np.random.normal(0, 1e-13),
                        muscle_cvs[0] + np.random.normal(0, 1e-13), 
                        muscle_cvs,
                        muscle_cvs[-1] + np.random.normal(0, 1e-13),
                        muscle_cvs[-1] + np.random.normal(0, 1e-13),
                        muscle_cvs[-1] + np.random.normal(0, 1e-13),
                        muscle_cvs[-1] + np.random.normal(0, 1e-13))).tolist()
#                muscle_cvs = np.vstack((muscle_cvs[0],muscle_cvs[0],muscle_cvs[0],muscle_cvs[0], muscle_cvs,  muscle_cvs[-1],  muscle_cvs[-1],  muscle_cvs[-1],  muscle_cvs[-1])).tolist()

                #muscle_cvs = cmds.getAttr('{}_{}_muscle_bspline.cv[*]'.format(muscle.eye, muscle_name))
                #muscle_cvs = rotatePoints(rotation_axis, rotation_degrees, muscle_cvs)
                muscle_cvs = (np.array(muscle_cvs) / scale).tolist()
                muscle_slices = getSlices('{}_{}_muscle_bspline'.format(muscle.eye, muscle_name), muscle)
                muscle_radii = (getSliceRadii(muscle_slices) / scale).tolist()
                muscle_centroids = (getCentroids(muscle_slices) / scale).tolist()
                radius = np.mean(muscle_radii)
            else:
                muscle_cvs = []
                muscle_slices = []
                muscle_radii = []
                muscle_centroids = []
            
            # Grab orbital points if applicable
            if muscle_name in orbital_names:
                orbital_points = getPositionsFromGroup('{}_{}_{}_orbital_points'.format(muscle.plane, muscle.eye, muscle_name))
                orbital_points = np.array(orbital_points)
#                orbital_points = (np.vstack((  orbital_points[0], orbital_points,  orbital_points[-1])) / scale).tolist()
                orbital_points = np.vstack((orbital_points[0] + np.random.normal(0, 1e-13),
                                        orbital_points[0] + np.random.normal(0, 1e-13),
                                        orbital_points[0] + np.random.normal(0, 1e-13),
                                        orbital_points[0] + np.random.normal(0, 1e-13), 
                                        orbital_points,
                                        orbital_points[-1] + np.random.normal(0, 1e-13),
                                        orbital_points[-1] + np.random.normal(0, 1e-13),
                                        orbital_points[-1] + np.random.normal(0, 1e-13),
                                        orbital_points[-1] + np.random.normal(0, 1e-13))).tolist()
                #orbital_points = rotatePoints(rotation_axis, rotation_degrees, orbital_points)
            else:
                orbital_points = []
            
            if not muscle_name == 'IO':
                # Grab muscle, tendon slices4
                tendon_slices = getSlices('{}_{}_tendon_bspline'.format(muscle.eye, muscle_name), muscle)
            else:
                io_muscle = muscles['Sagittal'][eye]['IO']
                tendon_slices = getSlices('{}_{}_tendon_bspline'.format(muscle.eye, muscle_name), io_muscle)
                temp_grp = cmds.group(em=True, name="{}_IO_muscle_bspline_starts".format(io_muscle.eye))
                for (z_index, points) in io_muscle.rawPoints.items():
                    pt = tuple(np.mean(points, axis=0))
                    loc = cmds.spaceLocator(p=pt)
                    cmds.parent(loc[0], temp_grp)
                    
            # Grab radii
            tendon_radii = getSliceRadii(tendon_slices) / scale
            tendon_radii = tendon_radii.tolist()
            print "RAdii: ", tendon_radii
            
            if muscle_name == "IO":
                radius = np.mean(tendon_radii)
            # Grab centroids
            tendon_centroids = getCentroids(tendon_slices) / scale
            tendon_centroids = tendon_centroids.tolist()
            
            muscle_data = {'muscle_centroids':muscle_centroids, 'muscle_bspline_cvs':muscle_cvs, 
                           'muscle_radii':muscle_radii, 'tendon_centroids':tendon_centroids,
                           'tendon_bspline_cvs':tendon_cvs, 'tendon_radii':tendon_radii,
                           'orbital_points':orbital_points, 'muscle_init_strain':0.5,
                           'tendon_init_strain':0.5, 'name':muscle_name}
            muscle_data = addBiomechanicalProperties(muscle_data, radius)
            data[eye]['muscles'][muscle_name] = muscle_data
            

        
    # Export
    with open(out_file_path, 'w') as out_file:
        json.dump(data, out_file)
    
def getComponentIndexList(componentList=[]):
	'''
	Return an list of integer component index values
	@param componentList: A list of component names. if empty will default to selection.
	@type componentList: list
	'''
	# Initialize return dictionary
	componentIndexList = {}
	
	# Check string input
	if type(componentList) == str or type(componentList) == unicode:
		componentList = [componentList]
	
	# Get selection if componentList is empty
	if not componentList: componentList = cmds.ls(sl=True,fl=True) or []
	if not componentList: return []
	
	# Get MSelectionList
	selList = OpenMaya.MSelectionList()
	for i in componentList: selList.add(str(i))
	
	# Iterate through selection list
	selPath = OpenMaya.MDagPath()
	componentObj = OpenMaya.MObject()
	componentSelList = OpenMaya.MSelectionList()
	for i in range(selList.length()):
		# Check for valid component selection
		selList.getDagPath(i,selPath,componentObj)
		if componentObj.isNull():
			# Clear component MSelectionList
			componentSelList.clear()
			# Get current object name
			objName = selPath.partialPathName()
			
			# Transform
			if selPath.apiType() == OpenMaya.MFn.kTransform:
				numShapesUtil = OpenMaya.MScriptUtil()
				numShapesUtil.createFromInt(0)
				numShapesPtr = numShapesUtil.asUintPtr()
				selPath.numberOfShapesDirectlyBelow(numShapesPtr)
				numShapes = OpenMaya.MScriptUtil(numShapesPtr).asUint()
				selPath.extendToShapeDirectlyBelow(numShapes-1)
			
			# Mesh
			if selPath.apiType() == OpenMaya.MFn.kMesh:
				meshFn = OpenMaya.MFnMesh(selPath.node())
				vtxCount = meshFn.numVertices()
				componentSelList.add(objName+'.vtx[0:'+str(vtxCount-1)+']')
			# Curve
			elif selPath.apiType() == OpenMaya.MFn.kNurbsCurve:
				curveFn = OpenMaya.MFnNurbsCurve(selPath.node())
				componentSelList.add(objName+'.cv[0:'+str(curveFn.numCVs()-1)+']')
			# Surface
			elif selPath.apiType() == OpenMaya.MFn.kNurbsSurface:
				surfaceFn = OpenMaya.MFnNurbsSurface(selPath.node())
				componentSelList.add(objName+'.cv[0:'+str(surfaceFn.numCVsInU()-1)+'][0:'+str(surfaceFn.numCVsInV()-1)+']')
			# Lattice
			elif selPath.apiType() == OpenMaya.MFn.kLattice:
				sDiv = cmds.getAttr(objName+'.sDivisions')
				tDiv = cmds.getAttr(objName+'.tDivisions')
				uDiv = cmds.getAttr(objName+'.uDivisions')
				componentSelList.add(objName+'.pt[0:'+str(sDiv - 1)+'][0:'+str(tDiv - 1)+'][0:'+str(uDiv - 1)+']')
			
			# Get object component MObject
			componentSelList.getDagPath(0,selPath,componentObj)
		
		# =======================
		# - Check Geometry Type -
		# =======================
		
		# MESH / NURBS CURVE
		if (selPath.apiType() == OpenMaya.MFn.kMesh) or (selPath.apiType() == OpenMaya.MFn.kNurbsCurve):
			indexList = OpenMaya.MIntArray()
			componentFn = OpenMaya.MFnSingleIndexedComponent(componentObj)
			componentFn.getElements(indexList)
			componentIndexList[selPath.partialPathName()] = list(indexList)
		
		# NURBS SURFACE
		if selPath.apiType() == OpenMaya.MFn.kNurbsSurface:
			indexListU = OpenMaya.MIntArray()
			indexListV = OpenMaya.MIntArray()
			componentFn = OpenMaya.MFnDoubleIndexedComponent(componentObj)
			componentFn.getElements(indexListU,indexListV)
			componentIndexList[selPath.partialPathName()] = zip(list(indexListU),list(indexListV))
		
		# LATTICE
		if selPath.apiType() == OpenMaya.MFn.kLattice:
			indexListS = OpenMaya.MIntArray()
			indexListT = OpenMaya.MIntArray()
			indexListU = OpenMaya.MIntArray()
			componentFn = OpenMaya.MFnTripleIndexedComponent(componentObj)
			componentFn.getElements(indexListS,indexListT,indexListU)
			componentIndexList[selPath.partialPathName()] = zip(list(indexListS),list(indexListT),list(indexListU))
	
	# Return Result
	return componentIndexList

def getSingleIndex(obj,index):
	'''
	Convert a 2 or 3 value index to a single value index.
	Returns the single element index of the given component.
	@param obj: Object parent of component.
	@type obj: str
	@param index: Multi element index of component.
	@type index: list
	'''
	
	# Mesh
	if cmds.objectType(obj) == 'mesh': return index
	
	# Nurbs Curve
	if cmds.objectType(obj) == 'nurbsCurve': return index
	
	# Nurbs Surface
	if cmds.objectType(obj) == 'nurbsSurface':
		# Get nurbsSurface function set
		surfList = OpenMaya.MSelectionList()
		surfObj = OpenMaya.MObject()
		OpenMaya.MGlobal.getSelectionListByName(obj,surfList)
		surfList.getDependNode(0,surfObj)
		surfFn = OpenMaya.MFnNurbsSurface(surfObj)
		# CV count in U an V directions
		numV = surfFn.numCVsInV()
		# Check for periodic surface
		if surfFn.formInV() == surfFn.kPeriodic:
			numV -= surfFn.degreeV()
		# Get Single Index
		return (index[0] * numV) + index[1]
	
	# Lattice
	elif cmds.objectType(obj) == 'lattice':
		sDiv = cmds.getAttr(obj+'.sDivisions')
		tDiv = cmds.getAttr(obj+'.tDivisions')
		return (index[0] + (index[1] * sDiv) + (index[2] * sdiv * tDiv) )
	
	# Return Result
	return None

def stringIndex(index,padding=2):
	'''
	Return the string equivalent for the specified iteger index.
	@param index: The index to get the string equivalent for
	@type index: int
	@param padding: The number of characters for the index string
	@type padding: int
	'''
	# Convert to String
	strInd = str(index)
	# Prepend Padding
	for i in range(padding-len(strInd)): strInd = '0'+strInd
	# Return Result
	return strInd

def getSurfacePoints(surface):
    
    controlPoints = getComponentIndexList(surface)[surface + "Shape"]
    spansU = cmds.getAttr("{}.spansU".format(surface))
    spansV = cmds.getAttr("{}.spansV".format(surface))
    # Create locators and connect to control points
    points = []
    for cv in controlPoints:
    	# Get control point world position
    	pos = cmds.pointPosition(surface+'.cv['+str(cv[0])+']['+str(cv[1])+']')
        points.append(pos)
    return points, spansU, spansV

        
def printMuscleLengths(muscles):

    muscle_names = ['MR', 'LR', 'SR', 'IR', 'SO', 'IO']

    # Iterate through each eye
    for eye, eye_muscles in muscles['Coronal'].items():
        
        print "Eye: ", eye
        
        for muscle_name in muscle_names:
            
            muscle = eye_muscles[muscle_name]
            
#            print muscle_name, "Tendon Length: ", cmds.arclen(muscle.tendonBspline), "\t Muscle Length: ", cmds.arclen(muscle.muscleBspline)
            print muscle_name, " Radius: ", muscle.radius

            
            
            

            
            
            
            
            
