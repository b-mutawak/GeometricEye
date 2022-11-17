# -*- coding: utf-8 -*-
"""
Created on Wed Jun 24 21:21:27 2020

@author: Bassam
"""


import numpy as np
from math import ceil, cos, sin, sqrt, pi, atan, degrees, asin
from parseInputPoints import parseDataFolder
from time import time
from os.path import exists
from collections import OrderedDict
import maya.cmds as cmds
from scipy import optimize
import json
import shape_utils
from shape_utils import closest_node




class EOS:
    """
        Data structure to hold major components of the extra ocular system.
    
    """
    
    def __init__(self):
        self.eyeball = None
        self.lr = None
        self.mr = None
        self.sr = None
        self.ir = None
        self.so = None
        self.io = None
        return
    

class muscleObject:
    """
        Data structure to hold different muscle objects.
        
        Plane  = Imaging plane and eye orientation.
            Ex: 'Sagittal Down'
        
        Eye    = Left or Right Eye
            Ex: 'OD'
        
        Points = Dictionary of 3D points. Each key corresponds to a Z-axis index, with the values being the 
        points that have that z-value
    
    """
    def __init__(self, plane, eye, muscle_name, sorted_points):
        self.plane = plane
        self.eye = eye
        self.sorted_points = sorted_points
        self.muscle_name = muscle_name
        self.rawPoints = []
        self.verts = []
        self.object = None
        self.object_list = []
        self.center = None
        self.curve = None
        self.radius = None
        self.bspline = None
        self.muscleBspline = None
        self.tendonBspline = None
        self.muscleInsertions = {}
        self.rawGroup = ''
        self.objectGroup = ''
        self.pulleyLocation = ''
        self.controlPoints = []
        self.radii = {}
        self.radius = 0
        self.PCSA = 0
        self.maxPassiveForce = 0
        self.maxActiveForce = 0
        self.maxTendonPassiveForce =0
        self.muscleDamping = 0
        self.tendonDamping =0
        self.initialStrain =0
        self.globalLayerPoints = []
        self.radii_points = {}
        self.centers = {}
        return
    
    def setObjects(self, objects):
        self.object_list = objects
        return
    
    def setPulleyLocation(self, pulley):
        self.pulleyLocation = pulley
        return

    def setNURBS(self, NURBS):
        self.object = NURBS
        return
    
    def setRadius(self, radius):
        self.radius = radius
        return
    
    def setCenter(self, center):
        self.center = center
        return
    
    def setMuscleVerts(self, verts):
        self.verts = verts
        return
    
    def setRawPoints(self, rawPoints):
        self.rawPoints = rawPoints
        return
    
    def setMuscleCurve(self, curve):
        self.curve = curve
        return
    
    def setMuscleBSpline(self, bspline):
        self.bspline = bspline
        return
    
    
    def setLocators(self, locators):
        self.locators = locators
        return
    
    def setMuscleInsertions(self, inserts):
        self.muscleInsertions = inserts
        return
    
    def setRawGroup(self, group):
        self.rawGroup = group
        return
    
    def setObjectGroup(self, group):
        self.objectGroup = group
        return
    
    def setControlPoints(self, cps):
        self.controlPoints = cps
        return
    
    def setMuscleBspline(self, bspline):
        self.muscleBspline = bspline
        return
    
    def setTendonBspline(self, bspline):
        self.tendonBspline = bspline
        return
    
    def setOriginPlane(self, originPlane):
        self.originPlane = originPlane
        return
    

def plotEye(muscles, showRaw=False):
    """
        Unhides all generated objects for each muscle.
    """
    
    # Go through each muscle
    for muscle in muscles:
        
        # Show object if available
        if not muscle.object is None:
            cmds.showHidden(muscle.object)
        
        # Show object_list if available
        if not len(muscle.object_list) == 0:
            
            # Iterate through each object in list
            for obj in muscle.object_list:
                cmds.showHidden(obj)
                
        # Show raw points if available
        if showRaw:
            
            # Iterate through each raw point
            for raw_point in muscle.rawPoints:
                cmds.showHidden(raw_point)
            
                
    return

def sortData(data):
    """
        Sorts incoming data into each MRI image slice. 
        
        Returns a dictionary for each muscle set, with keys being the z-index
    """
    
    toReturn = {}
    
    # Loop through each plane and muscles combo
    for eye_set, muscles in data.items():
        
        if "Sagittal" in eye_set:
            pass
        
        
        # Loop through each muscle
        for muscle, muscle_points in muscles.items():
            sorted_points = OrderedDict()
            
            # loop through each point
            for i, point in enumerate(np.array(muscle_points)):
                                    
                # Add the key if it doesn't exist
                if not str(point[2]) in sorted_points.keys():
                    sorted_points[str(point[2])] = np.zeros((0,3), dtype=np.float32)
                    
                sorted_points[str(point[2])] = np.concatenate((sorted_points[str(point[2])], np.array(point).reshape(1,3)), axis=0)
            
            if not eye_set in toReturn.keys():
                toReturn[eye_set] = {}
            
            # Add sorted points
            toReturn[eye_set][muscle] = sorted_points
            
            
    return toReturn



def combinePoints(data):
    
    x = np.zeros((0, 3), dtype=float)
    for z_index, points in data.items():
        x = np.append(x, points, axis=0)
    
    return x


def generateObjects(data):
    """
        Given a dictionary of points, this will populate an EOS object for each eye
    """
    
    # Muscles
    muscles = {}
    
    # Create OD and OS groups
    od_group = cmds.group(em=True, name='OD')
    os_group = cmds.group(em=True, name='OS')
    
    # Hide them
    cmds.hide(od_group)
    cmds.hide(os_group)
    ellipse_plot = {'OS':{}, 'OD':{}}
    
    # Create origin plane
    originPlane = shape_utils.createOriginPlane()
    
    
    # Loop through each plane/eye
    for point_set_name, point_sets in data.items():
        
        # Check which eye this is
        eye = 'OS' if 'OS' in point_set_name else 'OD' 
        
        # Plane
        plane = 'Sagittal' if 'Sagittal' in point_set_name else 'Coronal'
        
        # Create location for it
        if not plane in muscles.keys():
            muscles[plane] = {'OD':{}, 'OS':{}}
        
        # Create plane group
        parent_group = od_group if 'OD' in point_set_name else os_group
        plane_group = cmds.group(em=True, name=point_set_name, parent=parent_group)
        
        # Create holder for eye points in one of the planes
        eye_points = None
        
        # Iterate through each muscle in the set
        for muscle_name, points in point_sets.items():

            # Create muscle group
            muscle_group = cmds.group(em=True, name=muscle_name, parent=plane_group)
            
            # Create raw points group
            raw_group = cmds.group(em=True, name="{}_{}_{}_raw".format(plane, eye, muscle_name), parent=muscle_group)
            cmds.hide(raw_group)
            
            # Create objects group
            obj_group = cmds.group(em=True, name="{}_{}_{}_obj".format(plane, eye, muscle_name), parent=muscle_group)
            
            # Grab the points
            sorted_points = data[point_set_name][muscle_name]
            
            sagittal = True
            if "Coronal" in point_set_name:
                sagittal = False

            
            # Raw locator points
            locators = []
            
            # Loop through each z-index points
            for z_index, points in sorted_points.items():
                
                # Loop through each point
                for point in points:
                    
                    # Create locator
                    locator = cmds.spaceLocator(p=(point[0], point[1], point[2]))
                    
                    if sagittal:    
                        setRGBColor(locator[0], color=(0,0,0))
                    else:
                        setRGBColor(locator[0], color=(0,1,0))
                    
                    # Shrink
                    setLocalScale(locator[0], size=[0.25, 0.25, 0.25])
                    
                    # Group
                    cmds.parent(locator[0], raw_group)
                    
                    # Append
                    locators.append(locator[0])
            
            
            
            # Handle eyeball/socket separately
            if muscle_name == "Eyeball" and "Coronal" in point_set_name:
                
                # Create eye object
                eye_object = muscleObject(plane, eye, muscle_name, sorted_points)
                eye_object.objectGroup = obj_group
                eye_object.rawGroup = raw_group
                
                # Combine both point clouds
                cor_points = combinePoints(sorted_points)
                sag_points = combinePoints(data[point_set_name.replace("Coronal", "Sagittal")][muscle_name])
                sag_points = sag_points[sag_points[:, 2] < 21]
#                eye_points = np.append(sag_points, cor_points, axis=0)
                eye_points = cor_points
                
                # Now fit ellipsoid
                radii, eyeCenter = shape_utils.ellipsoid_fit_matlab(eye_points.T, 200)
                
                # Average the radii to get a sphere
                radii = np.ones((1, 3))*np.mean(radii)
                radii = radii.tolist()[0][0]
                eyeCenter = eyeCenter.tolist()[0]
                
                # Bound radius
                radii = max(radii, 11.6)
                radii = min(radii, 13)
                
                # Create eyeball
                eyeball = cmds.sphere(r=radii, n='{}_{}_Eyeball'.format(plane, eye), po=0)
                cmds.move(eyeCenter[0], eyeCenter[1], eyeCenter[2], eyeball[0])
                cmds.parent(eyeball[0], obj_group)
                
                # Set the nurbs object
                eye_object.setNURBS(eyeball)
                eye_object.setLocators(locators)
                eye_object.setRawPoints(sorted_points)
                eye_object.setCenter(eyeCenter)
                eye_object.setRadius(radii)

                # Add to list
                muscles[plane][eye][muscle_name] = eye_object

            
            elif muscle_name == "Socket" and "Coronal" in point_set_name:
                
                # Create socket object
                socket_object = muscleObject(plane, eye, muscle_name, sorted_points)
                socket_object.objectGroup = obj_group
                socket_object.rawGroup = raw_group
                curves = []
                
                # Go through each slice
                for key, values in sorted_points.items():

                    # Compute a least-squares fit of a circle to point cloud
                    xs = [val[0] for val in values]
                    ys = [val[1] for val in values]
                    xc, yc, R = shape_utils.circleFit(xs, ys)
                    
                    # Generate circle
                    circ = cmds.circle(center=[xc, yc, int(float(key))], radius=R)[0]
                    cmds.hide(circ)
                    cmds.parent(circ, obj_group)
                    curves.append(circ)
                
                # Generate a loft based on the circles generated
                socket_obj = cmds.loft(curves)[0]
                cmds.hide(socket_obj)
                cmds.parent(socket_obj, obj_group)
                
                # Set attributes
                socket_object.setObjects(curves)
                socket_object.setNURBS(socket_obj)
                socket_object.setRawPoints(sorted_points)
                socket_object.setLocators(locators)
                muscles[plane][eye][muscle_name] = socket_object
                
            else:
                
                # Create object
                other_obj = muscleObject(plane, eye, muscle_name, sorted_points)
                other_obj.objectGroup = obj_group
                other_obj.rawGroup = raw_group
                other_obj.setRawPoints(sorted_points)
                other_obj.setLocators(locators)
                other_obj.setOriginPlane(originPlane)
                muscles[plane][eye][muscle_name] = other_obj
                
    
    
    
    # Write ellipse plots
    with open(r"C:\Users\admin\Dropbox\Maya Model Generator\Project\Maya_New\Maya2020\bin\plug-ins\Model_Generator\ellipse_results.json", 'w') as fp:
        json.dump(ellipse_plot, fp)
                                
    return muscles
    
def intersection_of_lists(lst1, lst2):
    """
        Returns the elements that are part of both lists
    """
    intersect = []
    for val in lst1:
        if val in lst2:
            intersect.append(val)
            
    return intersect

def correctCoronalRotation(data):
    """
        Corrects any XY coronal rotation.
        Returns the corrected data points
    """
    
    toReturn = {}
    
    # Loop through each plane and muscles combo
    for eye_set, muscles in data.items():
        
        
        # Loop through each muscle
        for muscle, muscle_points in muscles.items():
            sorted_points = OrderedDict()
            
            # loop through each point
            for i, point in enumerate(np.array(muscle_points)):
                                    
                # Add the key if it doesn't exist
                if not str(point[2]) in sorted_points.keys():
                    sorted_points[str(point[2])] = np.zeros((0,3), dtype=np.float32)
                    
                sorted_points[str(point[2])] = np.concatenate((sorted_points[str(point[2])], np.array(point).reshape(1,3)), axis=0)
            
            if not eye_set in toReturn.keys():
                toReturn[eye_set] = {}
            
            # Add sorted points
            toReturn[eye_set][muscle] = sorted_points
            
    sorted_data = toReturn
    
    # Loop through each plane/eye
    for point_set_name, point_sets in sorted_data.items():
        
        if 'Sagittal' in point_set_name:
            continue
        
        # Grab LR, MR points
        mr_points = sorted_data[point_set_name]['MR']
        lr_points = sorted_data[point_set_name]['LR']
        
        # Grab the slice indexes, and the intersection of 
        common_keys = intersection_of_lists(mr_points.keys(), lr_points.keys())
        
        # Sort least to greatest to get anterior points first, remove the first element
        common_keys.sort()
        common_keys.pop(0)
        common_keys.pop(0)
        
        # How many slices to take into account
        num_slices_include = 3
        rotations = []
        for i in range(num_slices_include):
            
            # Grab current common slice and corresponding points
            current_slice_idx = common_keys[i]
            mr_rot_pts = np.array(mr_points[current_slice_idx])
            lr_rot_pts = np.array(lr_points[current_slice_idx])
            
            # Grab center of points
            mr_center = np.mean(mr_rot_pts, axis=0)
            lr_center = np.mean(lr_rot_pts, axis=0)
            
            # Compute angle of rotation
            if mr_center[1] > lr_center[1]:
                opposite_side = mr_center[1] - lr_center[1]
            else:
                opposite_side = lr_center[1] - mr_center[1]
            
            if mr_center[0] > lr_center[0]:
                adjacent = mr_center[0] - lr_center[0]
            else:
                adjacent = lr_center[0] - mr_center[0]
                
            hyp = sqrt((opposite_side**2) + (adjacent**2))
            angle = asin(opposite_side / hyp)
            rotations.append(angle)
        

        # Average rotations
        rotation_avg = degrees(np.mean(rotations)) * -1
        
        # Get rotation matrices
        rot_matrix = shape_utils.getRotationMatrix('z', rotation_avg)
        
        # Apply to all
        for key, points in data[point_set_name].items():
            
            data[point_set_name][key] = rot_matrix.dot(np.array(points).T).T
            
    
        
    return data
    


def computeBestTransform(r_set, t_set):
    """
        Computes the best rotation and translation transform to map
        the points in r_set to t_set.
        
        Returns the transform matrix that maps r_set -> t_set.
    """
    
    # Compute centroid of each set
    r_set_mean = np.mean(r_set, axis=0)
    t_set_mean = np.mean(t_set, axis=0)
    
    # Center the points of each set using centroid
    r_centered = r_set - r_set_mean
    t_centered = t_set - t_set_mean
    
    # Compute rotation matrix
    H = np.dot(r_centered.T, t_centered)
    U, S, Vt = np.linalg.svd(H)
    R = np.dot(Vt.T, U.T)
    
    # Length of each array
    m = r_set.shape[1]
    
    # special reflection case ???
    if np.linalg.det(R) < 0:
       Vt[m-1,:] *= -1
       R = np.dot(Vt.T, U.T)

    # translation
    t = t_set_mean.T - np.dot(R,r_set_mean.T)

    # homogeneous transformation
    T = np.identity(m+1)
    T[:m, :m] = R
    T[:m, m] = t

    return T, R, t
    

def sub_sample_points(points, perc):
    """
    
    Samples a set of points from the input point set. The number of
    sampled points is based on the passed in percentage value
    
    """
    
    # Compute new number of points
    points = np.array(points)
    num_sampled_points = int(ceil(points.shape[0] * perc))

    # Set of indices from old array
    rand_indices = np.arange(0, points.shape[0])

    # Shuffle array
    np.random.shuffle(rand_indices)
    
    # Get the new sampled indices
    return points[rand_indices[:num_sampled_points]]
    


# http://www.mathworks.com/matlabcentral/fileexchange/24693-ellipsoid-fit
# for arbitrary axes
def ellipsoid_fit(X):
    x = X[:, 0]
    y = X[:, 1]
    z = X[:, 2]
    D = np.array([x * x + y * y - 2 * z * z,
                 x * x + z * z - 2 * y * y,
                 2 * x * y,
                 2 * x * z,
                 2 * y * z,
                 2 * x,
                 2 * y,
                 2 * z,
                 1 - 0 * x])
    d2 = np.array(x * x + y * y + z * z).T # rhs for LLSQ
    u = np.linalg.solve(D.dot(D.T), D.dot(d2))
    a = np.array([u[0] + 1 * u[1] - 1])
    b = np.array([u[0] - 2 * u[1] - 1])
    c = np.array([u[1] - 2 * u[0] - 1])
    v = np.concatenate([a, b, c, u[2:]], axis=0).flatten()
    A = np.array([[v[0], v[3], v[4], v[6]],
                  [v[3], v[1], v[5], v[7]],
                  [v[4], v[5], v[2], v[8]],
                  [v[6], v[7], v[8], v[9]]])

    center = np.linalg.solve(- A[:3, :3], v[6:9])

    translation_matrix = np.eye(4)
    translation_matrix[3, :3] = center.T

    R = translation_matrix.dot(A).dot(translation_matrix.T)

    evals, evecs = np.linalg.eig(R[:3, :3] / -R[3, 3])
    evecs = evecs.T

    radii = np.sqrt(1. / np.abs(evals))
    radii *= np.sign(evals)

    return center, evecs, radii, v

def registerPoints(r_set, t_set, sub_sample_perc=0, std_dev_inc=0, num_iterations=0, exit_criteria=0):
    """
        Registers the registration set (r_set) to the target set (t_set) using ICP.
        Sub-samples the registration set by some percentage.
        Removes all distances not within std_dev_inc from consideration.
        Will loop for maximum of num_iterations.
        If current average distance is less than or equal to exit_criteria, we break.
    
        Returns the registered set.
    """
    
    # Convert to numpy arrays and copy
    r_set_copy = np.copy(r_set)
    t_set_copy = np.copy(t_set)
    
    # Sub sample points
    r_set_sampled = sub_sample_points(np.asarray(r_set), sub_sample_perc)
    
    # Reshape to add another dimension
    m = r_set_sampled.shape[1]
    r_set = np.ones((r_set_sampled.shape[0], r_set_sampled.shape[1] + 1))
    t_set = np.ones((t_set_copy.shape[0],  t_set_copy.shape[1] + 1))
    r_set_mod = np.ones((r_set_sampled.shape[0],  r_set_sampled.shape[1] + 1))
    final_r_set = np.ones((r_set_copy.shape[0], r_set_copy.shape[1] + 1))
    final_r_set[:, :m] = r_set_copy
    r_set[:, :m] = np.copy(r_set_sampled)
    t_set[:, :m] = np.copy(t_set_copy)
    r_set_mod[:, :m] = np.copy(r_set_sampled)
        
    # Holder
    closest_points = np.zeros((r_set.shape[0], r_set.shape[1]))
    distances = np.zeros((r_set.shape[0], 1))
    distances_2 = np.zeros((t_set.shape[0], 1))
    
    # Loop through num_iteration times
    for iteration in range(num_iterations):
     
        # Get the closest points
        for i, point in enumerate(r_set):
            closest_point, distance, _ = closest_node(point, t_set, m)
            closest_points[i] = closest_point
            distances[i] = distance
        
        # Compute standard deviation of distances
        std_dev = np.std(distances)
           
        # Remove all points not within std_dev_inc number of standard deviations
        good_indices = distances[:] <= std_dev * std_dev_inc
        source_points = closest_points[good_indices[:,0]]
        dest_points = r_set[good_indices[:,0]]       
        
        # Get transformation matrix
        T, _, _ = computeBestTransform(dest_points[:, :m], source_points[:, :m])
        
        # Apply transform 
        r_set = np.dot(T, r_set.T).T
        
        # Check error
        avg_dist = np.mean(distances)
        if avg_dist <= exit_criteria:
            break
        
    # Calculate final transform
    T, _, _ = computeBestTransform(r_set_mod[:,:m], r_set[:, :m])
    
    # Apply
    final_r_set = np.dot(T, final_r_set.T).T
    
    # Compute final errors and distances
    closest_points = np.zeros((final_r_set.shape[0], final_r_set.shape[1]))
    distances = np.zeros((final_r_set.shape[0], 1))
    
    for i, point in enumerate(final_r_set):
        closest_point, distance, _ = closest_node(point, t_set, m)
        closest_points[i] = closest_point
        distances[i] = distance
        
    
    return final_r_set, T, np.sqrt(np.mean(distances ** 2))

def setLocalScale(ctrl, size=[1, 1, 1]):
    
    cmds.setAttr(ctrl + ".localScaleX", size[0])
    cmds.setAttr(ctrl + ".localScaleY", size[1])
    cmds.setAttr(ctrl + ".localScaleZ", size[2])
    return


    
def setRGBColor(ctrl, color = (1,1,1)):
    
    rgb = ("R","G","B")
    
    cmds.setAttr(ctrl + ".overrideEnabled",1)
    cmds.setAttr(ctrl + ".overrideRGBColors",1)
    
    for channel, color in zip(rgb, color):
        
        cmds.setAttr(ctrl + ".overrideColor%s" %channel, color)
    
    return


def generateEye(folder_path, plot=False):
    
    # Make sure path exists
    if not exists(folder_path):
        return
    
    # The planes we care about
    eye_points = ['Sagittal OD Prim', 'Sagittal OS Prim', 'Coronal OD Prim', 'Coronal OS Prim']
    
    globe_on_index = 11
    # Parse the folder
    data = parseDataFolder(folder_path, 0)
    
    # Remove unneccessary data
    for key in data.keys():
        if not key in eye_points:
            data.pop(key)
    
    
    # Correct coronal rotation
    data = correctCoronalRotation(data)

    # Sort into respective slices
    sorted_data = sortData(data.copy())
        
    # Apply the change here. Sort by Z-value for all data. Make this a copy
    first_muscle = "Eyeball"
    second_muscle = "Socket"
    
    # Grab relevant information. For now, we only care about 
    # primary view.
    od_cor_prim_sr = np.copy(np.array(data['Coronal OD Prim'][first_muscle], dtype=np.float32))
    od_cor_prim_ir = np.copy(np.array(data['Coronal OD Prim'][second_muscle], dtype=np.float32))

    od_sag_prim_sr = np.copy(np.array(data['Sagittal OD Prim'][first_muscle], dtype=np.float32)) 
    od_sag_prim_ir = np.copy(np.array(data['Sagittal OD Prim'][second_muscle], dtype=np.float32))
    
    os_cor_prim_sr = np.copy(np.array(data['Coronal OS Prim'][first_muscle], dtype=np.float32))
    os_cor_prim_ir = np.copy(np.array(data['Coronal OS Prim'][second_muscle], dtype=np.float32))
    
    os_sag_prim_sr = np.copy(np.array(data['Sagittal OS Prim'][first_muscle], dtype=np.float32))
    os_sag_prim_ir = np.copy(np.array(data['Sagittal OS Prim'][second_muscle], dtype=np.float32))


    # Combine eyeball and sockets
    od_cor_prim_points = np.append(od_cor_prim_sr, od_cor_prim_ir, axis = 0)
    os_cor_prim_points = np.append(os_cor_prim_sr, os_cor_prim_ir, axis = 0)
    od_sag_prim_points = np.append(od_sag_prim_sr, od_sag_prim_ir, axis = 0)
    os_sag_prim_points = np.append(os_sag_prim_sr, os_sag_prim_ir, axis = 0)
    
    mirror_dim = 0
    
    # Resize them to match
    od_cor_prim_points[:, 0:2] *= 0.3125
    od_sag_prim_points[:, 0:2] *= 0.3125
    os_cor_prim_points[:, 0:2] *= 0.3125
    os_sag_prim_points[:, 0:2] *= 0.3125
    
    # Rotate sagittal
    od_rot = shape_utils.getRotationMatrix('y', 90)
    os_rot = shape_utils.getRotationMatrix('y', 90)
    od_sag_prim_points = od_rot.dot(od_sag_prim_points.T).T
    os_sag_prim_points = os_rot.dot(os_sag_prim_points.T).T
    
    # Mirror X for OD Sagittal
    od_sag_prim_points[:, mirror_dim] = -1 * od_sag_prim_points[:, mirror_dim]
        
    # Translate sagittal centroid to coronal centroid
    od_translation = np.mean(od_cor_prim_points, axis=0) - np.mean(od_sag_prim_points, axis=0)
    od_sag_prim_points += od_translation
    os_translation = np.mean(os_cor_prim_points, axis=0) - np.mean(os_sag_prim_points, axis=0)
    os_sag_prim_points += os_translation
    
    # Get start time
    start = time()
    
    # Register OD
    registered_od_sag_prim_points, od_T, rms_od = registerPoints(od_sag_prim_points, od_cor_prim_points,
                                                         sub_sample_perc=1.0, std_dev_inc=1, 
                                                         num_iterations=50, exit_criteria=2.5)
    
    # Register OS
    registered_os_sag_prim_points, os_T, rms_os = registerPoints(os_sag_prim_points, os_cor_prim_points,
                                                         sub_sample_perc=1.0, std_dev_inc=1, 
                                                         num_iterations=50, exit_criteria=2.5)
    
    print "Registration time: ", time() - start
    eye_points = ['Sagittal OD Prim', 'Sagittal OS Prim']
    RMS_Errors = {}
    numPoints = 0
    # Loop over each eye
    for eye in eye_points:
        if "OD" in eye:
            rot = od_rot
            transform = od_T
        else:
            rot = os_rot
            transform = os_T
            
        # Loop over other points, extend to 4D, apply transform, and return back to 3D
        for key in sorted_data[eye].keys():
            
            # Go through each z-index 
            for z_index, points in sorted_data[eye][key].items():
                numPoints += len(points)
                
                # Get attributes
                og_points = np.array(points,  dtype=np.float32)
                
                # Resize
                og_points[:, 0:2] *= 0.3125

                # Translate
                if "OD" in eye:
                    og_points = od_rot.dot(og_points.T).T
                    og_points[:,  mirror_dim] = -1 * og_points[:,  mirror_dim]
                    og_points += od_translation
                else:
                    og_points = os_rot.dot(og_points.T).T
                    og_points += os_translation
                
                # Initialize a 4D holder
                points = np.ones((og_points.shape[0], og_points.shape[1] + 1))
                
                # Copy over stuff
                points[:, :3] = og_points
                
                # Apply transform
                points = np.dot(transform, points.T).T
                
                points[:, 2] = points[:, 2] - (2*globe_on_index)

                
                sorted_data[eye][key][z_index] = np.copy(points[:, :3]).tolist()
            
            
                
    # Resize the coronal
    eye_points = ['Coronal OD Prim', 'Coronal OS Prim']

    # Loop over each eye
    for eye in eye_points:
            
        # Loop over other points, extend to 4D, apply transform, and return back to 3D
        for key in data[eye].keys():
            
            aggregated_points = np.zeros((0,3), dtype=np.float32)
            # Go through each z-index 
            for z_index, points in sorted_data[eye][key].items():
                numPoints += len(points)
                # Get points
                og_points = np.array(points,  dtype=np.float32)
                
                # Resize
                og_points[:, 0:2] *= 0.3125
                
                og_points[:, 2] = og_points[:, 2] - (2*globe_on_index)

                
                # Copy over points
                sorted_data[eye][key][z_index] = np.copy(og_points[:, :3]).tolist()
                
                aggregated_points = np.append(aggregated_points, sorted_data[eye][key][z_index], axis=0)
            
            # Compute RMS error
            other_points = sorted_data[eye.replace("Coronal", "Sagittal")][key]
            
            aggregated_other = np.zeros((0,3), dtype=np.float32)
            for z_index, points in other_points.items():
                aggregated_other = np.append(aggregated_other, points, axis=0)
                
            # Compute final errors and distances
            distances = np.zeros((aggregated_other.shape[0], 1))
            
            for i, point in enumerate(aggregated_other):
                closest_point, distance, _ = closest_node(point, aggregated_points, 3)
                distances[i] = distance
            
            if not eye in RMS_Errors.keys():
                RMS_Errors[eye] = OrderedDict({'LR':[], 'MR':[], 'SR':[], 'IR':[], 'IO':[], 'SO':[], 'Socket':[], 'Eyeball':[]})
                
            RMS_Errors[eye][key] = [np.sqrt(np.mean(distances ** 2)), np.max(distances)]
            
    
    
    # Compute average
    os_avg = 0
    os_ctr = 0
    od_avg = 0
    od_ctr = 0
    for eye in eye_points:

            
        for key in data[eye].keys():
            if "OD" in eye:
                od_avg += RMS_Errors[eye][key][0]
                od_ctr += 1
            else:
                os_avg += RMS_Errors[eye][key][0]
                os_ctr += 1
            
    print "OD AVG: ", od_avg / od_ctr
    print "OS AVG: ", os_avg / os_ctr
    print ""
    
    for key in RMS_Errors.keys():
        print "===", key, "==="
        for muscle in RMS_Errors[key].keys():
            print muscle, " : ", RMS_Errors[key][muscle][0], " ", RMS_Errors[key][muscle][1]
        
    # Correct muscle rotation
#    sorted_data = shape_utils.correctRotation(sorted_data)
    
    
    # Generate objects from everything
    muscles = generateObjects(sorted_data)
    
    # Center muscles on eye center
    muscles = shape_utils.globeCenter(muscles)
    
    # Create muscle insertion points
    muscles = shape_utils.generateInsertions(muscles)
    
    # Filter out points if needed
    muscles = shape_utils.filterPoints(muscles)
    
    # Interpolate points between muscle traces and insertion
    muscles = shape_utils.generateControlPoints(muscles, SOP=False)
    
    # Create bspline for muscles
    muscles = shape_utils.generateBSplines(muscles)

    # Extend
    muscles = shape_utils.extendMuscles(muscles)
    
    # Split
    muscles = shape_utils.splitMuscleTendon(muscles)
    
    # Extend to origin 
    muscles = shape_utils.extendToOrigin(muscles)

    
    # Create pulleys
    shape_utils.generatePulleyLocations(muscles)
    
    muscles = shape_utils.generateRadius(muscles)
    
    muscles = shape_utils.generateBiomechanicalProperties(muscles)
    
    muscles = shape_utils.generateOrbitalBspline(muscles)
    

    
    
    print muscles
    
    # Print muscle lengths
    shape_utils.printMuscleLengths(muscles)
    print "Took: ", (time() - start) * 1000
    print "For: ", numPoints, " points"
    return muscles

    
    
if __name__=="__main__":
    folder_path = r"D:\New folder\Dropbox\Maya Model Generator\9 SOP Datasets\Q6 Anonymized Stacks"
    registerEye(folder_path, plot=False)
    