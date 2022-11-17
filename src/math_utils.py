# -*- coding: utf-8 -*-
#!/usr/bin/env python -W ignore::DeprecationWarning
"""
    Functions to implement basic geoemtry calculations used throughout.
    
    Heavily relies on source code from: 
        - https://spacetelescope.github.io/spherical_geometry/_modules/spherical_geometry/great_circle_arc.html#interpolate
        -https://www.astropy.org/acknowledging.html
"""
import numpy as np
from numpy.core.umath_tests import inner1d
from math import ceil, cos, sin, sqrt, pi, atan, asin, atan2, radians, degrees

def length(A, B, degrees=True):
    r"""
    Returns the angular distance between two points (in vector space)
    on the unit sphere.

    Parameters
    ----------
    A, B : (*x*, *y*, *z*) triples or Nx3 arrays of triples
       The endpoints of the great circle arc, in vector space.

    degrees : bool, optional
        If `True` (default) the result is returned in decimal degrees,
        otherwise radians.

    Returns
    -------
    length : scalar or array of scalars
        The angular length of the great circle arc.

    Notes
    -----
    The length is computed using the following:

    .. math::

       \Delta = \arccos(A \cdot B)
    """

    A = np.asanyarray(A)
    B = np.asanyarray(B)

    A2 = A ** 2.0
    Al = np.sqrt(np.sum(A2, axis=-1))
    B2 = B ** 2.0
    Bl = np.sqrt(np.sum(B2, axis=-1))

    A = A / np.expand_dims(Al, 2)
    B = B / np.expand_dims(Bl, 2)

    dot = inner1d(A, B)
    dot = np.clip(dot, -1.0, 1.0)
    with np.errstate(invalid='ignore'):
        result = np.arccos(dot)

    if degrees:
        return np.rad2deg(result)
    else:
        return result
    
def interpolate(A, B, steps=50):
    r"""
    Interpolate along the great circle arc.

    Parameters
    ----------
    A, B : (*x*, *y*, *z*) triples or Nx3 arrays of triples
        The endpoints of the great circle arc.  It is assumed thats
        these points are already normalized.

    steps : int
        The number of interpolation steps

    Returns
    -------
    array : (*x*, *y*, *z*) triples
        The points interpolated along the great circle arc

    Notes
    -----

    This uses Slerp interpolation where *Ω* is the angle subtended by
    the arc, and *t* is the parameter 0 <= *t* <= 1.

    .. math::

        \frac{\sin((1 - t)\Omega)}{\sin \Omega}A + \frac{\sin(t \Omega)}{\sin \Omega}B
    """
    steps = int(max(steps, 2))
    t = np.linspace(0.0, 1.0, steps, endpoint=True).reshape((steps, 1))

    omega = length(A, B, degrees=False)
    if omega == 0.0:
        offsets = t
    else:
        sin_omega = np.sin(omega)
        offsets = np.sin(t * omega) / sin_omega

    return offsets[::-1] * A + offsets * B


def cart2Geo(r, x, y, z):
    """
        Computes lat/long given cartesian coordinates and radius
    """
    
    lon = atan2(y, x) * 360 / (2 * pi)
    lat = asin(z / r) * 360 / (2 * pi)
#    lon = degrees(lon)
#    lat = degrees(lat)
    return lat, lon

def geo2Cart(r, lat, lon):
    """
        Computes cartesian coordinates given radius, lat, lon
    """
    lat = radians(lat)
    lon = radians(lon)
    x = r * cos(lat) * cos(lon)
    y = r * cos(lat) * sin(lon)
    z = r * sin(lat)
    return x, y, z


def pointCart2Geo(r, pt):
    lat, lon = cart2Geo(r, pt[0], pt[1], pt[2])
    return (lat, lon)

def pointGeo2Cart(r, pt):
    x, y, z = geo2Cart(r, pt[0], pt[1])
    return (x, y, z)
    
def bearing(pt1, pt2):
    """
        Computes bearing between two geographic coordinates
    """
    pt1 = list(pt1)
    pt2 = list(pt2)
    for i in range(len(pt1)):
        pt1[i] = radians(pt1[i])
        pt2[i] = radians(pt2[i])
        
    delta_lon = pt2[1] - pt1[1]
    bearing = atan2(sin(delta_lon) * cos(pt2[0]), cos(pt1[0]) * sin(pt2[0]) - sin(pt1[0]) * cos(pt2[0]) * cos(delta_lon))
    return degrees(bearing)

    
    
def pointDist(pt1, pt2):
    """
        Computes euclidean distance between two points
    """
    
    dist = 0
    for i in range(len(pt1)):
        dist += (pt2[i] - pt1[i])**2
    
    return sqrt(dist)
    
    
def destinationPoint(pt, bearing, dist, radius):
    """
        Computes final destination point given starting point, bearing, distance,
        and radius.
    
    """
    lat, lon = radians(pt[0]), radians(pt[1])
    bearing = radians(bearing)
    angular_dist = dist / radius
    lat2 = asin(sin(lat) * cos(angular_dist) + cos(lat)*sin(angular_dist)*cos(bearing))
    lon2 = lon + atan2(sin(bearing) * sin(angular_dist) * cos(lat), cos(angular_dist) - sin(lat) * sin(lat2))
    lat2 = degrees(lat2)
    lon2 = degrees(lon2)
    return [lat2, lon2] #(lon2+540)%360-180)
    
    
    
def regr(X):
    """
        Fit a line to points.
    
    """
    y= np.average(X, axis=0)
    Xm = X-y
    u, s, v = np.linalg.svd((1./X.shape[0])*np.matmul(Xm.T,Xm))
        
        # Extra Credit: Going back
    z= np.matmul(u[:,0].T, Xm.T)
    c = np.array([z*n for n in u[:,0]])
    d = np.array(y.tolist()*c.shape[1]).reshape(c.shape[1],-1).T
    e = (c+d).T
    return u,s,v  
    
    
def generateParametricPoints(direction_vec, cps_mean, rng):
    """
        Generates points on a line given parametric equations
    """
    
    cps = []
    x1, y1, z1 = direction_vec
    xt, yt, zt = cps_mean
    for t in rng:
        x = x1 + xt*t
        y = y1 + yt*t
        z = z1 + zt*t
        cps.append([x,y,z])
    return cps

def closest_point_sphere(center, points, radius):
    """
        Finds the closest point to a sphere given radius r, center.
        Returns the closest point and the euclidean distance
    """
    deltas = np.array(points) - np.array(center)
    dist_2 = np.einsum('ij,ij->i', deltas, deltas) - radius
    min_index = np.argmin(dist_2)
    return [points[min_index], dist_2[min_index], min_index]

def findClosestPoints(direction_vec, cps_mean, rng, radius, center, neighbor_point_dist=5):
    """
        Finds the point closest to globe (pt1). Also returns a subsequent point in that direction.
        ** MAKE SURE THE SUPPLIED RANGE IS SUFFICIENTLY LARGE FOR THE GIVEN GLOBE **
    """
    
    # Generate points on this line
    parametric_points = generateParametricPoints(direction_vec,cps_mean, rng)
    
    # Find points closest to globe
    _, _, closest_index = closest_point_sphere(center, parametric_points, radius)
    
    # Grab that point, and also a neighboring point
    if closest_index + neighbor_point_dist > len(parametric_points):
        neighbor_point_dist = 1
        print "ERROR: NEIGHBOR_POINT_DIST IS SET TOO HIGH, DEFAULTING TO 1"
    
    pt1 = parametric_points[closest_index]
    pt2 = parametric_points[closest_index + neighbor_point_dist]
    
    # Project on globe
    pt1 = projectOnSphere(pt1, radius)
    pt2 = projectOnSphere(pt2, radius)
    
    return pt1, pt2, parametric_points
    

def projectOnSphere(pt, radius):
    """
        Projects a point onto a sphere
    """
    
    # Compute length of point vector
    length = pointDist(pt, [0,0,0])
    
    # Scale to length of sphere
    projected_point = (radius / length) * np.array(pt)
    
    return projected_point.tolist()
    
    
    
    
    
    
    