"""
    Set of functions for various curve manipulations.
    
    Taken from:
        https://github.com/bungnoid/glTools/blob/8ff0899de43784a18bd4543285655e68e28fb5e5/utils/curve.py
"""

import maya.mel as mm
import maya.cmds as cmds
import maya.OpenMaya as OpenMaya


# create color lookup dictionary
colorDict = {	'grey':0,
					'black':1,
					'dark grey':2,
					'light grey':3,
					'maroon':4,
					'dark blue':5,
					'blue':6,
					'dark green':7,
					'dark purple':8,
					'purple':9,
					'brown':10,
					'dark brown':11,
					'dull red':12,
					'red':13,
					'green':14,
					'dull blue':15,
					'white':16,
					'yellow':17,
					'light blue':18,
					'aqua':19,
					'pink':20,
					'pale orange':21,
					'pale yellow':22,
					'pale green':23,
					'pale brown':24,
					'dull yellow':25,
					'dull green':26,
					'dull aqua':27	}

def isCurve(curve):
	'''
	Check if the specified object is a nurbs curve or transform parent of a curve
	
	@param curve: Object to query
	@type curve: str
	
	@return: Boolean value indicating if the input object is a valid nurbs curve
	@returnType: bool
	'''
	# Check object exists
	if not cmds.objExists(curve): return False
	# Check shape
	if cmds.objectType(curve) == 'transform': curve = cmds.listRelatives(curve,s=True,ni=True,pa=True)[0]
	if cmds.objectType(curve) != 'nurbsCurve': return False
	
	# Return result
	return True

def getCurveFn(curve):
	'''
	Create an MFnNurbsCurve class object from the specified nurbs curve
	
	@param curve: Curve to create function class for
	@type curve: str
	
	@return: An MFnNurnrsCurve function class initialized with the input curve
	@returnType: MFnNurnrsCurve
	'''
	# Checks
	if not isCurve(curve): raise Exception('Object '+curve+' is not a valid curve!')
	
	# Get shape
	if cmds.objectType(curve) == 'transform':
		curve = cmds.listRelatives(curve,s=True,ni=True)[0]
	
	# Get MFnNurbsCurve
	curveSel = OpenMaya.MSelectionList()
	OpenMaya.MGlobal.getSelectionListByName(curve,curveSel)
	curvePath = OpenMaya.MDagPath()
	curveSel.getDagPath(0,curvePath)
	curveFn = OpenMaya.MFnNurbsCurve(curvePath)
	
	# Return result
	return curveFn

def getParamFromLength(curve,length):
	'''
	Get the U parameter fro the specified curve at a given length
	@param curve: Curve to get parameter from
	@type curve: str
	@param curve: Length along curve to sample parameter
	@type curve: float
	'''
	# Check curve
	if not isCurve(curve):
		raise Exception('Object "'+curve+'" is not a valid nurbs curve!')
	
	# Check length
	max_len = cmds.arclen(curve)
	if length > max_len:
		print('Input length is greater than actual curve length. Returning maximum U value!')
		return cmds.getAttr(curve+'.maxValue')
	
	# Get curve function set
	curveFn = getCurveFn(curve)
	
	# Get parameter from length
	param = curveFn.findParamFromLength(length)
	
	# Return result
	return param

def colorize(obj,color=None):
	'''
	Set override color for a specified object
	@param obj: The dag node to set the color for
	@type obj: str
	@param color: Dsiplay override colour. If None, colour is selected based on name prefix.
	@type color: str or int
	'''
	# Check object
	if not cmds.objExists(obj):
		raise Exception('Object "'+obj+'" does not exist!!')

	
	# Enable Overrides
	cmds.setAttr(obj+'.overrideEnabled',1)
	
	# Set Color
	if type(color) == str:
		if colorDict.has_key(color): cmds.setAttr(obj+'.overrideColor',colorDict[color])
		else: raise Exception('Color "'+color+'" is not recognised!!')
	elif type(color) == int:
		cmds.setAttr(obj+'.overrideColor',color)
	else:
		raise Exception('No valid color value supplied!!')