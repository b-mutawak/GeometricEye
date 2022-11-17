"""
    Surface manipulation functions.
    
    Taken from:
        -https://github.com/bungnoid/glTools/blob/8ff0899de43784a18bd4543285655e68e28fb5e5/utils/surface.py

"""




import maya.cmds as mc
import maya.OpenMaya as OpenMaya


def isSurface(surface):
    '''
    Check if the specified object is a nurbs surface or transform parent of a surface
    @param surface: Object to query
    @type surface: str
    '''
    # Check object exists
    if not mc.objExists(surface): return False
    # Check shape
    if mc.objectType(surface) == 'transform': surface = mc.listRelatives(surface,s=True,ni=True,pa=True)[0]
    if mc.objectType(surface) != 'nurbsSurface': return False
    
    # Return result
    return True

def getSurfaceFn(surface):
    '''
    Create an MFnNurbsSurface class object from the specified nurbs surface
    @param surface: Surface to create function class for
    @type surface: str
    '''
    # Checks
    if not isSurface(surface): raise Exception('Object '+surface+' is not a valid surface!')
    if mc.objectType(surface) == 'transform':
        surface = mc.listRelatives(surface,s=True,ni=True,pa=True)[0]
    
    # Get MFnNurbsSurface
    selection = OpenMaya.MSelectionList()
    OpenMaya.MGlobal.getSelectionListByName(surface,selection)
    surfacePath = OpenMaya.MDagPath()
    selection.getDagPath(0,surfacePath)
    surfaceFn = OpenMaya.MFnNurbsSurface()
    surfaceFn.setObject(surfacePath)
    
    # Return result
    return surfaceFn

def chordLength(surface,param=0.0,direction='u'):
    '''
    Return the length of a surface isoparm given a parameter and a direction
    @param surface: Surface to query closest point from
    @type surface: str
    @param param: The parameter on the surface to query length of
    @type param: float
    @param direction: Direction along the surface to measure length of
    @type direction: str
    '''
    # Check surface
    if not isSurface(surface): raise Exception('Object '+surface+' is not a valid surface!')
    # Duplicate surface curve
    curve = mc.duplicateCurve(surface+'.'+direction+'['+str(param)+']',ch=0,rn=0,local=0)
    # Measure curve length
    length = mc.arclen(curve[0])
    # Cleanup
    mc.delete(curve)
    # Return result
    return length

def closestPoint(surface,pos=(0,0,0)):
    '''
    Return closest point on surface to target position
    @param surface: Surface to query closest point from
    @type surface: str
    @param pos: Position to query surface from
    @type pos: tuple/list
    '''
    # Check surface
    if not isSurface(surface): raise Exception('Object '+surface+' is not a valid surface!')
    
    # Get point world position
    pt = OpenMaya.MPoint(pos[0],pos[1],pos[2],1.0)
    pt2 = OpenMaya.MPoint()
    
    # Get surface function set
    surfFn = getSurfaceFn(surface)
    
    # Get uCoord and vCoord pointer objects
    uCoord = OpenMaya.MScriptUtil()
    uCoord.createFromDouble(0.0)
    uCoordPtr = uCoord.asDoublePtr()
    vCoord = OpenMaya.MScriptUtil()
    vCoord.createFromDouble(0.0)
    vCoordPtr = vCoord.asDoublePtr()
    
    # get closest uCoord to edit point position
    surfFn.closestPoint(pt,uCoordPtr,vCoordPtr,True,0.00000000001,OpenMaya.MSpace.kWorld)

    return (OpenMaya.MScriptUtil(uCoordPtr).asDouble(),OpenMaya.MScriptUtil(vCoordPtr).asDouble())


def intersect(surface,source,direction):
	'''
	Return the intersection point on a specified nurbs surface given a source point and direction
	@param surface: Nurbs surface to perform intersection on
	@type surface: str
	@param source: Source point for the intersection ray
	@type source: list or tuple or str
	@param direction: Direction of the intersection ray intersection
	@type direction: list or tuple
	'''
	# Get surfaceFn
	surfaceFn = getSurfaceFn(surface)
	# Get source point
	source = OpenMaya.MPoint(source[0], source[1], source[2], 1.0)
	# Get direction vector
	direction = OpenMaya.MVector(direction[0],direction[1],direction[2])
	
	# Get uCoord and vCoord pointer objects
	uCoord = OpenMaya.MScriptUtil()
	uCoord.createFromDouble(0.0)
	uCoordPtr = uCoord.asDoublePtr()
	vCoord = OpenMaya.MScriptUtil()
	vCoord.createFromDouble(0.0)
	vCoordPtr = vCoord.asDoublePtr()
    
	# Calculate intersection
	hitPt = OpenMaya.MPoint()
	hit = surfaceFn.intersect(source,direction,uCoordPtr,vCoordPtr,hitPt,0.0001,OpenMaya.MSpace.kWorld,False,None,None)
	if not hit:
		print 'No intersection found!'
		hitPt = OpenMaya.MPoint.origin
	
	# Return intersection hit point
	return [hitPt[0],hitPt[1],hitPt[2]]




