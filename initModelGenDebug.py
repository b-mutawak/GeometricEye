import sys
import maya.OpenMaya as om
import maya.mel as mel
import maya.OpenMayaMPx as OpenMayaMPx
import os

kPluginCmdName = "__initDebug__"


# command
class scriptedCommand(OpenMayaMPx.MPxCommand):
    def __init__(self):
        OpenMayaMPx.MPxCommand.__init__(self)

    def doIt(self, argList):
        print("How do I manually source .mel scripts?")


# Creator
def cmdCreator():
    return OpenMayaMPx.asMPxPtr(scriptedCommand())


# Initialize the script plug-in
def initializePlugin(mobject):
    mplugin = OpenMayaMPx.MFnPlugin(mobject)
    try:
        mplugin.registerCommand(kPluginCmdName, cmdCreator)
    except:
        sys.stderr.write("Failed to register command: %s\n" % kPluginCmdName)
        raise
		
    path = "C:/Users/admin/Dropbox/Maya Model Generator/Project/Maya_New/Maya2020/bin/plug-ins/GeometricEye"
    if (path not in sys.path):
        upath = eval("u'"+path+"'")
        sys.path.append(upath)
        sys.stderr.write("Added Path\n")
	   
    for file in [file for file in os.listdir(path) if file.endswith(".py")]:
        mel.eval('source /GeometricEye/"' + file + '"')
		
		
# Uninitialize the script plug-in
def uninitializePlugin(mobject):
    mplugin = OpenMayaMPx.MFnPlugin(mobject)
    try:
        mplugin.deregisterCommand(kPluginCmdName)
    except:
        sys.stderr.write("Failed to unregister command: %s\n" % kPluginCmdName)
        raise