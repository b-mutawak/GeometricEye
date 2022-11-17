# GeometricEye

This library developed for Autodesk Maya generates 3D ocular models using MRI data. The following are instructions for installation, usage, and expected inputs.

## Installation

1. Install a version of [Autodesk Maya](https://www.autodesk.com/products/maya). This library was developed on a Windows 10 machine using Autodesk Maya 2020. While other versions and operating systems may work, there is no guarentee. Specifically, Maya 2020 uses the Python 2.7 interpretor while newer versions use Python 3.0. This alone will create incompatability issues (although they can be resolved). If Maya 2020 cannot be installed on your machine, I suggest looking into the source code and modifying the elements with incompatabilities (this will likely be src/read_roi.py or src/parseInputPoints.py)

2. After Maya is correctly installed (I recommend running the Maya program once for startup), two python libraries must be installed into the Maya python (mayapy) interpretor: numpy and scipy. To do this, we will use the mayapy.exe found in the install directory for Maya. In Windows, this would be *C:\Program Files\Autodesk\Maya<VersionNumber>\bin\mayapy.exe*. For other operating systems, visit the Maya documentation [here](https://knowledge.autodesk.com/support/maya/learn-explore/caas/CloudHelp/cloudhelp/2022/ENU/Maya-Scripting/files/GUID-72A245EC-CDB4-46AB-BEE0-4BBBF9791627-htm.html).

In a command prompt (or terminal), we install the two packages in the following way:
```console
  "C:\Users\admin\Dropbox\Maya Model Generator\Project\Maya_New\Maya2020\bin\mayapy.exe" -m pip install numpy scipy
```
  
Once completed, you should see scipy (v0.19.1) and numpy (v1.13.1) installed in the mayapy interpretor. You can check this by listing the packages installed in mayapy.
  
3. Once your packages are installed, the next step is to install the source files into Maya. 


