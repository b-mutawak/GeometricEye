# GeometricEye

This library developed for Autodesk Maya generates 3D ocular models using MRI data. The following are instructions for installation, usage, and expected inputs.

## Installation

1. Install a version of [Autodesk Maya](https://www.autodesk.com/products/maya). This library was developed on a Windows 10 machine using Autodesk Maya 2020. While other versions and operating systems may work, there is no guarentee. Specifically, Maya 2020 uses the Python 2.7 interpretor while newer versions use Python 3.0. This alone will create incompatability issues (although they can be resolved). If Maya 2020 cannot be installed on your machine, I suggest looking into the source code and modifying the elements with incompatabilities (this will likely be src/read_roi.py or src/parseInputPoints.py)

2. After Maya is correctly installed (I recommend running the Maya program once for startup), two python libraries must be installed into the Maya python (mayapy) interpretor: numpy and scipy. To do this, we will use the mayapy.exe found in the install directory for Maya. In Windows, this would be *C:\Program Files\Autodesk\Maya<VersionNumber>\bin\mayapy.exe*. For other operating systems, visit the Maya documentation [here](https://knowledge.autodesk.com/support/maya/learn-explore/caas/CloudHelp/cloudhelp/2022/ENU/Maya-Scripting/files/GUID-72A245EC-CDB4-46AB-BEE0-4BBBF9791627-htm.html).

* In a command prompt (or terminal), we install the two packages in the following way:
```console
  "C:\Users\admin\Dropbox\Maya Model Generator\Project\Maya_New\Maya2020\bin\mayapy.exe" -m pip install numpy scipy
```
  
* Once completed, you should see scipy (v0.19.1) and numpy (v1.13.1) installed in the mayapy interpretor. You can check this by listing the packages installed in mayapy.
  
3. Once your packages are installed, the next step is to install the source files into Maya. Navigate to the plug-ins folder of your Maya install (seperate from the above file path). On Windows, this can be: *C:\Path\To\Maya2020\bin\plug-ins*. Copy the GeometricEye folder and paste it into this plug-ins folder. Next, copy the "initModelGenDebug.py" file into the plug-ins folder. At this point, you should have both of these locations: *C:\Path\To\Maya2020\bin\plug-ins\initModelGenDebug.py* and *C:\Path\To\Maya2020\bin\plug-ins\GeometricEye*.

4. In your *initModelGenDebug.py* file, edit line 33 to be the FULL path to your GeometricEye folder. Make sure to replace all back-slashes (\) with forward slashes (/).

5. The next step is to load the plug-in in Maya. In Maya, go to Windows -> Settings/Preferences -> Plug-in Manager. If you copied all files correctly, you should be able to search for "initModelGenDebug.py" in the search bar and see a result. Make sure to check both "load" and "auto-load" options for this. Once completed, you have successfuly installed the GeometricEye plug-in. Note: you may see some warnings or errors in the console log but that does not impact the application. 

![image](https://user-images.githubusercontent.com/46249629/202595745-9d19138e-444a-40c4-baf7-8a1a045c8907.png)

## Usage

1. In Maya, open a script editor window using: Windows -> General Editors -> Script Editor. 

2. Click the "+" icon in the center of the window to create a new Script. Make sure to select "Python" as the language. 

3. In the Python script, paste the following code snippet:

```Python
import registration_utils

# Generate Eye
folder_path = r"C:\Users\admin\Dropbox\Maya Model Generator\9 SOP Datasets\Q9 Anonymized Stacks"
muscles = registration_utils.generateEye(folder_path)
```

4. Modify the "folder_path" variable to be the full path to the parent directory of all image planes / stacks (see expected inputs section). Once done, click the "Execute" button on the toolbar to execute the script. This may take some time to run. Once successfully completed, you should see several statistics printed in the console window and a model generated in the scene.

![image](https://user-images.githubusercontent.com/46249629/202597526-77ddbe73-aa66-4d1f-8c43-272004714f52.png)

5. Check the scene in Maya to ensure all components were properly generated. The eye models will be in the groups "OS" and "OD". All other items are not necessary. You may need to unhide the eye models for them to appear.

6. Once satisfied with the generated features, create a new Python script with the following code snippet to export 




