import sys
from read_roi import read_roi_zip
import os
import numpy as np
#import maya.cmds as cmds



def parseDataFolder(folder_path, globe_on_origin):
    """
        Recurses input folder_path looking for all ROI. 
        Reads ROI and returns the set of ROI for each muscle
    """
    
    # Make sure the path exists
    if not os.path.exists(folder_path):
        sys.stderr.write("Error, invalid input folder path")
        return
        
    # The set of expected gazes and planes
    eye_points = {'Coronal OD Prim':{}, 'Sagittal OD Prim':{},'Coronal OS Prim':{}, 'Sagittal OS Prim':{}}
    
    # Iterate through each set and try to find a folder that contains that data
    for eye_set in eye_points.keys():
        
        # Generate the expected path to folder
        set_path = os.path.join(folder_path, eye_set)
        if not os.path.exists(set_path):
            continue

    
        # Find each roi file
        rois = [os.path.join(set_path, file) for file in os.listdir(set_path) if file.endswith(".zip")]
        
        # Read each ROI set and parse the points
        for roi_set_path in rois:
        
            # Parse the zip file
            input_set_dict = read_roi_zip(roi_set_path)
            parsed_points = []
            
            # Iterate through each slice
            for slice_index, slice in enumerate(input_set_dict.keys()):
                data = input_set_dict[slice]
                
                # Iterate through each point
                for i in range(len(data['x'])):
                    
                    # Generate 3D points based on 2mm distance between each slice
                    parsed_points.append((data['x'][i], data['y'][i], 2 * (data['position']) - globe_on_origin))
                
            # Store the points
            muscle = os.path.basename(roi_set_path).split("_")[3].replace(".zip", "")
            
            # Mirror OS 
            if "OS" in eye_set:
                parsed_points = np.array(parsed_points, dtype=np.float)
                parsed_points = np.flipud(parsed_points)
                
            eye_points[eye_set][muscle] = parsed_points
            
            
    return eye_points
    
def generatePoints(eye_points):
    eye_set = eye_points['Coronal OD Prim']
    
    for muscle in eye_set.keys():
        
        points = eye_set[muscle]
        
        for point in points:
            locator = cmds.spaceLocator()
            cmds.move(point[0], point[1], point[2], locator[0])

    
    
if __name__=="__main__":
    folder_path = r"D:\New folder\Dropbox\Maya Model Generator\9 SOP Datasets\Q6 Anonymized Stacks"
    print(parseDataFolder(folder_path))

