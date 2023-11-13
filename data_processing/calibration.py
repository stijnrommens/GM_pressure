import os
import numpy as np
from nptdms import TdmsFile

import inputs
calibration_files, calibration_path = inputs.calibration_data() 

def calibration_fit(files, folder):
    ''' Fit a linear line to the calibration data. '''
    
    # Extract data
    data_to_fit = []
    for file in files:
        height = file[0]     # height, flow, etc.
        file_name = file[1] # -
        
        # Load file
        path = folder + file_name + '.tdms'
        loaded_file = TdmsFile(path)
        
        # Obtain voltage
        for group in loaded_file.groups():
            df = group.as_dataframe()
        voltage = df['Dev1/ai1']
        mean_voltage = np.mean(voltage) # V
        data_to_fit.append([height, mean_voltage])
    data_to_fit = np.array(data_to_fit)


    # Fit linear line over data
    x, y = data_to_fit[:,1], data_to_fit[:,0]
    degree = 1
    fit = np.polyfit(x, y, degree)
    return fit
fit = calibration_fit(calibration_files, calibration_path)
# print(fit)
