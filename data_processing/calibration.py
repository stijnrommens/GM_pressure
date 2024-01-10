import os
import openpyxl
import numpy as np
import pandas as pd
from nptdms import TdmsFile


# --- Experiment ---
excel_file = r'\\tudelft.net\student-homes\R\srommens\My Documents\GitHub\GM_pressure\data_processing\input_file.xlsx'
calibration_path = pd.read_excel(excel_file, sheet_name='Calibration', header=None, index_col=0, nrows=1)
calibration_path = calibration_path.to_numpy()[0][0]
calibration_files = pd.read_excel(excel_file, sheet_name='Calibration', header=2)
calibration_files = calibration_files.to_numpy()


# --- Calibration ---
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


def export_fit():
    # print(fit)
    return fit
