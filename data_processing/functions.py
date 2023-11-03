import os
import numpy as np
import pandas as pd
from nptdms import TdmsFile

def pressure_calibration(files_dict, rho, g):
    results = []
    for height in files_dict:
        path = r'u:\Bubble Column\Data\PXM419' + files_dict[height]
        file = TdmsFile(path)
        
        for group in file.groups():
            df = group.as_dataframe()
        voltage = df['Dev1/ai1']
        mean_voltage = sum(voltage)/len(voltage) # V
            
        pressure = rho*g*height/1000 /100 # mBar

        results.append([height, mean_voltage, pressure])
    results = np.array(results)
    return results

def pressure_measurement(files_dict, cal_results):
    results= []
    for flow in files_dict:
        path = r'u:\Bubble Column\Data\PXM419' + files_dict[flow]
        file = TdmsFile(path)
        
        for group in file.groups():
            df = group.as_dataframe()
        voltage = df['Dev1/ai1']
        mean_voltage = sum(voltage)/len(voltage) # V
        pressure = np.interp(mean_voltage, cal_results[:,1], cal_results[:,2]) # mBar
        results.append([flow, mean_voltage, pressure])
    results = np.array(results)
    return results

def fiber_measurement(files_dict):
    results = []
    for flow in files_dict:
        path = r'u:\Bubble Column\Data\A2 Fiber Probe\231101 - Flow variation in Water' + files_dict[flow]
        df = pd.read_csv(path, sep='\t', decimal=',') # Initialize dataframe

        df_valid = df[df.Valid == 1] # Make df with only valid bubbles

        velocity = df_valid['Veloc']
        mean_velocity = sum(velocity)/len(velocity) # m/s

        size = df_valid['Size']
        mean_size = sum(size)/len(size)/1000 # mm
        
        results.append([flow, mean_velocity, mean_size])
    results = np.array(results)
    return results

# %%
