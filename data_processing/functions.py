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

        pressure = rho*g*height # Pa

        results.append([height, mean_voltage, pressure])
    results = np.array(results)
    return results

def pressure_measurement(files_dict, cal_results, radius, liquid_height, sensor_height):
    results= []
    volume_L = np.pi * radius**2 * sensor_height # m3

    for flow in files_dict:
        path = r'u:\Bubble Column\Data\PXM419' + files_dict[flow]
        file = TdmsFile(path)
        for group in file.groups():
            df = group.as_dataframe()
        voltage = df['Dev1/ai1']
        mean_voltage = sum(voltage)/len(voltage) # V
        
        gas_height = np.interp(mean_voltage, cal_results[:,1], cal_results[:,0]) # m
        pressure = np.interp(mean_voltage, cal_results[:,1], cal_results[:,2]) # Pa
        
        volume_G =  np.pi * radius**2 * (gas_height-liquid_height) # m3
        hold_up = volume_G/volume_L # -
        
        results.append([flow, mean_voltage, pressure, gas_height, hold_up])
    results = np.array(results)
    return results

def fiber_measurement(files_dict, radius):
    results = []
    Ac = np.pi * radius**2 # m2
    for flow in files_dict:
        path = r'u:\Bubble Column\Data\A2 Fiber Probe\231101 - Flow variation in Water' + files_dict[flow] + '.evt'
        df = pd.read_csv(path, sep='\t', decimal=',') # Initialize dataframe

        df_valid = df[df.Valid == 1] # Make df with only valid bubbles

        velocity = df_valid['Veloc']
        mean_velocity = sum(velocity)/len(velocity) # m/s

        size = df_valid['Size']*1e-6
        mean_size = sum(size)/len(size) # m
        d32 = sum(size**3)/sum(size**2) # m
        
        path_stream = r'u:\Bubble Column\Data\A2 Fiber Probe\231101 - Flow variation in Water' + files_dict[flow] + '_stream.evt'
        df_stream = pd.read_csv(path_stream, sep='\t', decimal=',') # Initialize dataframe

        arrival = df_stream['Arrival']
        duration = df_stream['Duration']
        void_fraction = sum(duration)/arrival[len(arrival)-1]

        results.append([flow, mean_velocity, d32, void_fraction])
    results = np.array(results)
    return results

# %%
# import os
# import numpy as np
# import pandas as pd

# path_stream = r'u:\Bubble Column\Data\A2 Fiber Probe\231101 - Flow variation in Water\2023-11-01T150555_stream.evt'
# df_stream = pd.read_csv(path_stream, sep='\t', decimal=',') # Initialize dataframe
# arrival = df_stream['Arrival']
# duration = df_stream['Duration']

# lis=[]
# for i, ar in enumerate(arrival):
#     if i+1 == len(arrival):
#         break
#     else:
#         frac = (duration[i]/(arrival[i+1]-arrival[i]))
#         lis.append(frac)
#         # print(frac)

# print(lis)
# print(sum(lis)/len(lis))

# void_fraction = sum(duration)/arrival[len(arrival)-1]
# print(void_fraction)
# %%
