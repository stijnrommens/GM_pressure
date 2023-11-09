import os
import numpy as np
import pandas as pd
from nptdms import TdmsFile
import matplotlib.pyplot as plt

import inputs
import calibration

files = inputs.measurement_data()
calibration_files = inputs.calibration_data()

rho = 997 # kg/m3, density of water
g = 9.81  # m/s2, gravitational constant

radius = 192e-3/2    # m, inside radius column
Ac = np.pi*radius**2 # m2, cross-sectional area column

sensor_height = 563e-3+16e-3 # m, height pressure sensor from sparger
liquid_height = 224e-3       # m, water height from sensor with no flow


def fiber_probe(files):
    ''' From the measurements, extract fiber probe information. '''
    
    results = []
    for file in files:
        param = file[0]     # height, flow, etc.
        file_name = file[1] # -
        
        # Load 1st file
        path = r'u:\Bubble Column\Data\A2 Fiber Probe\231101 - Flow variation in Water' + file_name + '.evt'
        df = pd.read_csv(path, sep='\t', decimal=',')

        df_valid = df[df.Valid == 1] # Only valid bubbles

        # Obtain velocity and size
        velocity = df_valid['Veloc']
        mean_velocity = sum(velocity)/len(velocity) # m/s

        size = df_valid['Size']*1e-6
        d32 = sum(size**3)/sum(size**2) # m
        
        # Load 2nd file
        path_stream = r'u:\Bubble Column\Data\A2 Fiber Probe\231101 - Flow variation in Water' + file_name + '_stream.evt'
        df_stream = pd.read_csv(path_stream, sep='\t', decimal=',')

        # Obtain gas holdup
        arrival = df_stream['Arrival']
        duration = df_stream['Duration']
        void_fraction = sum(duration)/arrival[len(arrival)-1]

        results.append([param, mean_velocity, d32, void_fraction])
    return np.array(results)
fiber_probe_results = fiber_probe(files)


def pressure_sensor(files, fit):
    ''' From the measured mean voltage and the fitted line, obtain the gas holdup. '''
    
    results = []
    for file in files:
        param = file[0]     # height, flow, etc.
        file_name = file[2] # -
        
        # Load file
        path = r'u:\Bubble Column\Data\PXM419' + file_name + '.tdms'
        loaded_file = TdmsFile(path)
        
        # Obtain voltage
        for group in loaded_file.groups():
            df = group.as_dataframe()
        voltage = df['Dev1/ai1']
        mean_voltage = sum(voltage)/len(voltage) # V
        
        gas_height = fit[0]*mean_voltage + fit[1]
        volume_L = Ac * sensor_height # m3
        volume_G = Ac * (gas_height-liquid_height) # m3
        holdup = volume_G/volume_L # -
        
        results.append([param, mean_voltage, gas_height, holdup])
    return np.array(results)
fit = calibration.calibration_fit(calibration_files)
pressure_sensor_results = pressure_sensor(files, fit)
