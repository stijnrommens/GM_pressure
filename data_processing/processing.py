import os
import openpyxl
import numpy as np
import pandas as pd
from nptdms import TdmsFile
import matplotlib.pyplot as plt

import calibration

excel_file = r'\\tudelft.net\student-homes\R\srommens\My Documents\GitHub\GM_pressure\data_processing\input_file.xlsx'
fiber_probe_path = pd.read_excel(excel_file, sheet_name='Measurements', header=None, index_col=0, nrows=1)
fiber_probe_path = fiber_probe_path.to_numpy()[0][0]
pressure_sensor_path = pd.read_excel(excel_file, sheet_name='Measurements', header=0, index_col=0, nrows=1)
pressure_sensor_path = pressure_sensor_path.to_numpy()[0][0]
files = pd.read_excel(excel_file, sheet_name='Measurements', header=3)
files = files.to_numpy()

calibration_path = pd.read_excel(excel_file, sheet_name='Calibration', header=None, index_col=0, nrows=1)
calibration_path = calibration_path.to_numpy()[0][0]
calibration_files = pd.read_excel(excel_file, sheet_name='Calibration', header=2)
calibration_files = calibration_files.to_numpy()

rho = 997 # kg/m3, density of water
g = 9.81  # m/s2, gravitational constant

# radius = 192e-3/2    # m, inside radius column
radius = 150e-3/2    # m, inside radius column
Ac = np.pi*radius**2 # m2, cross-sectional area column

# sensor_height = 563e-3+16e-3 # m, height pressure sensor from sparger
sensor_height = 752e-3 # m, height pressure sensor from sparger
# liquid_height = 224e-3       # m, water height from sensor with no flow
liquid_height = 115e-3       # m, water height from sensor with no flow


def fiber_probe(files, folder):
    ''' From the measurements, extract fiber probe information. '''
    
    results = []
    for file in files:
        param = file[0]     # height, flow, etc.
        file_name = file[1] # -
        
        # Load 1st file
        path = folder + file_name + '.evt'
        df = pd.read_csv(path, sep='\t', decimal=',')

        df_valid = df[df.Valid == 1] # Only valid bubbles

        # Obtain velocity and size
        velocity = df_valid['Veloc']
        mean_velocity = np.mean(velocity) # m/s
        std_velocity = np.std(velocity)   # m/s

        size = df_valid['Size']*1e-6
        d32 = sum(size**3)/sum(size**2)   # m
        # d32 = np.mean(size)
        std_size = np.std(np.array(size)) # m

        # Load 2nd file
        path_stream = folder + file_name + '_stream.evt'
        df_stream = pd.read_csv(path_stream, sep='\t', decimal=',')

        # Obtain gas holdup
        arrival = df_stream['Arrival']
        duration = df_stream['Duration']
        void_fraction = sum(duration)/arrival[len(arrival)-1]

        results.append([param, mean_velocity, std_velocity, d32, std_size, void_fraction])
    return np.array(results)
fiber_probe_results = fiber_probe(files, fiber_probe_path)


def pressure_sensor(files, fit, folder):
    ''' From the measured mean voltage and the fitted line, obtain the gas holdup. '''
    
    results = []
    for file in files:
        param = file[0]     # height, flow, etc.
        file_name = file[2] # -
        
        # Load file
        path = folder + file_name + '.tdms'
        loaded_file = TdmsFile(path)
        
        # Obtain voltage
        for group in loaded_file.groups():
            df = group.as_dataframe()
        voltage = df['Dev1/ai1']
        mean_voltage = np.mean(voltage) # V
        
        gas_height = fit[0]*mean_voltage + fit[1]
        volume_L = Ac * sensor_height # m3
        volume_G = Ac * (gas_height-liquid_height) # m3
        holdup = volume_G/volume_L # -
        
        results.append([param, mean_voltage, gas_height, holdup])
    return np.array(results)
fit = calibration.calibration_fit(calibration_files, calibration_path)
pressure_sensor_results = pressure_sensor(files, fit, pressure_sensor_path)


# Plot velocity and size
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(   fiber_probe_results[:,0]/(Ac*60*1000), fiber_probe_results[:,1],     color="#69b3a2", lw=3)
ax1.errorbar(fiber_probe_results[:,0]/(Ac*60*1000), fiber_probe_results[:,1], yerr=fiber_probe_results[:,2], fmt='-o', color="#69b3a2")
ax2.plot(   fiber_probe_results[:,0]/(Ac*60*1000), fiber_probe_results[:,3]*1e3, color="#3399e6", lw=3)
ax2.errorbar(fiber_probe_results[:,0]/(Ac*60*1000), fiber_probe_results[:,3]*1e3, yerr=fiber_probe_results[:,4]*1e3, fmt='-o', color="#3399e6")
ax1.set_xlabel("Superficial gas velocity [m/s]")
ax1.set_ylabel("Bubble mean velocity [m/s]", color="#69b3a2", fontsize=14)
ax1.tick_params(axis="y", labelcolor="#69b3a2")
ax2.set_ylabel("d32 [mm]", color="#3399e6", fontsize=14)
ax2.tick_params(axis="y", labelcolor="#3399e6")
fig.suptitle("Bubble properties", fontsize=20)
plt.show()


# Plot gas holdup
fig, ax1 = plt.subplots()
ax1.plot(pressure_sensor_results[:,0]/(Ac*60*1000), pressure_sensor_results[:,-1], color="#3399e6", lw=3, label='Pressure sensor')
ax1.plot(fiber_probe_results[:,0]/(Ac*60*1000),     fiber_probe_results[:,-1],     color="#69b3a2", lw=3, label='Fiber probe')
ax1.set_xlabel("Superficial gas velocity [m/s]")
ax1.set_ylabel("Gas holdup [-]", fontsize=14)
ax1.tick_params(axis="y")
fig.suptitle("Gas holdup", fontsize=20)
plt.legend()
plt.show()