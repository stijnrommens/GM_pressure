import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

import inputs

rho = 997 # kg/m3, density of water
g = 9.81  # m/s2, gravitational constant

radius = 192e-3/2    # m, inside radius column
Ac = np.pi*radius**2 # m2, cross-sectional area column

sensor_height = 563e-3+16e-3 # m, height pressure sensor from sparger
liquid_height = 224e-3       # m, water height from sensor with no flow

f_probe_files = inputs.measurement_data()
f_probe_files = inputs.remove_column(f_probe_files, 2)


def data_processing(files):
    ''' From the measurements, extract information. '''
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
f_probe_results = data_processing(f_probe_files)


def plot_velocity_size(data):
    ''' Plot the measured bubble velocity and size. '''
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    ax1.plot(data[0], data[1], color="#69b3a2", lw=3)
    ax1.scatter(data[0], data[1], color="#69b3a2")
    ax2.plot(data[0], data[2], color="#3399e6", lw=3)
    ax2.scatter(data[0], data[2], color="#3399e6")
    ax1.set_xlabel("Superficial gas velocity [m/s]")
    ax1.set_ylabel("Bubble mean velocity [m/s]", color="#69b3a2", fontsize=14)
    ax1.tick_params(axis="y", labelcolor="#69b3a2")
    ax2.set_ylabel("d32 [mm]", color="#3399e6", fontsize=14)
    ax2.tick_params(axis="y", labelcolor="#3399e6")
    fig.suptitle("Bubble properties", fontsize=20)
    return plt.show()
f_probe_plot = [f_probe_results[:,0]/(Ac*60*1000), f_probe_results[:,1], f_probe_results[:,2]*1e3] # m/s, m/s, mm
plot = plot_velocity_size(f_probe_plot)

def plot_holdup(f_probe_data, p_sensor_data):
    ''' PLot the gas holdups from the fiber probe and the pressure sensor together. '''
    fig, ax1 = plt.subplots()
    ax1.plot(p_sensor_data[:,0], p_sensor_data[:,-1], color="#3399e6", lw=3, label='Pressure sensor')
    ax1.plot(f_probe_data[:,0], f_probe_data[:,-1], color="#69b3a2", lw=3, label='Fiber probe')
    ax1.set_xlabel("Gas flow [L/min]")
    ax1.set_ylabel("Gas holdup [-]", fontsize=14)
    ax1.tick_params(axis="y")
    fig.suptitle("Comparison", fontsize=20)
    plt.legend()
    return plt.show()
