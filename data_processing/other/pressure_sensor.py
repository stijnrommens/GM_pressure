import os
import numpy as np
import matplotlib.pyplot as plt
from nptdms import TdmsFile

import inputs

rho = 997 # kg/m3, density of water
g = 9.81  # m/s2, gravitational constant

radius = 192e-3/2    # m, inside radius column
Ac = np.pi*radius**2 # m2, cross-sectional area column

sensor_height = 563e-3+16e-3 # m, height pressure sensor from sparger
liquid_height = 224e-3       # m, water height from sensor with no flow

calibration_files = inputs.calibration_data() 
p_sensor_files = inputs.measurement_data()
p_sensor_files = inputs.remove_column(p_sensor_files, 1)


def get_voltage(files):
    ''' Get the mean voltage of the measurement file. '''
    results = []
    for file in files:
        param = file[0]     # height, flow, etc.
        file_name = file[1] # -
        
        # Load file
        path = r'u:\Bubble Column\Data\PXM419' + file_name + '.tdms'
        loaded_file = TdmsFile(path)
        
        # Obtain voltage
        for group in loaded_file.groups():
            df = group.as_dataframe()
        voltage = df['Dev1/ai1']
        mean_voltage = sum(voltage)/len(voltage) # V

        results.append([param, mean_voltage])
    return np.array(results)
calibration_voltage = get_voltage(calibration_files)
p_sensor_voltage = get_voltage(p_sensor_files)


def calibration_fit(data):
    ''' Fit the calibration data as a linear line. '''
    x, y = data[:,1], data[:,0]
    degree = 1
    fit = np.polyfit(x, y, degree)
    return fit
fit = calibration_fit(calibration_voltage)


def get_holdup(data, fit):
    ''' From the measured mean voltage and the fitted line, obtain the gas holdup. '''
    param = data[:,0]
    mean_voltage = data[:,1]
    gas_height = fit[0]*mean_voltage + fit[1]
    
    volume_L = Ac * sensor_height # m3
    volume_G = Ac * (gas_height-liquid_height) # m3
    holdup = volume_G/volume_L # -

    return param, mean_voltage, gas_height, holdup
p_sensor_holdup = get_holdup(p_sensor_voltage, fit)


def plot_voltage_height(calibration_data, measurement_data, fit):
    ''' Plot calibration line with measurements. '''
    calibration_voltage = np.linspace(0, 10, 11)
    calibration_height = fit[0]*calibration_voltage + fit[1]

    fig, ax1 = plt.subplots()
    ax1.plot(calibration_voltage, calibration_height, color="#69b3a2", lw=3, zorder=1)
    ax1.scatter(calibration_data[0], calibration_data[1], color="#69b3a2", lw=1, label='Calibration', zorder=2)
    ax1.scatter(measurement_data[0], measurement_data[1], color="#3399e6", lw=1, label='Measurement', zorder=2)
    ax1.set_xlabel("Voltage [V]")
    ax1.set_ylabel("Height [m]", fontsize=14)
    ax1.tick_params(axis="y")
    fig.suptitle("Pressure sensor calibration", fontsize=20)
    plt.legend()
    return plt.show()
calibration_plot = [calibration_voltage[:,1], calibration_voltage[:,0]]
p_sensor_plot = [p_sensor_holdup[1], p_sensor_holdup[2]]
plot = plot_voltage_height(calibration_plot, p_sensor_plot, fit)
