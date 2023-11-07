# Common libraries
import os
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

# Own files
import inputs
import functions

# Import files from inputs
fiber_meas_dict = inputs.fiber_meas()
pressure_cal_dict = inputs.pressure_cal() 
pressure_meas_dict = inputs.pressure_meas()
rho = inputs.constants()[0]
g = inputs.constants()[1]
radius = inputs.constants()[2]
Ac = inputs.constants()[5]
liquid_height = inputs.constants()[3]
sensor_height = inputs.constants()[4]

# Send inputs to functions
fiber_meas = functions.fiber_measurement(fiber_meas_dict, radius) # Fiber probe
pressure_cal = functions.pressure_calibration(pressure_cal_dict, rho, g) # Pressure sensor
pressure_meas = functions.pressure_measurement(pressure_meas_dict, pressure_cal, radius, liquid_height, sensor_height) # Pressure sensor
# print(fiber_meas[:,-1])
# print(pressure_meas[:,-1])

# Plot fiber probe
# fig, ax1 = plt.subplots()
# ax2 = ax1.twinx()
# ax1.plot(fiber_meas[:,0]/(Ac*60*1000), fiber_meas[:,1], color="#69b3a2", lw=3)
# ax1.scatter(fiber_meas[:,0]/(Ac*60*1000), fiber_meas[:,1], color="#69b3a2")
# ax2.plot(fiber_meas[:,0]/(Ac*60*1000), fiber_meas[:,2]*1e3, color="#3399e6", lw=3)
# ax2.scatter(fiber_meas[:,0]/(Ac*60*1000), fiber_meas[:,2]*1e3, color="#3399e6")
# ax1.set_xlabel("Superficial gas velocity [m/s]")
# ax1.set_ylabel("Bubble mean velocity [m/s]", color="#69b3a2", fontsize=14)
# ax1.tick_params(axis="y", labelcolor="#69b3a2")
# ax2.set_ylabel("d32 [mm]", color="#3399e6", fontsize=14)
# ax2.tick_params(axis="y", labelcolor="#3399e6")
# fig.suptitle("Bubble properties", fontsize=20)
# plt.show()

# Plot pressure sensor
# fig, ax1 = plt.subplots()
# ax1.plot(pressure_cal[:,1], pressure_cal[:,0], color="#69b3a2", lw=3)
# ax1.scatter(pressure_meas[:,1], pressure_meas[:,3], color="#3399e6", lw=3)
# ax1.set_xlabel("Voltage [V]")
# ax1.set_ylabel("Height [m]", fontsize=14)
# ax1.tick_params(axis="y")
# fig.suptitle("Pressure sensor calibration", fontsize=20)
# plt.show()

# Plot gas hold ups
# fig, ax1 = plt.subplots()
# ax1.plot(pressure_meas[:,0], pressure_meas[:,-1], color="#3399e6", lw=3, label='Pressure sensor')
# ax1.plot(fiber_meas[:,0], fiber_meas[:,-1], color="#69b3a2", lw=3, label='Fiber probe')
# ax1.set_xlabel("Gas flow [L/min]")
# ax1.set_ylabel("Gas holdup [-]", fontsize=14)
# ax1.tick_params(axis="y")
# fig.suptitle("Comparison", fontsize=20)
# plt.legend()
# plt.show()
