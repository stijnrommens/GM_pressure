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

# Send inputs to functions
fiber_meas = functions.fiber_measurement(fiber_meas_dict) # Fiber probe
pressure_cal = functions.pressure_calibration(pressure_cal_dict, rho, g) # Pressure sensor
pressure_meas = functions.pressure_measurement(pressure_meas_dict, pressure_cal) # Pressure sensor
print(fiber_meas)
print(pressure_meas)

# Plot fiber probe
fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(fiber_meas[:,0], fiber_meas[:,1], color="#69b3a2", lw=3)
ax2.plot(fiber_meas[:,0], fiber_meas[:,2], color="#3399e6", lw=3)
ax1.set_xlabel("Gas flow [L/min]")
ax1.set_ylabel("Bubble mean velocity [m/s]", color="#69b3a2", fontsize=14)
ax1.tick_params(axis="y", labelcolor="#69b3a2")
ax2.set_ylabel("Bubble mean size [mm]", color="#3399e6", fontsize=14)
ax2.tick_params(axis="y", labelcolor="#3399e6")
fig.suptitle("Bubble properties", fontsize=20)
plt.show()

# Plot pressure sensor
fig, ax1 = plt.subplots()
ax1.plot(pressure_cal[:,1], pressure_cal[:,2], color="#69b3a2", lw=3)
ax1.scatter(pressure_meas[:,1], pressure_meas[:,2], color="#3399e6", lw=3)
ax1.set_xlabel("Voltage [V]")
ax1.set_ylabel("Pressure [mbar]", fontsize=14)
ax1.tick_params(axis="y")
fig.suptitle("Pressure calibration", fontsize=20)
plt.show()
