import numpy as np
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# --- Processed Data ---
import processing
fiber_probe_results = processing.export_fiber_probe()
pressure_sensor_results = processing.export_pressure_sensor()

x_axis = processing.mass2strength(fiber_probe_results[:,0], fiber_probe_results[:,-2], 1)
x_label = 'Concentration [M]'
# x_axis = fiber_probe_results[:,0]
# x_label = 'Weight [grams]'

y_axis1 = fiber_probe_results[:,1] # m/s, velocity
errorbar_axis1 = np.array((fiber_probe_results[:,2], fiber_probe_results[:,3])) # m/s, velocity
y_label1 = 'Velocity [m/s]'

y_axis2 = fiber_probe_results[:,4]*1e3 # m, diameter
# y_axis2 = processing.diameter2percent(fiber_probe_results[:,4]) # %, percentage
errorbar_axis2 = np.array((fiber_probe_results[:,5], fiber_probe_results[:,6]))*1e3 # m, diameter
y_label2 = 'Diameter [mm]'

y_axis3 = fiber_probe_results[:,-1] # -, gas holdup
y_axis4 = pressure_sensor_results[:,-1] # -, gas holdup
y_label3 = 'Gas holdup [-]'


# --- Regressed Data ---
import ctrans_regression


# --- Plotting ---
def size_vel_plot(x, y1, y1_error, y2, y2_error):
    ''' Plot size and velocity from processed fiber probe data. '''
    
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    # y1 = np.hstack((y1[0], y1))
    # y1 = np.hstack((y1[0]-0.01, y1))
    # print(y1)
    # x = np.hstack((x[0]+0.001, x))
    # x = np.hstack((x[0]+0.005, x))
    # print(x)
    # ax1.scatter( x, y1, color="#69b3a2", lw=3)
    # ax1.errorbar(x, y1, yerr=y1_error, fmt='-o', color="#69b3a2", capsize=3)
    ax2.scatter( x, y2, color="#3399e6", lw=3)
    ax2.errorbar(x, y2, yerr=y2_error, fmt='-o', color="#3399e6", capsize=3)
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label1, color="#69b3a2", fontsize=14)
    ax1.tick_params(axis="y", labelcolor="#69b3a2")
    # ax1.set_xscale('log')
    ax2.set_xscale('log')
    ax2.set_ylabel(y_label2, color="#3399e6", fontsize=14)
    ax2.tick_params(axis="y", labelcolor="#3399e6")
    # ax1.set_ylim(0.15, 0.35)
    # ax2.set_ylim(1,2.5)
    # ax1.set_xlim(-1,3)
    # fig.suptitle("Bubble properties", fontsize=20)
    plt.show()
    return
size_vel_plot(x_axis, y_axis1, errorbar_axis1, y_axis2, errorbar_axis2)


def holdup_plot(x, y1, y2):
    ''' Plot gas holdup from processed fiber probe and pressure sensor data. '''
    
    fig, ax1 = plt.subplots()
    ax1.plot(x, y1, color="#3399e6", lw=3, label='Pressure sensor')
    ax1.plot(x, y2, color="#69b3a2", lw=3, label='Fiber probe')
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label3, fontsize=14)
    ax1.tick_params(axis="y")
    ax1.set_ylim(0, 0.1)
    # ax1.set_xlim(0 , 80)
    # fig.suptitle("Gas holdup", fontsize=20)
    plt.legend()
    plt.show()
    return 
# holdup_plot(x_axis, y_axis4, y_axis3)

