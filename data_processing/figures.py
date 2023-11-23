import numpy as np
import matplotlib.pyplot as plt

import processing

concentration = processing.mass2conc()

fiber_probe_results = processing.export_fiber_probe()
velocity_dispersion = np.array((fiber_probe_results[:,2], fiber_probe_results[:,3])) # m/s
size_dispersion = np.array((fiber_probe_results[:,5], fiber_probe_results[:,6])) # m
pressure_sensor_results = processing.export_pressure_sensor()

def size_vel_plot(x, y1, y1_error, y2, y2_error):
    ''' Plot size and velocity from processed fiber probe data. '''
    
    fig, ax1 = plt.subplots()
    ax2 = ax1.twinx()
    # ax1.plot(    x, y1, color="#69b3a2", lw=3)
    # ax1.errorbar(x, y1, yerr=y1_error, fmt='-o', color="#69b3a2", capsize=3)
    ax2.plot(    x, y2, color="#3399e6", lw=3)
    ax2.errorbar(x, y2, yerr=y2_error, fmt='-o', color="#3399e6", capsize=3)
    ax1.set_xlabel("Concentration [M]")
    ax1.set_ylabel("Velocity [m/s]", color="#69b3a2", fontsize=14)
    ax1.tick_params(axis="y", labelcolor="#69b3a2")
    # ax1.set_xscale('symlog')
    # ax2.set_xscale('symlog')
    # ax1.set_ylim(0.15, 0.45)
    ax2.set_ylabel("Diameter [mm]", color="#3399e6", fontsize=14)
    ax2.tick_params(axis="y", labelcolor="#3399e6")
    # ax2.set_ylim(0,3.5)
    fig.suptitle("Bubble properties", fontsize=20)
    # plt.plot(x, y2)
    return plt.show()
size_vel_plot(concentration, fiber_probe_results[:,1], velocity_dispersion, fiber_probe_results[:,4]*1e3, size_dispersion*1e3)


def holdup_plot(x, y1, y2):
    ''' Plot gas holdup from processed fiber probe and pressure sensor data. '''
    
    fig, ax1 = plt.subplots()
    ax1.plot(x, y1, color="#3399e6", lw=3, label='Pressure sensor')
    ax1.plot(x, y2, color="#69b3a2", lw=3, label='Fiber probe')
    ax1.set_xlabel("Concentration [M]")
    ax1.set_ylabel("Gas holdup [-]", fontsize=14)
    ax1.tick_params(axis="y")
    ax1.set_ylim(0, 0.1)
    fig.suptitle("Gas holdup", fontsize=20)
    plt.legend()
    return plt.show()
holdup_plot(concentration, pressure_sensor_results[:,-1], fiber_probe_results[:,-1])
