import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt


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

y_axis2 = fiber_probe_results[:,4]*1e3 # mm, diameter
# y_axis2 = processing.diameter2percent(fiber_probe_results[:,4]) # %, percentage
errorbar_axis2 = np.array((fiber_probe_results[:,5], fiber_probe_results[:,6]))*1e3 # mm, diameter
y_label2 = 'Diameter [mm]'


# --- Regression ---
def ctrans_regression(x, y):
    ''' Make a regression of the varied salt concentration measurements from the fiber probe. '''
    
    def func(x, a, b, c, d): 
        return (1 + np.exp(-b * (x - c)))/a + d # inverse logistic function
    popt, pcov = curve_fit(func, x, y)

    X = np.geomspace(1e-3, 0.5, num=101)
    # X = np.linspace(0, 0.5, num=401)
    Y = func(X, *popt)
    
    residuals = y - func(x, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum( (y-np.mean(y))**2 )
    r_squared = 1 - (ss_res/ss_tot)
    print(r_squared)
    return X, Y
regression_data = ctrans_regression(x_axis, y_axis2)


# --- Plotting ---
def size_vel_plot(x, y, y_error, X, Y):
    ''' Plot size and velocity from processed fiber probe data. '''
    
    fig, ax1 = plt.subplots()
    ax1.plot(X, Y, color="#3399e6", lw=3)
    ax1.scatter( x, y, color="#69b3a2", lw=3)
    # ax1.errorbar(x, y, yerr=y_error, fmt='o', color="#69b3a2", capsize=3) # color="#3399e6"
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label2, fontsize=14)#, color="#69b3a2"
    ax1.tick_params(axis="y")#, labelcolor="#69b3a2")
    ax1.set_xscale('log')
    # ax1.set_ylim(0.15, 0.35)
    # ax1.set_xlim(1e-3,0.5)
    # fig.suptitle("Bubble properties", fontsize=20)
    plt.show()
    return
size_vel_plot(x_axis, y_axis2, errorbar_axis2, regression_data[0], regression_data[1])


def export_measurements():
    return x_axis, y_axis1, errorbar_axis1, y_axis2, errorbar_axis2

def export_regression():
    return regression_data