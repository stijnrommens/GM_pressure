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

# errorbar_axis2 = y_axis2+errorbar_axis2
# y_axis2 = np.pi/6 * (y_axis2)**3 # mm3, volume
# errorbar_axis2 = np.pi/6 * (errorbar_axis2)**3 - y_axis2 # mm3, volume
# y_label2 = 'Volume [mm3]'
# y_axis2 = y_axis2/y_axis2[-1]

y_axis3 = fiber_probe_results[:,-1] # -, gas holdup
y_axis4 = pressure_sensor_results[:,-1] # -, gas holdup
y_label3 = 'Gas holdup [-]'


# --- Regression ---
def ctrans_regression(x, y):
    ''' Make a regression of the varied salt concentration measurements from the fiber probe. '''
    
    def func(x, a, b, c, d): 
        return a / (1 + np.exp(-b * (x - c))) + d
        # return (1 + np.exp(-b * (x - c)))/a + d # inverse logistic function
    popt, pcov = curve_fit(func, x, y)

    X = np.geomspace(1e-5, 1, num=201)
    # X = np.linspace(0, 0.5, num=401)
    Y = func(X, *popt)
    
    residuals = y - func(x, *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum( (y-np.mean(y))**2 )
    r_squared = 1 - (ss_res/ss_tot)
    print(f'R2 = {r_squared:.3f}')
    print(f'Fit = {popt[0]:.2f} / ( 1 + exp[-{popt[1]:.2f} * (x - {popt[2]:.2f})] ) + {popt[3]:.2f}')
    cct = Y[-1] + (Y[0] - Y[-1])*0.5 # % of original value (example: 95% inhibition = 5% of original)
    c_trans =  [np.interp( cct, Y[::-1],X[::-1]), cct ] # M, transition concentration

    print(f'Transition concentration = {1e3*c_trans[0]:.3f} mM')
    return X, Y, c_trans
regression_data = ctrans_regression(x_axis, y_axis2)


# --- Plotting ---
def size_vel_plot(x, y, y_error, X, Y, c_trans):
    ''' Plot size and velocity from processed fiber probe data. '''
    
    fig, ax1 = plt.subplots()
    ax1.plot(X, Y, color="#3399e6", lw=3, zorder=1)
    ax1.scatter(c_trans[0], c_trans[1], color="#8C0D07", lw=3, zorder=2)
    ax1.scatter( x, y, color="#69b3a2", lw=3, zorder=1)
    # ax1.errorbar(x, y, yerr=y_error, fmt='o', color="#69b3a2", capsize=3) # color="#3399e6"
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label2, fontsize=14)#, color="#69b3a2"
    ax1.tick_params(axis="y")#, labelcolor="#69b3a2")
    ax1.set_xscale('log')
    # ax1.set_ylim(1.7, 2.25)
    ax1.set_xlim(1e-4, 1)
    # fig.suptitle("Bubble properties", fontsize=20)
    plt.show()
    return
# size_vel_plot(x_axis, y_axis2, errorbar_axis2, 1, 1, 1)
size_vel_plot(x_axis, y_axis2, errorbar_axis2, regression_data[0], regression_data[1], regression_data[2])

def holdup_plot(x, y1, y2):
    ''' Plot gas holdup from processed fiber probe and pressure sensor data. '''
    
    fig, ax1 = plt.subplots()
    ax1.scatter(x, y1, color="#3399e6", lw=3, label='Pressure sensor')
    ax1.scatter(x, y2, color="#69b3a2", lw=3, label='Fiber probe')
    ax1.set_xlabel(x_label)
    ax1.set_ylabel(y_label3, fontsize=14)
    ax1.tick_params(axis="y")
    ax1.set_ylim(0, 0.1)
    # ax1.set_xlim(1e-4, 1)
    # fig.suptitle("Gas holdup", fontsize=20)
    plt.legend()
    plt.show()
    return 
holdup_plot(x_axis, y_axis4, y_axis3)


def export_measurements():
    return x_axis, y_axis1, errorbar_axis1, y_axis2, errorbar_axis2

def export_regression():
    return regression_data