import numpy as np
from scipy.optimize import curve_fit
from matplotlib import cm
import matplotlib.pyplot as plt
import matplotlib.animation as animation


# --- Experiment ---
x_data = [-70, -70, -70, -70, -60, -60, -60, -60, -50, -50, -40, -40, -30, -30, -20, -20, -10, -10,   
          0,   0,   0,   0,   0,   0,  10,  10,  20,  20, 30,  30,  40,  40,  50,  50,  50,  60,  60,  70, 70]
fiber_probe_data = [0.05893618, 0.04280574, 0.05917667, 0.0563055,  0.06051006, 0.06515897, 0.05066538, 0.05495068, 0.0586895,  0.06226933, 0.05995925, 0.06099276,
                    0.06328487, 0.06943013, 0.07171036, 0.0617046,  0.06777935, 0.06196412, 0.07065465, 0.07183665, 0.06626814, 0.07660479, 0.06426343, 0.06482662,
                    0.05882134, 0.06699567, 0.06094391, 0.0626539,  0.05863908, 0.07203397, 0.07286614, 0.07825246, 0.07175864, 0.06884522, 0.06762835, 0.07201519, 
                    0.0610767,  0.07974235, 0.08293386]

x_data2 = [-70, -70, -70, -70, -60, -60, -50, -50, -40, -40, -40, -30, -30, -30, -20, -20, -20, -10, -10,   
           0,   0,   0,   0,  10,  10,  20,  20,  30, 30,  30,  40,  40,  40,  50,  50,  50,  60,  60,  70,  70,  70]
fiber_probe_data2 = [0.05365972, 0.04664872, 0.05018213, 0.05420077, 0.04424736, 0.04858333, 0.04818651, 0.04613074, 0.05129754, 0.04516808, 0.0513079,  0.04698982,
                     0.047076,   0.04946231, 0.05252364, 0.04797888, 0.05436357, 0.04795493, 0.05073111, 0.04759599, 0.05260637, 0.05087095, 0.04961347, 0.04809385,
                     0.04278638, 0.05018265, 0.05551859, 0.04440523, 0.05574967, 0.04531603, 0.04536029, 0.04972584, 0.05009904, 0.04056823, 0.04388955, 0.04701655,
                     0.04760938, 0.04237517, 0.04852955, 0.05289236, 0.05150069]

x_x = x_data
y_x = np.zeros( len(x_x) )
z_x = fiber_probe_data

y_y = x_data2
x_y = np.zeros( len(y_y) )
z_y = fiber_probe_data2

x = np.hstack((x_x, x_y))
y = np.hstack((y_x, y_y))
z = np.hstack((z_x, z_y))


# --- Regression ---
def holdup_regression(x, y, z):
    ''' Make a regression of the radial gas holdup measurements from the fiber probe. '''
    
    def func(xy, a, b, c, d, e, f): 
        x, y = xy 
        f1 = a*x**2 + b*x + c
        f2 = d*y**2 + e*y + f
        return f1 + f2
    popt, pcov = curve_fit(func, (x, y), z)

    X, Y = np.meshgrid(np.arange(-75, 76, 1), np.arange(-75, 76, 1))
    Z = func((X, Y), *popt)
        
    radius = np.sqrt(X**2 + Y**2)
    X_masked = np.ma.masked_where(radius > 75, X)
    Y_masked = np.ma.masked_where(radius > 75, Y)
    Z_masked = np.ma.masked_where(radius > 75, Z)
    
    z_centre = func((0,0), *popt)
    residuals = z - func((x, y), *popt)
    ss_res = np.sum(residuals**2)
    ss_tot = np.sum( (z-np.mean(z))**2 )
    r_squared = 1 - (ss_res/ss_tot)
    return z_centre, r_squared, (X_masked, Y_masked, Z_masked)
regression_data = holdup_regression(x, y, z)


# --- Plotting ---
def regression_plot(measurement_data, regression_data):
    ''' Plot gas holdup regression from processed fiber probe. '''

    x, y, z = measurement_data
    X_masked, Y_masked, Z_masked = regression_data

    fig = plt.figure() 
    ax = fig.add_subplot(111, projection='3d') 
    ax.scatter(x, y, z, color='blue') 
    ax.set_xlabel('X') 
    ax.set_ylabel('Y') 
    ax.set_zlabel('Z') 
    surf = ax.plot_surface(X_masked, Y_masked, Z_masked, cmap=cm.coolwarm) 
    fig.colorbar(surf, shrink=0.5, aspect=10)

    # Rotate the axes and update
    for angle in range(0, 360+1):
        ax.view_init(30, angle)

        plt.draw()
        plt.pause(0.01)

    # Save animation
    # def update(angle):
    #     ax.view_init(30, angle)
    #     return fig, 
    # ani = animation.FuncAnimation(fig, update, interval=10, blit=True, repeat=True, frames=360)    
    # ani.save("TLI.gif", dpi=300, writer=animation.PillowWriter(fps=25))
    return
regression_plot((x,y,z), regression_data[-1])


def export_measurements():
    return x, y, z

def export_regression():
    return regression_data
