import numpy as np

# Constants
rho = 997 # kg/m3
g = 9.81 # m/s2
radius = 192e-3/2 # m
Ac = np.pi*radius**2 # m2
liquid_height = 224e-3 # m 
sensor_height = 563e-3+16e-3 # m

# Measurment files
fiber_meas_dict = { 10:r'\2023-11-01T134048',
                    20:r'\2023-11-01T134907',
                    30:r'\2023-11-01T135501',
                    40:r'\2023-11-01T142602',
                    50:r'\2023-11-01T143330',
                    60:r'\2023-11-01T144047',
                    70:r'\2023-11-01T144655',
                    80:r'\2023-11-01T145256',
                    90:r'\2023-11-01T145926',
                   100:r'\2023-11-01T150555',
                    10:r'\2023-11-01T152703',
                    40:r'\2023-11-01T151115'} # flow (L/min):file

pressure_cal_dict = { 0.0   :r'\2023-11-01_111637_pressure probe voltage.tdms',
                      0.106 :r'\2023-11-01_112949_pressure probe voltage.tdms',
                      0.195 :r'\2023-11-01_113547_pressure probe voltage.tdms',
                      0.302 :r'\2023-11-01_113855_pressure probe voltage.tdms',
                      0.3965:r'\2023-11-01_114332_pressure probe voltage.tdms',
                      0.484 :r'\2023-11-01_114833_pressure probe voltage.tdms'} # height (m):file

pressure_meas_dict = {10:r'\2023-11-01_134050_pressure probe voltage.tdms',
                    20:r'\2023-11-01_134909_pressure probe voltage.tdms',
                    30:r'\2023-11-01_135503_pressure probe voltage.tdms',
                    40:r'\2023-11-01_142604_pressure probe voltage.tdms',
                    50:r'\2023-11-01_143332_pressure probe voltage.tdms',
                    60:r'\2023-11-01_144049_pressure probe voltage.tdms',
                    70:r'\2023-11-01_144657_pressure probe voltage.tdms',
                    80:r'\2023-11-01_145257_pressure probe voltage.tdms',
                    90:r'\2023-11-01_145928_pressure probe voltage.tdms',
                   100:r'\2023-11-01_150557_pressure probe voltage.tdms',
                    10:r'\2023-11-01_152705_pressure probe voltage.tdms',
                    40:r'\2023-11-01_151117_pressure probe voltage.tdms'} # flow (L/min):file

def constants():
    return rho, g, radius, liquid_height, sensor_height, Ac

def fiber_meas():
    return fiber_meas_dict

def pressure_cal():
    return pressure_cal_dict

def pressure_meas():
    return pressure_meas_dict
