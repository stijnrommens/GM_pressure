# Constants
rho = 997 # kg/m3
g = 9.81 # m/s2

# Measurment files
fiber_meas_dict = { 10:r'\2023-11-01T134048.evt',
                    20:r'\2023-11-01T134907.evt',
                    30:r'\2023-11-01T135501.evt',
                    40:r'\2023-11-01T142602.evt',
                    50:r'\2023-11-01T143330.evt',
                    60:r'\2023-11-01T144047.evt',
                    70:r'\2023-11-01T144655.evt',
                    80:r'\2023-11-01T145256.evt',
                    90:r'\2023-11-01T145926.evt',
                   100:r'\2023-11-01T150555.evt',
                    10:r'\2023-11-01T152703.evt',
                    40:r'\2023-11-01T151115.evt'} # flow:file

pressure_cal_dict = {      0:r'\2023-11-01_111637_pressure probe voltage.tdms',
                      106:r'\2023-11-01_112949_pressure probe voltage.tdms',
                      195:r'\2023-11-01_113547_pressure probe voltage.tdms',
                      302:r'\2023-11-01_113855_pressure probe voltage.tdms',
                    396.5:r'\2023-11-01_114332_pressure probe voltage.tdms',
                      484:r'\2023-11-01_114833_pressure probe voltage.tdms'} # height:file

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
                    40:r'\2023-11-01_151117_pressure probe voltage.tdms'} # flow:file

def constants():
    return rho, g

def fiber_meas():
    return fiber_meas_dict

def pressure_cal():
    return pressure_cal_dict

def pressure_meas():
    return pressure_meas_dict
