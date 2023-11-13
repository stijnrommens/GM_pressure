# Measurment files [flow (L/min), fiber prober, pressure sensor]
# measurement = [ [10, r'\2023-11-01T134048', r'\2023-11-01_134050_pressure probe voltage'],
#                 [20, r'\2023-11-01T134907', r'\2023-11-01_134909_pressure probe voltage'],
#                 [30, r'\2023-11-01T135501', r'\2023-11-01_135503_pressure probe voltage'],
#                 [40, r'\2023-11-01T142602', r'\2023-11-01_142604_pressure probe voltage'],
#                 [50, r'\2023-11-01T143330', r'\2023-11-01_143332_pressure probe voltage'],
#                 [60, r'\2023-11-01T144047', r'\2023-11-01_144049_pressure probe voltage'],
#                 [70, r'\2023-11-01T144655', r'\2023-11-01_144657_pressure probe voltage'],
#                 [80, r'\2023-11-01T145256', r'\2023-11-01_145257_pressure probe voltage'],
#                 [90, r'\2023-11-01T145926', r'\2023-11-01_145928_pressure probe voltage'],
#                 [100, r'\2023-11-01T150555', r'\2023-11-01_150557_pressure probe voltage']]
#                 # [10, r'\2023-11-01T152703', r'\2023-11-01_152705_pressure probe voltage'],
#                 # [40, r'\2023-11-01T151115', r'\2023-11-01_151117_pressure probe voltage'] ]

measurement = [ [10, r'\2023-11-10T144823', r'\2023-11-10_144824_pressure probe voltage'],
                [20, r'\2023-11-10T145828', r'\2023-11-10_145830_pressure probe voltage'],
                [30, r'\2023-11-10T135726', r'\2023-11-10_135727_pressure probe voltage'],
                [40, r'\2023-11-10T140329', r'\2023-11-10_140330_pressure probe voltage'],
                [50, r'\2023-11-10T151231', r'\2023-11-10_151232_pressure probe voltage'],
                [60, r'\2023-11-10T144114', r'\2023-11-10_144115_pressure probe voltage'],
                [70, r'\2023-11-10T150718', r'\2023-11-10_150719_pressure probe voltage'] ]

# calibration = [ [0.000, r'\2023-11-10_153842_pressure probe voltage'],
#                 [0.103, r'\2023-11-10_152042_pressure probe voltage'],
#                 [0.115, r'\2023-11-10_132616_pressure probe voltage'],
#                 [0.2105, r'\2023-11-10_152425_pressure probe voltage'],
#                 [0.289, r'\2023-11-10_152609_pressure probe voltage'],
#                 [0.374, r'\2023-11-10_152748_pressure probe voltage'] ] # height (m):pressure sensor

calibration = [ [0.000, r'\2023-11-01_111637_pressure probe voltage'],
                [0.106, r'\2023-11-01_112949_pressure probe voltage'],
                [0.195, r'\2023-11-01_113547_pressure probe voltage'],
                [0.302, r'\2023-11-01_113855_pressure probe voltage'],
                [0.3965, r'\2023-11-01_114332_pressure probe voltage'],
                [0.484, r'\2023-11-01_114833_pressure probe voltage'] ] # height (m):pressure sensor

fiber_probe_path = r'u:\Bubble Column\Data\A2 Fiber Probe\231110 - Flow variation in Water (2)'
pressure_sensor_path = r'u:\Bubble Column\Data\PXM419\231110 - Flow variation in Water (2)'
calibration_path = r'u:\Bubble Column\Data\PXM419\231101 - Calibration'


def measurement_data():
    return measurement, fiber_probe_path, pressure_sensor_path


def calibration_data():
    return calibration, calibration_path

# %%
# import pandas as pd
# import openpyxl

# fiber_probe_path = pd.read_excel('input_file.xlsx', sheet_name='Measurements', header=None, index_col=0, nrows=1)
# fiber_probe_path = fiber_probe_path.to_numpy()[0][0]

# pressure_sensor_path = pd.read_excel('input_file.xlsx', sheet_name='Measurements', header=0, index_col=0, nrows=1)
# pressure_sensor_path = pressure_sensor_path.to_numpy()[0][0]

# measurement = pd.read_excel('input_file.xlsx', sheet_name='Measurements', header=3)
# measurement = measurement.to_numpy()
# print(measurement)

# cal = pd.read_excel('input_file.xlsx', sheet_name='Calibration')
# print(cal)
# %%
