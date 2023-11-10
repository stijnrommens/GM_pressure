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

calibration = [ [0.000, r'\2023-11-10_153842_pressure probe voltage'],
                [0.103, r'\2023-11-10_152042_pressure probe voltage'],
                [0.115, r'\2023-11-10_132616_pressure probe voltage'],
                [0.2105, r'\2023-11-10_152425_pressure probe voltage'],
                [0.289, r'\2023-11-10_152609_pressure probe voltage'],
                [0.374, r'\2023-11-10_152748_pressure probe voltage'] ] # height (m):pressure sensor


def measurement_data():
    return measurement


def calibration_data():
    return calibration
