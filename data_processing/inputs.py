# Measurment files [flow (L/min), fiber prober, pressure sensor]
measurement = [ [10, r'\2023-11-01T134048', r'\2023-11-01_134050_pressure probe voltage'],
                [20, r'\2023-11-01T134907', r'\2023-11-01_134909_pressure probe voltage'],
                [30, r'\2023-11-01T135501', r'\2023-11-01_135503_pressure probe voltage'],
                [40, r'\2023-11-01T142602', r'\2023-11-01_142604_pressure probe voltage'],
                [50, r'\2023-11-01T143330', r'\2023-11-01_143332_pressure probe voltage'],
                [60, r'\2023-11-01T144047', r'\2023-11-01_144049_pressure probe voltage'],
                [70, r'\2023-11-01T144655', r'\2023-11-01_144657_pressure probe voltage'],
                [80, r'\2023-11-01T145256', r'\2023-11-01_145257_pressure probe voltage'],
                [90, r'\2023-11-01T145926', r'\2023-11-01_145928_pressure probe voltage'],
                [100, r'\2023-11-01T150555', r'\2023-11-01_150557_pressure probe voltage']]
                # [10, r'\2023-11-01T152703', r'\2023-11-01_152705_pressure probe voltage'],
                # [40, r'\2023-11-01T151115', r'\2023-11-01_151117_pressure probe voltage'] ]

calibration = [ [0.000, r'\2023-11-01_111637_pressure probe voltage'],
                [0.106, r'\2023-11-01_112949_pressure probe voltage'],
                [0.195, r'\2023-11-01_113547_pressure probe voltage'],
                [0.302, r'\2023-11-01_113855_pressure probe voltage'],
                [0.3965, r'\2023-11-01_114332_pressure probe voltage'],
                [0.484, r'\2023-11-01_114833_pressure probe voltage'] ] # height (m):pressure sensor

def measurement_data():
    return measurement

def calibration_data():
    return calibration

def remove_column(matrix, col_to_remove):
    for row in matrix:
        del row[col_to_remove]
    return matrix
