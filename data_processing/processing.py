import os
import openpyxl
import numpy as np
import pandas as pd
from nptdms import TdmsFile

import calibration


# --- Experiment ---
excel_file = r'\\tudelft.net\student-homes\R\srommens\My Documents\GitHub\GM_pressure\data_processing\input_file.xlsx'
fiber_probe_path = pd.read_excel(excel_file, sheet_name='Transition_NaCl', header=None, index_col=0, nrows=1)
fiber_probe_path = fiber_probe_path.to_numpy()[0][0]
pressure_sensor_path = pd.read_excel(excel_file, sheet_name='Transition_NaCl', header=0, index_col=0, nrows=1)
pressure_sensor_path = pressure_sensor_path.to_numpy()[0][0]
files = pd.read_excel(excel_file, sheet_name='Transition_NaCl', header=3)
files = files.to_numpy()

# --- Constants ---
rho = 997 # kg/m3, density of water
g = 9.81 # m/s2, gravitational constant
radius = 150e-3/2 # m, inside radius column
sensor_height = 748e-3 # m, height pressure sensor from sparger
# mw_additive = 58.440 # g/mol, molar weight of NaCl 
mw_additive = 142.04 # g/mol, molar weight of Na2SO4 
# mw_additive = 53.491 # g/mol, molar weight of NH4Cl 
# mw_additive = 132.14 # g/mol, molar weight of (NH4)2SO4

Ac = np.pi*radius**2 # m2, cross-sectional area column


def fiber_probe(files, folder):
    ''' From the measurements, extract fiber probe information. '''
    
    results = []
    for i, file in enumerate(files):
        param = file[0]     # height, weight, etc.
        file_name = file[1] # -
        
        liquid_height   = file[3] # m, water height from sensor with no flow
        depleted_height = file[4] # m, correction for water loss in bottom port (0.104 mm/min)
        density_height  = file[5] # m, correction for density increase by salt
        volume_L = Ac*(sensor_height + liquid_height - depleted_height) # m3
        
        # Load the measurement files in dataframes + get void fraction
        path = folder + file_name + '.evt'
        df = pd.read_csv(path, sep='\t', decimal=',')
        path_stream = folder + file_name + '_stream.evt'
        df_stream = pd.read_csv(path_stream, sep='\t', decimal=',')
        arrival  = df_stream['Arrival']
        duration = df_stream['Duration']
        void_fraction = sum(duration)/arrival.iloc[-1] # -
        
        # Check if variable has multiple measurement files + add to dataframes
        prev_param_count = np.count_nonzero(files[:i+1,0] == param) # Check if same parameter has not been passed previously
        if prev_param_count > 1:
            continue
        tot_param_count = np.count_nonzero(files[:,0] == param) # Get all files with same parameter
        if tot_param_count > 1:
            for k, extra_file in enumerate(files[i+1:i+tot_param_count]):
                extra_param     = extra_file[0]
                extra_file_name = extra_file[1]
                extra_path = folder + extra_file_name + '.evt'
                extra_df   = pd.read_csv(extra_path, sep='\t', decimal=',')
                df = pd.concat([df, extra_df])
                extra_stream_path = folder + extra_file_name + '_stream.evt'
                extra_stream_df   = pd.read_csv(extra_stream_path, sep='\t', decimal=',')
                arrival  = extra_stream_df['Arrival']
                duration = extra_stream_df['Duration']
                void_fraction += sum(duration)/arrival.iloc[-1] # -
                
        df_valid = df[df.Valid == 1] # Only valid bubbles

        # Obtain velocity and size
        velocity = df_valid['Veloc'].sort_values()
        lower_velocity, median_velocity, upper_velocity = np.percentile(velocity, [25, 50, 75]) # m/s
        lower_velocity = abs(lower_velocity - median_velocity) # m/s
        upper_velocity = abs(upper_velocity - median_velocity) # m/s

        size = 1e-6*df_valid['Size'].sort_values()
        # d32 = sum(size**3)/sum(size**2)
        lower_size, median_size, upper_size = np.percentile(size, [25, 50, 75]) # m
        lower_size = abs(lower_size - median_size) # m
        upper_size = abs(upper_size - median_size) # m

        # Obtain average gas holdup
        void_fraction /= tot_param_count # -
        
        results.append([param, median_velocity, lower_velocity, upper_velocity, median_size, lower_size, upper_size, volume_L, void_fraction])
    return np.array(results)
fiber_probe_results = fiber_probe(files, fiber_probe_path)


def pressure_sensor(files, folder, fit):
    ''' From the measured mean voltage and the fitted line, obtain the gas holdup. '''
    
    results = []
    for i, file in enumerate(files):
        param = file[0]     # height, flow, etc.
        file_name = file[2] # -
        
        liquid_height   = file[3] # m, water height from sensor with no flow
        depleted_height = file[4] # m, correction for water loss in bottom port (0.104 mm/min)
        density_height  = file[5] # m, correction for density increase by salt
        
        # Load file
        path = folder + file_name + '.tdms'
        loaded_file = TdmsFile(path)
        
        # Obtain holdup from voltage
        for group in loaded_file.groups():
            df = group.as_dataframe()
        voltage = df['Dev1/ai1']
        mean_voltage = np.mean(voltage) # V
        
        gas_height = fit[0]*mean_voltage + fit[1] # m
        volume_L = Ac * sensor_height + depleted_height + density_height # m3
        volume_G = Ac * (gas_height - liquid_height) # m3
        holdup   = volume_G/volume_L # -
        
        # Check if variable has multiple measurement files + add to holdup
        prev_param_count = np.count_nonzero(files[:i+1,0] == param) # Check if same parameter has not been passed previously
        if prev_param_count > 1:
            continue
        tot_param_count = np.count_nonzero(files[:,0] == param) # Get all files with same parameter
        if tot_param_count > 1:
            for k, extra_file in enumerate(files[i+1:i+tot_param_count]):
                extra_param     = extra_file[0]
                extra_file_name = extra_file[2]
                liquid_height   = extra_file[3] # m, water height from sensor with no flow
                depleted_height = extra_file[4] # m, correction for water loss in bottom port (0.104 mm/min)
                density_height  = extra_file[5] # m, correction for density increase by salt
                extra_path = folder + extra_file_name + '.tdms'
                extra_loaded_file = TdmsFile(extra_path)
                for group in extra_loaded_file.groups():
                    df = group.as_dataframe()
                voltage = df['Dev1/ai1']
                mean_voltage = np.mean(voltage) # V
                gas_height   = fit[0]*mean_voltage + fit[1]
                volume_L = Ac * sensor_height + depleted_height + density_height # m3
                volume_G = Ac * (gas_height - liquid_height) # m3
                holdup  += volume_G/volume_L # -
        
        # Obtain average gas holdup
        holdup /= tot_param_count # -
        
        results.append([param, holdup])
    return np.array(results)
fit = calibration.export_fit()
pressure_sensor_results = pressure_sensor(files, pressure_sensor_path, fit)


def mass2strength(param, volume_L, valence):
    ''' Run this function to convert mass (an expimental variable) -> ionic strength. Can also be used to obtain concentration (use valence=1). '''
    conc = param / mw_additive / (1000*volume_L) # M, convert added weight to concentration
    return conc * valence**2 # M, concentration to ionic strength

def export_fiber_probe():
    return fiber_probe_results

def export_pressure_sensor():
    return pressure_sensor_results
