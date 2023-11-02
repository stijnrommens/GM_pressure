import os
import input_file
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

files_dict = input_file.fiber() # Import files from fiber probe

results = []
for flow in files_dict:
    path = r'u:\Bubble Column\Data\A2 Fiber Probe\231101 - Flow variation in Water' + files_dict[flow]
    df = pd.read_csv(path, sep='\t', decimal=',') # Initialize dataframe

    df_valid = df[df.Valid == 1] # Make df with only valid bubbles

    velocity = df_valid['Veloc']
    mean_velocity = sum(velocity)/len(velocity) # m/s

    size = df_valid['Size']
    mean_size = sum(size)/len(size)/1000 # mm
    
    results.append([flow, mean_velocity, mean_size])
results = np.array(results)

fig, ax1 = plt.subplots()
ax2 = ax1.twinx()
ax1.plot(results[:,0], results[:,1], color="#69b3a2", lw=3)
ax2.plot(results[:,0], results[:,2], color="#3399e6", lw=3)
ax1.set_xlabel("Gas flow [L/min]")
ax1.set_ylabel("Bubble mean velocity [m/s]", color="#69b3a2", fontsize=14)
ax1.tick_params(axis="y", labelcolor="#69b3a2")
ax2.set_ylabel("Bubble mean size [mm]", color="#3399e6", fontsize=14)
ax2.tick_params(axis="y", labelcolor="#3399e6")
fig.suptitle("Bubble properties", fontsize=20)
plt.show()
