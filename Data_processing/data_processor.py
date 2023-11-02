import os
import numpy as np
import pandas as pd

path = r'u:\Bubble Column\Data\A2 Fiber Probe\231101 - Flow variation in Water\2023-11-01T134048.evt'
df = pd.read_csv(path, sep='\t')


valid = df['Valid']
df = df.drop()
# for i, val in enumerate(valid):
#     print(val)
    # if val != "1":
    #     df = df.drop(val)