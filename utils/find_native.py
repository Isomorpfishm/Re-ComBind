import numpy as np
import pandas as pd

col_names = ["col1", "col2", "col3", "col4", "col5", 
             "col6", "col7", "col8", "col9", "col10", 
             "col11", "col12", "col13", "col14", "col15", 
             "col16", "col17", "col18", "col19", "col20"]
             
path = r'/home/yuyang/temp02/vina/'
filename = path + 'Trial01_rmsd_results.csv'
filesave = path + 'Trial01_rmsd_bool.txt'
             
df = pd.read_csv(filename, header=None, sep=",", names=col_names)

data = []

for i in range(len(df)):
    for j in range(len(df.columns)):
        if df.iloc[i,j] > 15.0 :
            data.append(0)
        elif np.isnan(df.iloc[i,j]):
            data.append(1)
        elif df.iloc[i,j] <= 2.0:
            data.append(2)
        else:
            data.append(3)
            
df2 = pd.DataFrame(data)

np.savetxt(filesave, df2.values, fmt='%d')

arr = df2.values.copy()
arr.resize(228, 20)

df3 = pd.DataFrame(arr)

print(df3.head())

#df3.to_csv("Trial01_rmsd_bool.csv", header=False, index=False)
