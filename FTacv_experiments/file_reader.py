import numpy as np
import os
import matplotlib.pyplot as plt
dir_path = os.path.dirname(os.path.realpath(__file__))
data_path="/Experiment_data"
folder="/Black"
Method ="Improved"
path=dir_path+data_path+folder
files= os.listdir(path)
for data in files:
    if Method in data:
        x=np.fromfile((path+"/"+data))
plt.plot(np.diff(results[:,0]))
plt.show()
