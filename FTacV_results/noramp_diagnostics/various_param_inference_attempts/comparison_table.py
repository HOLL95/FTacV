import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import matplotlib.ticker as ticker
def latex_table_plotter(data1, data2):
    num_rows=len(data1)
    num_cols=len(data1[0])

    for i in range(0, num_rows):
        string1=""
        for j in range(0, num_cols):
            string1=string1+ str(round(data2[i][j],3))+ " & "
        print string1[:-2] + "\\\\"


noises=["0.005", "0.01", "0.02"]
files=[".cm", ".lowkf"]
titles=['$E^0$', '$k_0$', '$C_{dl}$', '$Ru$', '$\sigma$']
true_params=[-0.4, 1e1, 0.000134,20.0]
pos=[0,1,3,2,4]
param_cols=len(titles)
rows=len(noises)
cols=len(files)
summary_vals=np.zeros((rows, cols))
cm_param_means=np.zeros((rows, param_cols))
#nr_param_means=np.zeros((row, param_cols))
cm_param_sds=np.zeros((rows, param_cols))
#nr_param_sds=np.zeros((row, param_cols))
for i in range(0, rows):
    print noises[i]+files[0],noises[i]+files[1]
    cm_chains=np.load(noises[i]+files[0])
    nr_chains=np.load(noises[i]+files[1])
    for k in range(0, param_cols):
        cm_param=cm_chains[:, 5000:, k]
        nr_param=nr_chains[:, 5000:, pos[k]]
        #cm_err=(np.mean(cm_param)/true_params[k])
        #nr_err=np.mean(nr_param)/true_params[k]
        #cm_param_means[i,k]=nr_err/cm_err
        cm_param_sds[i,k]=np.std(nr_param)/np.std(cm_param)
#latex_table_plotter(cm_param_means)
latex_table_plotter(cm_param_means, cm_param_sds)
