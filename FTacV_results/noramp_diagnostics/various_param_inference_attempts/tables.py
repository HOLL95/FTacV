import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import matplotlib.ticker as ticker
from decimal import Decimal
files=os.listdir('.')
#f=open(filename, "r")
titles=['$E^0$', '$k_0$', '$C_{dl}$', '$Ru$', '$\\sigma$']
pos=[0,1,2,3, 4]
true_params=[-0.4, 1e4, 0.000134,20.0]
noise_vals=[0.01, 0.02, 0.005, 0]
#fig, ax=plt.subplots(2,3)
counter=-1
extension=".cm"
chain_select=1
j=0
filename_list=[]
number_list=[]
def latex_table_plotter(data):
    num_rows=len(data)
    for i in range(0, num_rows):
        string=""
        for j in range(0, len(data[i,:])):
            string=string+str(data[i,j] + ' & ')
        string = string[:-2] + "\\\\"
        print string


for j in range(0, len(files)):
    if files[j][-len(extension):-1]== extension[:-1]:
        filename_list.append(files[j])
        labels=files[j][:files[j].index(extension[1])]
        number_list.append(labels[:-1])
number_list=[float(x) for x in number_list]
number_list=np.asarray(number_list)
idx=np.argsort(number_list)
print number_list
number_list=number_list[idx]
#number_list=np.flip(number_list)
print number_list
filename_list=np.asarray(filename_list)
filename_list=filename_list[idx]
#filename_list=np.delete(filename_list, [1,3])
#number_list=np.delete(number_list, [1,3])
#filename_list=np.flip(filename_list)
table_data=np.chararray((len(number_list)+1, len(titles)+1), itemsize=30)
table_data[0,0]="Noise\\%"
table_data[0, 1:]=titles
table_data[1:, 0]=number_list*100
for j in range(1, len(filename_list)+1):
        print filename_list[j-1]
        print number_list[j-1]
        chains=np.load(filename_list[j-1])
        plt.hist(chains[0, 5000:, 1], alpha=0.5, label=str(filename_list[j-1]))
        for i in range(1, len(titles)+1):
            mean=np.mean(chains[:, 5000:, pos[i-1]])
            std=np.std(chains[:, 5000:, pos[i-1]])
            if titles[i-1]=='$\\sigma$' :
                pc=number_list[j-1]
                x=np.median((chains[:, 5000:, pos[i-1]]))/pc
                mean=mean/x
                std=std/x
            if (abs(mean) <0.005)==1 :
                table_data[j, i] = "{:.2E}".format(Decimal(str(mean)))
            else:
                table_data[j, i] = str(round(mean,3))
            if  (std <0.01)==1:
                table_data[j, i]=table_data[j, i] + ' ('+"{:.2E}".format(Decimal(str(std))) +')'
            else:
                table_data[j, i]=table_data[j, i] + ' ('+str(round(std,3)) +')'
plt.legend()
latex_table_plotter(table_data)
plt.show()
