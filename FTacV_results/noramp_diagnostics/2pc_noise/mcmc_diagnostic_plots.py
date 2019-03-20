import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import scipy.stats as stats
from matplotlib import gridspec
import matplotlib.ticker as ticker
files=os.listdir('.')
fig=plt.figure(0)
lists=[]
extension=".1kk"
chain_select=0
file_list=[]
value_list=[]
def sort_by_idx(list_for_sort, indexes):
    sorted_list=[list_for_sort[x] for x in indexes]
    return sorted_list
for j in range(0, len(files)):
    if files[j][-len(extension):-1]== extension[:-1]:
        file_list.append(files[j])
        value_list.append(files[j][:files[j].index("_")])
print file_list
int_values=[float(x) for x in value_list]
sort_index=np.argsort(int_values)
sorted_files=sort_by_idx(file_list, sort_index)
sorted_values=sort_by_idx(value_list, sort_index)
sorted_uniques= np.sort(list(set(int_values)))
k=0
temp=""
param_select=1
num_plots=len(sorted_uniques)
mean_values=np.zeros(num_plots)
std_values=np.zeros(num_plots)
labels=[]
stds=[]
k_chain=np.array([])
test_stat=np.zeros(num_plots)
ax1=plt.subplot2grid((2,num_plots), (0,0), colspan=num_plots)
for i in range(0, num_plots):
    lists.append(plt.subplot2grid((2,num_plots), (1,i)))
#plt.subplot(
for i in range(0, num_plots):
    while (k<len(sorted_values)) and (str(int(sorted_uniques[i]))==sorted_values[k]):
        print sorted_uniques[i], sorted_values[k]
        chains=np.load(sorted_files[k])
        for j in range(0, 3):
            k_chain=np.append(k_chain, chains[j, 2000:, param_select])
        k+=1
        stds.append(np.std(chains[:, 2000:, param_select]))
    lists[i].hist(k_chain)
    lists[i].set_yticklabels([])
    lists[i].set_xticklabels([])
    lists[i].axvline(10000, color="black", linestyle="--")
    mean_values[i]=np.mean(k_chain)
    std_values[i]=np.mean(stds)
    k_chain=np.array([])
ax1.errorbar(range(len(mean_values)), mean_values, xerr= None, yerr=std_values)
ax1.axhline(10000, color="black", linestyle="--")
ax1.set_xlabel('Input frequency')
ax1.set_ylabel('Inferred k0 mean')
labels=[str(sorted_uniques[x])+" pi" for x in range(0, num_plots)]
ax1.set_xticks(range(0, num_plots))
ax1.set_xticklabels(labels)
plt.show()
