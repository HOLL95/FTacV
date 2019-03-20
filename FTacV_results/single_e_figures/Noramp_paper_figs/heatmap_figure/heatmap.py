import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import scipy.stats as stats
from matplotlib.ticker import FormatStrFormatter
import seaborn as sns
files=os.listdir('.')
freq_list=[]
point_list=[]
print [np.int(x) for x in freq_list]
num_files=len(files)
extension=".0.hm"
for i in range(0, num_files):
    if files[i][files[i].index('.'):]==extension:
        chains=np.load(files[i])
        freq_number=int(files[i][:files[i].index("_")])
        point_number=int(files[i][files[i].index("_")+1:files[i].index(".")])
        if (freq_number in freq_list) == False:
            freq_list.append(freq_number)
        if (point_number in point_list) ==False:
            point_list.append(point_number)
freq_list=np.sort(freq_list)
freq_list=np.flip(freq_list)
point_list=np.sort(point_list)
num_freqs=len(freq_list)
num_points=len(point_list)
mean_array=np.zeros((num_freqs, num_points))
#mean_array=[mean_array]*num_points
#fig, ax=plt.subplots(num_freqs, num_points)
std_array=np.zeros((num_freqs, num_points))
for i in range(0, num_files):
    if files[i][files[i].index('.'):]==extension:
        chains=np.load(files[i])
        thinned_chain=chains[:, 5000:, 1]
        freq_number=int(files[i][:files[i].index("_")])
        point_number=int(files[i][files[i].index("_")+1:files[i].index(".")])
        freq_idx=np.where(freq_list==freq_number)[0][0]
        point_idx=np.where(point_list==point_number)[0][0]
        print freq_idx, point_idx,files[i]
        mean_array[freq_idx, point_idx]=abs(np.mean(thinned_chain)-0)
        std_array[freq_idx, point_idx]=np.std(thinned_chain)
        #ax[freq_idx, point_idx].hist(thinned_chain)
#plt.show()
freq_list=[str(x)+'$\pi$' for x in freq_list]
#plt.subplot(1,2,1)
#ax=sns.heatmap(mean_array, cmap='viridis', xticklabels=point_list, yticklabels=freq_list)
ax=sns.heatmap(std_array, cmap='plasma_r',xticklabels=point_list, yticklabels=freq_list,  annot=mean_array,  fmt='.3g',  cbar_kws={'label': 'standard deviation'})
plt.xlabel('Number of points')
plt.ylabel('Input frequency')
plt.show()
