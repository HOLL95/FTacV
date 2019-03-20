import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import matplotlib.ticker as ticker
def chain_appender(chains, param):
    new_chain=chains[0, 5000:, param]
    for i in range(1, len(chains)):
        new_chain=np.append(new_chain, chains[i, 5000:, param])
    return new_chain
files=os.listdir('.')


#f=open(filename, "r")
titles=['$E^0$', '$k_0$', '$C_{dl}$      $x10^4$', '$Ru$', '$\sigma$']
pos=[0,1,3,2, 4]
true_params=[-0.4, 1e1, 0.000134*1e4,20.0]
noise_vals=[0.01, 0.02, 0.005, 0]
counter=-1
extension=".lowkf"
chain_select=1
j=0
fig=plt.figure(num=None, figsize=(12,9), dpi=120, facecolor='w', edgecolor='k')
filename_list=[]
number_list=[]
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
number_list=np.flip(number_list)
print number_list
filename_list=np.asarray(filename_list)
filename_list=filename_list[idx]
filename_list=np.flip(filename_list)
#filename_list=np.delete(filename_list, [1,3])
#number_list=np.delete(number_list, [1,3])
#filename_list=np.flip(filename_list)
for j in range(0, len(filename_list)):
        chains=np.load(filename_list[j])
        print filename_list[j]
        counter+=1
        for i in range(0, len(titles)):
            print pos[i]
            axes=plt.subplot(2,3,pos[i]+1)
            labels=number_list[j]
            if titles[i]=="$\sigma$":
                chains=chain_appender(chains, i)
                mean_noise=np.mean(chains)
                pc=labels
                x=mean_noise/pc
                chains=np.divide(chains, x)
                axes.hist(chains, alpha=0.4,bins=20, label=str(labels), edgecolor='black')
                axes.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0e'))
            else:
                plot_chain=chain_appender(chains, i)
                if titles[i]=='$Ru$':
                    plot_chain=plot_chain*1e4
                axes.hist(plot_chain, alpha=0.4,bins=20,label=str(labels), stacked=True, edgecolor='black')
                lb, ub = axes.get_xlim( )
                #axes.set_xticks([lb,true_params[pos[i]], ub])
            axes.set_xlabel(titles[pos[i]])
            axes.set_ylabel('frequency')
            if titles[i]!='$\sigma$':
                axes.axvline(true_params[pos[i]],color="black", linestyle="--")
            axes.legend()
            #axes.set_xlim((0.9*true_params[i], 1.1*true_params[i]))
axes=plt.subplot(2,3,i+2)
axes.axis('off')
plt.subplots_adjust(left=0.08, bottom=0.09, right=0.95, top=0.92, wspace=0.30)
plt.show()
