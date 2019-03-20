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
titles=['$E^0$', '$k_0$', '$C_{dl}$      $x10^4$', '$Ru$','$\Gamma$' ]
true_params=[-0.4, 1e1, 0.000134*1e4,20.0]
noise_vals=[0.01, 0.02, 0.005, 0]
pos=[0,1,3,2,4]
counter=-1
extension=".lowkf"
chain_select=1
j=0
fig=plt.figure(num=None, figsize=(12,9), dpi=120, facecolor='w', edgecolor='k')
filename_list=[]
number_list=[]

chains=np.load('results_so_far.txt')
for i in range(0, len(titles)):
    print pos[i]
    axes=plt.subplot(2,3,pos[i]+1)
    if titles[i]=="$\sigma$":
        chains=chain_appender(chains, i)
        mean_noise=np.mean(chains)
        x=mean_noise/pc
        chains=np.divide(chains, x)
        axes.hist(chains, alpha=0.4,bins=20, edgecolor='black')
        axes.xaxis.set_major_formatter(ticker.FormatStrFormatter('%0.0e'))
    else:
        plot_chain=chain_appender(chains, i)
        if titles[pos[i]]=='$C_{dl}$      $x10^4$':
            plot_chain=plot_chain*1e4
        axes.hist(plot_chain, alpha=0.4,bins=20, stacked=True, edgecolor='black')
        lb, ub = axes.get_xlim( )
        axes.set_xticks(np.linspace(lb, ub, 4))
    axes.set_xlabel(titles[pos[i]])
    axes.set_ylabel('frequency')
    #if titles[i]!='$\sigma$':
    #    axes.axvline(true_params[pos[i]],color="black", linestyle="--")
    #axes.legend()
            #axes.set_xlim((0.9*true_params[i], 1.1*true_params[i]))
axes=plt.subplot(2,3,i+2)
axes.axis('off')
plt.subplots_adjust(left=0.08, bottom=0.09, right=0.95, top=0.92, wspace=0.30)
plt.show()
