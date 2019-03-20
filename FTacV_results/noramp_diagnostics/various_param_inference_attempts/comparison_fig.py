import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import matplotlib.ticker as ticker
files=os.listdir('.')
#f=open(filename, "r")
titles=['$E^0$', '$k_0$', '$C_{dl}$      $x10^4$','$Ru$']

pos=[0,1,3,2]
true_params=[-0.4, 1e1, 0.000134*1e4,20.0]
noise_vals=[0.01, 0.02, 0.005, 0]
fig=plt.figure(num=None, figsize=(8,8), dpi=120, facecolor='w', edgecolor='k')
counter=-1
file_1="0.005.cm"
file_2="0.005.lowkf"
experiment=["Ramp-free","Ramped"]
chain_select=0
files=[file_2, file_1]
def chain_appender(chains, param):
    new_chain=chains[0, 5000:, param]
    for i in range(1, len(chains)):
        new_chain=np.append(new_chain, chains[i, 5000:, param])
    return new_chain
for j in range(0, len(files)):
        chains=np.load(files[j])
        print files[j]
        counter+=1
        for i in range(0, len(titles)):
            axes=plt.subplot(2,2,pos[i]+1)
            labels=experiment[j]
            if j==0:
                plot_chain=chain_appender(chains, i)
                if titles[pos[i]]=='$C_{dl}$      $x10^4$':
                    axes.hist(1e4*plot_chain, alpha=0.4,bins=20,label=str(labels),  edgecolor='black')#edgecolor='black'
                else:
                    axes.hist(plot_chain, alpha=0.4,bins=20,label=str(labels),  edgecolor='black')

            else:
                plot_chain=chain_appender(chains, pos[i])
                if titles[pos[i]]=='$C_{dl}$      $x10^4$':
                    axes.hist(1e4*plot_chain, alpha=0.4,bins=20,label=str(labels),  edgecolor='black')
                else:
                    axes.hist(plot_chain, alpha=0.4,bins=20,label=str(labels),  edgecolor='black')
            axes.set_xlabel(titles[pos[i]])
            axes.set_ylabel('frequency')
            if titles[i]!='$\sigma$':
                axes.axvline(true_params[pos[i]],color="black", linestyle="--")
            axes.legend()
            #axes.set_xlim((0.9*true_params[i], 1.1*true_params[i]))
#axes=plt.subplot(2,3,i+2)
#axes.axis('off')
plt.subplots_adjust(left=0.08, bottom=0.09, right=0.95, top=0.92)
plt.show()
