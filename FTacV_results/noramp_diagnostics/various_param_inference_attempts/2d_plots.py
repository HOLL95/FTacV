import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import scipy.stats as stats
from matplotlib.ticker import FormatStrFormatter
#chains_100=np.load("100_k.highk")
#chains_10=np.load("0.005.tsn")
chains_5000_0_005=np.load("0.005.hifk")
#chains_5000_2=np.load("5000_k_2")
titles=['E_0', 'k_0', 'Cdl', 'Ru']
burn_in_val=5000
true_params=[-0.4, 5000, 0.000133985,20]
n_param=len(titles)
chains=chains_5000_0_005
chain_select=2
fig_size=(12,12)
fig, ax=plt.subplots(n_param, n_param)
def plot_kde_1d(x, ax):
    """ Creates a 1d histogram and an estimate of the PDF using KDE. """
    xmin = np.min(x)
    xmax = np.max(x)
    x1 = np.linspace(xmin, xmax, 100)
    x2 = np.linspace(xmin, xmax, 50)
    ax.hist(x, bins=x2)
def plot_kde_2d(x, y, ax):
    # Get minimum and maximum values
    xmin, xmax = np.min(x), np.max(x)
    ymin, ymax = np.min(y), np.max(y)

    # Plot values
    values = np.vstack([x, y])
    ax.set_xlim(xmin, xmax)
    ax.set_ylim(ymin, ymax)
    ax.imshow(np.rot90(values), cmap=plt.cm.Blues, extent=[xmin, xmax, ymin, ymax])

    # Create grid
    xx, yy = np.mgrid[xmin:xmax:100j, ymin:ymax:100j]
    positions = np.vstack([xx.ravel(), yy.ravel()])

    # Get kernel density estimate and plot contours
    kernel = stats.gaussian_kde(values)
    f = np.reshape(kernel(positions).T, xx.shape)
    ax.contourf(xx, yy, f, cmap='Blues')
    ax.contour(xx, yy, f, colors='k')

    # Fix aspect ratio
    print ((xmax - xmin)/ (ymax - ymin))
    ax.set_aspect((xmax - xmin)/ (ymax - ymin))
for i in range(0,n_param):
    for j in range(0, n_param):
        if i==j:
            axes=ax[i,j]
            if titles[i]=='Cdl':
                axes.xaxis.set_major_formatter(FormatStrFormatter('%.1e'))
            plot_kde_1d(chains[chain_select, burn_in_val:, j], ax=axes)
            axes.axvline(true_params[i],color="black", linestyle="--")
            ax[i, j].set_ylabel("frequency")
        elif i<j:
            ax[i,j].axis('off')
        else:
            axes=ax[i,j]
            if titles[j]=='Cdl':
                axes.xaxis.set_major_formatter(FormatStrFormatter('%.2e'))
            if titles[i]=='Cdl':
                axes.yaxis.set_major_formatter(FormatStrFormatter('%.2e'))
            plot_kde_2d(chains[chain_select, burn_in_val:, j], chains[chain_select, burn_in_val:, i], ax=axes)
            axes.axvline(true_params[j],color="black", linestyle="--")
            axes.axhline(true_params[i],color="black", linestyle="--")
        if i!=0:
            ax[i, 0].set_ylabel(titles[i])
        if j!=n_param:
            ax[-1, i].set_xlabel(titles[i])

plt.show()
