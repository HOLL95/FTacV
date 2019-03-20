import matplotlib.pyplot as plt
import numpy as np
import pints
import pints.plot
import os
import matplotlib.ticker as ticker
names=["time_series", "fourier","harmonics"]
extensions=["_r.fig", "_nr.fig"]
counter=-1
def time_series_plotter(times, data, ax):
    ax.plot(times, np.multiply(data, 1000))
    ax.set_xlabel('Time(s)')
    ax.set_ylabel('Current(mA)')
    return ax
def fourier_plotter(locs, labels, frequencies, data, ax):
    ax.semilogy(frequencies, data)
    ax.set_xticks(locs)
    ax.set_xticklabels(labels)
    ax.set_xlabel('Frequency(Hz)')
    ax.set_ylabel('Amplitude')
    return ax
def harmonics_plotter(times, harmonics, ax, harm_no):
    harmonics=np.multiply(harmonics, 1e9)
    print harm_no
    ax.plot(times, abs(np.fft.ifft(harmonics)))
    ax.yaxis.set_label_position("right")
    ax.set_ylabel(str(harm_no), rotation=0)
    ax.set_xlabel('Time(s)')
    if harm_no!=7:
        ax.set_xticks([])
        ax.set_xlabel("")
non_harm_list=np.zeros(14)
width=6
height=4
fig=plt.figure()
time_list=[]
letters=['A', 'B', 'C', 'D', 'E', 'F']
letter_counter=-1
for i  in range(0, 2):

    for j  in range(0, 2):
        letter_counter+=1
        ax=(plt.subplot2grid((14, 13), (i*height+(1*i), j*width+(1*j)), colspan=width, rowspan=height))
        ax.text(-0.1, 1.15, letters[letter_counter], transform=ax.transAxes,
            fontsize=16, fontweight='bold', va='top', ha='right')
        if i==0:
            filename=names[i]+extensions[j]
            file=np.load(filename)
            time_list.append(file[0])
            time_series_plotter(file[0], file[1], ax)
        else:
            filename=names[i]+extensions[j]
            file=np.load(filename)
            fourier_plotter(file[0], file[1], file[2], file[3], ax)
harmonic_range=range(4,9)
for j in range(0, 2):
    filename=names[2]+extensions[j]
    file=np.load(filename)
    counter=harmonic_range[0]
    for i in range(10, 14):
        ax=plt.subplot2grid((14, 13), (i, j*width+(1*j)), colspan=width, rowspan=1)
        if i == 10:
            letter_counter+=1
            ax.text(-0.1, 1.15, letters[letter_counter], transform=ax.transAxes,
                fontsize=16, fontweight='bold', va='top', ha='right')
        harmonics_plotter(time_list[j], file[counter, :], ax, counter)
        counter +=1
plt.subplots_adjust(left=0.08, bottom=0.09, right=0.95, top=0.92)
fig.text(0.02, 0.21, 'Current(nA)', ha='center', va='center', rotation='vertical')
fig.text(0.5, 0.21, 'Current(nA)', ha='center', va='center', rotation='vertical')
plt.show()
"""
fig, axes=plt.subplots(3,2)
for i in range(0, 2):
    idx=0
    filename=names[idx]+extensions[i]
    file=np.load(filename)
    time_series_plotter(file[0], file[1], axes[idx, i])
for i in range(0, 2):
    idx=1
    filename=names[idx]+extensions[i-2]
    file=np.load(filename)
    fourier_plotter(file[0], file[1], file[2], file[3], axes[idx, i])
axes[2,0]=fig
plt.show()
"""

"""

                plt.subplot(self.num_harmonics, 1, i+1)
                frame1 = plt.gca()
                if i<self.num_harmonics -1:
                    frame1.xaxis.set_major_formatter(plt.NullFormatter())
                frame1.yaxis.set_label_position("right")
                frame1.set_ylabel(str(self.harmonic_range[i]), rotation=0)
                plt.plot(abs(np.fft.ifft(harmonics[i,:])))
        if harmonical==True:
            frame1.set_xlabel('Time(s)')
            fig.text(0.04, 0.5, 'Current(A)', va='center', rotation='vertical')
            fig.text(0.94, 0.5, 'Harmonic', va='center', rotation=270)
            plt.show()
        #print labels
"""
