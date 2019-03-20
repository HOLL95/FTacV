import numpy as np
import os
import matplotlib.pyplot as plt
import isolver_noramp
import math
import numpy as np
import matplotlib.pyplot as plt
from single_e_class_noramp  import single_electron
from harmonics_plotter import harmonics
import pints
import pints.plot
import copy
import time
params_for_opt=[]

dir_path = os.path.dirname(os.path.realpath(__file__))
data_path="/Experiment_data"
folder="/Black"
Method ="N_"
type="current"
path=dir_path+data_path+folder
files= os.listdir(path)
for data in files:
    if (Method in data)  and (type in data):
        results=np.loadtxt(path+"/"+data)
length_list=[1e4, 2e4, 3e4]
dec_list=[8, 16, 32, 64]
repeat_num=5
desired_length=int(length_list[0])
dec_amount=dec_list[2]
current_results=results[0::dec_amount, 1]
time_results=results[0::dec_amount, 0]
de=300e-3
estart=260e-3-de
ereverse=estart+2*de
param_list={
    'E_start': estart, #(starting dc voltage - V)
    'E_reverse': ereverse,    #  (reverse dc voltage - V)
    'omega':8.94,#8.88480830076,  #    (frequency Hz)
    'd_E': 300e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 10.36e-3,   #       (scan rate s^-1)
    'area': 0.03, #(electrode surface area cm^2)
    'Ru': 200.0,  #     (uncompensated resistance ohms)
    'Cdl': 0.0000134, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 6.5e-12,          # (surface coverage per unit area)
    'k_0': 1000.0, #(reaction rate s-1)
    'k0_std': 0.0,
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    'phase' : 3*(math.pi/2),
    'time_end':1000,
    'num_peaks': 50
}
param_list['E_0']=(param_list['E_reverse']-param_list['E_start'])/2
harmonic_range=np.arange(4,9,1)
noramp_fit=single_electron(param_list, params_for_opt, harmonic_range, 0.5)
time_results=time_results[:desired_length]/noramp_fit.nd_param.c_T0
current_results=current_results[:desired_length]/noramp_fit.nd_param.c_I0
noramp_fit.time_vec=time_results
signal_length=len(current_results)
noramp_fit.num_points=signal_length
frequencies=np.fft.fftfreq(signal_length, noramp_fit.time_vec[1]-noramp_fit.time_vec[0])
frequencies=frequencies[np.where(frequencies>0)]
noramp_fit.frequencies=frequencies
last_point= (harmonic_range[-1]*noramp_fit.nd_param.omega)+(noramp_fit.nd_param.omega*0.5)
plot_frequencies=frequencies[np.where(frequencies<last_point)]
noramp_fit.test_frequencies=plot_frequencies
harm_class=harmonics(harmonic_range, noramp_fit.nd_param.omega*noramp_fit.nd_param.c_T0, 0.5)
data_harmonics=harm_class.generate_harmonics(time_results, current_results)
noramp_fit.optim_list=['E_0', 'k_0', 'Ru', 'Cdl','gamma']
noramp_fit.label="MCMC"
chains=np.load("results_so_far.txt")
means=[np.mean(chains[:, 5000:, x]) for x in np.arange(0, 5)]
print means
#fig, ax=plt.subplots(harm_class.num_harmonics,1)
for i in range(0, 1):#len(param_guesses)):
    time_series=noramp_fit.simulate(means, time_results, "no", "timeseries", "no")
    exp_harmonics=harm_class.generate_harmonics(time_results, time_series)
    #for j in range(0, len(exp_harmonics)):
#        harm_class.comparison_harmonics_plot(time_results, data_harmonics[j,:], exp_harmonics[j,:], ax[j])
    harm_class.plot_harmonics(time_results, data_harmonics, exp_harmonics)
    #plt.plot(time_results, time_series)
    #plt.plot(time_results, current_results)
    #plt.show()




"""
param_guesses=np.array([
[2.40666629e-01, 1.33659203e+04, 8.32138327e+00, 6.48986896e-0,  6.99137283e-10 ],
 [2.89424312e-01, 1.42829102e+04, 1.16728294e+00, 2.95794563e-02, 1.54470074e-10 ],
 [2.37426876e-01, 1.40219177e+04, 1.19240444e+01, 2.29657002e-01, 1.42238592e-09 ],
 [2.41984331e-01, 1.03875610e+00, 2.04067494e+01, 5.02952291e-03, 3.85177644e-10 ],
[2.34130799e-01, 4.02039142e+01, 2.46522095e+02, 3.91211938e-04, 2.59101097e-10 ],
 [2.77009135e-01, 1.01677918e+04, 1.22365956e+01, 1.17780590e+00, 1.21045676e-09 ],
 [2.85687934e-01, 1.39493817e+04, 5.61819373e-01, 4.73712354e-02, 1.46875670e-10 ],
 [2.85658802e-01, 1.31128704e+04, 5.60503182e-01, 4.73852493e-02, 1.46876400e-10 ],
 [2.74746341e-01, 1.57993170e+03, 1.46569159e+02, 3.42081253e-03, 6.72772614e-10],
 [2.84937019e-01, 1.24812926e+04, 1.51167967e+02, 1.58554929e-04, 1.47501166e-10]
])
#plt.show()
"""
