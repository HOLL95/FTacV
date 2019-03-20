import numpy as np
import os
import matplotlib.pyplot as plt
import isolver_ramped
import math
import numpy as np
import matplotlib.pyplot as plt
from single_e_class_ramped  import single_electron
from harmonics_plotter import harmonics
import pints
import pints.plot
import copy
import time
params_for_opt=[]

dir_path = os.path.dirname(os.path.realpath(__file__))
data_path="/Experiment_data"
folder="/Black"
Method ="O_"
type="current"
type2="voltage"
path=dir_path+data_path+folder
files= os.listdir(path)
for data in files:
    if (Method in data)  and (type in data):
        results=np.loadtxt(path+"/"+data)
    elif (Method in data)  and (type2 in data):
        voltages=np.loadtxt(path+"/"+data)
dec_amount=1
current_results=results[0::dec_amount, 1]
time_results=results[0::dec_amount, 0]
de=300e-3
estart=260e-3-de
ereverse=estart+2*de
param_list={
    'E_start': -180e-3, #(starting dc voltage - V)
    'E_reverse': 620e-3,    #  (reverse dc voltage - V)
    'omega':8.9389,#8.88480830076,  #    (frequency Hz)
    'd_E': 150e-3,   #(ac voltage amplitude - V) freq_range[j],#
    'v': 29.8e-3,   #       (scan rate s^-1)
    'area': 0.03, #(electrode surface area cm^2)
    'Ru': 1.09154850e+01 ,  #     (uncompensated resistance ohms)
    'Cdl': 0.0000134, #(capacitance parameters)
    'CdlE1': 0,#0.000653657774506,
    'CdlE2': 0,#0.000245772700637,
    'CdlE3': 0,#1.10053945995e-06,
    'gamma': 1.73291079e-11,          # (surface coverage per unit area)
    'k_0': 3.33567800e+03, #(reaction rate s-1)
    'k0_std': 0.0,
    'alpha': 0.5,
    'sampling_freq' : (1.0/200),
    'phase' : 0,
    'time_end':1000,
    'num_peaks': 50
    }


param_list['E_0']=0.23471918314326964
harmonic_range=np.arange(4,8,1)
ramp_fit=single_electron(param_list, params_for_opt, harmonic_range, 0.1)
ramp_fit.label="cmaes"
ramp_fit.voltages=voltages[0::dec_amount, 1]/ramp_fit.nd_param.c_E0
time_results=time_results/ramp_fit.nd_param.c_T0
print time_results
current_results=current_results/ramp_fit.nd_param.c_I0
plt.plot(time_results, current_results, label="data")
ramp_fit.time_vec=time_results
signal_length=len(current_results)
ramp_fit.num_points=signal_length
frequencies=np.fft.fftfreq(signal_length, ramp_fit.time_vec[1]-ramp_fit.time_vec[0])
frequencies=frequencies[np.where(frequencies>0)]
ramp_fit.frequencies=frequencies
last_point= (harmonic_range[-1]*ramp_fit.nd_param.omega)+(ramp_fit.nd_param.omega*0.5)
plot_frequencies=frequencies[np.where(frequencies<last_point)]
ramp_fit.test_frequencies=plot_frequencies
harm_class=harmonics(harmonic_range, ramp_fit.nd_param.omega*ramp_fit.nd_param.c_T0, 0.05)
likelihood_func=ramp_fit.kaiser_filter(current_results)
data_harmonics=harm_class.generate_harmonics(time_results, current_results)
ramp_fit.pass_extra_data(current_results, likelihood_func)
param_boundaries=[[param_list['E_start'],1, 0, 0, 1.0e-12, 0.98*param_list['omega'], 0.4], \
                    [param_list['E_reverse'], 1e4,500, 0.1, 9e-11, 1.02*param_list['omega'], 0.6]]# #
ramp_fit.define_boundaries(param_boundaries)
ramp_fit.optim_list=['E_0', 'k_0', 'Ru', 'Cdl','gamma']
current_guess=np.array([0.73906793, 0.33350115, 0.10915485, 0.00086753, 0.85942673])
current_guess=ramp_fit.change_norm_group(current_guess, "un_norm")
current_guess[0]=param_list['E_0']
test=ramp_fit.simulate(current_guess,frequencies, "no", "timeseries", "no" )
test_harmonics=harm_class.generate_harmonics(time_results, test)
plt.plot(time_results, test, label="numerical")
plt.legend()
plt.show()
harm_class.plot_harmonics(time_results, test_harmonics, data_harmonics)

ramp_fit.optim_list=['E_0', 'k_0', 'Ru', 'Cdl','gamma', 'omega', 'alpha']
dummy_times=np.linspace(0, 1, len(likelihood_func))
#cmaes_problem=pints.SingleOutputProblem(ramp_fit, time_results, current_results)
cmaes_problem=pints.SingleOutputProblem(ramp_fit, dummy_times, likelihood_func)
plt.plot(test)
plt.show()
score = pints.SumOfSquaresError(cmaes_problem)
CMAES_boundaries=pints.RectangularBoundaries([np.zeros(len(param_boundaries[0]))], [np.ones(len(param_boundaries[0]))])
x0=[0.89709251, 0.00615058, 0.99999, 0.00111652, 0.71079056, 0.56844907, 0.09457175]


print
for i in range(0, 1):
    found_parameters, found_value=pints.optimise(
                                                score,
                                                x0,
                                                boundaries=CMAES_boundaries,
                                                method=pints.CMAES
                                                )
    x0=found_parameters
    print found_parameters
cmaes_results=ramp_fit.change_norm_group(found_parameters, "un_norm")
print cmaes_results
cmaes_time_prediction=ramp_fit.simulate(found_parameters,frequencies, "optimise", "timeseries", "yes" )
cmaes_prediction=ramp_fit.simulate(found_parameters,frequencies, "optimise", "fourier", "yes" )

test_harmonics=harm_class.generate_harmonics(time_results,cmaes_time_prediction)
harm_class.plot_harmonics(time_results, test_harmonics, data_harmonics)
#error=np.std(np.subtract(cmaes_prediction, current_results))
error=np.std(np.subtract(cmaes_prediction, likelihood_func))
mcmc_problem=pints.SingleOutputProblem(ramp_fit, dummy_times, likelihood_func)
#mcmc_problem=pints.SingleOutputProblem(ramp_fit, time_results, current_results)
updated_lb=np.append(cmaes_results*0.5, [0])#found_parameters[3]*0.97,
updated_ub=np.append(cmaes_results*1.5, [max(likelihood_func)])#found_parameters[3]*1.03,
updated_boundaries=[updated_lb, updated_ub]
updated_boundaries=np.sort(updated_boundaries, 0)
ramp_fit.define_boundaries(updated_boundaries)
log_liklihood=pints.UnknownNoiseLogLikelihood(mcmc_problem)
log_prior=pints.UniformLogPrior(updated_boundaries[0],
                                updated_boundaries[1])
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
mcmc_parameters=np.append(cmaes_results, error)
xs=[mcmc_parameters,
    mcmc_parameters,
    mcmc_parameters
    ]
ramp_fit.label="MCMC"
mcmc = pints.MCMCSampling(log_posterior, 3, xs,method=pints.AdaptiveCovarianceMCMC)
mcmc.set_max_iterations(10000)
chains=mcmc.run()
pints.plot.trace(chains)
plt.show()
