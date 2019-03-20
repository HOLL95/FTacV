import numpy as np
import os
import matplotlib.pyplot as plt
import isolver
import math
import numpy as np
import matplotlib.pyplot as plt
from single_e_class  import single_electron
import pints
import pints.plot
import copy
import time
params_for_opt=[]
param_list={
    'E_start': 260e-3, #(starting dc voltage - V)
    'E_reverse': 260e-3+(300e-3*2),    #  (reverse dc voltage - V)
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
param_list['E_0']=(param_list['E_start']+param_list['E_reverse'])/2
harmonic_range=np.arange(4,10,1)
noramp_fit=single_electron(param_list, params_for_opt, harmonic_range, 0.5)
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
desired_length=int(1e4)
current_results=results[0::64, 1]
time_results=results[0::64, 0]
time_results=time_results[:desired_length]/noramp_fit.nd_param.c_T0
current_results=current_results[:desired_length]/noramp_fit.nd_param.c_I0
signal_length=len(current_results)
noramp_fit.num_points=signal_length
noramp_fit.time_vec=time_results
frequencies=np.fft.fftfreq(signal_length, noramp_fit.time_vec[1]-noramp_fit.time_vec[0])
frequencies=frequencies[np.where(frequencies>0)]
noramp_fit.frequencies=frequencies
last_point= (harmonic_range[-1]*noramp_fit.nd_param.omega)+(noramp_fit.nd_param.omega*0.5)
plot_frequencies=frequencies[np.where(frequencies<last_point)]
noramp_fit.test_frequencies=plot_frequencies
likelihood_func=noramp_fit.kaiser_filter(current_results)
noramp_fit.pass_extra_data(current_results, likelihood_func)
test=noramp_fit.simulate([],frequencies, "no", "fourier", "yes" )
param_boundaries=[[param_list['E_start'],1, 0, 0, 1.0e-12], \
                    [param_list['E_reverse'], 1e4,100*param_list['omega'], 10, 9.5e-10]]# #
noramp_fit.define_boundaries(param_boundaries)
noramp_fit.label="cmaes"
noramp_fit.optim_list=['E_0', 'k_0', 'Ru', 'Cdl','gamma']
dummy_times=np.linspace(0, 1, len(likelihood_func))
#cmaes_problem=pints.SingleOutputProblem(noramp_fit, time_results, current_results)
cmaes_problem=pints.SingleOutputProblem(noramp_fit, dummy_times, likelihood_func)
score = pints.SumOfSquaresError(cmaes_problem)
CMAES_boundaries=pints.RectangularBoundaries([np.zeros(len(param_boundaries[0]))], [np.ones(len(param_boundaries[0]))])
x0=abs(np.random.rand(len(param_boundaries[0])))
for i in range(0, 1):
    found_parameters, found_value=pints.optimise(
                                                score,
                                                x0,
                                                boundaries=CMAES_boundaries,
                                                method=pints.CMAES
                                                )
    x0=found_parameters
cmaes_params=noramp_fit.change_norm_group(found_parameters, "un_norm")
cmaes_result=noramp_fit.simulate(found_parameters,frequencies, "optimise", "fourier", "yes" )
#noramp_fit.simulate(found_parameters,frequencies, "optimise", "fourier", "yes" )
noramp_fit.label="MCMC"
error=np.std(np.subtract(cmaes_result, likelihood_func))
mcmc_problem=pints.SingleOutputProblem(noramp_fit, dummy_times, likelihood_func)
updated_lb=np.append(cmaes_params*0.5, [0])#found_parameters[3]*0.97,
updated_ub=np.append(cmaes_params*1.5, [max(likelihood_func)])#found_parameters[3]*1.03,
updated_boundaries=[updated_lb, updated_ub]
updated_boundaries=np.sort(updated_boundaries, 0)
noramp_fit.define_boundaries(updated_boundaries)
print cmaes_params
print updated_boundaries
log_liklihood=pints.UnknownNoiseLogLikelihood(mcmc_problem)
log_prior=pints.UniformLogPrior(updated_boundaries[0],
                                updated_boundaries[1])
log_posterior=pints.LogPosterior(log_liklihood, log_prior)
mcmc_parameters=np.append(cmaes_params, error)
xs=[mcmc_parameters,
    mcmc_parameters,
    mcmc_parameters
    ]
mcmc = pints.MCMCSampling(log_posterior, 3, xs,method=pints.AdaptiveCovarianceMCMC)
mcmc.set_max_iterations(10000)
chains=mcmc.run()
pints.plot.trace(chains)
plt.show()
