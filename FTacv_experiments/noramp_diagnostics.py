
import isolver
import math
import numpy as np
import matplotlib.pyplot as plt
from single_e_class  import single_electron
import pints
import pints.plot
import copy
import time
freq_values=[4, 16, 32, 64, 128]
freq_values=np.sort(freq_values)
freq_range=np.multiply(freq_values , math.pi)
#points=[1e3, 4e3, 7e3, 1e4, 2e4]
noises=[0.05, 0.01, 0.02]
points=[5e5]
for lcv_1 in range(0, len(noises)):
    for lcv_2 in range(0, len(points)):
        print points[lcv_2]
        filename="Convergence_results_"
        points_range=int(points[lcv_2])
        params_for_opt=['E_0', 'k_0', 'Ru', 'Cdl']
        true_k=1e1
        true_omega=4*math.pi
        param_list={
            'E_start': -0.85, #(starting dc voltage - V)
            'E_reverse': -0.1,    #  (reverse dc voltage - V)
            'omega':true_omega,#8.88480830076,  #    (frequency Hz)
            'd_E': 150e-3,   #(ac voltage amplitude - V) freq_range[j],#
            'v': 24.97e-3,   #       (scan rate s^-1)
            'area': 0.03, #(electrode surface area cm^2)
            'Ru': 20.0,  #     (uncompensated resistance ohms)
            'Cdl': 0.000134, #(capacitance parameters)
            'CdlE1': 0,#0.000653657774506,
            'CdlE2': 0,#0.000245772700637,
            'CdlE3': 0,#1.10053945995e-06,
            'gamma': 6.5e-12,          # (surface coverage per unit area)
            'E_0': -0.4,      #       (reversible potential V)
            'E0_std':0.0312279186927,# (reversible potential dispersion)
            'k_0': true_k, #(reaction rate s-1)
            'k0_std': 0.0,
            'alpha': 0.5,
            'sampling_freq' : (1.0/200),
            'phase' : 3*(math.pi/2),
            'time_end':1000,
            'num_peaks': 50
        }
        harmonic_range=np.arange(1,13,1)
        noramp=single_electron(param_list, params_for_opt, harmonic_range, 0.5)
        noramp.label="cmaes"
        noramp.nd_param.time_end=(noramp.nd_param.num_peaks/noramp.nd_param.nd_omega)*2*math.pi
        noramp.nd_param.time_end=(2*(noramp.nd_param.E_reverse-noramp.nd_param.E_start)/noramp.nd_param.v)
        noramp.times(points_range)
        frequencies=np.fft.fftfreq(int(points_range), noramp.time_vec[1]-noramp.time_vec[0])
        frequencies=frequencies[np.where(frequencies>0)]
        noramp.frequencies=frequencies
        last_point= (harmonic_range[-1]*noramp.nd_param.omega)+(noramp.nd_param.omega*0.5)
        plot_frequncies=frequencies[np.where(frequencies<last_point)]
        noramp.test_frequencies=plot_frequncies
        true_params=[-0.4, true_k, 6.5e-12,param_list['Cdl'], param_list['CdlE1'],param_list['CdlE2'],param_list['CdlE3']]
        param_boundaries=[ [-0.6, 1, 1.5e-12,-1,0,-0.1,-0.1], \
                            [-0.2, 1e4, 1.5e-11,1,1,0.1,0.1]]# #
        noramp.optim_list=[]
        synthetic_data=noramp.simulate([], frequencies, "nah", "timeseries", "no")
        noise_val=noises[lcv_1]
        noise_max=max(synthetic_data)*noise_val
        noise=np.random.normal(0,noise_max, len(synthetic_data))
        noramp.define_boundaries(param_boundaries)
        noisy_data=np.add(synthetic_data, noise)
        #plt.plot(noisy_data)
        #plt.show()
        noisy_data_filtered=noramp.kaiser_filter(noisy_data)
        true_data_filtered=noramp.kaiser_filter(synthetic_data)
        error=np.subtract(noisy_data_filtered, true_data_filtered)
        noramp.pass_extra_data(noisy_data, noisy_data_filtered)
        noramp.optim_list=['E_0', 'k_0', 'Ru', 'Cdl']
        found_parameters=[param_list[x]for x in noramp.optim_list]
        mcmc_problem=pints.SingleOutputProblem(noramp, np.linspace(0, 1, len(noisy_data_filtered)), noisy_data_filtered)
        noramp.label="MCMC"
        updated_lb=np.append(np.multiply([param_list[x]for x in noramp.optim_list], 0.5), [0*noise_max])#found_parameters[3]*0.97,
        updated_ub=np.append(np.multiply([param_list[x]for x in noramp.optim_list], 1.5), [max(noisy_data_filtered)])#found_parameters[3]*1.03,
        updated_boundaries=[updated_lb, updated_ub]
        updated_boundaries=np.sort(updated_boundaries, 0)
        noramp.define_boundaries(updated_boundaries)
        print updated_boundaries
        log_liklihood=pints.UnknownNoiseLogLikelihood(mcmc_problem)
        log_prior=pints.UniformLogPrior(updated_boundaries[0],
                                        updated_boundaries[1])
        log_posterior=pints.LogPosterior(log_liklihood, log_prior)
        found_parameters=np.append(found_parameters, np.std(error)*1.5)
        print found_parameters
        xs=[found_parameters,
            found_parameters,
            found_parameters
            ]
        mcmc = pints.MCMCSampling(log_posterior, 3, xs,method=pints.AdaptiveCovarianceMCMC)
        mcmc.set_max_iterations(10000)
        chains=mcmc.run()
        #pints.plot.trace(chains)
        #plt.show()
        filename=str(noises[lcv_1])+".cm"
        f=open(filename, "w")
        np.save(f, chains)
        f.close()


"""
save=raw_input("Do you want to save this data? Y/N?")
if save =="Y":
    filename=raw_input("filename?")
    f=open(filename, "w")
    np.save(f, chains)
    f.close()
"""
"""
Cdl_range=np.array([0, 10, 100])#[4*pi, 16*pi,  64*pi, 128*pi, 256*math.pi]
#fig, ax=plt.subplots(noramp.num_harmonics,1)
for i in range(0,len(Cdl_range)):
    #plt.subplot(2,2,i+1)
    synthetic_data=noramp.simulate([Cdl_range[i]], frequencies, "nah", "timeseries")
    #plt.plot(synthetic_data)
    #plt.plot(caps, alpha=0.5)
    #noisy_data=np.add(synthetic_data,noise)
    harmonics=noramp.kaiser_filter(noisy_data, False)
    #print len(harmonics)

#plt.show()
#plt.plot(synthetic_data)
#plt.show()
plt.legend()
plt.show()
"""
