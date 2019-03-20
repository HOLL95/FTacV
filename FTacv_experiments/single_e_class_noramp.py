import isolver_noramp
import math
import numpy as np
import matplotlib.pyplot as plt
from params_class import params
import copy
import time
class single_electron:
    def __init__(self, dim_paramater_dictionary, optim_list, harmonic_range, filter_val):
        key_list=dim_paramater_dictionary.keys()
        self.nd_param=params(dim_paramater_dictionary)
        for i in range(0, len(key_list)):
            self.nd_param.non_dimensionalise(key_list[i], dim_paramater_dictionary[key_list[i]])
        nd_dict=self.nd_param.__dict__
        self.optim_list=optim_list
        self.harmonic_range=harmonic_range
        self.num_harmonics=len(harmonic_range)
        self.filter_val=filter_val
    def define_boundaries(self, boundaries):
        self.boundaries=boundaries
    def normalise(self, norm, boundaries):
        return  (norm-boundaries[0])/(boundaries[1]-boundaries[0])
    def un_normalise(self, norm, boundaries):
        return (norm*(boundaries[1]-boundaries[0]))+boundaries[0]
    def dimensionalise_data(self,data, times):
        data=np.multiply(data, self.nd_param.c_I0)
        times=np.multiply(times,self.nd_param.c_T0)
        return times, data
    def nondim_data(self, data, times):
        data=np.divide(data, self.nd_param.c_I0)
        times=np.divide(times,self.nd_param.c_T0)
        return dtimes, data
    def n_outputs(self):
        return 1
    def n_parameters(self):
        return len(self.optim_list)
    def pass_extra_data(self, time_series, fourier):
        self.secret_data_time_series=time_series
        self.secret_data_fourier=fourier
    def kaiser_filter(self, time_series, harmonical=False):
        frequencies=self.frequencies
        L=len(time_series)
        window=np.hanning(L)
        time_series=np.multiply(time_series, window)
        #f=np.fft.fftfreq(len(time_series), self.time_vec[1]-self.time_vec[0])
        Y=np.fft.fft(time_series)
        #plt.plot(f,Y)

        #Y_pow=np.power(copy.deepcopy(Y[0:len(frequencies)]),2)
        top_hat=copy.deepcopy(Y[0:len(frequencies)])
        inverse_top_hat=np.flip(copy.deepcopy(Y[-len(frequencies):-1]))
        results=np.zeros(len(frequencies), dtype=complex)
        inverse_results=np.zeros(len(frequencies), dtype=complex)
        #top_hat=np.zeros(len(frequencies))
        #inverse_top_hat=np.zeros(len(frequencies))
        xlocs=np.zeros(self.num_harmonics)
        labels=[""]*self.num_harmonics
        if harmonical==True:
            harmonics=np.zeros((self.num_harmonics, len(time_series)))
            #fig=plt.figure(num=None, figsize=(8,5), dpi=120, facecolor='w', edgecolor='k')
        for i in range(0, self.num_harmonics):
            true_harm=self.harmonic_range[i]*self.nd_param.omega*self.nd_param.c_T0
            #print true_harm
            #plt.axvline(true_harm, color="black", linestyle="--" )
            #xlocs[i]=true_harm
            #labels[i]=str(self.harmonic_range[i])+"f"
            #print true_harm, true_harm+(self.nd_param.omega*self.filter_val), true_harm-(self.nd_param.omega*self.filter_val)
            filter_bit=top_hat[np.where((frequencies<(true_harm+(true_harm*self.filter_val))) & (frequencies>true_harm-(true_harm*self.filter_val)))]
            #filter_bit=np.multiply(filter_bit, np.kaiser(len(filter_bit), 14))
            filter_bit_inverse=top_hat[np.where((frequencies<(true_harm+(true_harm*self.filter_val))) & (frequencies>true_harm-(true_harm*self.filter_val)))]
            #filter_bit_inverse=np.multiply(filter_bit, np.kaiser(len(filter_bit), 14))
            results[np.where((frequencies<(true_harm+(true_harm*self.filter_val))) & (frequencies>true_harm-(true_harm*self.filter_val)))]=filter_bit
            inverse_results[np.where((frequencies<(true_harm+(true_harm*self.filter_val))) & (frequencies>true_harm-(true_harm*self.filter_val)))]=filter_bit_inverse
            if harmonical==True:
                harmonics[i,0:len(filter_bit)]=filter_bit
        #if harmonical==True:
            #f=open('harmonics_nr.fig', 'w')
            #np.save(f, harmonics)
            #f.close()
        #print 'saved'
        #f=open('fourier_nr.fig', 'w')
        #np.save(f, [xlocs, labels, self.test_frequencies, np.abs(np.power(Y[:len(self.test_frequencies)],2))])
        #f.close()
        #print 'saved'
        #plt.plot(filter_bit)
        #plt.show()
        #if harmonical==True:
        #    return harmonics
        #plt.plot(frequencies, results)
        #plt.show()
        last_harm=(self.harmonic_range[-1]*true_harm)+(true_harm*self.filter_val)
        last_point=np.where(frequencies<last_harm)
        #print last_point
        results=results[:last_point[0][-1]]
        inverse_results=inverse_results[:last_point[0][-1]]
        real_hat=np.append(np.real(results), np.real(np.flip(inverse_results)))
        imag_hat=np.append(np.imag(results), np.imag(np.flip(inverse_results)))
        top_hat=np.append(real_hat, imag_hat)
        return abs(top_hat)
    def times(self, num_points):
        self.num_points=num_points
        #self.time_vec=np.arange(0, self.nd_param.time_end, self.nd_param.sampling_freq)
        self.time_vec=np.linspace(0, self.nd_param.time_end, num_points)
    def change_norm_group(self, param_list, method):
        normed_params=copy.deepcopy(param_list)
        if method=="un_norm":
            for i in range(0,len(param_list)):
                normed_params[i]=self.un_normalise(normed_params[i], [self.boundaries[0][i],self.boundaries[1][i]])
        elif method=="norm":
            for i in range(0,len(param_list)):
                normed_params[i]=self.normalise(normed_params[i], [self.boundaries[0][i],self.boundaries[1][i]])
        return normed_params
    def simulate(self,parameters, frequencies, flag='optimise', flag2='fourier', test="no"):

        if flag=='optimise' and self.label=="cmaes":
            normed_params=self.change_norm_group(parameters, "un_norm")
        else:
            normed_params=copy.deepcopy(parameters)
        for i in range(0, len(self.optim_list)):
                self.nd_param.non_dimensionalise(self.optim_list[i], normed_params[i])
        time_series=isolver_noramp.e_surface(self.nd_param.Cdl,self.nd_param.CdlE1,self.nd_param.CdlE2,self.nd_param.CdlE3,self.nd_param.nd_omega,self.nd_param.v, self.nd_param.alpha , \
                                    self.nd_param.E_start,  self.nd_param.E_reverse,  self.nd_param.d_E,  self.nd_param.Ru,200, self.time_vec,  self.nd_param.gamma, \
                                     self.nd_param.E_0, self.nd_param.k_0, self.nd_param.phase, math.pi, self.num_points)
        if flag2=='fourier':
            filtered=self.kaiser_filter(time_series)
            if test=="yes":
                plt.plot(self.secret_data_fourier)
                plt.plot(filtered, alpha=0.7)
                plt.show()
                #plt.plot(frequencies, filtered)
                #plt.plot(frequencies, self.secret_data)
                #plt.show()
                #original_harmonics=self.kaiser_filter(self.secret_data_time_series, frequencies, harmonical=True)
                #predicted_harmonics=self.kaiser_filter(time_series, frequencies, harmonical=True)
                #self.harmonics_plotter(original_harmonics, predicted_harmonics)

            return filtered
        elif flag2=='timeseries':
            if test=="yes":
                plt.plot(self.time_vec, time_series)
                plt.plot(self.time_vec, self.secret_data_time_series, alpha=0.7)
                plt.show()
            return time_series
