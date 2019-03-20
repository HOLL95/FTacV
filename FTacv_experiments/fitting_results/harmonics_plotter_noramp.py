import  matplotlib.pyplot as plt
import pints.plot
import numpy as np

class harmonics:
    def __init__(self, harmonics, input_frequency, filter_val):
        self.harmonics=harmonics
        self.num_harmonics=len(harmonics)
        self.input_frequency=input_frequency
        self.filter_val=filter_val
    def generate_harmonics(self, times, data):
        L=len(data)
        window=np.hanning(L)
        time_series=np.multiply(time_series, window)
        f=np.fft.fftfreq(len(time_series), times[1]-times[0])
        Y=np.fft.fft(time_series)
        #plt.plot(self.test_frequencies, Y[:len(self.test_frequencies)])
        #plt.plot(f, Y)
        #plt.show()
        #Y_pow=np.power(copy.deepcopy(Y[0:len(frequencies)]),2)
        last_harm=(self.harmonics[-1]*self.input_frequency)
        frequencies=f[np.where((f>0) & (f<(last_harm+(0.5*self.input_frequency))))]
        top_hat=copy.deepcopy(Y[0:len(frequencies)])
        plt.plot(top_hat)
        plt.show()
        harmonics=np.zeros((self.num_harmonics, len(time_series)))
        for i in range(0, self.num_harmonics):
            true_harm=self.harmonics[i]*self.input_frequency
            filter_bit=top_hat[np.where((frequencies<(true_harm+(self.nd_param.omega*self.filter_val))) & (frequencies>true_harm-(self.nd_param.omega*self.filter_val)))]
            harmonics[i,0:len(filter_bit)]=filter_bit
