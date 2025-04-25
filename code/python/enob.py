
from file_reader import file_reader
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt
import math
import os



class EnobPlot:
    def __init__(self, filename = "", freq = 0, ):
        self.filename = filename
        self.lsb = 1
        self.n_bits = 12
        self.v_fs = 0.6
        self.v_max = 0.6
        self.v_min = -0.6
        self.lsb = (self.v_max - self.v_min) / 2 ** 12

        

    
    def show_data(self, x_data,y_data, label):
        plt.scatter(x_data,y_data, label = label)
        
        

    def fit_enob(self, y_fit, y_data, verbose = False):
        sigma_q = ((self.lsb * self.lsb)/12) ** (0.5)
        nad = self.variance(y_fit, y_data) ** 0.5

        if verbose:
            print(f"sigma_q: {sigma_q}")
            print(f"nad: {nad}")

        enob = self.n_bits - math.log(nad / sigma_q, 2)
        return enob
    
    def variance(self, fit, data, verbose = False):
        sum = 0

        for i in range(len(fit)):
            sum += (fit[i] - data[i]) ** 2  #MSE
            if verbose:
                print(f"fit {fit[i]} - data {data[i]}")
                print (f"{sum}")
        
        return sum * (self.lsb ** 2) / len(fit)
    

    def sine_function(self, x, A, F, P, C):     # x = argument, A = Amplitude, F = Frequency, P = phase, C = constant
      s = A * np.sin( 2 * np.pi * F * x * 6.25 + P) + C       ## 6.25 = sampling period [ns] if removed it breaks everything
      return s
    
    def sine_fit(self, x, y, p0, bounds, verbose):
        try:
            parameters, covariance = curve_fit(self.sine_function, x, y, p0=p0, bounds = bounds)
            A, F, P, C = parameters

            if verbose:
                print("- amplitude = {:.3f} ADC".format(A))
                print("- frequency = {:.9f} MHz".format(F))
                print("- phase     = {:.3f}".format(P))
                print("- offset    = {:.3f} ADC".format(C))

        except Exception as e:
            print("Fitting error - Unable to fit data")
            parameters = [-1, -1, -1, -1]
            covariance = -1

        return parameters, covariance



    def do_fit(self, y_data, freq=0, verbose = False, graph = False):

        

        x_data          = np.arange(0, len(y_data), 1)
        amplitude_guess = ((max(y_data) - min(y_data)) / 2)
        offset_guess    = ((max(y_data) + min(y_data) )/ 2)
        phase_guess     = 3.14
        frequence_guess = self.freq 
        guess           = [amplitude_guess, frequence_guess,phase_guess, offset_guess]
        lower_bounds = [amplitude_guess * 0.95, self.freq * 0.999, -10, offset_guess - 10]
        upper_bounds = [amplitude_guess * 1.05, self.freq * 1.001, +10, offset_guess + 10]
        param_bounds = (lower_bounds, upper_bounds)

        if verbose:
            print("The guesses are as following:")
            print(f"Amplitude guess:{amplitude_guess}")
            print(f"Offset guess:{offset_guess}")
            print(f"frequency guess {self.freq}")
            print("phase guess")

        y_guess = self.sine_function(x_data, amplitude_guess ,frequence_guess , phase_guess, offset_guess)

        p, c = self.sine_fit(x_data, y_data, guess, param_bounds, verbose=verbose)
        A, F, P, C = p

        y_fit = self.sine_function(x_data, A, F, P, C)

        if graph:
            self.show_data(x_data, y_fit, label="Fit")
            self.show_data(x_data, y_data, label = "Dati")
            self.show_data(x_data, y_guess, label = "Guess")

            plt.legend()
            plt.show()

        return 0, y_fit, y_data, p
        

    def enob(self, filepath = " ", adc = "h",freq = 0, graph = False, verbose = False):
        self.freq = freq * 0.001
        h_data, l_data = file_reader(filepath)
    
        if adc == "l":
            if verbose: print("Ricavando dati da ADCL")
            y_data = l_data
        else:
            if verbose: print("Ricavando dati da ADCH")
            y_data = h_data

        
        result, y_fit, y_data, p = self.do_fit(y_data, self.freq, verbose = verbose, graph = graph)
        enob = self.fit_enob(y_fit, y_data)
        return enob
    
    def enob_plot(self, path):
        txt_files = []
        for root, dirs, files in os.walk(path):  # Walk through directory and subdirectories
            for file in files:
                if file.endswith('.txt'):
                    full_path = os.path.join(root, file)
                    txt_files.append(full_path)

        enob_values = []
        freqs = []

        for file in txt_files:
            print(f"File: {file}")
            filename = os.path.basename(file)
            freq_str = filename.replace('.txt', '')

            try:
                freq = float(freq_str)  # Convert filename to float
            except ValueError:
                print(f"Skipping file {file}, filename does not represent a valid frequency.")
                continue

            freqs.append(freq)
            enob_val = self.enob(file, adc="l", freq=freq)
            enob_values.append(enob_val)

        plt.scatter(freqs, enob_values)
        plt.xlabel("Frequency")
        plt.ylabel("ENOB")
        plt.title("ENOB vs Frequency")
        plt.grid(True)
        plt.show()

        return enob_values
        
    
        




