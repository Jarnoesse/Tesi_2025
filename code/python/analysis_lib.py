import sys
import math
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit
from pathlib import Path
#from termcolor import cprint


def lin_func(x, params):
   a = params[0]
   b = params[1]
   return a*x + b


########################################################################################################################
########################################################################################################################
##
##
class pll_an :

   def pll_analysis(self, pll_dict, outfile, save_mode, quiet=True, show=False):     ## function to analyze and display results from PLL scan
      outfile_txt = outfile[:-3] + "dat"
      fout = open(outfile_txt, "w")
      for val in pll_dict.values():
         fout.write("{}\n".format(val))
      fout.close()

      if save_mode > 0:
         fig = plt.figure(figsize=(8,5))
         plt.plot(pll_dict.values(), 'bo--')
         plt.title("PLL Manual Lock Scan")
         plt.xlabel("# code")
         plt.ylabel("NO PLL-Lock counts")
         plt.grid()
      regs = [reg for reg, cnt in pll_dict.items() if cnt == 0]
      if not(quiet):
         print(regs)
      self.n_lock = len(regs)
      if self.n_lock > 0:        ## check if PLL locked for at least one setting
         ind = [i for i, j in enumerate(list(pll_dict.values())) if j == 0]
         if (len(ind) > 0) and (( np.max(ind) - np.min(ind) + 1 )) == self.n_lock:     ## check if settings with PLL locked are successive
            r = regs[int(self.n_lock/2)]
            reg_0f = r[0]
            reg_10 = r[1:]
            self.reg0f = hex(int(r[0 ], 2))
            self.reg10 = hex(int(r[1:], 2))
            if save_mode > 0:
               plt.plot(ind, [0 for i in range(len(ind))], 'ro', label='PLL lock range = {}'.format(self.n_lock))
               center = int( ( np.max(ind)+0.5 + np.min(ind)+0.5 ) / 2)
               plt.plot([center], [0], 'ks', markersize=8,  label='PLL lock center = {}'.format(center))
               plt.plot([], [], ' ', label="reg0f = {}".format(self.reg0f))
               plt.plot([], [], ' ', label="reg10 = {}".format(self.reg10))
            ret = 0
         else:
            ret = 1
         if save_mode > 0:
            plt.legend()
            if save_mode == 2:
               fig.savefig(outfile)       ## save PLL scan plot
      else:
         ret = 1
      if save_mode > 0:
         if show:
            plt.show()
      
      return ret

   def get_pllReg(self):
      return [self.reg0f, self.reg10]
      
   def get_pllRange(self):
      return self.n_lock
##
##
########################################################################################################################
########################################################################################################################



########################################################################################################################
########################################################################################################################
##
##
class plt_an :

   def check_catiaTP(self, x, y, outfile, save_mode, show=False):
      ret = 0
      slope, const = -1, -1
      #y_th = [4 + 4*xx for xx in x]
      #if y_th != y:
      #   print("\n\nCATIA test-pulse scan points not ok")
      #   print("TP duration measured: {}".format(y))
      #   print("TP duration expected: {}\n\n".format(y_th))
      #   ret = 1
      slope, const = np.polyfit(x, y, 1)
      if ( abs(slope-4) > 0.1 ) or ( abs(const-4) > 1 ):
         print("\n\nCATIA test-pulse scan fit not ok: slope = {0:.3f}, const = {1:.3f}".format(slope, const))
         #print("Expected slope = 4, const = 4\n\n")
         ret = 1
      x_fit = np.arange(np.min(x), np.max(x), 0.1)
      y_fit = [slope*xx + const for xx in x_fit]

      outfile_txt = outfile[:-3] + "dat"
      fout = open(outfile_txt, "w")
      for i in range(len(x)):
         fout.write("{}, {}\n".format(x[i], y[i]))
      fout.close()

      if save_mode > 0:
         fig = plt.figure(figsize=(8,5))
         plt.plot(x, y, 'bo', label="data")
         plt.plot(x_fit, y_fit, 'r--', label="fit: y = {0:.1f} x + {1:.1f}".format(slope, const))
         plt.title("CATIA TP length")
         plt.xlabel("TP reg")
         plt.ylabel("TP duration [clk cycles]")
         plt.grid()
         plt.legend()
         fig.tight_layout()
         if save_mode == 2:
            fig.savefig(outfile)
         if show:
            plt.show()
   
      return ret, slope, const


   def final_result_plot(self, tests_passed, test_time):
      if tests_passed:
         txt = "PASS"
         color = "green"
      else:
         txt = "FAIL"
         color = "red"
      fig = plt.figure(figsize=(5,2))
      ax = fig.add_subplot()
      ax.text(0.35, 0.5, txt, fontsize=40, bbox={'facecolor': color, 'alpha': 0.8, 'pad': 200})
      ax.text(0.2, 0.05, "Test Time = {} s".format(test_time), fontsize=20, color="k")
      plt.axis('off')

##
##
########################################################################################################################
########################################################################################################################



########################################################################################################################
########################################################################################################################
##
##
class sine_wave :
   def __init__(self, filename):
      self.Nbits = 12
      self.LSB = 1
      self.V_fs =  0.6
      self.vmax =  0.6
      self.vmin = -0.6
      self.lsb = (self.vmax - self.vmin) / 2**12
      self.fn = filename
      ## Open decoded file
      f = open(self.fn, "r")
      self.adc_h = []
      self.adc_l = []

      for line in f:
         l = (line.strip()).split(", ")
         self.adc_h.append(int(l[0], 16))
         self.adc_l.append(int(l[1], 16))

   def plot_codes(self, adcH, adcL):
      fig = plt.figure(figsize=(12, 8))
      binning = range(2**12)
      plt.subplot(2, 1, 1)
      nH, codesH, __ = plt.hist(adcH, bins=binning, color="blue")
      plt.title('ADCH')
      plt.xlabel('ADC code')
      plt.ylabel('N')
      plt.grid()
      plt.subplot(2, 1, 2)
      nL, codesL, __  = plt.hist(adcL, bins=binning, color="red")
      plt.title('ADCL')
      plt.xlabel('ADC code')
      plt.ylabel('N')
      plt.grid()
      fig.tight_layout()
   
   def sine_pdf(self, n, A):
      ## Sinewave Probability Density Function
      ## n = ADC code
      ## A = sinewave amplitude
      va = ( n - 0.5 ) * self.lsb;
      vb = ( n + 0.5 ) * self.lsb;
      r = 1/math.pi * ( math.asin( (vb + self.vmin)/A ) - math.asin( (va + self.vmin)/A ) )
      return r
   
   
   def dnl(self, yfit, ydata):
      
      ## Build histogram for ADC codes measured counts
      code_meas, n_meas = np.unique(ydata, return_counts = True)
      zip_iterator = zip(code_meas, n_meas)
      Nmeas = dict(zip_iterator)

      Nsamples = len(ydata)                  # Total number of samples 
      offset = 0.0
      A      = self.V_fs / ( math.sin( Nsamples / ( Nsamples + Nmeas[0] + Nmeas[4095] ) * math.pi/2 ) )      #A = 2.45/4
      print("\nNumber of samples = {}".format(Nsamples))
      print("Sinewave amplitude = {:.3f} V".format(A))
      
      ## Build histogram for ADC codes theory counts
      code_th = np.arange(0, 4097, 1)
      func = np.vectorize(self.sine_pdf)
      n_th = func(code_th, A) * Nsamples
      zip_iterator = zip(code_th, n_th)
      Nth = dict(zip_iterator)
      start = min(code_th)
      stop  = max(code_th)

      sum_low = 0
      sum_high = 0
      for k in Nmeas.keys():
         if k < 2**11:
            sum_low += Nmeas[k]
         else:
            sum_high += Nmeas[k]
      print("Number of samples between [   0,2047]: {}".format(sum_low))
      print("Number of samples between [2048,4095]: {}".format(sum_high))

      ## Build histograms for cumulative counts
      H_meas = dict()
      H_th = dict()
      sum_meas = 0
      sum_th = 0
      for i in range(0, 4096):
         try:
            sum_meas += Nmeas[i]
         except:
            Nmeas[i] = 0
            print("\n********************************************************")
            print("Missing code: {0}\t(expected {1:.0f} counts)".format(i, Nth[i]))
            print("********************************************************")
         H_meas[i] = sum_meas
         sum_th += Nth[i]
         H_th[i] = sum_th
      Nmeas_sort = dict(sorted(Nmeas.items()))
      
      fig = plt.figure(figsize=(12, 9))
      plt.plot(H_meas.keys(), H_meas.values(), '-b', label='measured')
      plt.plot(H_th.keys(),   H_th.values(),   '-r', label='theory')
      plt.title('Cumulative counts')
      plt.xlabel('ADC code')
      plt.ylabel('Cumulative counts')
      plt.legend()
      plt.grid()
      fig.tight_layout()
      
      ## Evaluate transition levels for measured samples and theory
      H_tran_meas = dict()
      H_tran_th = dict()
      for i in range(1, 4096):
         H_tran_meas[i] = offset - A * math.cos( ( math.pi / Nsamples ) * H_meas[i-1] )
         H_tran_th[i] = (i - 0.5) * self.lsb + self.vmin

      ## Evaluate gain and offset corrections to measured transition levels
      slope_meas, const_meas = np.polyfit(list(H_tran_meas.keys()), list(H_tran_meas.values()), 1)
      slope_th, const_th     = np.polyfit(list(H_tran_th.keys()), list(H_tran_th.values()), 1)
      Gain = slope_meas / slope_th
      Voff = const_meas - const_th
      print("Gain correction = {:.5f}".format(Gain)) 
      print("Offset correction = {:.5f} V".format(Voff))

      ## Apply gain and offset corrections to measured transition levels
      H_tran_corr = dict()
      for i in range(1, 4096):
         H_tran_corr[i] = ( H_tran_meas[i] - const_meas ) / Gain + const_th

      fig = plt.figure(figsize=(12, 9))
      plt.plot(H_tran_meas.keys(), H_tran_meas.values(), '-b',  label='measured')
      plt.plot(H_tran_th.keys(),   H_tran_th.values(),   '-r',  label='theory')
      plt.plot(H_tran_corr.keys(), H_tran_corr.values(), '--k', label='corrected')
      plt.title('Transition Levels')
      plt.xlabel('ADC code')
      plt.ylabel('Transition voltages [V]')
      plt.legend()
      plt.grid()
      fig.tight_layout()
		
      ## Evaluate DNL
      dnl_dict = dict()
      for i in range(1, 4095):
         dnl_dict[i] = ( H_tran_corr[i+1] - H_tran_corr[i] ) / self.lsb - 1
      self.plot_dnl(Nmeas_sort.keys(), Nmeas_sort.values(), code_th, n_th, dnl_dict)
         
      return dnl_dict, H_tran_corr, H_tran_th

      
   def plot_dnl(self, xdata, ydata, xfit, yfit, dnl_dict):
      fig = plt.figure(figsize=(12, 9))
      plt.subplot(2, 1, 1)
      plt.plot(xdata, ydata, '-b', label='measured')
      plt.plot(xfit,  yfit,  '-r', label='theory')
      plt.title('Sine-wave ADC codes Histogram')
      plt.xlabel('ADC code')
      plt.ylabel('N')
      plt.legend()
      plt.grid()
      plt.subplot(2, 1, 2)
      plt.plot(dnl_dict.keys(), dnl_dict.values(), '-k')
      plt.title('Differential Non Linearity')
      plt.xlabel('ADC code')
      plt.ylabel('DNL [LSB]')
      plt.grid()
      fig.tight_layout()
      plt.show()
      
   def inl(self, Htran_meas, Htran_th):
      ## Evaluate INL
      inl_dict = dict()
      for i in range(1, 4095):
         epsilon = Htran_meas[i] - Htran_th[i]
         inl_dict[i] = epsilon / self.lsb
      fig = plt.figure(figsize=(12, 9))
      plt.plot(inl_dict.keys(), inl_dict.values(), '-k')
      plt.title('Integral Non Linearity')
      plt.xlabel('ADC code')
      plt.ylabel('INL [LSB]')
      plt.grid()
      fig.tight_layout()
      plt.show()



   def plot_sine(self, zoom=False, show=True, sat=False):
      if zoom:
         self.fig = plt.figure(figsize=(28, 11), dpi=80)
         plt.subplot(2, 2, 1)
      else:
         self.fig = plt.figure(figsize=(14, 9), dpi=80)
         plt.subplot(2, 1, 1)
      if sat:
         plt.title('ADCH (ADC saturation test)')
      else:
         plt.title('ADCH')
      plt.ylabel('ADC code')
      plt.xlabel('Sample #')
      plt.grid()
      plt.plot(self.adc_h, color="dodgerblue",  markeredgecolor = 'black',
               marker = ".", markersize = 8, linestyle = 'dotted')

      if zoom:
         plt.subplot(2, 2, 3)
      else:
         plt.subplot(2, 1, 2)
      if sat:
         plt.title('ADCL (ADC saturation test)')
      else:
         plt.title('ADCL')
      plt.ylabel('ADC code')
      plt.xlabel('Sample #')
      plt.grid()
      plt.plot(self.adc_l, color="orange", markeredgecolor = 'black',
               marker = ".", markersize = 8, linestyle = 'dotted')
      if zoom:
         N = int(len(self.adc_h)/2)
         MAX = max(self.adc_h[N-100 : N+100])
         ind_max = self.adc_h[N-100 : N+100].index(MAX)
         xmin = (N-100) + ind_max - 100
         xmax = (N-100) + ind_max + 450
         plt.subplot(2, 2, 2)
         if sat:
            plt.title('ADCH zoom (ADC saturation test)')
         else:
            plt.title('ADCH zoom')
         plt.ylabel('ADC code')
         plt.xlabel('Sample #')
         plt.grid()
         plt.plot(self.adc_h, color="dodgerblue",  markeredgecolor = 'black',
                  marker = ".", markersize = 8, linestyle = 'dotted')
         plt.xlim(xmin, xmax)
         plt.subplot(2, 2, 4)
         if sat:
            plt.title('ADCL zoom (ADC saturation test)')
         else:
            plt.title('ADCL zoom')
         plt.ylabel('ADC code')
         plt.xlabel('Sample #')
         plt.grid()
         plt.plot(self.adc_l, color="orange", markeredgecolor = 'black',
                  marker = ".", markersize = 8, linestyle = 'dotted')
         plt.xlim(xmin, xmax)

      self.fig.tight_layout()
      if show:
         plt.show(block=True)
      
      
   def save_plot(self, outf="dummy", quiet=True):
      if outf == "dummy":
         outfile = self.fn.split(".")[0] + "_sine.png"
      else:
         outfile = outf.split(".")[0] + "_sine.png"
      if not(quiet):
         print("Saving figure as {}".format(outfile))
      self.fig.savefig(outfile)



   def sine_func(self, x, A, F, P, C):     # x = argument, A = Amplitude, F = Frequency, P = phase, C = constant
      s = A * np.sin( 2 * np.pi * F / 1000 * x * 6.25 + P) + C       ## 6.25 = sampling period [ns]
      return s
	   
   def var(self, fit, data):
      s = 0
      i = 0
      while i < len(fit):
         s = s + ( (fit[i] - data[i]) ** 2 )
         i += 1
      r = s * (self.LSB ** 2) / len(fit)
      return r
	   
   def res_list(self, fit, data):
      res_list = []
      i = 0
      while i < len(fit): 
         res_list.append(fit[i] - data[i])
         i += 1
      return res_list
	   
   def do_fit(self, freq=0, adc="H", quiet=False):
      if freq == 0:
         ydata = self.adc_h
         temp = self.fn.split("_")
         self.in_freq = float(temp[-2] + "." + temp[-1][ : -4])
         if self.in_freq > 60:
            self.in_freq = self.in_freq / 1000     ## transform kHz into MHz for 500kHz input frequency
         print("Input frequency from filename: {} MHz".format(self.in_freq))
      else:
         self.in_freq = freq
         if adc == "H":
            if not(quiet):
               print("\nFitting sine wave of ADCH")
            ydata = self.adc_h[0 : 1000]
         else:
            if not(quiet):
               print("\nFitting sine wave of ADCL")
            ydata = self.adc_l[0 : 1000]
         
      xdata = np.arange(0, len(ydata), 1)
      amplitude_guess = (( max(ydata) - min(ydata) ) / 2)
      offset_guess    = (( max(ydata) + min(ydata) ) / 2)
      guess = [ amplitude_guess, self.in_freq, 0.1, offset_guess ] 
      lower_bounds = [amplitude_guess * 0.95, self.in_freq * 0.999, -10, offset_guess - 10]
      upper_bounds = [amplitude_guess * 1.05, self.in_freq * 1.001, +10, offset_guess + 10]
      param_bounds = (lower_bounds, upper_bounds)
      p, c = self.sine_fit(xdata, ydata, guess, param_bounds, quiet=quiet)
      A, F, P, C = p
      if A < 1500:
         if not(quiet):
            print("Fit Error: sine waveform amplitude is too low: {}".format(A))
         return 1, p
      elif abs(C - 2048) > 10:
         if not(quiet):
            print("Fit Error: sine waveform offset is not centered: {}".format(C))
         return 1, p
      else:
         xfit = np.arange(min(xdata), max(xdata)+1, 1)
         yfit = self.sine_func(xfit, A, F, P, C)
         return 0, yfit, ydata, p


   def enob_from_fit(self, yfit, ydata):
      sigma_Q = ((self.LSB * self.LSB)/12) ** (0.5)
      NAD     = self.var(yfit, ydata) ** (0.5)
      enob    = self.Nbits - math.log( NAD/sigma_Q , 2 )
      return enob


   def do_fft(self, x, data, realFFT=False):
      delta_t = x[1] - x[0]      ## time bin
      nPts = len(data)
      #print("Number of samples = {}".format(nPts))
      if realFFT:
         ft = np.fft.rfft(data)
         freq = np.fft.rfftfreq(len(x)) / delta_t           # Frequency bin of each component
      else:
         ft = np.fft.fft(data)
         freq = np.fft.fftfreq(len(x)) / delta_t            # Frequency bin of each component

      #print("FFT freq bin = {} MHz".format(freq[1]-freq[0]))
      ampl = abs(ft)             # Amplitude of each frequency is given by the absolute value
      return freq, ampl   
		
   def enob_from_fft(self, ADC="ADCH", is_plot=False, Ncut=1, quiet=True):
      ## Spectral leakage is larger at higher frequencies: to avoid any spectral leakage around the fundamental, 
      ## often several bins around the fundamental are ignored. This is set by the "Ncut" parameter
      Fs = 160
      realFFT = False
      if ADC == "ADCH":
         if not(quiet):
            print("\nEvaluating ENOB from FFT for data from ADCH")
         adc = self.adc_h
      else:
         if not(quiet):
            print("\nEvaluating ENOB from FFT for data from ADCL")
         adc = self.adc_l
      
      Nsamples = len(adc)
      ns = 2**14
      Ndiv = int(len(adc) / ns)
      if not(quiet):
         print("Data has {} samples".format(Nsamples))
         print("Data will be split in {} parts, with {} samples each".format(Ndiv, ns))
      a = []
      for i in range(Ndiv):
         y0 = adc[i * ns : (i+1) * ns]      ## split data samples in parts, with ns=16384 samples each
         t0 = np.arange(0, len(y0) * 1/Fs - 1/Fs, 1/Fs)
         freq, ampl = self.do_fft(t0, y0, realFFT)
         a.append(ampl)
      fft_average = np.average(np.array(a), axis=0)   ## take average to reduce floor noise
      enob, enob2, maxFreq = self.fft_enob(freq, fft_average, Ncut, quiet=quiet)
      if is_plot:
         self.plot_fft(freq, fft_average, ADC, '-b')
      
      return enob, enob2, maxFreq
      
      
   def fft_enob(self, freq, ampl, Ncut=1, quiet=True):
      nBins = len(freq)
      #print("Number of bins = {}".format(nBins))
      ## Remove DC component
      a0 = ampl[1: ]
      f0 = freq[1: ]

      maxInd   = np.argmax(a0[ : 2**13 - 1])    ## main frequency index (should be input signal frequency)
      maxInd2  = 2**14 - maxInd - 2             ## main frequency index mirrored in the range [-Fs/2, 0]
      maxFreq  = f0[maxInd]                     ## main frequency
      maxFreq2 = f0[maxInd2]                    ## main frequency mirrored in the range [-Fs/2, 0]
      maxVal   = a0[maxInd]                     ## value of main frequency
      maxVal2  = a0[maxInd2]                    ## value of mirrored main frequency 
      if not(quiet):
         print("Main freq = {} MHz".format(maxFreq))
      Ps  = a0[maxInd]**2                       ## signal power level
      a_rms = ( (2 * Ps)**0.5 ) / nBins         ## rms signal amplitude	
      sqSum = 0                                 ## noise (floor) + distortion (harmonic) power level (N+D) 			
      for b in range(len(a0)):		
         if ( abs(b - maxInd) < Ncut ) or ( abs(b - maxInd2) < Ncut ):
            pass
         else:
            sqSum += a0[b]**2

      NAD = ( sqSum**0.5 ) / ( ( nBins * (nBins-3) )**0.5 )
      enob = 12 - math.log2( NAD / ( 1 / (12)**0.5 ) )      
      SINAD = 20 * math.log10(a_rms / NAD)
      enob2 = (SINAD - 1.76) / 6.02 						## ENOB from SINAD
      #print("NAD   = {}".format(NAD))
      #print("SINAD = {} dB".format(SINAD))
      return enob, enob2, maxFreq


   def plot_fft(self, freq, ampl, title, style, ms=8):
      fig = plt.figure(figsize=(12, 8))
      maxInd = np.argmax(ampl[ 1 : int(len(ampl)/2) ]) + 1
      maxVal = ampl[maxInd]
      Fmax   = freq[maxInd]
      ## Shift negative frequency values to range above Nyquist frequency
      for i, f in enumerate(freq):
         if f < 0:
            freq[i] = f + 160
      f2 = []
      a2 = []
      for i, ai in enumerate(ampl):
         if ai > 0:
            f2.append(freq[i])
            a2.append(20 * math.log10(ai / maxVal))
      plt.plot(f2, a2, style, markersize=ms, label="average")
      plt.xlim([-1,80.001])
      plt.xlabel('frequency [MHz]', fontsize=14)
      plt.ylabel(r"$ 20 \cdot$ log $_{10} $ $ (\frac{A(f)}{A_0}) $  [dB]", fontsize=14)
      plt.title('{}, Fin = {} MHz'.format(title, Fmax), fontsize=16)
      plt.tick_params(axis='both', which='major', labelsize=12)
      plt.grid()


   def check_enob(self, enob, adc="ADCH", enob_limit=8, quiet=False):
      if enob < enob_limit:
         if not(quiet):
            print("Error: ENOB of {} too low!!!".format(adc))
         return 1
      else:
         return 0

   def sine_fit(self, x, y, p0, bounds, quiet=False):
      try:
         param, cov = curve_fit(self.sine_func, x, y, p0 = p0, bounds = bounds)
         A, F, P, C = param
      except:
         print("Unable to fit data")
         param = [-1, -1, -1, -1]
         cov = -1
      if not(quiet):
         print("- amplitude = {:.3f} ADC".format(A))
         print("- frequency = {:.9f} MHz".format(F))
         print("- phase     = {:.3f}".format(P))
         print("- offset    = {:.3f} ADC".format(C))
      return param, cov

   def check_sat(self, quiet=False):
      retH, retL = 0, 0

      ADCH_min = np.min(self.adc_h) 
      ADCH_max = np.max(self.adc_h)
      ADCL_min = np.min(self.adc_l)
      ADCL_max = np.max(self.adc_l)

      if ADCH_min > 0:
         retH = 1
         if not(quiet):
            print("Error: ADCH saturation MIN not 0: {}".format(ADCH_min))
      if ADCH_max < 4095:
         retH = 1
         if not(quiet):
            print("Error: ADCH saturation MAX not 0: {}".format(ADCH_max))
      if ADCL_min > 0:
         retL = 1
         if not(quiet):
            print("Error: ADCL saturation MIN not 0: {}".format(ADCL_min))
      if ADCL_max < 4095:
         retL = 1
         if not(quiet):
            print("Error: ADCL saturation MAX not 0: {}".format(ADCL_max))

      return retH, retL

   def delete_rawData(self, path, filename, quiet=False):
      path_to_remove = Path(path + filename)
      if not(quiet):
         print("Removing file {}".format(path_to_remove))
      Path.unlink(path_to_remove, missing_ok=False)

   def get_input_freq(self):
      return self.in_freq
      
   def plot_enob(self, freq, enobFIT, enobFFT):
      fig = plt.figure(figsize=(10, 6))
      plt.plot(freq, enobFIT, 'bo', label="FIT")
      plt.plot(freq, enobFFT, 'ro', label="FFT")
      plt.ylim([8, 11])
      plt.xlabel("freq [ MHz ]")
      plt.ylabel("ENOB [ bit ]")
      plt.title("ENOB vs freq")
      plt.grid()
      plt.legend()
      plt.show()
##
##
########################################################################################################################
########################################################################################################################





########################################################################################################################
########################################################################################################################
##
##
class pulse_signal :
   def __init__(self, filename, style="raw"):
      self.fn = filename                                          ## path + filename
      lines = open(self.fn).readlines()                           ## read full file
      #countlines = len(lines)                                    ## check lenght of file
      #print("number of lines in file: ", countlines)
      if style == "raw":                                             ## 4 columns with binary format data
         data = [l.split(", ") for l in lines]                       ## generate matrix data
         data_list = list(map(lambda *temp: list(temp), *data))      ## transpose matrix
         data_list[3] = [ d.strip() for d in data_list[3] ]          ## remove "\n" after data for dout_3 column
         self.dout_0 = data_list[0]

      elif style == "hex":                                              ## one column with hex data
         self.dout_0 = [bin(int(l, 16))[2:].zfill(32) for l in lines]   ## convert hex to binary and keep 32-bit format
      elif style == "hex_pd":                                           ## more column, first is hex data
         data = [l.split(" ")[0] for l in lines]                        ## take first column
         self.dout_0 = [bin(int(d, 16))[2:].zfill(32) for d in data]    ## convert hex to binary and keep 32-bit format

      # Bit-pattern for 32-bit words header
      self.align_pattern   = "01011010010110100101101001011010"   ## 5a5a5a5a
      self.base_pattern    = "01"                                 ## Baseline data format    (full)
      self.base_pattern_1  = "10"                                 ## Baseline data format 1  (not full)
      self.sign_pattern    = "001010"                             ## Signal data format      (full)
      self.sign_pattern_1  = "001011"                             ## Signal data format 1    (not full)
      self.frame_pattern   = "1101"                               ## Frame delimeter --> inserted after 50 data words (signal or baseline)
      self.idle_pattern    = "1110"                               ## Idle word       --> inserted to keep sync with receiver (when only baseline words)



   def crc12(self, data_in, crc_in):
       ## Function to calculate CRC12 for a sequence of 32-bit data words (data_in)
       ## crc_in is the CRC12 evaluated in the previous step; it is equal to 12'b0 for the first step

       bit_0  = data_in[30] ^data_in[29] ^data_in[26] ^data_in[25] ^data_in[24] ^data_in[23] ^data_in[22] ^data_in[17] ^data_in[16] ^data_in[15] ^data_in[14] ^data_in[13] ^data_in[12] ^data_in[11] ^data_in[8] ^data_in[7] ^data_in[6] ^data_in[5] ^data_in[4] ^data_in[3] ^data_in[2] ^data_in[1] ^data_in[0] ^crc_in[2] ^crc_in[3] ^crc_in[4] ^crc_in[5] ^crc_in[6] ^crc_in[9] ^crc_in[10]
       bit_1  = data_in[31] ^data_in[29] ^data_in[27] ^data_in[22] ^data_in[18] ^data_in[11] ^data_in[9] ^data_in[0] ^crc_in[2] ^crc_in[7] ^crc_in[9] ^crc_in[11]
       bit_2  = data_in[29] ^data_in[28] ^data_in[26] ^data_in[25] ^data_in[24] ^data_in[22] ^data_in[19] ^data_in[17] ^data_in[16] ^data_in[15] ^data_in[14] ^data_in[13] ^data_in[11] ^data_in[10] ^data_in[8] ^data_in[7] ^data_in[6] ^data_in[5] ^data_in[4] ^data_in[3] ^data_in[2] ^data_in[0] ^crc_in[2] ^crc_in[4] ^crc_in[5] ^crc_in[6] ^crc_in[8] ^crc_in[9]
       bit_3  = data_in[27] ^data_in[24] ^data_in[22] ^data_in[20] ^data_in[18] ^data_in[13] ^data_in[9] ^data_in[2] ^data_in[0] ^crc_in[0] ^crc_in[2] ^crc_in[4] ^crc_in[7]
       bit_4  = data_in[28] ^data_in[25] ^data_in[23] ^data_in[21] ^data_in[19] ^data_in[14] ^data_in[10] ^data_in[3] ^data_in[1] ^crc_in[1] ^crc_in[3] ^crc_in[5] ^crc_in[8]
       bit_5  = data_in[29] ^data_in[26] ^data_in[24] ^data_in[22] ^data_in[20] ^data_in[15] ^data_in[11] ^data_in[4] ^data_in[2] ^crc_in[0] ^crc_in[2] ^crc_in[4] ^crc_in[6] ^crc_in[9]
       bit_6  = data_in[30] ^data_in[27] ^data_in[25] ^data_in[23] ^data_in[21] ^data_in[16] ^data_in[12] ^data_in[5] ^data_in[3] ^crc_in[1] ^crc_in[3] ^crc_in[5] ^crc_in[7] ^crc_in[10]
       bit_7  = data_in[31] ^data_in[28] ^data_in[26] ^data_in[24] ^data_in[22] ^data_in[17] ^data_in[13] ^data_in[6] ^data_in[4] ^crc_in[2] ^crc_in[4] ^crc_in[6] ^crc_in[8] ^crc_in[11]
       bit_8  = data_in[29] ^data_in[27] ^data_in[25] ^data_in[23] ^data_in[18] ^data_in[14] ^data_in[7] ^data_in[5] ^crc_in[3] ^crc_in[5] ^crc_in[7] ^crc_in[9]
       bit_9  = data_in[30] ^data_in[28] ^data_in[26] ^data_in[24] ^data_in[19] ^data_in[15] ^data_in[8] ^data_in[6] ^crc_in[4] ^crc_in[6] ^crc_in[8] ^crc_in[10]
       bit_10 = data_in[31] ^data_in[29] ^data_in[27] ^data_in[25] ^data_in[20] ^data_in[16] ^data_in[9] ^data_in[7] ^crc_in[0] ^crc_in[5] ^crc_in[7] ^crc_in[9] ^crc_in[11]
       bit_11 = data_in[29] ^data_in[28] ^data_in[25] ^data_in[24] ^data_in[23] ^data_in[22] ^data_in[21] ^data_in[16] ^data_in[15] ^data_in[14] ^data_in[13] ^data_in[12] ^data_in[11] ^data_in[10] ^data_in[7] ^data_in[6] ^data_in[5] ^data_in[4] ^data_in[3] ^data_in[2] ^data_in[1] ^data_in[0] ^crc_in[1] ^crc_in[2] ^crc_in[3] ^crc_in[4] ^crc_in[5] ^crc_in[8] ^crc_in[9]

       crc = [bit_11, bit_10, bit_9, bit_8, bit_7, bit_6, bit_5, bit_4, bit_3, bit_2, bit_1, bit_0]
       crc_out = crc[::-1]

       return crc_out



   def decode_data(self):

      ret = 0
      
      default_crc = [0 for i in range(12)]
      new_crc = default_crc

      is_first_frame = True
      frame_number = -99
      Ntot = 0
      ndatawords = 0
      Nsamples = 0
      frame_nsamples = 0

      base_count = 0
      base_full_count = 0
      base_incomplete_count = 0
      signal_count = 0
      signal_full_count = 0
      signal_incomplete_count = 0
      frame_delimeter_count = 0
      idle_count = 0
      error_count = 0

      self.samples = []
      types = []

      self.time_base = []
      self.time_g10 = []
      self.time_g1 = []
      self.samples_base = []
      self.samples_g10 = []
      self.samples_g1 = []

      frame_rollover = 0
      frame_num_list = []
      frame_size_list = []


      for d0 in self.dout_0:
      
          if d0 == self.align_pattern:
              print("\n\n************************************")
              print("Error in decoding data!!")
              print("Found align pattern")
              print("************************************\n\n")
              ret = 1
              break

          else:
          
              if d0[0:2] == self.base_pattern:
                  is_data_word = True
                  data_word = d0
                  ndatawords += 1
                  base_full_count += 1
                  Nsamples += 5
                  sign_out = [ int(d, 2) for d in [ d0[26:32], d0[20:26], d0[14:20], d0[8:14], d0[2:8] ] ]
                  for s in sign_out:
                      self.samples.append(s)
                      self.samples_base.append(s)
                      self.time_base.append(Ntot)
                      Ntot += 1
                      base_count += 1
                      types.append(0)

              elif d0[0:2] == self.base_pattern_1:
                  is_data_word = True
                  data_word = d0
                  ndatawords += 1
                  base_incomplete_count += 1
                  nbase = int ( d0[2:8], 2 )      ## number of baseline samples contained in this data word
                  Nsamples += nbase
                  sign_out = []
                  for n in range(nbase):
                      sign_out.append( int( d0[26-n*6 : 32-n*6], 2 ) )
                  for s in sign_out:
                      self.samples.append(s)
                      types.append(1)
                      self.samples_base.append(s)
                      self.time_base.append(Ntot)
                      Ntot += 1
                      base_count += 1

              elif d0[0:6] == self.sign_pattern:
                  is_data_word = True
                  data_word = d0
                  ndatawords += 1
                  signal_full_count += 1
                  Nsamples += 2
                  for n in range(2):
                      s = ( int( d0[ 20-n*13 : 32-n*13 ], 2 ) )
                      self.samples.append(s)
                      gain_sel = d0[19-n*13]
                      if gain_sel == "1":
                          types.append(2)
                          self.samples_g1.append(s)
                          self.time_g1.append(Ntot)
                      elif gain_sel == "0":
                          types.append(3)
                          self.samples_g10.append(s)
                          self.time_g10.append(Ntot)
                      Ntot += 1
                      signal_count += 1

              elif d0[0:6] == self.sign_pattern_1:
                  is_data_word = True
                  data_word = d0
                  ndatawords += 1
                  signal_incomplete_count += 1
                  Nsamples += 1
                  s = int ( d0[20:32], 2 )
                  self.samples.append(s)
                  gain_sel = d0[19]
                  if gain_sel == "1":
                      types.append(2)
                      self.samples_g1.append(s)
                      self.time_g1.append(Ntot)
                  elif gain_sel == "0":
                      types.append(3)
                      self.samples_g10.append(s)
                      self.time_g10.append(Ntot)
                  Ntot += 1
                  signal_count += 1

              elif d0[0:4] == self.idle_pattern:
                  is_data_word = False
                  idle_count += 1
                  if d0[4:32] != "1010101010101010101010101010":
                      print("\n\n************************************")
                      print("Error in decoding data!!")
                      print("Error in IDLE pattern: {}".format( d0[4:32] ) )
                      ret = 1
                      print("************************************\n\n")
                      break
                  
              elif d0[0:4] == self.frame_pattern:
                  is_data_word = False
                  frame_delimeter_count += 1
                  frame_num = int( d0[24:32], 2 )
                  crc_12 = int( d0[12:24], 2 )
                  frame_nsamples = int( d0[4:12], 2 )

                  if frame_num == 0:
                      frame_rollover +=1

                  frame_number = frame_num + 2**8 * frame_rollover
                  if (not(is_first_frame)) and (frame_number != frame_number_old + 1):
                      print("\n\n************************************")
                      print("Error in decoding data!!")
                      print("Error on frame number: {} (new) vs {} (old)".format(frame_number, frame_number_old))
                      print("************************************\n\n")
                      ret = 2
                      break
                     
                  frame_num_list.append(frame_number)
                  frame_size_list.append(frame_nsamples)
                  frame_number_old = frame_number

                  if not(is_first_frame):
                      if (ndatawords != 50) and (ndatawords != 51):
                          print("\n\n************************************")
                          print("Error in decoding data!!")
                          print("Error on number of data words inside the frame (n.{}): {} vs 50".format(frame_number, ndatawords))
                          print("************************************\n\n")
                          ret = 2
                          break
                          
                      if Nsamples != frame_nsamples:
                          print("\n\n************************************")
                          print("Error in decoding data!!")
                          print("Error on number of samples inside the frame: {} vs {} ".format(Nsamples, frame_nsamples))
                          print("************************************\n\n")
                          ret = 2 
                          break
          
                      #############################################################################################
                      ##
                      final_crc = int(''.join([str(b) for b in new_crc[::-1]]), 2)
                      #print( "\nFrame Number = {}".format(frame_number) )
                      #print( "Frame # samples = {}".format(frame_nsamples) )
                      #print("calc CRC12 = {0:b}".format(final_crc))
                      #print("data CRC12 = {0:b}".format(crc_12) )
                      if final_crc != crc_12:
                          print("\n\n************************************")
                          print("Error in decoding data!!")
                          print("Error on CRC12: {} vs {}".format(final_crc, crc_12))
                          print("************************************\n\n")
                          ret = 2
                          break
                      ##
                      #############################################################################################

                  if is_first_frame:
                      is_first_frame = False

                  ndatawords = 0
                  Nsamples = 0
                  new_crc = default_crc

              else:
                  is_data_word = False
                  error_count += 1
                  print("\n\n************************************")
                  print("Error in decoding data!!")
                  print("Unknown word: {}".format(d0))
                  print("************************************\n\n")
                  ret = 2
                  break

              if is_data_word:        ## update CRC12
                  input_bitstring = data_word
                  rev_bitstring = input_bitstring[::-1]           ## invert bits order
                  data_in = [int(v,2) for v in rev_bitstring]     ## convert string to bit list
                  old_crc = new_crc
                  new_crc = self.crc12(data_in, old_crc)

      if len(self.samples) != len(types):
          print("Error on lists length")
          ret = 1

      return ret
      
      
   def check_bsln(self, quiet=True):
      MIN = np.min(self.samples)
      MAX = np.max(self.samples)
      MEAN = np.mean(self.samples)
      STD = np.std(self.samples)
      if not(quiet):
         print("")
         print("MIN baseline value = {}".format(MIN))
         print("MAX baseline value = {}".format(MAX))
         print("MEAN baseline value = {:.1f}".format(MEAN))
         print("STD baseline value = {:.1f}".format(STD))
         print("")
      return MEAN

   def bsln_fit(self, bsln_sub, bsln_list, target_bsln):
      try:
         slope, const = np.polyfit(bsln_sub, bsln_list, 1)
         baseline_sub = int( ( target_bsln - const ) / slope )
      except:
         print("Unable to fit data")
         baseline_sub = -1    ## np.max(bsln_sub)
      return baseline_sub

   def bsln_txt(self, bsln_sub, bsln_list, outfile):
      fout = open(outfile, "w")
      for i in range(len(bsln_sub)):
         fout.write("{}, {}\n".format(bsln_sub[i], round(bsln_list[i], 2)))
      fout.close()
      
   def bsln_plot(self, bsln_sub, bsln_list, sub_value, outfile, save_fig=False):
      fig, ax = plt.subplots(1, figsize=(8,5))
      ax.plot(bsln_sub, bsln_list, "bo")
      ax.plot([sub_value, sub_value], [0, np.max(bsln_list)], "k--", label="reg={}".format(hex(sub_value)))
      ax.set_title('ADCH Baseline Subtraction', fontsize = 16)
      ax.set_xlabel('baseline subtraction value [digits]', fontsize = 14)
      ax.set_ylabel('baseline mean [ADC]', fontsize = 14)
      ax.grid()
      ax.legend()
      fig.tight_layout()
      if save_fig:
         fig.savefig(outfile)


   def check_samples(self, pulse_type=None, quiet=True):
      N_tot  = len(self.samples)
      N_bsln = len(self.samples_base)
      N_G10  = len(self.samples_g10)
      N_G1   = len(self.samples_g1)
      bsln_mean = np.mean(self.samples_base)
      if not(quiet):
         print("\nTotal number of samples = {}".format(N_tot))
         print("- Number of Baseline samples = {} ({:.2f}%)".format(N_bsln, N_bsln/N_tot*100))
         print("- Number of Gain 10 samples = {} ({:.2f}%)".format(N_G10, N_G10/N_tot*100))
         print("- Number of Gain  1 samples = {} ({:.2f}%)".format(N_G1, N_G1/N_tot*100))
         print("\nBaseline mean value = {:.1f} ADC\n".format(bsln_mean))
      if N_tot > 0:
         """
         pulse frequency is 4.5 kHz --> 0.222 ms
         acquisition window is about 40k samples --> 0.25 ms
         sometimes we can have 2 separated pulses in the same acquired window
         """
         index_max = [self.samples.index(i) for i in sorted(self.samples, reverse=True)][:4]
         index_sorted = sorted(index_max)
         index_diff = [index_sorted[i] - index_sorted[i-1] for i in range(1, len(index_sorted))]
         if not(quiet):
            print("MAX at: {}".format(index_max))
         if np.all(np.array(index_diff)==1):
            min_ratio = 0.999
         else:
            min_ratio = 0.998
            if not(quiet):
               print("Found two separate peaks at {} and {}\n".format(index_max[0], index_max[1]))
         
         if N_bsln/N_tot < min_ratio:
            if not(quiet):
               print("Error: too many signals!!!\n")
            return 1, N_tot, N_bsln, N_G10, N_G1
         elif ((pulse_type=="normal") and (N_G1>0)):
            if not(quiet):
               print("Error: found G1 data!!!\n")
            return 1, N_tot, N_bsln, N_G10, N_G1
         elif ((pulse_type=="normal") and (N_G10==0)):
            if not(quiet):
               print("Error: pulse with G10 data not found!!!\n")
            return 1, N_tot, N_bsln, N_G10, N_G1
         elif ((pulse_type=="saturation") and (N_G10>0)):
            if not(quiet):
               print("Error: found G10 data!!!\n")
            return 1, N_tot, N_bsln, N_G10, N_G1
         elif ((pulse_type=="saturation") and (N_G1==0)):
            if not(quiet):
               print("Error: pulse with G1 data not found!!!\n")
            return 1, N_tot, N_bsln, N_G10, N_G1
         elif (pulse_type!="normal") and (pulse_type!="saturation"):
            if not(quiet):
               print("Warning: no pulse type ('normal' or 'saturation') provided\n")
            return 1, N_tot, N_bsln, N_G10, N_G1
         else:
            return 0, N_tot, N_bsln, N_G10, N_G1
      else:
         return 1, N_tot, N_bsln, N_G10, N_G1

   def pulse_func(self, t, bsln, scale, T, t0):
      return bsln + scale * ( 1.0 / (np.exp( (-t + t0) / T ) + 1) )

   def plot_pulse(self, save_mode, pulse_type, isFit=True, show=True, isZoom=True, reduce_data=True, quiet=False):
      ret = 0
      param = [-1, -1, -1, -1]
      max_val = max(self.samples[30:-30])         ## find ADC max: FC 24.11.2023 --> avoid to select max at the edges of the samples
      if max_val == 0:
         print("Error: no signal in data")
         ret = 1
      else:
         if save_mode > 0:
            self.fig, ax = plt.subplots(1, figsize=(10,6))

         max_index = self.samples.index(max_val)   ## find index of max point
         start_index = 0                           ## find first point to be used in the fit
         for i in range(max_index):
            if self.samples[max_index-i] == 0:
               start_index = max_index - i - 3     ## take 4 baseline points (ADC=0) before pulse signal
               break
         if isFit:
            param = self.pulse_fit(start_index, max_index, max_val, quiet=quiet)
            if param[0] == -1:      ## --> FC 21.11.2023
               ret = 1
            else:
               if save_mode > 0:
                  xfit = np.arange(start_index, max_index+0.2, 0.1)
                  yfit = self.pulse_func(xfit, *param)
                  ax.plot(xfit, yfit, "b--", label = "fit")
         if save_mode > 0:
            ax.plot(self.samples, "-k", label="all")
            ax.plot(self.time_base, self.samples_base, "og", label="baseline")
            ax.plot(self.time_g10, self.samples_g10, "or", label="gain10")
            ax.plot(self.time_g1, self.samples_g1, "ob", label="gain1")
            ax.set_title('Signal (Pulse type: {})'.format(pulse_type), fontsize = 16)
            ax.set_xlabel('samples', fontsize = 14)
            ax.set_ylabel('ADC', fontsize = 14)
            if isZoom:
               ax.set_xlim(max_index-15, max_index+15)   ## set plot default display range
            ax.grid()
            ax.legend()
            if show:
               plt.show()

         if reduce_data:
            self.reduce_data(max_index, saved_samples=50, quiet=quiet)

      return ret, max_val, *param


   def reduce_data(self, max_index, saved_samples=50, quiet=False):
      outfile = self.fn[:-3] + "reduced.dat"
      fout = open(outfile, "w")
      N = int(saved_samples/2)
      for i in range(saved_samples):
         fout.write("{}, {}\n".format(max_index-N+i, self.samples[max_index-N+i]))
      if len(self.time_g10) > 0:
         fout.write("\n\nG10\n")
         for i in range(len(self.time_g10)):
            fout.write("{}, {}\n".format(self.time_g10[i], self.samples_g10[i]))
      if len(self.time_g1) > 0:
         fout.write("\n\nG1\n")
         for i in range(len(self.time_g1)):
            fout.write("{}, {}\n".format(self.time_g1[i], self.samples_g1[i]))
      fout.close()
      path_to_remove = Path(self.fn)
      if not(quiet):
         print("Removing file {}".format(path_to_remove))
      Path.unlink(path_to_remove, missing_ok=False)

   
   def delete_rawData(self, path, filename, quiet=False):
      path_to_remove = Path(path + filename)
      if not(quiet):
         print("Removing file {}".format(path_to_remove))
      Path.unlink(path_to_remove, missing_ok=False)


   def pulse_fit(self, start_index, max_index, max_val, quiet=False):
      xdata = np.arange(start_index, max_index+1)
      ydata = self.samples[start_index : max_index+1]
      guess = [ 0.1, max_val, 0.5, max_index-3 ] 
      lower_bounds = [0, max_val * 0.98, 0.01, max_index-9]
      upper_bounds = [100, max_val * 1.02, 10,   max_index]
      param_bounds = (lower_bounds, upper_bounds)
      try:
         param, cov = curve_fit(self.pulse_func, xdata, ydata, p0 = guess, bounds = param_bounds)
         B, A, T, t0 = param
         if not(quiet):
            print("\nFit results:")
            print("- Baseline = {:.2f} ADC".format(B))
            print("- MAX      = {:.2f} ADC".format(A))
            print("- Width    = {:.2f} ".format(T))
            print("- t0       = {:.1f} samples".format(t0))
         yfit = self.pulse_func(xdata, *param)
         RMSE = self.rms(ydata, yfit)
         if not(quiet):
            print("\nRoot-mean-square error = {:.2f}\n".format(RMSE))
         if (RMSE > 600) or (A < 2500):
            print("Bad results from fit")
            param = [-1, -1, -1, -1]         ## --> FC 21.11.2023
      except:
         print("Unable to fit data")
         param = [-1, -1, -1, -1]         ## --> FC 20.02.2024
      return param
      
   def rms(self, ydata, yfit):
      err = yfit - ydata
      SE = np.square(err)        # squared errors
      MSE = np.mean(SE)          # mean squared errors
      RMSE = np.sqrt(MSE)        # Root Mean Squared Error
      return RMSE
      
   def save_plot(self, outf=None, gain="dummy", quiet=True):
      if outf == None:
         fout = self.fn.split(".")[0]
      else:
         fout = outf.split(".")[0]

      if gain == "normal":
         outfile = fout + "_pulseG10.png"
      elif gain == "saturation":
         outfile = fout + "_pulseG1.png"
      else:
         outfile = fout + "_pulse.png"
      
      if not(quiet):
         print("Saving figure as {}".format(outfile))
      self.fig.savefig(outfile)
##
##
########################################################################################################################
########################################################################################################################











