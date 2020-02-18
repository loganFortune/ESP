# -*- coding: utf-8 -*-
"""
Author : LF - SM -TD

Date : 01/02/2020
"""

import nolds.examples
import nolds
from edf import edf_data
import matplotlib.pyplot as plt
import os
import PyGnuplot as gp

def test_nolds_package_service_accuracy():
    nolds.examples.plot_lyap("tent")
    nolds.examples.plot_lyap("logistic")

def test_tisean_package_service_accuracy():
    
    # Input Data
    print("We are using the amplitude.dat file given by Tisean website.")
    datafile = open("./Tisean_test_accuracy/amplitude.dat","r")
    data = []
    for ligne in datafile:
        data.append(float(ligne))   
    datafile.close()
    
    # Histogram
    print("Histogram :")
    gp.c('plot "< histogram -b50 ./Tisean_test_accuracy/amplitude.dat" with boxes')
    gp.p('./Tisean_test_accuracy/histogram.ps')
    
    # Correlation
    stream = os.popen('corr ./Tisean_test_accuracy/amplitude.dat')
    output = stream.read()
    print("[Corr] :")
    print("The first two lines contain: 1. the average and 2. the standard deviation of the data. The following lines are the autocorrelations (first column: delay, second column: autocorrelation). ")
    print(output)
    

def compute_lyapunov_exponent():
    
    time_series_epylepsie = edf_data('./Epileptic_Seizure_Data/chb01_03.edf')
        
    L_window = 2000
    
    lyap_exp = []
    time_debug = []
    
    # NOLDS Method 
    for i in range (0, int(len(time_series_epylepsie.sigbufs[0]) / L_window)):
        lyap_exp.append(nolds.lyap_r(time_series_epylepsie.sigbufs[0][i*L_window:(i+1)*L_window]))
        time_debug.append(i*L_window/time_series_epylepsie.sample_frequency)
        print(i*L_window/int(len(time_series_epylepsie.sigbufs[0])), "%")
    """"""
    # Result
    # naming the x axis 
    plt.xlabel('Time [sec]') 
    # naming the y axis 
    plt.ylabel('Maximum Lyapunov Exponent') 
    plt.plot(time_series_epylepsie.sigbufs[0])
    plt.show()

def compute_lyapunov_exponent_tisean():
        time_series_epylepsie = edf_data('./Epileptic_Seizure_Data/chb01_03.edf')
        
        L_window = 5000
    
        lyap_exp = []
        time_debug = []
    
        for i in range (0, int(len(time_series_epylepsie.sigbufs[0]) / L_window)):
            lyap_exp.append(time_series_epylepsie.compute_lyap_tisean_wrapper(time_series_epylepsie.sigbufs[0][i*L_window:(i+1)*L_window], 'Linux', 'lyap_r'))
            time_debug.append(i*L_window/time_series_epylepsie.sample_frequency)
        # Result
        # naming the x axis 
        plt.xlabel('Time [sec]') 
        # naming the y axis 
        plt.ylabel('Maximum Lyapunov Exponent') 
        plt.plot(time_debug, lyap_exp)
        plt.show()

test_tisean_package_service_accuracy()