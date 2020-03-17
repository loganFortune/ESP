#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author : LF - SM -TD

Date : 22/02/2020
"""

# Python Library
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import PyGnuplot as gp
import os

# MY LIBRARY
from edf import edf_data


def ESP_data_analysis():
    
    data = edf_data('./Epileptic_Seizure_Data/chb07_01.edf')
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[0]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0,len(data.sigbufs[0]))
    time = time/data.sample_frequency
    
    """
    # Data Visualization
    plt.title('EEG Recording : A Seizure begins at 2996 seconds.')
    plt.xlabel('seconds')
    plt.plot(time[700000:900000], data.sigbufs[0][700000:900000])
    #plt.plot(time[0:100], data.sigbufs[0][0:100])
    plt.show()
    """
    
    """
    # Data Comparison
    fig, axs = plt.subplots(6, sharex=True, sharey=True, gridspec_kw={'hspace': 0})
    fig.suptitle('Same Recording From Different Parts Of The Skull')
    axs[0].plot(data.sigbufs[0])
    axs[1].plot(data.sigbufs[18])
    axs[2].plot(data.sigbufs[16])
    axs[3].plot(data.sigbufs[17])
    axs[4].plot(data.sigbufs[12])
    axs[5].plot(data.sigbufs[22])

    # Hide x labels and tick labels for all but bottom plot.
    for ax in axs:
        ax.label_outer()
    """
    """
    # Spectrum
    freqs, psd = signal.welch(data.sigbufs[0], data.sample_frequency, signal.windows.hann(8192), 8192)

    plt.figure(figsize=(5, 4))
    plt.semilogx(freqs, psd)
    plt.title('PSD: power spectral density')
    plt.xlabel('Frequency')
    plt.ylabel('Power')
    plt.tight_layout()
    plt.show()
    """
    
    """
    means = []
    variances = []
    window_length = 8192
    
    for i in range (0, int(len(data.sigbufs[0])/window_length)):
        means.append(np.mean(data.sigbufs[0][i*window_length:(i+1)*window_length]))
        variances.append(np.cov(data.sigbufs[0][i*window_length:(i+1)*window_length]))
        
    plt.title('Running Variances Over The Entire Time Series')
    plt.plot(variances)
    plt.show()
    """
    
    # Spectrogram
    powerSpectrum, frequenciesFound, time, imageAxis = plt.specgram(data.sigbufs[0], 2046, data.sample_frequency)
   
    plt.title('Spectrogram')
    plt.xlabel('Time')
    plt.ylabel('Frequency')
    plt.show()   
    
def ESP_phase_space():
    #delay -d14 -m2 -o delay_output.dat temp_data.dat
    data = edf_data('./Epileptic_Seizure_Data/chb01_01.edf')
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[0]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0,len(data.sigbufs[0]))
    time = time/data.sample_frequency

    f = open("temp_data.dat", "w")
    for i in range (0, len(data.sigbufs[0])):
        f.write(str(data.sigbufs[0][i]))
        f.write("\n")
    f.close()
    
    os.system('delay -d14 -m2 -o delay_output.dat temp_data.dat')
    
    f = open("delay_output.dat", "r")
    j = 0
    x=[]
    y=[]
    lecture = f.readlines()
    N_lignes = len(lecture)
    
    for j in range(0, int(N_lignes/2)):
        for i in range(0, len(lecture[j])):
           if(lecture[j][i]==" "):
                c = lecture[j][0:i]
                c2 = lecture[j][i+1:]
                x.append(float(c))
                y.append(float(c2))
                break

    x_f=[x[i] for i in range(0,len(x),10)]
    y_f=[y[i] for i in range(0,len(y),10)]
    plt.figure()
    plt.plot(x_f,y_f)
    plt.show()

    


ESP_phase_space()