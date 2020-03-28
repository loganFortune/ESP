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
from scipy.stats import linregress

# MY LIBRARY
from edf import edf_data


def ESP_data_analysis(patient, channel):
    
    data = edf_data(patient)
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel]))
    time = time / data.sample_frequency

    # Data Visualization
    plt.title('EEG Recording : A Seizure begins at 2996 seconds.')
    plt.xlabel('seconds')
    plt.plot(time[700000:900000], data.sigbufs[channel][700000:900000])
    #plt.plot(data.sigbufs[channel])
    plt.show()

    """
    # Data Comparison. The channel provided is not used !
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
    freqs, psd = signal.welch(data.sigbufs[channel], data.sample_frequency, signal.windows.hann(8192), 8192)

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
    
    for i in range (0, int(len(data.sigbufs[channel])/window_length)):
        means.append(np.mean(data.sigbufs[channel][i*window_length:(i+1)*window_length]))
        variances.append(np.cov(data.sigbufs[channel][i*window_length:(i+1)*window_length]))
        
    plt.title('Running Variances Over The Entire Time Series')
    plt.plot(variances)
    plt.show()
    """

    """
    # Spectrogram
    powerSpectrum, frequenciesFound, time, imageAxis = plt.specgram(data.sigbufs[channel], 2046, data.sample_frequency)
   
    plt.title('Spectrogram')
    plt.xlabel('Time')
    plt.ylabel('Frequency')
    plt.show()   
    """

def visualization_phase_space(patient, channel, useowndelay = False, delay=5):

    # Write Data in a file : temp_data_viz.dat
    data = edf_data(patient)
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel]))
    time = time / data.sample_frequency

    f = open("temp_data_viz.dat", "w")
    for i in range(0, len(data.sigbufs[channel])):
        f.write(str(data.sigbufs[channel][i]))
        f.write("\n")
    f.close()
    
    # Launch Autocorrelation Measurement 
    
    os.system('autocor temp_data_viz.dat -p -o temp_data_viz.dat_co')
    
    # Result Analysis
    f = open("temp_data_viz.dat_co", "r")
    
    i_corr = 0
    timelag_autocorr = 0
    timelaginfzerofirst = False
    
    x_corr = []
    y_corr = []
    lecture = f.readlines()
    N_lignes = len(lecture)
    blankspacefound = False
    
    for j in range(0, 100):
        for i in range(0, len(lecture[j])):
            if (lecture[j][i] != " ") and (blankspacefound == False):
               i_corr = i
               blankspacefound = True
               continue
            if (lecture[j][i] == " ") and (blankspacefound == True):
                c = lecture[j][i_corr:i]
                c2 = lecture[j][i+1:]
                x_corr.append(float(c))
                y_corr.append(float(c2))
                if (float(c2) < 0 and timelaginfzerofirst == False):
                    timelaginfzerofirst = True
                    timelag_autocorr = j-1
                blankspacefound = False
                break
    f.close()

    print(" Autocor --> Time Lag :", timelag_autocorr)
    
    if (useowndelay==True):
        command1 = 'delay -d'+str(delay)+' -m2 -o delay_output.dat temp_data_viz.dat'
    else:
        command1 = 'delay -d'+str(timelag_autocorr)+' -m2 -o delay_output.dat temp_data_viz.dat'
    
    os.system(command1)
    
    f = open("delay_output.dat", "r")

    x = []
    y = []
    lecture = f.readlines()
    N_lignes = len(lecture)

    for j in range(0, int(N_lignes / 2)):
        for i in range(0, len(lecture[j])):
            if (lecture[j][i] == " "):
                c = lecture[j][0:i]
                c2 = lecture[j][i + 1:]
                x.append(float(c))
                y.append(float(c2))
                break
    f.close()
    
    x_f = [x[i] for i in range(0, len(x), 20)]
    y_f = [y[i] for i in range(0, len(y), 20)]
    plt.figure(1)
    plt.plot(x_f, y_f)
    plt.figure(2)
    plt.plot(x_corr, y_corr)
    plt.show()


def mutual_info_phase_space():
    
    channel = 0

    data = edf_data('./Epileptic_Seizure_Data/chb01_03.edf')
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel]))
    time = time / data.sample_frequency

    plt.figure()

    firstminimumfromallchannels = []

    for channel_i in range(0, int(data.numberofchannels)):

        f = open("temp_data.dat", "w")

        for i in range(0, len(data.sigbufs[channel_i])):
            f.write(str(data.sigbufs[channel_i][i]))
            f.write("\n")
        f.close()

        os.system('mutual -D50 -b100 -o time_lag.mut temp_data.dat')

        f = open("time_lag.mut", "r")
        j = 0
        mut = []
        lecture = f.readlines()
        N_lignes = len(lecture)

        for j in range(1, int(N_lignes)):
            for i in range(0, len(lecture[j])):
                if (lecture[j][i] == " "):
                    c = lecture[j][i + 1:]
                    mut.append(float(c))
                    break
        f.close()

        minimumexistence = False

        for i in range(1, len(mut) - 1):
            if (mut[i - 1] > mut[i]) and (mut[i] < mut[i + 1]):
                firstminimumfromallchannels.append(i)
                minimumexistence = True
                break

        if (minimumexistence == False):
            firstminimumfromallchannels.append(-1)  # -1 means no minimum found

        plt.plot(mut)

    print("\n Minimums for all channels : ")
    print(firstminimumfromallchannels)
    print("\n")
    plt.show()

def space_time_separation_plot():
    
    channel = 0

    data = edf_data('./Epileptic_Seizure_Data/chb01_03.edf')
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel]))
    time = time / data.sample_frequency

    f = open("temp_data_stp.dat", "w")
    for i in range(0, len(data.sigbufs[channel])):
        f.write(str(data.sigbufs[channel][i]))
        f.write("\n")
    f.close()

    gp.c('set terminal png size 1000,900')
    gp.c('set output \'./stp_graph.png\'')
    gp.c('plot "< stp temp_data_stp.dat -x1000 -%0.01 -d39 -m20 -t500"')
    
    print("\n You can see the result with the image stp_graph.png placed at the source of the project.")
    
def false_nearest_phase_space():
    
    channel = 0

    data = edf_data('./Epileptic_Seizure_Data/chb01_03.edf')
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel]))
    time = time / data.sample_frequency

    f = open("temp_data_fnn.dat", "w")
    for i in range(0, len(data.sigbufs[channel])):
        f.write(str(data.sigbufs[channel][i]))
        f.write("\n")
    f.close()

    # Be careful : time lag and embedding dimension must be inserted in the command below ! 
    os.system('false_nearest temp_data_fnn.dat -x1000 -m3 -M1,20 -d32 -t100 -o output_data_chb01_03.fnn')
    
    f = open("output_data_chb01_03.fnn", "r")
    x = []
    y = []
    lecture = f.readlines()
    N_lignes = len(lecture)
    first = False

    for j in range(0, int(N_lignes)):
        for i in range(0, len(lecture[j])):
            if lecture[j][i] == " " and first != True:
                first_i = i
                first = True
                continue
            if lecture[j][i] == " " and first == True:
                dim = lecture[j][0:first_i]
                ratio = lecture[j][first_i + 1: i]
                x.append(int(dim))
                y.append(float(ratio)*100.)
                first = False
            if i == len(lecture[j]) - 1:
                first = False
    print(x)
    print(y)
    f.close()
    plt.figure()
    plt.title('Find the Embedding dimension with the False Neighbours method \n Channel 01_03')
    plt.xlabel('Dimension')
    plt.ylabel('Fraction of false neighbors (%)')
    plt.xticks(np.arange(3, 13, step=1))
    plt.plot(x, y)
    plt.savefig('false_nearest_chb01_03.png')
    
def Lyapunov_exponent():
    init = 4000
    end = 6000
    channel = 0

    data = edf_data('./Epileptic_Seizure_Data/chb01_03.edf')
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel][init:end]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel][init:end]))
    time = time / data.sample_frequency

    f = open("temp_data_lyap.dat", "w")
    for i in range(init, end):
        f.write(str(data.sigbufs[channel][i]))
        f.write("\n")
    f.close()

    os.system('lyap_r temp_data_lyap.dat -m12 -d32 -t100 -s500 -o lyap_output.ros')

    f = open("lyap_output.ros", "r")
    lyap = []
    lecture = f.readlines()
    N_lignes = len(lecture)

    for j in range(1, int(N_lignes)):
        for i in range(0, len(lecture[j])):
            if lecture[j][i] == " ":
                c = lecture[j][i + 1:]
                lyap.append(float(c))
                break
    f.close()

    lenLyap = 300
    y = np.linspace(0, len(lyap[0:lenLyap]), len(lyap[0:lenLyap]), False, False, np.dtype('int16'))

    slope, intercept, r_value, p_value, std_err = linregress(y, lyap[0:lenLyap])
    x = np.arange(0, 500, step=1)
    
    print("Lyapunov Exponent : ", slope)
    plt.figure()
    plt.plot(lyap)
    plt.plot(x, intercept + slope*x, 'r', label='fitted line')
    plt.show()

def Lyapunov_exponent_dynamic():
    
    windowlength = 10000
    channel = 0
    
    data = edf_data('./Epileptic_Seizure_Data/chb01_03.edf')
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel]))
    time = time / data.sample_frequency    
    
    itwindow = 0
    
    LyapDyn = []
    
    while ((itwindow+windowlength) <  len(data.sigbufs[channel])):
        
        f = open("temp_data_lyap.dat", "w")
        for i in range(itwindow, itwindow+windowlength):
            f.write(str(data.sigbufs[channel][i]))
            f.write("\n")
        f.close()
        
        os.system('lyap_r temp_data_lyap.dat -m12 -d32 -t100 -s500 -o lyap_output.ros')
    
        print("\n Analysing... ", itwindow*100/len(data.sigbufs[channel]), "%")
        
        f = open("lyap_output.ros", "r")
        lyap = []
        lecture = f.readlines()
        N_lignes = len(lecture)
    
        for j in range(1, int(N_lignes)):
            for i in range(0, len(lecture[j])):
                if lecture[j][i] == " ":
                    c = lecture[j][i + 1:]
                    lyap.append(float(c))
                    break
        f.close()
        
        lenLyap = 300
        y = np.linspace(0, len(lyap[0:lenLyap]), len(lyap[0:lenLyap]), False, False, np.dtype('int16'))
    
        slope, intercept, r_value, p_value, std_err = linregress(y, lyap[0:lenLyap])
        LyapDyn.append(slope)
        itwindow=itwindow+windowlength
    
    plt.figure()
    plt.plot(LyapDyn)


patient = './Epileptic_Seizure_Data/chb01_03.edf'
channel = 19
delay = 50
#ESP_data_analysis(patient, channel)
visualization_phase_space(patient, channel)
#space_time_separation_plot()
#Lyapunov_exponent_dynamic()
#Lyapunov_exponent()
#mutual_info_phase_space()
#false_nearest_phase_space()
