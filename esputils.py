#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author : LF - SM -TD

Date : 19/04/2020
"""

# Python Library
import matplotlib.pyplot as plt
from scipy import signal
import numpy as np
import PyGnuplot as gp
import os
from scipy.stats import linregress
from statsmodels.tsa.stattools import adfuller
from statistics import mean, variance, stdev

# MY LIBRARY
from edf import edf_data

def esp_data_analysis(patient, channel):
    
    data = edf_data(patient)
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel]))
    time = time / data.sample_frequency
    
    
    # Data Visualization
    x = np.linspace(0, 3600, len(data.sigbufs[channel]), False, False, np.dtype('int16'))
    print(" Max value :", max(data.sigbufs[channel]))
    plt.xlabel('seconds')
    plt.plot(x, data.sigbufs[channel])
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
# Space-time separation plot is useful to get the right value for the Theiler window even if three times the time lag could be a valid value.
def space_time_separation_plot(patient, channel):
    data = edf_data(patient)
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel]))
    time = time / data.sample_frequency

    f = open("temp_data.dat", "w")
    for i in range(0, len(data.sigbufs[channel])):
        f.write(str(data.sigbufs[channel][i]))
        f.write("\n")
    f.close()

    gp.c('set terminal png size 1000,900')
    gp.c('set output \'./stp_graph.png\'')
    gp.c('plot "< stp temp_data.dat -x1000 -%0.01 -d39 -m20 -t500"')

    print("\n You can see the result with the image stp_graph.png placed at the source of the project.")

# Get the right dimension value with correlation dimension !
def correlation_dimension(patient, channel, delay, maximumdimension, theilerwindow, timeinitsec = 0, timeendsec = 10):
    
    data = edf_data(patient)
    
    if timeinitsec != 0:
        timeinitsec = timeinitsec * data.sample_frequency
    if timeendsec == -1:
        timeendsec = len(data.sigbufs[channel]) - 1
    elif timeendsec != -1:
        timeendsec = timeendsec * data.sample_frequency
        
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel]))
    time = time / data.sample_frequency

    print(" Data used between ", timeinitsec, "and ", timeendsec)

    f = open("temp_data.dat", "w")
    for i in range(int(timeinitsec), int(timeendsec)):
        f.write(str(data.sigbufs[channel][i]))
        f.write("\n")
    f.close()
    
    command = 'd2 temp_data.dat -d' + str(delay) + ' -M 1,' + str(maximumdimension) + ' -t'+str(theilerwindow)+' -o corrdimension'
    os.system(command)
    command2 = 'av-d2 corrdimension.d2 -o'
    os.system(command2)
    
    gp.c('set terminal png size 800,700')
    gp.c('set output \'./corrdimensionC2.png\'')
    gp.c('set logscale') # a = 2.0, b = 0, 10stre000 iterations for the Henon map, 20 iterations for the Lyapunov Exponent
    gp.c('plot \'corrdimension.c2\' with lines')
    gp.c('reset session')
    gp.c('set output \'./corrdimensiond2.png\'')
    gp.c('set logscale x') # a = 2.0, b = 0, 10stre000 iterations for the Henon map, 20 iterations for the Lyapunov Exponent
    gp.c('plot \'corrdimension.d2.av\' with lines')

        
# Get the right window length to compute the Lyapunov exponent
def find_optimal_window_for_near_stationary(patient, channel):
    
    data = edf_data(patient)
    
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel]))
    time = time / data.sample_frequency
    
    minimumlengthwindow = 3072 # You must take in account the representation space and the vectors length associated.
    maximumlengthwindow = 7680 # Same thing !
    precisionstep = 100 # not too small if you don't want to be too long !
    
    xitwindow = np.arange(minimumlengthwindow, maximumlengthwindow, step=precisionstep)
    meanVarMean = []
    for itwindow in range(minimumlengthwindow, maximumlengthwindow, precisionstep):
        meanresult = []
        i = 0
        while i+itwindow < len(data.sigbufs[channel]):
            meanresult.append(mean(data.sigbufs[channel][i:i+itwindow]))
            i = i + itwindow
        meanVarMean.append(variance(meanresult))
        
    plt.figure()
    plt.plot(xitwindow, meanVarMean)
    plt.show()
    
def single_Window_Lyapunov_exponent(patient, channel, dimension=5, delay=5, theilerwindow=100, initsec=10,
                                    windowlengthsec=23):
    data = edf_data(patient)

    end = (initsec + windowlengthsec) * data.sample_frequency
    initsec = initsec * data.sample_frequency

    assert(end < len(data.sigbufs[channel]))
    
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel][initsec:end]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel][initsec:end]))
    time = time / data.sample_frequency

    print(" Data used between ", initsec, " and ", end)

    f = open("temp_data.dat", "w")
    for i in range(int(initsec), int(end)):
        f.write(str(data.sigbufs[channel][i]))
        f.write("\n")
    f.close()

    rparam = (int(end)-int(initsec))/10
    
    print(" You're actually computing the maximum Lyapunov exponent with dimension/delay/theilerwindow:", dimension,"/", delay, "/", theilerwindow)
    command = 'lyap_r temp_data.dat -m' + str(dimension) + ' -d' + str(delay) + ' -t' + str(theilerwindow) + ' -r' + str(rparam) + ' -s500 -o lyap_output.ros'
    os.system(command)

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

    if len(lyap) == 0:
        print(" There was a problem during execution ! Be careful that the window length window is sufficiently large to operate with the dimension.")
        return

    lenLyap1 = 100
    lenLyap2 = 100
    y1 = np.linspace(0, lenLyap1, lenLyap1, False, False, np.dtype('int16'))
    y2 = np.linspace(0, lenLyap2, lenLyap2, False, False, np.dtype('int16'))
    
    slope1, intercept1, r_value1, p_value1, std_err1 = linregress(y1, lyap[0:lenLyap1])
    slope2, intercept2, r_value2, p_value2, std_err2 = linregress(y2, lyap[len(lyap)-lenLyap2:len(lyap)])
    
    intersect12point = int((intercept2-intercept1)/slope1)
    yfinal = np.linspace(0, intersect12point, intersect12point, False, False, np.dtype('int16'))
    slopef, interceptf, r_valuef, p_valuef, std_errf = linregress(yfinal, lyap[0:intersect12point])
    
    x = np.arange(0, len(lyap), step=1)

    print(" Lyapunov Exponent : ", slope1)
    plt.figure()
    plt.plot(lyap)
    plt.plot(x, interceptf + slopef * x, 'r', label='fitted line')
    plt.show()

def show_lyapunovexponent(patient, channelstouse):
    
    lyapfromchannels = []
    
    print(" Don't forget to modify the name of the file (where Lyapunov exponents are stored) manually ! ")
        
    for i in range(0, len(channelstouse)):

        DIR = os.path.abspath(os.path.dirname(__file__))
        filesname = os.path.join(DIR, 'lyap_data/lyap_data_04/lyap_output_channel_' + str(channelstouse[i]) + '.ros')

        f = open(filesname, "r")
        lyap = []
        lecture = f.readlines()
        N_lignes = len(lecture)
    
        for j in range(1, int(N_lignes)):
            for i in range(0, len(lecture[j])):
                    c = lecture[j]
                    lyap.append(float(c))
                    break
        f.close()
        lyapfromchannels.append(lyap)
    
    lyapfromchannelsMean = []
    for t in range (0, len(lyapfromchannels[0])):
        allvalue = []
        for i in range(0, len(lyapfromchannels)):
            allvalue.append(lyapfromchannels[i][t])
        lyapfromchannelsMean.append(mean(allvalue))
    
    x = np.linspace(0, 3600, len(lyapfromchannels[0]), False, False, np.dtype('int16'))
    
    plt.figure()
    plt.plot(x,lyapfromchannelsMean)
    plt.show()