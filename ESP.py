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


def visualization_phase_space(patient, channel, useowndelay=False, delay=5, timeinitsec=0, timeendsec=-1):
    
    # Write Data in a file : temp_data_viz.dat

    data = edf_data(patient)
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel]))
    time = time / data.sample_frequency

    if timeinitsec != 0:
        timeinitsec = timeinitsec * data.sample_frequency
    if timeendsec == -1:
        timeendsec = len(data.sigbufs[channel]) - 1
    elif timeendsec != -1:
        timeendsec = timeendsec * data.sample_frequency

    print(" Data used between ", timeinitsec, "and ", timeendsec)

    f = open("temp_data_viz.dat", "w")
    for i in range(int(timeinitsec), int(timeendsec)):
        f.write(str(data.sigbufs[channel][i]))
        f.write("\n")
    f.close()

    # Launch Autocorrelation Measurement 

    os.system('autocor temp_data_viz.dat -p -o temp_data_viz.dat_co')

    # Autocor Results Analysis

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
            if (lecture[j][i] != " ") and not blankspacefound:
                i_corr = i
                blankspacefound = True
                continue
            if (lecture[j][i] == " ") and blankspacefound:
                c = lecture[j][i_corr:i]
                c2 = lecture[j][i + 1:]
                x_corr.append(float(c))
                y_corr.append(float(c2))
                if float(c2) < 0 and not timelaginfzerofirst:
                    timelaginfzerofirst = True
                    timelag_autocorr = j - 1
                blankspacefound = False
                break
    f.close()

    print(" Autocor --> Time Lag :", timelag_autocorr)

    # Phase Space 2D-Time Lag representation 

    if useowndelay:
        command1 = 'delay -d' + str(delay) + ' -m2 -o delay_output.dat temp_data_viz.dat'
    else:
        command1 = 'delay -d' + str(timelag_autocorr) + ' -m2 -o delay_output.dat temp_data_viz.dat'

    os.system(command1)

    f = open("delay_output.dat", "r")

    x = []
    y = []
    lecture = f.readlines()
    N_lignes = len(lecture)

    for j in range(0, int(N_lignes / 2)):
        for i in range(0, len(lecture[j])):
            if lecture[j][i] == " ":
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


def autocor_mutual_info_phase_space(patient, timeinitsec=0, timeendsec=-1):
    
    data = edf_data(patient)
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[0]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[0]))
    time = time / data.sample_frequency

    if timeinitsec != 0:
        timeinitsec = timeinitsec * data.sample_frequency
    if timeendsec == -1:
        timeendsec = len(data.sigbufs[channel]) - 1
    elif timeendsec != -1:
        timeendsec = timeendsec * data.sample_frequency

    print(" Data used between ", timeinitsec, "and ", timeendsec)
    print(" Mutual Information & Autocor Algorithms for patient ", patient, " ...")

    plt.figure()

    firstminimumfromallchannels = []
    timelag_autocorr = []

    # Autocor and Mutual Info Computation for all channels
    for channel_i in range(0, int(data.numberofchannels)):

        # Store data from the specific channel
        f = open("temp_data_viz.dat", "w")
        for i in range(int(timeinitsec), int(timeendsec)):
            f.write(str(data.sigbufs[channel_i][i]))
            f.write("\n")
        f.close()

        # Launch System Command
        os.system('mutual -D100 -b100 -o delay_output.dat temp_data_viz.dat')
        os.system('autocor temp_data_viz.dat -p -o temp_data_viz.dat_co')

        # Mutual Info Analysis
        f = open("delay_output.dat", "r")
        j = 0
        mut = []
        lecture = f.readlines()
        N_lignes = len(lecture)

        for j in range(1, int(N_lignes)):
            for i in range(0, len(lecture[j])):
                if lecture[j][i] == " ":
                    c = lecture[j][i + 1:]
                    mut.append(float(c))
                    break
        f.close()

        minimumexistence = False

        for i in range(1, len(mut) - 1):
            if (mut[i - 1] > mut[i]) and (mut[i] < mut[i + 1]) and (minimumexistence == False):
                firstminimumfromallchannels.append(i)
                minimumexistence = True

        if not minimumexistence:
            firstminimumfromallchannels.append(-1)  # -1 means no minimum found

        plt.plot(mut)

        # Autocorr Analysis
        f = open("temp_data_viz.dat_co", "r")

        i_corr = 0
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
                    c2 = lecture[j][i + 1:]
                    x_corr.append(float(c))
                    y_corr.append(float(c2))
                    if float(c2) < 0 and timelaginfzerofirst == False:
                        timelaginfzerofirst = True
                        timelag_autocorr.append(j - 1)
                    blankspacefound = False
                    break
        f.close()

    print("\n Mutual First Minimums for all channels : ")
    print(firstminimumfromallchannels)
    print("\n")
    print(" Autocor - Index Crossing Zero :")
    print(timelag_autocorr)
    print("\n")

    plt.show()

    print(" For more precisions, use the \'visualization_phase_space\' function.")


def space_time_separation_plot(patient, channel):
    data = edf_data(patient)
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


def false_nearest_phase_space(patient, channel, maximumdimension=5, delay=5, theilerwindow=100, timeinitsec=0,
                              timeendsec=-1):
    data = edf_data(patient)
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel]))
    time = time / data.sample_frequency

    if timeinitsec != 0:
        timeinitsec = timeinitsec * data.sample_frequency
    if timeendsec == -1:
        timeendsec = len(data.sigbufs[channel]) - 1
    elif timeendsec != -1:
        timeendsec = timeendsec * data.sample_frequency

    print(" Data used between ", timeinitsec, "and ", timeendsec)

    f = open("temp_data_fnn.dat", "w")
    for i in range(int(timeinitsec), int(timeendsec)):
        f.write(str(data.sigbufs[channel][i]))
        f.write("\n")
    f.close()

    # Be careful : time lag and embedding dimension must be inserted in the command below ! 
    command = 'false_nearest temp_data_fnn.dat -x1000 -m3 -M1,' + str(maximumdimension) + ' -d' + str(
        delay) + ' -t' + str(theilerwindow) + ' -o output_data_FNN.fnn'
    os.system(command)

    f = open("output_data_FNN.fnn", "r")
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
                y.append(float(ratio) * 100.)
                first = False
            if i == len(lecture[j]) - 1:
                first = False

    print(x)
    print(y)
    f.close()
    plt.figure()
    titlespec = 'Find the Embedding dimension with the False Neighbours method \n Patient/Experience:' + patient
    plt.title(titlespec)
    plt.xlabel('Dimension')
    plt.ylabel('Fraction of false neighbors (%)')
    plt.xticks(np.arange(3, maximumdimension + 1, step=1))
    plt.plot(x, y)
    plt.savefig('false_nearest_chb01_03.png')

def find_optimal_window_for_near_stationary(patient, channel):
    
    data = edf_data(patient)
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel]))
    time = time / data.sample_frequency
    
    xitwindow = np.arange(3072, 7680, step=100)
    meanVarMean = []
    for itwindow in range(3072,7680, 100):
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

    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel][initsec:end]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel][initsec:end]))
    time = time / data.sample_frequency

    print(" Data used between ", initsec, " and ", end)

    f = open("temp_data_lyap.dat", "w")
    for i in range(int(initsec), int(end)):
        f.write(str(data.sigbufs[channel][i]))
        f.write("\n")
    f.close()

    print(" You're actually computing the maximum Lyapunov exponent with dimension/delay/theilerwindow:", dimension,"/", delay, "/", theilerwindow)
    command = 'lyap_r temp_data_lyap.dat -m' + str(dimension) + ' -d' + str(delay) + ' -t' + str(theilerwindow) + ' -s500 -o lyap_output.ros'
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

    lenLyap = 300
    y = np.linspace(0, len(lyap[0:lenLyap]), len(lyap[0:lenLyap]), False, False, np.dtype('int16'))

    slope, intercept, r_value, p_value, std_err = linregress(y, lyap[0:lenLyap])
    x = np.arange(0, 500, step=1)

    print(" Lyapunov Exponent : ", slope)
    plt.figure()
    plt.plot(lyap)
    plt.plot(x, intercept + slope * x, 'r', label='fitted line')
    plt.show()


def dynamic_Lyapunov_exponent(patient, channel, dimension=5, delay=5, theilerwindow=100, windowlength=25, timeendsec=-1):

    data = edf_data(patient)
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel]))
    time = time / data.sample_frequency

    timeendsec = timeendsec * data.sample_frequency
    windowlength = windowlength*data.sample_frequency

    itwindow = 0

    LyapDyn = []

    while (itwindow + windowlength) < timeendsec:

        f = open("temp_data_lyap.dat", "w")
        for i in range(itwindow, itwindow + windowlength):
            f.write(str(data.sigbufs[channel][i]))
            f.write("\n")
        f.close()

        output_file = 'lyap_output_channel.ros'
        command = 'lyap_r temp_data_lyap.dat -m' + str(dimension) + ' -d' + str(delay) + ' -t' + str(theilerwindow) + ' -s500 -o' + str(output_file)
        os.system(command)

        print("\n Analysing... ", itwindow * 100 / len(data.sigbufs[channel]), "%")

        f = open(output_file, "r")
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
        itwindow += int(windowlength/2)

    DIR = os.path.abspath(os.path.dirname(__file__))
    output_file = os.path.join(DIR, 'lyap_data\lyap_output_channel_' + str(channel) + '.ros')

    f = open(output_file, "w")
    for i in range(0, len(LyapDyn)):
        f.write(str(LyapDyn[i]))
        f.write("\n")
    f.close()
    
    plt.figure()
    plt.plot(LyapDyn)
    plt.savefig('dynamic_lyap_result_channel'+str(channel)+'.png')

def index():
    
    channelstouse = [0, 5, 7]
    lyapfromchannels = []
    
    for i in range(0, len(channelstouse)):

        DIR = os.path.abspath(os.path.dirname(__file__))
        filesname = os.path.join(DIR, 'lyap_data\lyap_output_channel_' + str(channel) + '.ros')

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
    
    
    N = 30
       
    Tindexfinal = []
    
    for t in range (0, len(lyapfromchannels[0])-N-1):
        Tij = []
        for i in range(0, len(lyapfromchannels)):  # electrode i
            for j in range(i+1, len(lyapfromchannels)):  # electrode j
                diffwindow = []
                for k in range(t, t+N+1):
                    diffwindow.append(abs(lyapfromchannels[i][k]-lyapfromchannels[j][k]))
                Tij.append(mean(diffwindow)*np.sqrt(N)/np.std(diffwindow))  # check axis=0 or 1
        Tindexfinal.append(mean(Tij))
    
    """
    plt.figure(dpi=120)
    x = np.linspace(0, 3200, len(lyapfromchannels[0]), False, False, np.dtype('int16'))
    plt.plot(x, lyapfromchannels[1])
    plt.plot(x, lyapfromchannels[0])
    plt.plot(x, lyapfromchannels[2])
    plt.show()
    """
    plt.figure(dpi=80)
    plt.plot(Tindexfinal)
    plt.show()
        
"""
    
    EXPERIMENTATIONS / TESTS
    
channel 0 : dimension->13 time lag -> 34
channel 5 : dimension->12 time lag -> 27
channel 7 : dimension->12 time lag -> 35

"""

patient = './Epileptic_Seizure_Data/chb01_03.edf'
theilerwindow = 200
initsec = 200
windowlengthsec = 100

#  Seizure
start_time = 2996
end_time = 3036
channel = 7
delay = 28
maximumdimension = 33
dimension = 30
# autocor_mutual_info_phase_space(patient, start_time, end_time)
# visualization_phase_space(patient, channel, False, delay, start_time, end_time)
# false_nearest_phase_space(patient, channel, maximumdimension, delay, theilerwindow, start_time, end_time)
# single_Window_Lyapunov_exponent(patient, channel, dimension, delay, theilerwindow, initsec, windowlengthsec)

# Preictal
start_time = 0
end_time = 2995
channel = 5
delay = 49
maximumdimension = 33
# autocor_mutual_info_phase_space(patient, start_time, end_time)
# visualization_phase_space(patient, channel, False, delay, start_time, end_time)
# false_nearest_phase_space(patient, channel, maximumdimension, delay, theilerwindow, start_time, end_time)

#  Postictal
start_time = 3037
end_time = -1
channel = 5
delay = 34
maximumdimension = 33
dimension = 30
#autocor_mutual_info_phase_space(patient, start_time, end_time)
#visualization_phase_space(patient, channel, False, delay, start_time, end_time)
#false_nearest_phase_space(patient, channel, maximumdimension, delay, theilerwindow, start_time, end_time)
#single_Window_Lyapunov_exponent(patient, channel, dimension, delay, theilerwindow, initsec, windowlengthsec)

# autocor_mutual_info_phase_space(patient, 2996, 3036)
# visualization_phase_space(patient, channel, False, delay, 2996, 3036)
# false_nearest_phase_space(patient, channel, maximumdimension, delay, theilerwindow, 2996, 3036)
# space_time_separation_plot(patient, channel)
# single_Window_Lyapunov_exponent(patient, channel, dimension, delay, theilerwindow, initsec, windowlengthsec)
# dynamic_Lyapunov_exponent(patient, channel, dimension, delay, theilerwindow, windowlength)

# Research - time window
channel = 0 # 5 7
dimension = 13
delay = 34
end_time = 3200

#esp_data_analysis(patient, channel)
dynamic_Lyapunov_exponent(patient, channel, dimension, delay, theilerwindow=3*delay, windowlength=24, timeendsec = end_time)
index()

