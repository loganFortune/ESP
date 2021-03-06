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
from esputils import esp_data_analysis, space_time_separation_plot, single_Window_Lyapunov_exponent, find_optimal_window_for_near_stationary, show_lyapunovexponent

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

    f = open("temp_data.dat", "w")
    for i in range(int(timeinitsec), int(timeendsec)):
        f.write(str(data.sigbufs[channel][i]))
        f.write("\n")
    f.close()

    # Launch Autocorrelation Measurement

    os.system('autocor temp_data.dat -p -o temp_data_viz.dat_co')

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

    for j in range(0, N_lignes):
        for i in range(0, len(lecture[j])):
            if (lecture[j][i] != " ") and not blankspacefound:
                i_corr = i
                blankspacefound = True
                continue
            if (lecture[j][i] == " ") and blankspacefound and not timelaginfzerofirst:
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

    if timelag_autocorr == 0 or timelag_autocorr < 0:
        print(" The representation is not possible because of the timelag_autocoor value.")
        return

    # Phase Space 2D-Time Lag representation

    if useowndelay:
        command1 = 'delay -d' + str(delay) + ' -m2 -o delay_output.dat temp_data.dat'
    else:
        command1 = 'delay -d' + str(timelag_autocorr) + ' -m2 -o delay_output.dat temp_data.dat'

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
    titlespec = 'Visual Inspection \n Patient/Experience:' + patient
    plt.title(titlespec)
    plt.plot(x_f, y_f)
    plt.figure(2)
    titlespec = 'Autocorrelation \n Patient/Experience:' + patient
    plt.title(titlespec)
    plt.plot(x_corr, y_corr)
    plt.axhline(y=0, color='r')
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
        f = open("temp_data.dat", "w")
        for i in range(int(timeinitsec), int(timeendsec)):
            f.write(str(data.sigbufs[channel_i][i]))
            f.write("\n")
        f.close()

        # Launch System Command
        os.system('mutual -D100 -b100 -o delay_output.dat temp_data.dat')
        os.system('autocor temp_data.dat -p -o temp_data_viz.dat_co')

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

        for i in range(1, len(mut)-1):
            if (mut[i - 1] > mut[i]) and (mut[i] < mut[i + 1]) and (minimumexistence == False):
                firstminimumfromallchannels.append(i)
                minimumexistence = True

        if not minimumexistence:
            firstminimumfromallchannels.append(-1)  # -1 means no minimum found

        titlespec = 'Mutual Information \n Patient/Experience:' + patient
        plt.title(titlespec)
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

        for j in range(0, N_lignes):
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

        if timelaginfzerofirst == False:
            timelag_autocorr.append(-1)

    print("\n Mutual First Minimums for all channels : ")
    print(firstminimumfromallchannels)
    print("\n")
    print(" Autocor - Index Crossing Zero :")
    print(timelag_autocorr)
    print("\n")

    plt.show()

    print(" For more precisions, use the \'visualization_phase_space\' function.")

# Get the right dimension value !
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

    f = open("temp_data.dat", "w")
    for i in range(int(timeinitsec), int(timeendsec)):
        f.write(str(data.sigbufs[channel][i]))
        f.write("\n")
    f.close()

    # Be careful : time lag and embedding dimension must be inserted in the command below !

    command = 'false_nearest temp_data.dat -x1000 -m3 -M1,' + str(maximumdimension) + ' -d' + str(
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

    dimensionfinal10 = 0
    dimensionfinal5 = 0
    for i in range(0, len(y)):
        if (y[i] < 10 and dimensionfinal10 == 0):
            dimensionfinal10 = i
        elif (y[i] < 1):
            dimensionfinal5 = i
            break

    print(" FNN under 10%: ", dimensionfinal10)
    print(" FNN under 1%: ", dimensionfinal5)

    print(y)
    f.close()
    plt.figure(dpi=200)
    titlespec = 'False Neighbours method \n Patient/Experience:' + patient
    plt.title(titlespec)
    plt.xlabel('Dimension')
    plt.ylabel('Fraction of false neighbors (%)')
    plt.xticks(np.arange(3, maximumdimension + 1, step=2))
    plt.xticks(fontsize=12)
    plt.plot(x, y)
    plt.axhline(y=0, color='r')
    plt.axhline(y=10, color='r')
    plt.savefig('false_nearest_chb01_03.png')


def dynamic_Lyapunov_exponent(patient, channel, dimension=5, delay=5, theilerwindow=200, windowlength=25, timeendsec=-1):

    data = edf_data(patient)
    print(" Data Analysis:")
    print(" Data length (number of points) :", len(data.sigbufs[channel]))
    print(" Data Sample Frequency (Hz) :", data.sample_frequency)
    time = range(0, len(data.sigbufs[channel]))
    time = time / data.sample_frequency

    if timeendsec == -1:
        timeendsec = len(data.sigbufs[channel]) - 1
    elif timeendsec != -1:
        timeendsec = timeendsec * data.sample_frequency

    assert(windowlength > 0)

    windowlength = windowlength*data.sample_frequency

    itwindow = 0

    LyapDyn = []

    while (itwindow + windowlength) < timeendsec:

        f = open("temp_data.dat", "w")
        for i in range(itwindow, itwindow + windowlength):
            f.write(str(data.sigbufs[channel][i]))
            f.write("\n")
        f.close()

        output_file = 'lyap_output_channel.ros'
        command = 'lyap_r temp_data.dat -m' + str(dimension) + ' -d' + str(delay) + ' -t' + str(theilerwindow) + ' -s500 -o' + str(output_file)
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

        lenLyap1 = 200
        lenLyap2 = 100
        y1 = np.linspace(0, lenLyap1, lenLyap1, False, False, np.dtype('int16'))
        y2 = np.linspace(0, lenLyap2, lenLyap2, False, False, np.dtype('int16'))

        slope1, intercept1, r_value1, p_value1, std_err1 = linregress(y1, lyap[0:lenLyap1])
        slope2, intercept2, r_value2, p_value2, std_err2 = linregress(y2, lyap[len(lyap)-lenLyap2:len(lyap)])

        intersect12point = int((intercept2-intercept1)/slope1)

        if(intersect12point > len(lyap)):
            print(" One Lyapunov exponent is not consistent. Thus, we will use the continuity of the Lyapunov exponent through time.")
            slopef= LyapDyn[-1]
        else:
            yfinal = np.linspace(0, intersect12point, intersect12point, False, False, np.dtype('int16'))
            slopef, interceptf, r_valuef, p_valuef, std_errf = linregress(yfinal, lyap[0:intersect12point])

        LyapDyn.append(slopef)
        itwindow += int(windowlength/2)

    DIR = os.path.abspath(os.path.dirname(__file__))
    output_file = os.path.join(DIR, 'lyap_data/lyap_output_channel_' + str(channel) + '.ros')

    f = open(output_file, "w")
    for i in range(0, len(LyapDyn)):
        f.write(str(LyapDyn[i]))
        f.write("\n")
    f.close()

    plt.figure()
    plt.plot(LyapDyn)
    plt.savefig('dynamic_lyap_result_channel'+str(channel)+'.png')

def show_lyapunovexponent(patient, channelstouse):

    lyapfromchannels = []

    print(" Don't forget to modify the name of the file (where Lyapunov exponents are stored) manually ! ")

    for i in range(0, len(channelstouse)):

        DIR = os.path.abspath(os.path.dirname(__file__))
        filesname = os.path.join(DIR, 'lyap_data_04/lyap_output_channel_' + str(channelstouse[i]) + '.ros')

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

def index(channelstouse, timeendsec):

    if timeendsec == -1:
        timeendsec = 3600

    lyapfromchannels = []

    print(" Don't forget to modify the name of the file (where Lyapunov exponents are stored) manually ! ")

    for i in range(0, len(channelstouse)):

        DIR = os.path.abspath(os.path.dirname(__file__))
        filesname = os.path.join(DIR, 'lyap_data_01/lyap_output_channel_' + str(channelstouse[i]) + '.ros')

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
                Tij.append(mean(diffwindow)*np.sqrt(N)/np.std(diffwindow))
        Tindexfinal.append(mean(Tij))

    """
    plt.figure(dpi=120)
    x = np.linspace(0, 3600, len(lyapfromchannels[0]), False, False, np.dtype('int16'))
    plt.plot(x, lyapfromchannels[1])
    #plt.plot(x, lyapfromchannels[0])
    plt.plot(x, lyapfromchannels[2])
    #plt.plot(x, lyapfromchannels[3])
    plt.plot(x, lyapfromchannels[4])
    plt.plot(x, lyapfromchannels[5])
    #plt.plot(x, lyapfromchannels[6])
    #plt.plot(x, lyapfromchannels[7])
    plt.show()
    """


    x = np.linspace(0, timeendsec, len(lyapfromchannels[0])-N-1, False, False, np.dtype('int16'))
    print(len(Tindexfinal))
    plt.figure(dpi=80)
    plt.plot(x, Tindexfinal)
    plt.show()


"""
    EXPERIMENTATIONS / TESTS

The following is written as below :
    
    "channel x : dimension -> a (FNN under 10%) time lag -> b (mutual;autocor) = mean(mutual;autocor)"
    
    Experiments that must be used for tests :
        chb01_01, chb01_02, chb01_03 (seizure occured), chb01_04 (seizure occured)
    
    Results for chb01_03 :

    Seizure Start Time: 2996 seconds
    Seizure End Time: 3036 seconds
    
    channel 0 : dimension->13 (10) time lag -> 34 (5;34) = 34 L.F
    channel 5 : dimension->13 (11) time lag -> 27 (26;28) = 27 L.F
    channel 7 : dimension->12 (11) time lag -> 33 (31;35) = 33 L.F
    channel 10 : dimension->13 (11) time lag -> 29 (26;32) = 29 S.M
    channel 11 : dimension->15 (12) time lag -> 36 (34;38) = 36 S.M
    channel 15 : dimension->13 (11) time lag -> 30 (28;32) = 30 T.D
    channel 16 : dimension->13 (10) time lag -> 24 (22;27) = 24 T.D
    channel 20 : dimension->13 (11) time lag -> 23 (20;27) = 23 T.D
    
    Results for chb01_04 :

    Seizure Start Time: 1467 seconds
    Seizure End Time: 1494 seconds

    ! Lyapunov Exponent computation until 2500 sec !    
    
    channel 0 : dimension->15 (10) time lag -> 26 (23;31) = 26 L.F
    channel 4 : dimension->15 (9) time lag -> 29 (29;29) = 29 L.F
    (X) channel 5 : dimension->13 (11) time lag -> 23 (18;29) = 23 L.F (not very good)
    channel 8 : dimension->15 (9) time lag -> 26 (23;30) = 26 L.F
    channel 10 : dimension->13 (11) time lag -> 29 (26;32) = 29 S.M
    channel 11 : dimension->13 (10) time lag -> 36 (34;38) = 36 S.M
    channel 15 : dimension->13 (9) time lag -> 31 (30;32) = 31 T.D
    channel 16 : dimension->13 (8) time lag -> 26 (21;32) = 26 T.D (not very good)
    channel 19 : dimension->13 (10) time lag -> 28 (26;30) = 28 T.D
    
    (X) channel 20 : dimension->X (X) time lag -> X (X;X) = X (not good at all)
    
    Results for chb01_01 :
         
    channel 0 : dimension->13 time lag -> 30 L.F
    channel 5 : dimension->13 time lag -> 30 L.F
    channel 7 : dimension->13 time lag -> 30 L.F
    channel 10 : dimension->13 time lag -> 30 L.F
    channel 11 : dimension->15 time lag -> 30 L.F
    channel 15 : dimension->13 time lag -> 30 T.D
    channel 16 : dimension->13 time lag -> 30 T.D
    channel 20 : dimension->13 time lag -> 30 T.D
    
"""

patient = './Epileptic_Seizure_Data/chb01_03.edf'
channel = 5

start_time = 2996 #seizure
end_time = 3036 # seizure
initsec = 200 # start_time
timeendsec = -1 # end_time

delay = 28
maximumdimension = 40
dimension = 13
theilerwindow = 200

windowlengthsec = 800

# channelstouse = [0, 5, 7, 10, 11, 15, 16, 20]
# channelstouse = [0, 5, 7, 10, 11, 15, 16, 20] #, 10, 11, 15, 16, 19]
# channelstouse = [0,4,8,11,15,16,19]

############## Analysis #######################################################

autocor_mutual_info_phase_space(patient, start_time, end_time)
visualization_phase_space(patient, channel, False, delay, start_time, end_time)
#false_nearest_phase_space(patient, channel, maximumdimension, delay, theilerwindow, start_time, end_time)
# correlation_dimension(patient, channel, delay, maximumdimension, theilerwindow, timeinitsec = start_time, timeendsec = end_time)
# single_Window_Lyapunov_exponent(patient, channel, dimension, delay, theilerwindow, initsec, windowlengthsec)
# space_time_separation_plot(patient, channel)
# find_optimal_window_for_near_stationary(patient,channel)

############## Prediction #####################################################

# show_lyapunovexponent(patient, channelstouse)
# esp_data_analysis(patient, channel)
# single_Window_Lyapunov_exponent(patient, channel, dimension, delay, 200, 2996, 23)
# dynamic_Lyapunov_exponent(patient, channel, dimension, delay, theilerwindow=3*delay, windowlength=23, timeendsec=timeendsec)
# index(channelstouse, timeendsec)