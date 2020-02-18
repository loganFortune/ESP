# -*- coding: utf-8 -*-
"""
Author : LF - SM -TD

Date : 01/02/2020
"""

# Library Python Basics
import os
import matplotlib.pyplot as plt

#Library Chaos Libraries
import nolds.examples
import nolds
import PyGnuplot as gp

def test_nolds_package_service_accuracy():
    print("You are testing the 'nolds' package...")
    nolds.examples.plot_lyap("tent")
    nolds.examples.plot_lyap("logistic")

def test_tisean_package_service_PartOne():
    
    print("\n You are testing the Tisean Package (Part1) ... \n")
    
    # Input Data
    print("We are testing the Henon map.")
    gp.c('set yrange [-1:1]')
    gp.c('plot \'< henon -B0 -A2. -l100\' using 0:1 with linespoints') # a = 2.0, b=0 , 100 iterations
    gp.p('./Tisean_test_accuracy/part1_henonmap.ps')
    
    # Bifurcation Diagram
    print("Bifurcation Henon Map.")
    gp.c('set autoscale') # Important line because if you have plotted something before the next plot might be impacted !
    gp.c('plot \'< henon -B0 -A2. -l50000 | delay\' with dots') # a = 2.0, b = 0, 50000 iterations
    gp.p('./Tisean_test_accuracy/part1_henonmap_bifurcation.ps')
    
    # Histogram
    print("Histogram.")
    os.system("henon -l10000 -B0 -A2. -o ./Tisean_test_accuracy/henon.dat") # a = 2.0, b = 0, 10000 iterations
    gp.c('set autoscale')
    gp.c('plot "< histogram -b50 ./Tisean_test_accuracy/henon.dat" with boxes') # 50 bins
    gp.p('./Tisean_test_accuracy/histogram.ps')
    
    # Lyapunov
    print("Lyapunov.")
    os.system('henon -l10000 -B0 -A2. | lyap_r -s20 -o ./Tisean_test_accuracy/lyap_r_henon.dat') # a = 2.0, b = 0, 10000 iterations for the Henon map, 20 iterations for the Lyapunov Exponent
    gp.c('set autoscale')
    gp.c('plot \'./Tisean_test_accuracy/lyap_r_henon.dat\' with lines, x*log(2.)-8, -log(2.), \'lyap_r.dat\' with lines ')
    gp.p('./Tisean_test_accuracy/Lyapunov_henon.ps')
    
    print("Finished. See Tisean_test_accuracy directory !")
    
def test_tisean_package_service_PartTwo():
    
    print("\n You are testing the Tisean Package (Part2) ... \n")
    
    # Input Data
    print("We are using the data given by the Tisean Package : amplitude.dat")
    datafile = open("./Tisean_test_accuracy/amplitude.dat","r")
    data = []
    for ligne in datafile:
        data.append(float(ligne))   
    datafile.close()
    plt.plot(data)
    plt.show()
    
    # Histogram
    print("Histogram.")
    gp.c('set autoscale')
    gp.c('plot "< histogram -b50 ./Tisean_test_accuracy/amplitude.dat" with boxes')
    gp.p('./Tisean_test_accuracy/histogram_amplitude.ps')
    
    # Correlation
    stream = os.popen('corr -D10 ./Tisean_test_accuracy/amplitude.dat') # Only 10 correlations values
    output = stream.read()
    print("[Corr] :")
    print("The first two lines contain: 1. the average and 2. the standard deviation of the data. The following lines are the autocorrelations (first column: delay, second column: autocorrelation). ")
    print(output)
    
test_tisean_package_service_PartOne()
test_tisean_package_service_PartTwo()