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
    gp.c('set terminal png size 400,300')
    gp.c('set output \'./Tisean_test_accuracy/part1_henonmap.png\'')
    gp.c('set yrange [-1:1]')
    gp.c('plot \'< henon -B0 -A2. -l100\' using 0:1 with linespoints') # a = 2.0, b=0 , 100 iterations

<<<<<<< HEAD
=======
    
>>>>>>> 96858bce892e991d82c3216ad3c134a9dffc3e4f
    # Bifurcation Diagram
    print("Bifurcation Henon Map.")
    gp.c('reset session')
    gp.c('set terminal png size 400,300')
    gp.c('set output \'./Tisean_test_accuracy/part1_henonmap_bifurcation.png\'')
<<<<<<< HEAD
    gp.c('load \'./Tisean_test_accuracy/henon.gnu\'')
=======
    gp.c('reset session')
    gp.c('load \'./Tisean_test_accuracy/henon.gnu\'')
    gp.c('set terminal png size 400,300')
    gp.c('set output \'./Tisean_test_accuracy/part1_henonmap_delay.png\'')
>>>>>>> 96858bce892e991d82c3216ad3c134a9dffc3e4f
    gp.c('reset session') # Important line because if you have plotted something before the next plot might be impacted !
    gp.c('set terminal png size 400,300')
    gp.c('set output \'./Tisean_test_accuracy/part1_henonmap_delay.png\'')
    gp.c('plot \'< henon -B0 -A2. -l50000 | delay\' with dots') # a = 2.0, b = 0, 50000 iterations
<<<<<<< HEAD
=======

>>>>>>> 96858bce892e991d82c3216ad3c134a9dffc3e4f
    
    # Histogram
    print("Histogram.")
    gp.c('set terminal png size 400,300')
    gp.c('set output \'./Tisean_test_accuracy/histogram.png\'')
    gp.c('plot "< histogram -b50 ./Tisean_test_accuracy/henon.dat" with boxes')  # 50 bins
    gp.c('set terminal png size 400,300')
    stream = os.system("henon -l10000 -B0 -A2. -o ./Tisean_test_accuracy/henon.dat") # a = 2.0, b = 0, 10000 iterations
    gp.c('reset session')
<<<<<<< HEAD
    gp.c('set terminal png size 400,300')
    gp.c('set output \'./Tisean_test_accuracy/part1_histogram.png\'')
    gp.c('plot "< histogram -b50 ./Tisean_test_accuracy/henon.dat" with boxes') # 50 bins
=======
>>>>>>> 96858bce892e991d82c3216ad3c134a9dffc3e4f
    
    # Lyapunov
    print("Lyapunov.")
    gp.c('set terminal png size 400,300')
    gp.c('set output \'./Tisean_test_accuracy/Lyapunov_henon.png\'')
    stream = os.system('henon -l10000 -B0 -A2. | lyap_r -s20 -o ./Tisean_test_accuracy/lyap_r_henon.dat') # a = 2.0, b = 0, 10stre000 iterations for the Henon map, 20 iterations for the Lyapunov Exponent
    gp.c('reset session')
<<<<<<< HEAD
    gp.c('set terminal png size 400,300')
    gp.c('set output \'./Tisean_test_accuracy/Lyapunov_henon.png\'')
    gp.c('plot \'./Tisean_test_accuracy/lyap_r_henon.dat\' with lines, x*log(2.)-8, -log(2.), \'lyap_r.dat\' with lines ')
    
=======
    gp.c('plot \'./Tisean_test_accuracy/lyap_r_henon.dat\' with lines, x*log(2.)-8, -log(2.) with lines ')
>>>>>>> 96858bce892e991d82c3216ad3c134a9dffc3e4f
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
<<<<<<< HEAD
    gp.c('reset session')
    gp.c('set terminal png size 400,300')
    gp.c('set output \'./Tisean_test_accuracy/histogram_amplitude.png\'')
    gp.c('plot "< histogram -b50 ./Tisean_test_accuracy/amplitude.dat" with boxes')
    
=======
    gp.c('set terminal png size 400,300')
    gp.c('set output \'./Tisean_test_accuracy/histogram_amplitude.png\'')
    gp.c('plot "< histogram -b50 ./Tisean_test_accuracy/amplitude.dat" with boxes')

>>>>>>> 96858bce892e991d82c3216ad3c134a9dffc3e4f
    # Correlation
    stream = os.popen('corr -D10 ./Tisean_test_accuracy/amplitude.dat') # Only 10 correlations values : we need to see it crossing zero...
    output = stream.read()
    print("[Corr] :")
    print("The first two lines contain: 1. the average and 2. the standard deviation of the data. The following lines are the autocorrelations (first column: delay, second column: autocorrelation). ")
    print(output)
    stream.close()

    """
    #Spectrum
    print("Spectrum")
    gp.c('reset session')
    gp.c('set terminal png size 400,300')
    gp.c('set output \'./Tisean_test_accuracy/part2_amplitude_spectrum.png\'')
    os.system("spectrum ./Tisean_test_accuracy/amplitude.dat -o ./Tisean_test_accuracy/amplitude.dat_sp ")
    gp.c('set logscale y')
    gp.c('plot \'amplitude.dat_sp\'')
<<<<<<< HEAD
    
=======
    """
>>>>>>> 96858bce892e991d82c3216ad3c134a9dffc3e4f
test_tisean_package_service_PartOne()
test_tisean_package_service_PartTwo()