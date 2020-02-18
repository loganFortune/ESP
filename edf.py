#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author : LF - SM -TD

Date : 01/02/2020
"""

# Libraries
import pyedflib
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import os

class edf_data :
    
    def __init__(self, file_path):
        self.file_path = file_path
        self.file_read_pyedflib = pyedflib.EdfReader(self.file_path) # './Epileptic_Seizure_Data/'
        self.sigbufs = np.zeros((self.file_read_pyedflib.signals_in_file , self.file_read_pyedflib.getNSamples()[0]))
        self.sample_frequency = self.file_read_pyedflib.getSampleFrequencies()[0]
        
        for i in np.arange(self.file_read_pyedflib.signals_in_file):
            self.sigbufs[i, :] = self.file_read_pyedflib.readSignal(i)

        self.file_read_pyedflib.close()
                
    def compute_lyap_tisean_wrapper(self, data, your_os='linux', lyap_method='lyap_r'):
        if(your_os=='linux' or your_os=='Linux'):
            #print("You are working with a Linux operating system.")
            if(lyap_method=='lyap_r'):
                #print("You choose the algorithm of Rosenstein et al.")
                datafile =open("temp_lyap_data.dat","w")
                for i in range(0, len(data)):
                    datafile.write(str(data[i]))
                    datafile.write('\n')
                datafile.close()
                os.system('lyap_r -s500 -o temp_lyap_result.dat temp_lyap_data.dat')
                datafile = open("temp_lyap_result.dat")
                id_max = 0
                for ligne in datafile:
                    id_max += 1
                    if(id_max==501):
                        lyap_max= float(ligne[3::])
                datafile.close()
                return lyap_max
            elif(lyap_method=='lyap_k'):
                print("You choose the algorithm of Kantz.")
            else:
                print("This algorithm is not compatible with TISEAN package.")
        elif(your_os=='windows' or your_os=='Windows'):
            print(" The Windows Wrapper has not been implemented yet... ")
        else :
            print(" This operating system is not supported with the TISEAN package ")
        

