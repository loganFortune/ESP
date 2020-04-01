#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Author : LF - SM -TD

Date : 01/02/2020
"""

# Libraries
import pyedflib
import numpy as np
import os

class edf_data :
    
    def __init__(self, file_path):
        self.file_path = file_path
        self.file_read_pyedflib = pyedflib.EdfReader(self.file_path) # './Epileptic_Seizure_Data/name_of_the_file.edf'
        self.sigbufs = np.zeros((self.file_read_pyedflib.signals_in_file , self.file_read_pyedflib.getNSamples()[0]))
        self.sample_frequency = self.file_read_pyedflib.getSampleFrequencies()[0]
        self.numberofchannels = self.file_read_pyedflib.signals_in_file
        
        for i in np.arange(self.file_read_pyedflib.signals_in_file):
            self.sigbufs[i, :] = self.file_read_pyedflib.readSignal(i)

        self.file_read_pyedflib.close()

    def detect_seizure(self, data):
        print("next step")
        

