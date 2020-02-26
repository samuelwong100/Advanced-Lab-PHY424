# -*- coding: utf-8 -*-
"""
File Name: compton_analysis.py
Purpose: 
Author: Samuel Wong
"""
import numpy as np
import matplotlib.pyplot as plt
from my_odr_fit import my_odr_fit

def load_and_plot(file,title):
    data = np.loadtxt(file)
    channel = data[:,0]
    counts = data[:,1]
    plt.figure()
    plt.title(title)
    plt.scatter(channel,counts)
    plt.show()
    return channel, counts

load_and_plot("data/aluminum#1.txt","aluminum#1")
load_and_plot("data/noise.txt","noise")