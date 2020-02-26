# -*- coding: utf-8 -*-
"""
File Name: compton_analysis.py
Purpose: 
Author: Samuel Wong
"""
import numpy as np
import matplotlib.pyplot as plt
from my_odr_fit import my_odr_fit

def load(file):
    data = np.loadtxt(file)
    channel = data[:,0]
    counts = data[:,1]
    return channel, counts

def load_adjusted(file):
    channel, counts = load(file)
    counts = counts - noise
    return counts

def plot(channel,counts,title):
    plt.figure()
    plt.title(title)
    plt.scatter(channel,counts)
    plt.savefig(title+".png")
    plt.show()
    
def load_and_plot(file,title):
    channel, counts = load(file)
    plot(channel,counts,title)
    return channel, counts

def load_adjusted_and_plot(file,title):
    counts = load_adjusted(file)
    plot(channel,counts,title)
    return counts

#load and plot the noise counts
channel, noise = load_and_plot("data/noise.txt","noise")
#load and plot counts for different absorbers with noise subtracted
alum1 = load_adjusted_and_plot("data/aluminum#1.txt","aluminum#1")
                               

