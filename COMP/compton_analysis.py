# -*- coding: utf-8 -*-
"""
File Name: compton_analysis.py
Purpose: 
Author: Samuel Wong
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def load(file):
    data = np.loadtxt(file)
    channel = data[:,0]
    counts = data[:,1]
    return channel, counts

def load_adjusted(file):
    channel, counts = load(file)
    counts = np.abs(counts - noise)
    return counts

def plot(channel,counts,title):
    plt.figure()
    plt.title(title)
    plt.scatter(channel,counts)
    plt.savefig(title+".png")
    plt.show()
    
def plot_with_error_bar(channel,counts,error,title):
    plt.figure()
    plt.title(title)
    plt.errorbar(x=channel,y=counts,yerr=error,ecolor='r',fmt='o',markersize=1)
    plt.ylabel('Counts (error=$\sqrt{N}$)')
    plt.xlabel('Channels')
    plt.savefig(title+".png")
    plt.show()
    
def load_and_plot(file,title):
    channel, counts = load(file)
    plot(channel,counts,title)
    return channel, counts

def load_adjusted_and_plot(file,title):
    counts = load_adjusted(file)
    counts[0:410] = 0
    counts[600:] = 0
    uncertainty = np.sqrt(counts)
    plot_with_error_bar(channel,counts,uncertainty,title)
    return counts, uncertainty

def gaussian(x,a,b,c):
    return a*np.exp(-((x-b)**2)/(2*c**2))
    
    
#load and plot the noise counts
channel, noise = load_and_plot("data/noise.txt","noise")
#load and plot counts for different absorbers with noise subtracted
alum1, d_alum1 = load_adjusted_and_plot("data/aluminum#1.txt","aluminum#1")
                         

popt, pcov = curve_fit(gaussian, xdata=channel, ydata=alum1, sigma=d_alum1, 
                       p0=(150,540,20), maxfev=20000)                

