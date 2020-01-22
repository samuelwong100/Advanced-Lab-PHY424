# -*- coding: utf-8 -*-
"""
File Name: data_analysis_exercise.py
Purpose: data analysis exercise
Author: Samuel Wong
"""
import numpy as np

scale=np.array([2,2,1,1,2,1,2,2,2,1,1,1]) #exact
#all units below are in cm
#sigma are the uncertainty (1 std deviation)
#the measured y
y_meas=np.array([2.65,3.55,2.00,3.50,3.60,3.30,2.60,3.35,3.70,1.40,0.08,0.25])
sigma_y_meas=np.array([0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.02,0.05])
x_meas=np.array([0.60,0.90,2.70,2.55,1.00,1.25,1.65,1.10,0.80,1.00,0.90,2.90])
sigma_x_meas=np.array([0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.05,0.02,0.02])
calibration_x=2.00 #25 unit
sigma_cal_x = 0.02
calibration_y=1.90 #25 unit
sigma_cal_y=0.02

#convert to (x,y) unit using the formula:
# unit = measured*scale*(calibration units/calibration)
x = x_meas*scale*(25/calibration_x)
y = y_meas*scale*(25/calibration_y)
