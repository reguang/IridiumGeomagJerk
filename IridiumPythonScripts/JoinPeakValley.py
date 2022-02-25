#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 26 13:52:55 2021

@author: rangapp1
"""
##############################################################################
##############################################################################
# Import Necessary Modules

from pylab import *
from glob import glob
import os
import numpy as np
from scipy.io import savemat, loadmat

##############################################################################
##############################################################################
# Load File

inputs = np.sort(glob('MasterJerk*_Btheta.mat'))
print(inputs)

lat_array_peak = []
lon_array_peak = []
date_peak_a = []
month_peak_a = []
year_peak_a = []
halfwidth_peak_a = []
#halfwidth_mean_peak_a =[]
mag_peak_a = []

for filename in inputs:
    #print(filename)        
    data = loadmat(filename)
    locals().update(data)
    
    if filename[10] == 'P':
        if len(date_p) > 0:
            for i in range(len(date_p)):
            
                lat = lat_p.flatten()[i]
                lon = lon_p.flatten()[i]
            
                lat_array_peak.append(lat)
                lon_array_peak.append(lon)
                
                date = date_p.flatten()[i]
                month = month_p.flatten()[i]
                year = year_p.flatten()[i]
                hw = hw_p.flatten()[i]
                #hw_mean = hwm_p.flatten()[i]
                mag = mag_p.flatten()[i]
                
                date_peak_a.append(date)
                month_peak_a.append(month) 
                year_peak_a.append(year) 
                halfwidth_peak_a.append(hw)
                #halfwidth_mean_peak_a.append(hw_mean)
                mag_peak_a.append(mag)
         
    if filename[10] == 'V':
        if len(date_v) > 0:
            for i in range(len(date_v)):
            
                lat = lat_v.flatten()[i]
                lon = lon_v.flatten()[i]
                
                lat_array_peak.append(lat)
                lon_array_peak.append(lon)    
                
                date = date_v.flatten()[i]
                month = month_v.flatten()[i]
                year = year_v.flatten()[i]
                hw = hw_v.flatten()[i]
                if hw > 0.:
                    print(hw, 'lat & long',lat, lon, date)
                #hw_mean = hwm_v.flatten()[i]
                mag = mag_v.flatten()[i]
                
                date_peak_a.append(date)
                month_peak_a.append(month) 
                year_peak_a.append(year) 
                halfwidth_peak_a.append(hw)
                #halfwidth_mean_peak_a.append(hw_mean)
                mag_peak_a.append(mag)
                        
#savemat('MasterJerk_Btheta.mat', {'lat_p':lat_array_peak,'lon_p':lon_array_peak,'date_p':date_peak_a,'month_p':month_peak_a,'year_p':year_peak_a,'hw_p':halfwidth_peak_a,'hwm_p':halfwidth_mean_peak_a,'mag_p':mag_peak_a})
savemat('MasterJerk_Btheta.mat', {'lat_p':lat_array_peak,'lon_p':lon_array_peak,'date_p':date_peak_a,'month_p':month_peak_a,'year_p':year_peak_a,'hw_p':halfwidth_peak_a,'mag_p':mag_peak_a})
