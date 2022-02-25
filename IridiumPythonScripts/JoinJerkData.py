#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Sep 12 16:26:52 2019

@author: rangapp1
"""

##############################################################################
##############################################################################
# Import Necessary Modules

from pylab import *
from glob import glob
import os
import numpy as np
import scipy as sp
from scipy.io import savemat, loadmat

##############################################################################
##############################################################################
# Load File

def getLatLon(filename):
    lat = int32(filename.split('_')[2])
    lon = int32((filename.split('_')[3]).split('.')[0])
    return lat, lon

inputs = np.sort(glob('JerkTimes_Bphi_*_*.mat'))

lat_array_peak = []
lon_array_peak = []
date_peak_a = []
month_peak_a = []
year_peak_a = []
halfwidth_peak_a = []
#halfwidth_mean_peak_a =[]
mag_peak_a = []

lat_array_valley = []
lon_array_valley = []
date_valley_a = []
month_valley_a = []
year_valley_a = []
halfwidth_valley_a = []
#halfwidth_mean_valley_a =[]
mag_valley_a = []


for filename in inputs:
    print(filename)        
    data = loadmat(filename)
    locals().update(data)
    lat, lon = getLatLon(filename)
    
    if len(date_peak) > 0:
        for i in range(len(date_peak)):
        
            lat_array_peak.append(lat)
            lon_array_peak.append(lon)
            
            date = date_peak.flatten()[i]
            month = month_peak.flatten()[i]
            year = year_peak.flatten()[i]
            hw = halfwidth_peak.flatten()[i]
            #hw_mean = halfwidth_mean_peak.flatten()[i]
            mag = mag_peak.flatten()[i]
            
            date_peak_a.append(date)
            month_peak_a.append(month) 
            year_peak_a.append(year) 
            halfwidth_peak_a.append(hw)
            #halfwidth_mean_peak_a.append(hw_mean)
            mag_peak_a.append(mag)
    else:
        print('lat', lat, 'lon', lon, 'without peak')
        
    if len(date_valley) > 0:
        for i in range(len(date_valley)):
        
            lat_array_valley.append(lat)
            lon_array_valley.append(lon)    
            
            date = date_valley.flatten()[i]
            month = month_valley.flatten()[i]
            year = year_valley.flatten()[i]
            hw = halfwidth_valley.flatten()[i]
            #hw_mean = halfwidth_mean_valleys.flatten()[i]
            mag = mag_valley.flatten()[i]
            
            date_valley_a.append(date)
            month_valley_a.append(month)
            year_valley_a.append(year)
            halfwidth_valley_a.append(hw)
            #halfwidth_mean_valley_a.append(hw_mean)
            mag_valley_a.append(mag)
            
    else:
        print('lat', lat, 'lon', lon, 'without valley')
 
#np.savetxt('MasterJerkPeak_Bphi.txt', (lat_array_peak,lon_array_peak,date_peak_a,month_peak_a,year_peak_a,halfwidth_peak_a,mag_peak_a), delimiter=' ')
#np.savetxt('MasterJerkValley_Bphi.txt', (lat_array_valley,lon_array_valley,date_valley_a,month_valley_a,year_valley_a,halfwidth_valley_a,mag_valley_a), delimiter=' ')
#savemat('MasterJerkValley_Bphi.mat' , {'lat_v':lat_array_valley,'lon_v':lon_array_valley,'date_v':date_valley_a,'month_v':month_valley_a,'year_v':year_valley_a,'hw_v':halfwidth_valley_a,'hwm_v':halfwidth_mean_valley_a,'mag_v':mag_valley_a})
#savemat('MasterJerkPeak_Bphi.mat' , {'lat_p':lat_array_peak,'lon_p':lon_array_peak,'date_p':date_peak_a,'month_p':month_peak_a,'year_p':year_peak_a,'hw_p':halfwidth_peak_a,'hwm_p':halfwidth_mean_peak_a,'mag_p':mag_peak_a})
savemat('MasterJerkValley_Bphi.mat' , {'lat_v':lat_array_valley,'lon_v':lon_array_valley,'date_v':date_valley_a,'month_v':month_valley_a,'year_v':year_valley_a,'hw_v':halfwidth_valley_a,'mag_v':mag_valley_a})
savemat('MasterJerkPeak_Bphi.mat' , {'lat_p':lat_array_peak,'lon_p':lon_array_peak,'date_p':date_peak_a,'month_p':month_peak_a,'year_p':year_peak_a,'hw_p':halfwidth_peak_a,'mag_p':mag_peak_a})
