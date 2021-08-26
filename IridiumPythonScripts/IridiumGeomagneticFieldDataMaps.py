#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Apr 20 22:03:51 2019

@author: Regupathi (Regu) Angappan (he/him)
"""

##############################################################################
##############################################################################
#Import Necessary Modules

from pylab import *
from libgauss import *
from glob import glob
import os
import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from datetime import datetime
from scipy.io import savemat, loadmat
import csv 
from datetime import datetime, time

##############################################################################
##############################################################################

#############
# Load File #
#############

# Lets get all the dates we need #
input1 = ('qh_SHA_raw_L00_M00.csv')

with open(input1) as csvfile:
        
    # Read CSV file #        
    readCSV = csv.reader(csvfile, delimiter=',')
    
    # Skip firs two header lines #
    next(readCSV)
    next(readCSV)
    
    date = []
    time = []
    # Assign array to relevant data #        
    for row in readCSV:
            
        # Extracting date & time
        date_d = row[0]
        date.append(date_d)
        time_d = row[1]
        time.append(time_d)

inputs = np.sort(glob('*.csv'))

for k in range(len(inputs)):
     
    data = genfromtxt(inputs[k],delimiter=',',skip_header=2)
  
    # Figure out l & m that data file corresponds to        
    split1 = inputs[k].split('.')[0]
    split2 = split1.split('_')
    
    l_series_loc = split2[3]
    l_series = int32(l_series_loc[1:3])
    
    m_series_loc = split2[4]
    m_series = int32(m_series_loc[1:3])
    
    
        
    glm_tmp = data[:,2::2] 
    hlm_tmp = data[:,3::2] 

    if m_series == 0:       # Comment this to produce files with the m=0 coefficients inlcuded
        glm_tmp[:,-1] = 0.  # Set only m=0 in phi to 0 
        #glm_tmp[...] = 0.  # Set m=0 in all components to 0

    if l_series == 0. and m_series ==0.:
        glm_tmp[:,0] = 0.

    if k == 0:
        glm = glm_tmp
        hlm = hlm_tmp
    else:
        glm = hstack([glm,glm_tmp])
        hlm = hstack([hlm,hlm_tmp])
        
    

for a in range(len(date)):
    aa = date[a]
    tt = time[a]
    tt = tt[0:2]      
    bb_res = 9. # Grid size
    theta_cut = (180/bb_res) 
    theta_cut = int32(theta_cut)
    phi_cut = (360/bb_res) 
    phi_cut = int32(phi_cut)
   
    theta = np.linspace(0.,np.pi,theta_cut)    # Polar angle from 0 to pi
    phi = np.linspace(0.,2*np.pi,phi_cut)      # Azimuthal angle from 0 to 2pi
    
    p2D,th2D = get_grid(phi,theta)
    
    # Plotting the spherical harmonic fit using the derived Gauss coeficients
    lmax = 13 # Change as necessary
    mmax = 13 # Change as necessary
    r = 1.
    lArr, mArr, idx = gen_array(lmax)
    Br = getB(lmax,glm[a,0::3],hlm[a,0::3],idx,r,p2D,th2D)
    Btheta = getB(lmax,glm[a,1::3],hlm[a,1::3],idx,r,p2D,th2D)
    Bphi = getB(lmax,glm[a,2::3],hlm[a,2::3],idx,r,p2D,th2D)
    
    theta = pi/2 - theta
    phi = phi - pi
    
    plt.figure()
    plt.pcolor(phi*(180/pi), theta*(180/pi),transpose(Br),cmap='seismic')
    xlim(-180,180)
    ylim(-90,90)
    plt.xlabel('longitude(degree)')
    plt.ylabel('latitude(degree)')
    plt.title(r"$B_{r}$ on %s" %aa)
    colorbar()
    clim(-400,400)
    plt.savefig('delta_Br_%s.png'%aa, dpi=400)
    plt.show()
    
    plt.figure()
    plt.pcolor(phi*(180/pi), theta*(180/pi),transpose(Btheta),cmap='seismic')
    xlim(-180,180)
    ylim(-90,90)
    plt.xlabel('longitude(degree)')
    plt.ylabel('latitude(degree)')
    plt.title(r"$B_{\theta}$ on %s" %aa)
    colorbar()
    clim(-400,400)
    plt.savefig('delta_Btheta_%s.png'%aa, dpi=400)
    plt.show()
    
    plt.figure()
    plt.pcolor(phi*(180/pi), theta*(180/pi),transpose(Bphi),cmap='seismic')
    xlim(-180,180)
    ylim(-90,90)
    plt.xlabel('longitude(degree)')
    plt.ylabel('latitude(degree)')
    plt.title(r"$B_{\phi}$ on %s" %aa)
    colorbar()
    clim(-400,400)
    plt.savefig('delta_Bphi_%s.png'%aa, dpi=400)
    plt.show()

    # Save new filtered binned data points
    savemat('filteredbin_%s-%s.mat'%(aa,tt) ,{'br_bin':Br,
                                              'btheta_bin':Btheta,
                                              'bphi_bin':Bphi,
                                              'theta':theta,
                                              'phi':phi})