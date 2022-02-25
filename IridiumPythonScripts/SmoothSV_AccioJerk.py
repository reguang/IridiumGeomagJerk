#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Aug 22 14:30:28 2019

@author: Regupathi (Regu) Angappan (he/him)
"""

##############################################################################
##############################################################################
# Import Necessary Modules

from pylab import *
from glob import glob
import numpy as np
from scipy import signal
import matplotlib.pyplot as plt
import matplotlib as mpl
from datetime import datetime
from scipy.io import savemat, loadmat
from numpy.polynomial.polynomial import polyfit
import sys
from datetime import datetime, time, timedelta
from dateutil.relativedelta import *

##############################################################################
##############################################################################
# Load File

inputs = np.sort(glob('grad_timeseries_*.mat'))

# Create directories to store output

dirname1 = 'JerkTimeSeries'
dirname2 = 'JerkData'

if not os.path.exists(dirname1):
    os.mkdir(dirname1)

if not os.path.exists(dirname2):
    os.mkdir(dirname2)

##############################################################################
# Define necessary functions

##########################
#Peak selection criteria #
##########################
def discard_peaks(t,b,peaks,valleys,jahr):
    
    pIdx = array([])
    
    tpeaks = t[peaks]
    
    for k in range(len(peaks)-1):
        dpeaks = t[peaks[k+1]] - t[peaks[k]]
    
        if dpeaks > jahr:
            pIdx = append(pIdx, peaks[k])
        else:
            tmpIdx = argmax([b[peaks[k]], b[peaks[k+1]]])
            pIdx = append(pIdx, peaks[k+tmpIdx])

    else:
        pIdx = peaks

    vIdx = array([])
    pIdx2 = array([])
    
    if len(pIdx) > 0:

        # Get new array of only peaks with valleys in between
        for k in range(len(pIdx)-1):
            idx2 = where((valleys > pIdx[k]) & (valleys < pIdx[k+1]))[0]
            if len(idx2) > 0:
                pIdx2 = append(pIdx2,[pIdx[k],pIdx[k+1]])
            else:
                tmpIdx = argmax([b[pIdx[k]], b[pIdx[k+1]]])
                pIdx2 = append(pIdx2,pIdx[k+tmpIdx])
                
        if len(pIdx2) > 0:
            pIdx2 = unique(sort(pIdx2))
            pIdx = int32(pIdx2)
        
        
        # Use remaining peaks as before to detect valleys
        for k in range(len(pIdx)-1):
            idx2 = where((valleys > pIdx[k]) & (valleys < pIdx[k+1]))[0]
            vIdx = append(vIdx,valleys[idx2])

            if len(idx2) > 1:
                tmpIdx = argmin(b[valleys[idx2]])
                idx2 = idx2[tmpIdx]
                vIdx = append(vIdx,valleys[idx2])
                vIdx = unique(sort(vIdx))

        idx3 = where(valleys < pIdx[0])[0]

        if len(idx3) > 1:
            tmpIdx = argmin(b[valleys[idx3]])
            idx3 = idx3[tmpIdx]
        vIdx = int32(append(vIdx,valleys[idx3]))

        idx4 = where(valleys > pIdx[-1])[0]

        if len(idx4) > 1:

            tmpIdx = argmin(b[valleys[idx4]])
            idx4 = idx4[tmpIdx]

        vIdx = int32(append(vIdx,valleys[idx4]))

    else:
        idx5 = where(valleys)
        idx5 = argmin(b[valleys[idx5]])
        vIdx = int32(append(vIdx,valleys[idx5]))

    vIdx = array(vIdx).flatten()
    return pIdx,vIdx

#######################
#Date time conversion #
#######################
def jerk_date(dt_seconds):
    start_date = datetime.strptime('2010-01-01 00', '%Y-%m-%d %H')
    jerk_date_time = (start_date + timedelta(seconds = dt_seconds)).isoformat(' ')
    jerk_date = jerk_date_time.split(' ')[0]
    m_dj = jerk_date.split('-')[1]
    y_dj = (jerk_date.split('_')[0])[2:4]
    return jerk_date, m_dj, y_dj

##########################
#Getting Time Normalized #
##########################

def get_normalizedyear(t):
    for i in range(len(t)):
        yt = t[i]
        t1 = datetime(2010, 1, 1, 00, 00, 00)
        t2 = datetime(yt, 1, 1, 00, 00, 00) 
        delta_t_norm = ((t2-t1).total_seconds() - min_t)/(max_t-min_t)
        year_norm.append(delta_t_norm)
    return year_norm

###############################
#Getting  Normalized Jerk Time#
###############################

def GetJerkLine(yyyy, mm):
   t1 = datetime(2010, 1, 1, 00, 00, 00) 
   t2 = datetime(yyyy, mm, 1, 00, 00, 00)
   delta_t_norm = ((t2-t1).total_seconds() - min_t)/(max_t-min_t)
   return delta_t_norm

############
#Curvature #
############
def fwhm(t,f,idx):
    t1 = t[idx]
    dt = 6./72.
    tmask = abs(t - t1) <= dt
    tFit = t[tmask]
    fFit = f[tmask]
    
    p = polyfit(tFit,fFit,2)
    
    mag = p[0]+p[1]*(tFit)+p[2]*(tFit**2)
    
    fitmax = -p[1]**2/(4*p[2]) + p[0]
    fitT = -p[1]/(2*p[2])
    
    t6months = fitT + dt
    t6prior = fitT - dt

    magt6 = p[0]+p[1]*(t6months)+p[2]*(t6months**2)
    magt6prior = p[0]+p[1]*(t6prior)+p[2]*(t6prior**2)
    
    mag_last = mag[-1]
    mag_first = mag[0]
    mag_ave = (mag_last+mag_first)/2.

    hw = (fitmax - magt6)
    hw_ave = (f[idx] - mag_ave)
    return hw 

##############################################################################
# Analayze data

t = [2010, 2011, 2012, 2013, 2014, 2015, 2016]
min_t = 43200
max_t = 185911200
year_norm = []
year_norm = get_normalizedyear(t)

for filename in inputs:
    print(filename)        
    data = loadmat(filename)
    lat = filename.split('_')[3]
    long = (filename.split('_')[4]).split('.')[0]
    
    Br_grad = data['Br_grad'].flatten()
    Btheta_grad = data['Btheta_grad'].flatten()
    Bphi_grad = data['Bphi_grad'].flatten()
    t_norm = data['time'].flatten()
    
    min_t = 43200
    max_t = 185911200

##############################################################################    
    # For Br
    
    Br_grad2 = Br_grad[~isnan(Br_grad)]
    t_norm2 = t_norm[~isnan(Br_grad)]
    
    if len(Br_grad2) !=0. :
    
        
        Br_grad_smooth = signal.savgol_filter(Br_grad2,21,2)
        Br_grad_smooth = signal.savgol_filter(Br_grad_smooth,21,2)
       
        peaks = signal.argrelextrema(Br_grad_smooth,np.greater, order =6)[0]
        valleys = signal.argrelextrema(Br_grad_smooth,np.less, order =6)[0]
        
        plt.figure(figsize=(16,9))
        plt.plot(t_norm, Br_grad, 'o-' ,ms=10 , alpha=0.3)
        plt.plot(t_norm2, Br_grad_smooth,'r-')
        xlabel('Time, yr', fontsize=50)
        ylabel(r"$\dot{B_{r}}$, nT/yr", fontsize=50)
        xticks(year_norm, ['2010','2011','2012','2013','2014','2015','2016'], fontsize=30)
        yticks(fontsize=30)
        plt.xlim(year_norm[1], year_norm[-2])
        plt.title(r'VGO at (lat.:%.1f$^{\degree}$,long.:%.1f$^{\degree}$)'%(lat_c,long_c), fontsize=50)
        
        jahr = 1./6.
        if (len(peaks) > 1) or (len(valleys) > 1):
            peakIdx, valleyIdx = discard_peaks(t_norm2,Br_grad_smooth,peaks,valleys,jahr)
        else:
            peakIdx = peaks
            valleyIdx = valleys
                
        ####################
        # Get date of jerk #
        ####################
        
        peakIdx = np.sort(peakIdx)                                        
        valleyIdx = np.sort(valleyIdx)                                    
        
        date_jerk_valleys = []
        month_jerk_valleys =[]
        year_jerk_valleys =[]
        curve_valleys = []
        Br_jerkmag_valley = []
    
        date_jerk_peaks = []
        month_jerk_peaks =[]
        year_jerk_peaks =[]
        curve_peaks = []
        Br_jerkmag_peak = []
        RemIdxv = []
        RemIdxp = []        
        if size(valleyIdx) != 0. :
            
            ###########################
            # Get halfwidth of valley #
            ###########################
            
            curve_valleys = []
            RemIdxv = []
            
            for j in range(len(valleyIdx)):
                curve = fwhm(t_norm2,Br_grad_smooth,valleyIdx[j])
                curve_valleys.append(curve)
                if curve > 0.:
                    RemIdxv.append(j)
            
            curve_valleys = np.delete(curve_valleys,RemIdxv)
            valleyIdx = np.delete(valleyIdx,RemIdxv)
            
                
        if size(peakIdx) != 0. :
            
            #########################
            # Get halfwidth of peak #
            #########################
            
            curve_peaks = []
            RemIdxp = []
            
            for j in range(len(peakIdx)):
                curve = fwhm(t_norm2,Br_grad_smooth,peakIdx[j])
                curve_peaks.append(curve)
                if curve < 0.: 
                    RemIdxp.append(j) 
            
            curve_peaks = np.delete(curve_peaks,RemIdxp)
            peakIdx = np.delete(peakIdx,RemIdxp)
            
        if (len(RemIdxv) > 0) or (len(RemIdxp) > 0):
            
            #print('!!SELECTION NEEDS CHANGING!!')
                
            if (len(peakIdx) > 0) or (len(valleyIdx) > 0):
                peakIdx, valleyIdx = discard_peaks(t_norm2,Br_grad_smooth, peakIdx, valleyIdx,jahr)
            
        ########################
        # Get timing of valley #
        ########################
        
        t_valleys_r = t_norm2[valleyIdx]
    
        t_valleys_r2 = t_valleys_r*(max_t-min_t)+min_t
        
        date_jerk_valleys = []
        month_jerk_valleys =[]
        year_jerk_valleys =[]
        
        for j in range(len(t_valleys_r2)):                             
            dt_seconds = t_valleys_r2[j]                               
            date_jerk, month_jerk, year_jerk = jerk_date(dt_seconds)                        
            date_jerk_valleys.append(date_jerk)
            month_jerk_valleys.append(month_jerk)
            year_jerk_valleys.append(year_jerk)

        
        ###########################
        # Get magnitude of valley #
        ###########################
        
        Br_jerkmag_valley = []
        
        for j in range(len(valleyIdx)):
            a = valleyIdx[j]
            B_mag_valley = Br_grad_smooth[a]
            Br_jerkmag_valley.append(B_mag_valley)
        
        ######################
        # Get timing of peak #
        ######################
        
        t_peaks_r = t_norm2[peakIdx]
    
        t_peaks_r2 = t_peaks_r*(max_t-min_t)+min_t
        
        date_jerk_peaks = []                                          
        month_jerk_peaks =[]
        year_jerk_peaks =[]
        
        for j in range(len(t_peaks_r2)):                             
            dt_seconds = t_peaks_r2[j]                               
            date_jerk, month_jerk, year_jerk = jerk_date(dt_seconds)                        
            date_jerk_peaks.append(date_jerk)
            month_jerk_peaks.append(month_jerk)
            year_jerk_peaks.append(year_jerk)                       
            
        
        #########################
        # Get magnitude of peak #
        #########################
        
        Br_jerkmag_peak = []        
        
        for j in range(len(peakIdx)):
            a = peakIdx[j]
            B_mag_peak = Br_grad_smooth[a]
            Br_jerkmag_peak.append(B_mag_peak)
                
        ######################
        # Plotting continued #
        ######################
        for i in range(len(peakIdx)):
            idx = peakIdx[i]
            if t_norm2[idx] > (1./6.) and t_norm2[idx] < (1 - (1./6.)):
                plot(t_norm2[idx],Br_grad_smooth[idx],'*y',ms=30)
            else:
                continue
            
        for i in range(len(valleyIdx)):
            idx = valleyIdx[i]
            if t_norm2[idx] > (1./6.) and t_norm2[idx] < (1 - (1./6.)):
                plot(t_norm2[idx],Br_grad_smooth[idx],'*m',ms=30)
            else:
                continue
            
        for j in range(len(valleyIdx)):
            valleyIdx = np.sort(valleyIdx)
            a = valleyIdx[j]
            if t_norm2[a] > (1./6.) and t_norm2[a] < (1 - (1./6.)):
                plt.text(t_norm2[a]-0.045, Br_grad_smooth[a]-4, '%s/%s'%(month_jerk_valleys[j],year_jerk_valleys[j]), fontsize=30, fontweight='bold', bbox=dict(facecolor='m', ec='m', alpha=0.4))
            else:
                continue
            
        for j in range(len(peakIdx)):
            a = peakIdx[j]
            if t_norm2[a] > (1./6.) and t_norm2[a] < (1 - (1./6.)):
                plt.text(t_norm2[a]-0.045, Br_grad_smooth[a]+3, '%s/%s'%(month_jerk_peaks[j],year_jerk_peaks[j]), fontsize=30, fontweight='bold', bbox=dict(facecolor='y', ec='y', alpha=0.4))
            else:
                continue
        
        plt.tight_layout()
        plt.savefig('JerkTimeSeries/Jerks_Br_lat_%s_long_%s.pdf'%(lat,long), dpi=600, bbox_inches="tight")
        savemat('JerkData/JerkTimes_Br_%s_%s.mat' %(lat,long) , {'date_peak':date_jerk_peaks, 'month_peak':month_jerk_peaks, 'year_peak':year_jerk_peaks, 'halfwidth_peak':curve_peaks, 'mag_peak':Br_jerkmag_peak, 'date_valley':date_jerk_valleys, 'month_valley':month_jerk_valleys, 'year_valley':year_jerk_valleys, 'halfwidth_valley':curve_valleys, 'mag_valley':Br_jerkmag_valley })    
        if size(valleyIdx) != 0. :
            del date_jerk_valleys, month_jerk_valleys, year_jerk_valleys, curve_valleys, Br_jerkmag_valley
        if size(peakIdx) != 0. :
            del date_jerk_peaks, month_jerk_peaks, year_jerk_peaks, curve_peaks, Br_jerkmag_peak
            
    else:
        print('br grad length 0 for',lat,'&', long)
##############################################################################    
    # For Btheta
    
    Btheta_grad2 = Btheta_grad[~isnan(Btheta_grad)]
    t_norm2 = t_norm[~isnan(Btheta_grad)]
    
    if len(Btheta_grad2) !=0. :
    
        Btheta_grad_smooth = signal.savgol_filter(Btheta_grad2,21,2)
        Btheta_grad_smooth = signal.savgol_filter(Btheta_grad_smooth,21,2)
       
        peaks = signal.argrelextrema(Btheta_grad_smooth,np.greater, order =6)[0]
        valleys = signal.argrelextrema(Btheta_grad_smooth,np.less, order =6)[0]
        
        
        plt.figure(figsize=(16,9))
        plt.plot(t_norm, Btheta_grad, 'o-' , ms=10, alpha=0.3)
        plt.plot(t_norm2, Btheta_grad_smooth,'r-')
        #plt.xlim(2011, 2015)
        xlabel('Time, yr', fontsize=50)
        ylabel(r"$\dot{B_{\theta}}$, nT/yr", fontsize=50)
        xticks(year_norm, ['2010','2011','2012','2013','2014','2015','2016'], fontsize=30)
        yticks(fontsize=30)
        #plt.xlim(1/6, 1-(1/6))
        plt.xlim(year_norm[1], year_norm[-2])
        plt.title(r'VGO at (lat.:%.1f$^{\degree}$,long.:%.1f$^{\degree}$)'%(lat_c,long_c), fontsize=50)
        
        jahr = 1./6.
        if (len(peaks) > 1) or (len(valleys) > 1):
            peakIdx, valleyIdx = discard_peaks(t_norm2,Btheta_grad_smooth,peaks,valleys,jahr)
        else:
            peakIdx = peaks
            valleyIdx = valleys
        
        ####################
        # Get date of jerk #
        ####################
        
        peakIdx = np.sort(peakIdx)                                        
        valleyIdx = np.sort(valleyIdx)                                    
    
        date_jerk_valleys = []
        month_jerk_valleys =[]
        year_jerk_valleys =[]
        curve_valleys = []
        Btheta_jerkmag_valley = []
    
        date_jerk_peaks = []
        month_jerk_peaks =[]
        year_jerk_peaks =[]
        curve_peaks = []
        Btheta_jerkmag_peak = []
        RemIdxv =[]
        RemIdxp = []
        if size(valleyIdx) != 0. :
            
            ###########################
            # Get halfwidth of valley #
            ###########################
            
            curve_valleys = []
            RemIdxv = []
            
            for j in range(len(valleyIdx)):
                curve = fwhm(t_norm2,Btheta_grad_smooth,valleyIdx[j])
                curve_valleys.append(curve)
                if curve > 0.:
                    RemIdxv.append(j)
            
            curve_valleys = np.delete(curve_valleys,RemIdxv)
            valleyIdx = np.delete(valleyIdx,RemIdxv)
            
                
        if size(peakIdx) != 0. :
            
            #########################
            # Get halfwidth of peak #
            #########################
            
            curve_peaks = []
            RemIdxp = []
            
            for j in range(len(peakIdx)):
                curve = fwhm(t_norm2,Btheta_grad_smooth,peakIdx[j])
                curve_peaks.append(curve)
                if curve < 0.: 
                    RemIdxp.append(j) 
            
            curve_peaks = np.delete(curve_peaks,RemIdxp)
            peakIdx = np.delete(peakIdx,RemIdxp)
            
        if (len(RemIdxv) > 0) or (len(RemIdxp) > 0):
            
            #print('!!SELECTION NEEDS CHANGING!!')
        
            if (len(peakIdx) > 0) or (len(valleyIdx) > 0):
                peakIdx, valleyIdx = discard_peaks(t_norm2,Btheta_grad_smooth, peakIdx, valleyIdx,jahr)
            
        ########################
        # Get timing of valley #
        ########################
        
        t_valleys_theta = t_norm2[valleyIdx]
    
        t_valleys_theta2 = t_valleys_theta*(max_t-min_t)+min_t
        
        date_jerk_valleys = []
        month_jerk_valleys =[]
        year_jerk_valleys =[]
        
        for j in range(len(t_valleys_theta2)):                             
            dt_seconds = t_valleys_theta2[j]                               
            date_jerk, month_jerk, year_jerk = jerk_date(dt_seconds)                        
            date_jerk_valleys.append(date_jerk)
            month_jerk_valleys.append(month_jerk)
            year_jerk_valleys.append(year_jerk)
                        

        ###########################
        # Get magnitude of valley #
        ###########################
        
        Btheta_jerkmag_valley = []
        
        for j in range(len(valleyIdx)):
            a = valleyIdx[j]
            B_mag_valley = Btheta_grad_smooth[a]
            Btheta_jerkmag_valley.append(B_mag_valley)
        
        ######################
        # Get timing of peak #
        ######################
        
        t_peaks_theta = t_norm2[peakIdx]
    
        t_peaks_theta2 = t_peaks_theta*(max_t-min_t)+min_t
        
        date_jerk_peaks = []                                          
        month_jerk_peaks =[]
        year_jerk_peaks =[]
        
        for j in range(len(t_peaks_theta2)):                             
            dt_seconds = t_peaks_theta2[j]                               
            date_jerk, month_jerk, year_jerk = jerk_date(dt_seconds)                        
            date_jerk_peaks.append(date_jerk)
            month_jerk_peaks.append(month_jerk)
            year_jerk_peaks.append(year_jerk)                       
            
        
        #########################
        # Get magnitude of peak #
        #########################
        
        Btheta_jerkmag_peak = []        
        
        for j in range(len(peakIdx)):
            a = peakIdx[j]
            B_mag_peak = Btheta_grad_smooth[a]
            Btheta_jerkmag_peak.append(B_mag_peak)
    
        ######################
        # Plotting continued #
        ######################
        for i in range(len(peakIdx)):
            idx = peakIdx[i]
            if t_norm2[idx] > (1./6.) and t_norm2[idx] < (1 - (1./6.)):
                plot(t_norm2[idx],Btheta_grad_smooth[idx],'*y',ms=30)
            else:
                continue
            
        for i in range(len(valleyIdx)):
            idx = valleyIdx[i]
            if t_norm2[idx] > (1./6.) and t_norm2[idx] < (1 - (1./6.)):
                plot(t_norm2[idx],Btheta_grad_smooth[idx],'*m',ms=30)
            else:
                continue
    
        for j in range(len(valleyIdx)):
            valleyIdx = np.sort(valleyIdx)
            a = valleyIdx[j]
            if t_norm2[idx] > (1./6.) and t_norm2[idx] < (1 - (1./6.)):
                plt.text(t_norm2[a]-0.045, Btheta_grad_smooth[a]-3, '%s/%s'%(month_jerk_valleys[j],year_jerk_valleys[j]), fontsize=30, fontweight='bold', bbox=dict(facecolor='m', ec='m', alpha=0.4))
            else:
                continue
        
        for j in range(len(peakIdx)):
            a = peakIdx[j]
            if t_norm2[a] > (1./6.) and t_norm2[a] < (1 - (1./6.)):
                plt.text(t_norm2[a]-0.045, Btheta_grad_smooth[a]+3, '%s/%s'%(month_jerk_peaks[j],year_jerk_peaks[j]), fontsize=30, fontweight='bold', bbox=dict(facecolor='y', ec='y', alpha=0.4))
            else:
                continue
        
        plt.tight_layout()
        plt.savefig('JerkTimeSeries/Jerks_Btheta_lat_%s_long_%s.pdf'%(lat,long), dpi=600, bbox_inches="tight")
        savemat('JerkData/JerkTimes_Btheta_%s_%s.mat' %(lat,long) , {'date_peak':date_jerk_peaks, 'month_peak':month_jerk_peaks, 'year_peak':year_jerk_peaks, 'halfwidth_peak':curve_peaks, 'mag_peak':Btheta_jerkmag_peak, 'date_valley':date_jerk_valleys, 'month_valley':month_jerk_valleys, 'year_valley':year_jerk_valleys, 'halfwidth_valley':curve_valleys, 'mag_valley':Btheta_jerkmag_valley })
        if size(valleyIdx) != 0. :
            del date_jerk_valleys, month_jerk_valleys, year_jerk_valleys, curve_valleys, Btheta_jerkmag_valley
        if size(peakIdx) != 0. :
            del date_jerk_peaks, month_jerk_peaks, year_jerk_peaks, curve_peaks, Btheta_jerkmag_peak
            
    else:
        print('btheta grad length 0 for',lat,'&', long)
##############################################################################    
    # # For Bphi
    
    Bphi_grad2 = Bphi_grad[~isnan(Bphi_grad)]
    t_norm2 = t_norm[~isnan(Bphi_grad)]
    
    if len(Bphi_grad2) !=0. :
    
        Bphi_grad_smooth = signal.savgol_filter(Bphi_grad2,21,2)
        Bphi_grad_smooth = signal.savgol_filter(Bphi_grad_smooth,21,2)
       
        peaks = signal.argrelextrema(Bphi_grad_smooth,np.greater, order =6)[0]
        valleys = signal.argrelextrema(Bphi_grad_smooth,np.less, order =6)[0]
        
        
        plt.figure(figsize=(16,9))
        plt.plot(t_norm, Bphi_grad, 'o-' , ms=10, alpha=0.3)
        plt.plot(t_norm2, Bphi_grad_smooth,'r-')
        xlabel('Time, yr', fontsize=50)
        ylabel(r"$\dot{B_{\phi}}$, nT/yr", fontsize=50)
        xticks(year_norm, ['2010','2011','2012','2013','2014','2015','2016'], fontsize=30)
        yticks(fontsize=30)
        plt.xlim(year_norm[1], year_norm[-2])
        plt.title(r'VGO at (lat.:%.1f$^{\degree}$,long.:%.1f$^{\degree}$)'%(lat_c,long_c), fontsize=50)
        
        jahr = 1./6.
        if (len(peaks) > 1) or (len(valleys) > 1):
            peakIdx, valleyIdx = discard_peaks(t_norm2,Bphi_grad_smooth,peaks,valleys,jahr)
        else:
            peakIdx = peaks
            valleyIdx = valleys
        
        ####################
        # Get date of jerk #
        ####################
        
        peakIdx = np.sort(peakIdx)                                        
        valleyIdx = np.sort(valleyIdx)                         
    
        date_jerk_valleys = []
        month_jerk_valleys =[]
        year_jerk_valleys =[]
        curve_valleys = []
        Bphi_jerkmag_valley = []
    
        date_jerk_peaks = []
        month_jerk_peaks =[]
        year_jerk_peaks =[]
        curve_peaks = []
        Bphi_jerkmag_peak = []
        RemIdxv = []
        RemIdxp = []
        if size(valleyIdx) != 0. :
            
            ###########################
            # Get halfwidth of valley #
            ###########################
            
            curve_valleys = []
            RemIdxv = []
            
            for j in range(len(valleyIdx)):
                curve = fwhm(t_norm2,Bphi_grad_smooth,valleyIdx[j])
                curve_valleys.append(curve)
                if curve > 0.:
                    RemIdxv.append(j)
            
            curve_valleys = np.delete(curve_valleys,RemIdxv)
            valleyIdx = np.delete(valleyIdx,RemIdxv)
            
                
        if size(peakIdx) != 0. :
            
            #########################
            # Get halfwidth of peak #
            #########################
            
            curve_peaks = []
            RemIdxp = []
            
            for j in range(len(peakIdx)):
                curve = fwhm(t_norm2,Bphi_grad_smooth,peakIdx[j])
                curve_peaks.append(curve)
                if curve < 0.: 
                    RemIdxp.append(j) 
            
            curve_peaks = np.delete(curve_peaks,RemIdxp)
            peakIdx = np.delete(peakIdx,RemIdxp)
            
        if (len(RemIdxv) > 0) or (len(RemIdxp) > 0):
            
            #print('!!SELECTION NEEDS CHANGING!!')
        
            if (len(peakIdx) > 0) or (len(valleyIdx) > 0):
                peakIdx, valleyIdx = discard_peaks(t_norm2,Bphi_grad_smooth, peakIdx, valleyIdx,jahr)
            
        ########################
        # Get timing of valley #
        ########################
        
        t_valleys_phi = t_norm2[valleyIdx]
    
        t_valleys_phi2 = t_valleys_phi*(max_t-min_t)+min_t
        
        date_jerk_valleys = []
        month_jerk_valleys =[]
        year_jerk_valleys =[]
        
        for j in range(len(t_valleys_phi2)):                             
            dt_seconds = t_valleys_phi2[j]                               
            date_jerk, month_jerk, year_jerk = jerk_date(dt_seconds)                        
            date_jerk_valleys.append(date_jerk)
            month_jerk_valleys.append(month_jerk)
            year_jerk_valleys.append(year_jerk)
                        

        ###########################
        # Get magnitude of valley #
        ###########################
        
        Bphi_jerkmag_valley = []
        
        for j in range(len(valleyIdx)):
            a = valleyIdx[j]
            B_mag_valley = Bphi_grad_smooth[a]
            Bphi_jerkmag_valley.append(B_mag_valley)
        
        ######################
        # Get timing of peak #
        ######################
        
        t_peaks_phi = t_norm2[peakIdx]
    
        t_peaks_phi2 = t_peaks_phi*(max_t-min_t)+min_t
        
        date_jerk_peaks = []                                          
        month_jerk_peaks =[]
        year_jerk_peaks =[]
        
        for j in range(len(t_peaks_phi2)):                             
            dt_seconds = t_peaks_phi2[j]                               
            date_jerk, month_jerk, year_jerk = jerk_date(dt_seconds)                        
            date_jerk_peaks.append(date_jerk)
            month_jerk_peaks.append(month_jerk)
            year_jerk_peaks.append(year_jerk)                       
            
        
        #########################
        # Get magnitude of peak #
        #########################
        
        Bphi_jerkmag_peak = []        
        
        for j in range(len(peakIdx)):
            a = peakIdx[j]
            B_mag_peak = Bphi_grad_smooth[a]
            Bphi_jerkmag_peak.append(B_mag_peak)
        
        ######################
        # Plotting continued #
        ######################
        for i in range(len(peakIdx)):
            idx = peakIdx[i]
            if t_norm2[idx] > (1./6.) and t_norm2[idx] < (1 - (1./6.)):
                plot(t_norm2[idx],Bphi_grad_smooth[idx],'*y',ms=30)
            else:
                continue
            
        for i in range(len(valleyIdx)):
            idx = valleyIdx[i]
            if t_norm2[idx] > (1./6.) and t_norm2[idx] < (1 - (1./6.)):
                plot(t_norm2[idx],Bphi_grad_smooth[idx],'*m',ms=30)
            else:
                continue
            
        for j in range(len(valleyIdx)):
            valleyIdx = np.sort(valleyIdx)
            a = valleyIdx[j]
            if lat=='-27' and long=='18':
                if t_norm2[idx] > (1./6.) and t_norm2[idx] < (1 - (1./6.)):
                    plt.text(t_norm2[a]-0.045, Bphi_grad_smooth[a]-0.5, '%s/%s'%(month_jerk_valleys[j],year_jerk_valleys[j]), fontsize=30, fontweight='bold', bbox=dict(facecolor='m', ec='m', alpha=0.4))
                else:
                    continue                
            else:
                if t_norm2[idx] > (1./6.) and t_norm2[idx] < (1 - (1./6.)):
                    plt.text(t_norm2[a]-0.045, Bphi_grad_smooth[a]-2, '%s/%s'%(month_jerk_valleys[j],year_jerk_valleys[j]), fontsize=30, fontweight='bold', bbox=dict(facecolor='m', ec='m', alpha=0.4))
                else:
                    continue
            
        for j in range(len(peakIdx)):
            a = peakIdx[j]
            if lat=='-27' and long=='18':
                if t_norm2[a] > (1./6.) and t_norm2[a] < (1 - (1./6.)):
                    plt.text(t_norm2[a]-0.045, Bphi_grad_smooth[a]+1, '%s/%s'%(month_jerk_peaks[j],year_jerk_peaks[j]), fontsize=30, fontweight='bold', bbox=dict(facecolor='y', ec='y', alpha=0.4))
                else:
                    continue
            else:
                if t_norm2[a] > (1./6.) and t_norm2[a] < (1 - (1./6.)):
                    plt.text(t_norm2[a]-0.045, Bphi_grad_smooth[a]+2, '%s/%s'%(month_jerk_peaks[j],year_jerk_peaks[j]), fontsize=30, fontweight='bold', bbox=dict(facecolor='y', ec='y', alpha=0.4))
                else:
                    continue
        
        plt.tight_layout()
        plt.savefig('JerkTimeSeries/Jerks_Bphi_lat_%s_long_%s.pdf'%(lat,long), dpi=600, bbox_inches="tight")
        savemat('JerkData/JerkTimes_Bphi_%s_%s.mat' %(lat,long) , {'date_peak':date_jerk_peaks, 'month_peak':month_jerk_peaks, 'year_peak':year_jerk_peaks, 'halfwidth_peak':curve_peaks, 'mag_peak':Bphi_jerkmag_peak, 'date_valley':date_jerk_valleys, 'month_valley':month_jerk_valleys, 'year_valley':year_jerk_valleys, 'halfwidth_valley':curve_valleys, 'mag_valley':Bphi_jerkmag_valley })
        if size(valleyIdx) != 0. :
            del date_jerk_valleys, month_jerk_valleys, year_jerk_valleys, curve_valleys, Bphi_jerkmag_valley
        if size(peakIdx) != 0. :
            del date_jerk_peaks, month_jerk_peaks, year_jerk_peaks, curve_peaks, Bphi_jerkmag_peak
            
    else:
        print('bphi grad length 0 for',lat,'&', long)
############################################################################### 