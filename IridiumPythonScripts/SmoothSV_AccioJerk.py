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
#inputs = ['grad_timeseries_b_18_-162.mat','grad_timeseries_b_-18_-54.mat'] # Uncomment to test specific VGOs

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
    tpeaks = t[peaks]
    dpeaks = diff(tpeaks)
    if len(dpeaks) > 0:
        dpeaks = append(dpeaks,dpeaks[-1])
        pIdx = peaks[where(dpeaks > jahr)]
    else:
        pIdx = peaks
    
    vIdx = array([])
    
    if len(pIdx) > 0:
        for k in range(len(pIdx)-1):
            idx2 = where((valleys > pIdx[k]) & (valleys < pIdx[k+1]))

            if len(idx2[0]) > 1:
                idx2 = argmin(b[valleys[idx2]])

            vIdx = append(vIdx,valleys[idx2])
            
        idx3 = where(valleys < pIdx[0])
        if len(idx3[0]) > 1:
            idx3 = argmin(b[valleys[idx3]])
        vIdx = int32(append(vIdx,valleys[idx3]))
        
        idx4 = where(valleys > pIdx[-1])
        if len(idx4[0]) > 1:
            idx4 = argmin(b[valleys[idx4]])
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
    return hw, hw_ave  

##############################################################################
# Analayze data

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

        peaks = signal.argrelextrema(Br_grad_smooth,np.greater)[0]
        valleys = signal.argrelextrema(Br_grad_smooth,np.less)[0]
        
        plt.figure()
        plt.plot(t_norm, Br_grad, 'o-' , alpha=0.3)
        plt.plot(t_norm2, Br_grad_smooth,'r-')
        xlabel('Time, yr')
        ylabel(r"$B_{r}$ slope/secular variation, nT/yr")
        xticks(arange(0.,1.1,1./6.), ['2010','2011','2012','2013','2014','2015','2016'])
        plt.title(r'Secular Variation of $B_{r}$ for latitude %s & longitude %s'%(lat,long))
        
        jahr = 1./6.
        if len(peaks) > 1:
            peakIdx, valleyIdx = discard_peaks(t_norm2,Br_grad_smooth,peaks,valleys,jahr)
        else:
            peakIdx = peaks
            valleyIdx = valleys
        
        ####################
        # Get date of jerk #
        ####################
        
        paekIdx = np.sort(peakIdx)                                        
        valleyIdx = np.sort(valleyIdx)                                    
        
        date_jerk_valleys = []
        month_jerk_valleys =[]
        year_jerk_valleys =[]
        curve_valleys = []
        curve_mean_valleys = []
        Br_jerkmag_valley = []
    
        date_jerk_peaks = []
        month_jerk_peaks =[]
        year_jerk_peaks =[]
        curve_peaks = []
        curve_mean_peaks = []
        Br_jerkmag_peak = []
        
        if size(valleyIdx) != 0. :
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
                            
            #########################
            # Get halfwidth of jerk #
            #########################
            
            curve_valleys = []
            curve_mean_valleys = []
            RemIdx = []
            
            for j in range(len(valleyIdx)):
                curve,curve_mean = fwhm(t_norm2,Br_grad_smooth,valleyIdx[j])
                curve_valleys.append(curve)
                curve_mean_valleys.append(curve_mean)
                if curve > 0.:                                                 
                    RemIdx.append(j)
            
            curve_valleys = np.delete(curve_valleys,RemIdx)
            curve_mean_valleys = np.delete(curve_mean_valleys,RemIdx)
            date_jerk_valleys = np.delete(date_jerk_valleys,RemIdx)
            month_jerk_valleys = np.delete(month_jerk_valleys,RemIdx)
            year_jerk_valleys = np.delete(year_jerk_valleys,RemIdx)
            valleyIdx = np.delete(valleyIdx,RemIdx) 
            
            #########################
            # Get magnitude of jerk #
            #########################
            
            Br_jerkmag_valley = []
            
            for j in range(len(valleyIdx)):
                a = valleyIdx[j]
                B_mag_valley = Br_grad_smooth[a]
                Br_jerkmag_valley.append(B_mag_valley)

            
        if size(peakIdx) != 0. :
            
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
            # Get halfwidth of jerk #
            #########################
            
            curve_peaks = []
            curve_mean_peaks = []
            RemIdx = []
            for j in range(len(peakIdx)):
                curve, curve_mean = fwhm(t_norm2,Br_grad_smooth,peakIdx[j])
                curve_peaks.append(curve)
                curve_mean_peaks.append(curve_mean)
                if curve < 0.:                                                 
                    RemIdx.append(j)
            
            curve_peaks = np.delete(curve_peaks,RemIdx)
            curve_mean_peaks = np.delete(curve_mean_peaks,RemIdx)
            date_jerk_peaks = np.delete(date_jerk_peaks,RemIdx)
            month_jerk_peaks = np.delete(month_jerk_peaks,RemIdx)
            year_jerk_peaks = np.delete(year_jerk_peaks,RemIdx)
            peakIdx = np.delete(peakIdx,RemIdx)                                
            
            if len(peakIdx) > 0:
                peakIdx, valleyIdx = discard_peaks(t_norm2,Br_grad_smooth,peakIdx, valleyIdx,jahr)

            #########################
            # Get magnitude of jerk #
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
            plot(t_norm2[idx],Br_grad_smooth[idx],'*y',ms=15)
                
        for i in range(len(valleyIdx)):
            idx = valleyIdx[i]
            plot(t_norm2[idx],Br_grad_smooth[idx],'*m',ms=15)
            
        for j in range(len(valleyIdx)):
            valleyIdx = np.sort(valleyIdx)
            a = valleyIdx[j]
            plt.text(t_norm2[a]-0.045, Br_grad_smooth[a]-15, '%s/%s'%(month_jerk_valleys[j],year_jerk_valleys[j]), fontsize=12, fontweight='bold', bbox=dict(facecolor='m', ec='m', alpha=0.4))

        for j in range(len(peakIdx)):
            a = peakIdx[j]    
            plt.text(t_norm2[a]-0.045, Br_grad_smooth[a]+15, '%s/%s'%(month_jerk_peaks[j],year_jerk_peaks[j]), fontsize=12, fontweight='bold', bbox=dict(facecolor='y', ec='y', alpha=0.4))
        
        plt.savefig('JerkTimeSeries/Jerks_Br_lat_%s_long_%s.png'%(lat,long), dpi=400)
        savemat('JerkData/JerkTimes_Br_%s_%s.mat' %(lat,long) , {'date_peak':date_jerk_peaks, 'month_peak':month_jerk_peaks, 'year_peak':year_jerk_peaks, 'halfwidth_peak':curve_peaks, 'mag_peak':Br_jerkmag_peak, 'date_valley':date_jerk_valleys, 'month_valley':month_jerk_valleys, 'year_valley':year_jerk_valleys, 'halfwidth_valley':curve_valleys, 'mag_valley':Br_jerkmag_valley, 'halfwidth_mean_peak': curve_mean_peaks, 'halfwidth_mean_valleys': curve_mean_valleys })    
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

        peaks = signal.argrelextrema(Btheta_grad_smooth,np.greater)[0]
        valleys = signal.argrelextrema(Btheta_grad_smooth,np.less)[0]
        
        plt.figure()
        plt.plot(t_norm, Btheta_grad, 'o-' , alpha=0.3)
        plt.plot(t_norm2, Btheta_grad_smooth,'r-')
        xlabel('Time')
        ylabel(r"$B_{\theta}$ slope/secular variation, nT/yr")
        xticks(arange(0.,1.1,1./6.), ['2010','2011','2012','2013','2014','2015','2016'])
        plt.title(r'Secular Variation of $B_{\theta}$ for latitude %s & longitude %s'%(lat,long))#, x=0.5, y=-0.25)
        
        jahr = 1./6.
        if len(peaks) > 1:
            peakIdx, valleyIdx = discard_peaks(t_norm2,Btheta_grad_smooth,peaks,valleys,jahr)
        else:
            peakIdx = peaks
            valleyIdx = valleys

        ####################
        # Get date of jerk #
        ####################
        
        paekIdx = np.sort(peakIdx)                                        
        valleyIdx = np.sort(valleyIdx)                                    
    
        date_jerk_valleys = []
        month_jerk_valleys =[]
        year_jerk_valleys =[]
        curve_valleys = []
        curve_mean_valleys = []
        Btheta_jerkmag_valley = []
    
        date_jerk_peaks = []
        month_jerk_peaks =[]
        year_jerk_peaks =[]
        curve_peaks = []
        curve_mean_peaks = []
        Btheta_jerkmag_peak = []
    
        if size(valleyIdx) != 0. :
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
                
            #########################
            # Get halfwidth of jerk #
            #########################
            
            curve_valleys = []
            curve_mean_valleys = []
            RemIdx = []

            for j in range(len(valleyIdx)):
                curve, curve_mean = fwhm(t_norm2,Btheta_grad_smooth,valleyIdx[j])
                curve_valleys.append(curve)
                curve_mean_valleys.append(curve_mean)
                if curve > 0.:                                                 
                    RemIdx.append(j)
                
            curve_valleys = np.delete(curve_valleys,RemIdx)
            curve_mean_valleys = np.delete(curve_mean_valleys,RemIdx)
            date_jerk_valleys = np.delete(date_jerk_valleys,RemIdx)
            month_jerk_valleys = np.delete(month_jerk_valleys,RemIdx)
            year_jerk_valleys = np.delete(year_jerk_valleys,RemIdx)
            valleyIdx = np.delete(valleyIdx,RemIdx)                            

            #########################
            # Get magnitude of jerk #
            #########################
            
            Btheta_jerkmag_valley = []
            
            for j in range(len(valleyIdx)):
                a = valleyIdx[j]
                B_mag_valley = Btheta_grad_smooth[a]
                Btheta_jerkmag_valley.append(B_mag_valley)
            
        if size(peakIdx) != 0. :
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
            # Get halfwidth of jerk #
            #########################
            
            curve_peaks = []
            curve_mean_peaks = []
            RemIdx = []
            
            for j in range(len(peakIdx)):
                curve,curve_mean= fwhm(t_norm2,Btheta_grad_smooth,peakIdx[j])
                curve_peaks.append(curve)
                curve_mean_peaks.append(curve_mean)
                if curve < 0.:                                                
                    RemIdx.append(j)
                
            curve_peaks = np.delete(curve_peaks,RemIdx)
            curve_mean_peaks = np.delete(curve_mean_peaks,RemIdx)
            date_jerk_peaks = np.delete(date_jerk_peaks,RemIdx)
            month_jerk_peaks = np.delete(month_jerk_peaks,RemIdx)
            year_jerk_peaks = np.delete(year_jerk_peaks,RemIdx)
            peakIdx = np.delete(peakIdx,RemIdx)                                
            
            if len(peakIdx) > 0:
                peakIdx, valleyIdx = discard_peaks(t_norm2,Btheta_grad_smooth, peakIdx, valleyIdx,jahr)
            
            #########################
            # Get magnitude of jerk #
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
            plot(t_norm2[idx],Btheta_grad_smooth[idx],'*y',ms=15)
        
        for i in range(len(valleyIdx)):
            idx = valleyIdx[i]
            plot(t_norm2[idx],Btheta_grad_smooth[idx],'*m',ms=15)
    
        for j in range(len(valleyIdx)):
            valleyIdx = np.sort(valleyIdx)
            a = valleyIdx[j]
            plt.text(t_norm2[a]-0.045, Btheta_grad_smooth[a]-15, '%s/%s'%(month_jerk_valleys[j],year_jerk_valleys[j]), fontsize=12, fontweight='bold', bbox=dict(facecolor='m', ec='m', alpha=0.4))
        
        for j in range(len(peakIdx)):
            a = peakIdx[j]    
            plt.text(t_norm2[a]-0.045, Btheta_grad_smooth[a]+15, '%s/%s'%(month_jerk_peaks[j],year_jerk_peaks[j]), fontsize=12, fontweight='bold', bbox=dict(facecolor='y', ec='y', alpha=0.4))
                
        plt.savefig('JerkTimeSeries/Jerks_Btheta_lat_%s_long_%s.png'%(lat,long), dpi=400)
        savemat('JerkData/JerkTimes_Btheta_%s_%s.mat' %(lat,long) , {'date_peak':date_jerk_peaks, 'month_peak':month_jerk_peaks, 'year_peak':year_jerk_peaks, 'halfwidth_peak':curve_peaks, 'mag_peak':Btheta_jerkmag_peak, 'date_valley':date_jerk_valleys, 'month_valley':month_jerk_valleys, 'year_valley':year_jerk_valleys, 'halfwidth_valley':curve_valleys, 'mag_valley':Btheta_jerkmag_valley, 'halfwidth_mean_peak': curve_mean_peaks, 'halfwidth_mean_valleys': curve_mean_valleys })
    
        if size(valleyIdx) != 0. :
            del date_jerk_valleys, month_jerk_valleys, year_jerk_valleys, curve_valleys, Btheta_jerkmag_valley
        if size(peakIdx) != 0. :
            del date_jerk_peaks, month_jerk_peaks, year_jerk_peaks, curve_peaks, Btheta_jerkmag_peak
            
    else:
        print('btheta grad length 0 for',lat,'&', long)
##############################################################################    
    # For Bphi
    
    Bphi_grad2 = Bphi_grad[~isnan(Bphi_grad)]
    t_norm2 = t_norm[~isnan(Bphi_grad)]
    
    if len(Bphi_grad2) !=0. :
    
        Bphi_grad_smooth = signal.savgol_filter(Bphi_grad2,21,2)
        Bphi_grad_smooth = signal.savgol_filter(Bphi_grad_smooth,21,2)

        peaks = signal.argrelextrema(Bphi_grad_smooth,np.greater)[0]
        valleys = signal.argrelextrema(Bphi_grad_smooth,np.less)[0]
        
        plt.figure()
        plt.plot(t_norm, Bphi_grad, 'o-' , alpha=0.3)
        plt.plot(t_norm2, Bphi_grad_smooth,'r-')
        xlabel('Time')
        ylabel(r"$B_{\phi}$ slope/secular variation, nT/yr")
        xticks(arange(0.,1.1,1./6.), ['2010','2011','2012','2013','2014','2015','2016'])
        plt.title(r'Secular Variation of $B_{\phi}$ for latitude %s & longitude %s'%(lat,long))
        
        jahr = 1./6.
        if len(peaks) > 1:
            peakIdx, valleyIdx = discard_peaks(t_norm2,Bphi_grad_smooth,peaks,valleys,jahr)
        else:
            peakIdx = peaks
            valleyIdx = valleys

        ####################
        # Get date of jerk #
        ####################
        
        paekIdx = np.sort(peakIdx)                                        
        valleyIdx = np.sort(valleyIdx)                                    
    
        date_jerk_valleys = []
        month_jerk_valleys =[]
        year_jerk_valleys =[]
        curve_valleys = []
        curve_mean_valleys = []
        Bphi_jerkmag_valley = []
    
        date_jerk_peaks = []
        month_jerk_peaks =[]
        year_jerk_peaks =[]
        curve_peaks = []
        cureve_mean_peaks = []
        Bphi_jerkmag_peak = []
    
        if size(valleyIdx) != 0. :
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

            #########################
            # Get halfwidth of jerk #
            #########################
            
            curve_valleys = []
            curve_mean_valleys = []
            RemIdx = []
            
            for j in range(len(valleyIdx)):
                curve, curve_mean = fwhm(t_norm2,Bphi_grad_smooth,valleyIdx[j])
                curve_valleys.append(curve)
                curve_mean_valleys.append(curve_mean)
                if curve > 0.:                                                 
                    RemIdx.append(j)
                
            curve_valleys = np.delete(curve_valleys,RemIdx)
            curve_mean_valleys = np.delete(curve_mean_valleys,RemIdx)
            date_jerk_valleys = np.delete(date_jerk_valleys,RemIdx)
            month_jerk_valleys = np.delete(month_jerk_valleys,RemIdx)
            year_jerk_valleys = np.delete(year_jerk_valleys,RemIdx)
            valleyIdx = np.delete(valleyIdx,RemIdx)                            
    
            #########################
            # Get magnitude of jerk #
            #########################
            
            Bphi_jerkmag_valley = []
            
            for j in range(len(valleyIdx)):
                a = valleyIdx[j]
                B_mag_valley = Bphi_grad_smooth[a]
                Bphi_jerkmag_valley.append(B_mag_valley)    
                   
        if size(peakIdx) != 0. :
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
            # Get halfwidth of jerk #
            #########################
            
            curve_peaks = []
            curve_mean_peaks = []
            RemIdx = []
            
            for j in range(len(peakIdx)):
                curve, curve_mean = fwhm(t_norm2,Bphi_grad_smooth,peakIdx[j])
                curve_peaks.append(curve)
                curve_mean_peaks.append(curve_mean)
                if curve < 0.:                                                 
                    RemIdx.append(j)
            curve_peaks = np.delete(curve_peaks,RemIdx)
            curve_mean_peaks = np.delete(curve_mean_peaks,RemIdx)
            date_jerk_peaks = np.delete(date_jerk_peaks,RemIdx)
            month_jerk_peaks = np.delete(month_jerk_peaks,RemIdx)
            year_jerk_peaks = np.delete(year_jerk_peaks,RemIdx)
            peakIdx = np.delete(peakIdx,RemIdx)                                
            
            if len(peakIdx) > 0:
                peakIdx, valleyIdx = discard_peaks(t_norm2,Bphi_grad_smooth, peakIdx, valleyIdx,jahr)
            
            #########################
            # Get magnitude of jerk #
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
            plot(t_norm2[idx],Bphi_grad_smooth[idx],'*y',ms=15)
        
        for i in range(len(valleyIdx)):
            idx = valleyIdx[i]
            plot(t_norm2[idx],Bphi_grad_smooth[idx],'*m',ms=15)
        
        for j in range(len(valleyIdx)):
            valleyIdx = np.sort(valleyIdx)
            a = valleyIdx[j]
            plt.text(t_norm2[a]-0.045, Bphi_grad_smooth[a]-15, '%s/%s'%(month_jerk_valleys[j],year_jerk_valleys[j]), fontsize=12, fontweight='bold', bbox=dict(facecolor='m', ec='m', alpha=0.4))

    
        for j in range(len(peakIdx)):
            a = peakIdx[j]
            plt.text(t_norm2[a]-0.045, Bphi_grad_smooth[a]+15, '%s/%s'%(month_jerk_peaks[j],year_jerk_peaks[j]), fontsize=12, fontweight='bold', bbox=dict(facecolor='y', ec='y', alpha=0.4))
        
        plt.savefig('JerkTimeSeries/Jerks_Bphi_lat_%s_long_%s.png'%(lat,long), dpi=400)
        savemat('JerkData/JerkTimes_Bphi_%s_%s.mat' %(lat,long) , {'date_peak':date_jerk_peaks, 'month_peak':month_jerk_peaks, 'year_peak':year_jerk_peaks, 'halfwidth_peak':curve_peaks, 'mag_peak':Bphi_jerkmag_peak, 'date_valley':date_jerk_valleys, 'month_valley':month_jerk_valleys, 'year_valley':year_jerk_valleys, 'halfwidth_valley':curve_valleys, 'mag_valley':Bphi_jerkmag_valley, 'halfwidth_mean_peak': curve_mean_peaks, 'halfwidth_mean_valleys': curve_mean_valleys })
        
        if size(valleyIdx) != 0. :
            del date_jerk_valleys, month_jerk_valleys, year_jerk_valleys, curve_valleys, Bphi_jerkmag_valley
        if size(peakIdx) != 0. :
            del date_jerk_peaks, month_jerk_peaks, year_jerk_peaks, curve_peaks, Bphi_jerkmag_peak
            
    else:
        print('bphi grad length 0 for',lat,'&', long)
###############################################################################   