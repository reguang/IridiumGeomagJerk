#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jul 18 17:08:45 2019

@author: Regupathi (Regu) Angappan (he/him)
"""
##############################################################################
##############################################################################
#Import Necessary Modules

import matplotlib
from pylab import *
from libgauss import *
from glob import glob
import os
import numpy as np
from datetime import datetime
from scipy.io import savemat, loadmat
from numpy.polynomial.polynomial import polyfit
import sys
from datetime import datetime, time, timedelta
from dateutil.relativedelta import *

##############################################################################
##############################################################################
#############
# Load File #
#############

inputs = np.sort(glob('filteredbin*.mat'))

bb_res = 9
lon_ideal = np.arange(-180, 180+bb_res, bb_res)
lat_ideal = np.arange(-90, 90+bb_res, bb_res)
sc = 1.5

year_norm = [-0.00023242300987797793, 0.1694363742010459, 0.3391051714119698, 0.5092388146426496, 0.6789076118535735, 0.8485764090644974, 1.0182452062754213]

# Let's make a directory to store all the outputs! 
dirname = 'TimeSeriesPlots'

if not os.path.exists(dirname):
    os.mkdir(dirname)

# Time to analyze the data! 
for long in range(len(lon_ideal)-1):
    
    for lat in range(len(lat_ideal)-1):
        
        long_oc = lon_ideal[long]    #Longitude of concern - that is being plotted/analyzed
        lat_oc = lat_ideal[lat]      #Latitude of concern - that is being plotted/analyzed
        
        Br_plot = []
        Btheta_plot = []
        Bphi_plot =[]
        date_v = []
        t_s_plot = []
        date1_a =[]
        
        ######################################################################
        for filename in inputs:
            
            dat = loadmat(filename)
            title = filename.split('.')[0]
            date = title.split('_')[1]
            
            year = date.split('-')[0]
            month = date.split('-')[1]
            day = date.split('-')[2]
            time = date.split('-')[3]
            
            YY = year[2:4]
            MM = month
            DD = day
            HH = time
            
            ##################################################################
            #######################
            #Date time conversion #
            #######################
            dt = '20'+YY+'-'+MM+'-'+DD+' '+HH
            
            def date_diff_in_Seconds(dt2, dt1):
                timedelta_r = dt2 - dt1 
                return timedelta_r.days * 24 * 3600 + timedelta_r.seconds
            
            date1 = datetime.strptime(dt, '%Y-%m-%d %H')
            date1_a.append(date1)
            
            # Reference date
            date2 = datetime.strptime('2010-01-01 00', '%Y-%m-%d %H')
            
            t_s = date_diff_in_Seconds(date1, date2)
            t_s_plot.append(t_s)
            
            ##################################################################
            Br = dat['br_bin']
            Btheta = dat['btheta_bin']
            Bphi = dat['bphi_bin']

            Br1 = Br[long, lat]
            Br_plot.append(Br1)
            Btheta1 = Btheta[long, lat]
            Btheta_plot.append(Btheta1)
            Bphi1 = Bphi[long, lat]
            Bphi_plot.append(Bphi1)
            
            Br_plot_a = np.asarray(Br_plot)
            Btheta_plot_a = np.asarray(Btheta_plot)
            Bphi_plot_a = np.asarray(Bphi_plot)
            
            ##################################################################
        
        # Plot the resulting overall time series
        t_s_plot_a = np.asarray(t_s_plot)
        t_s_plot_norm = (t_s_plot_a - min(t_s_plot_a))/(max(t_s_plot_a) - min(t_s_plot_a))
        
        plt.figure(1)
        #plt.plot(t_s_plot_norm, Br_plot_a,'o-',label='data')
        plt.plot(t_s_plot_norm, Br_plot_a,'o-')
        #legend()
        xlabel('Time')
        ylabel(r"$B_{r}$, nT")
        xticks(year_norm, ['2010','2011','2012','2013','2014','2015','2016'])
        plt.title(r"$B_{r}$ for latitude %d & longitude %d"%(lat_oc,long_oc))
        
        
        plt.figure(2)
        #plt.plot(t_s_plot_norm, Btheta_plot_a,'o-',label='data')
        plt.plot(t_s_plot_norm, Btheta_plot_a,'o-')
        #legend()
        xlabel('Time')
        ylabel(r"$B_{\theta}$, nT")
#       xticks(arange(0.,1.1,0.16666), ['2010','2011','2012','2013','2014','2015','2016'])
        xticks(year_norm, ['2010','2011','2012','2013','2014','2015','2016'])
        plt.title(r"$B_{\theta}$ for latitude %d & longitude %d"%(lat_oc,long_oc))
        
        
        plt.figure(3)
        #plt.plot(t_s_plot_norm, Bphi_plot_a,'o-',label='data')
        plt.plot(t_s_plot_norm, Bphi_plot_a,'o-')
        #legend()
        xlabel('Time')
        ylabel(r"$B_{\phi}$, nT")
        xticks(year_norm, ['2010','2011','2012','2013','2014','2015','2016'])
        plt.title(r"$B_{\phi}$ for latitude %d & longitude %d"%(lat_oc,long_oc))
        
        # Calculating the polynomial fit for the 18 month windows, sliding by 1 month
        date_end_d = date1_a[-1]                          # Last date in data
        date_end_f = date_end_d+relativedelta(months=-3)  # Last date in fit
        
        def nearest(items, pivot):
                return min(items, key=lambda x: abs(x - pivot))
        
        move_d = date1_a[0]
        start_previous = 0
        
        Br_fit_slope = []
        Btheta_fit_slope = []
        Bphi_fit_slope = []
        t_slope_a = []
        t_start_a =[]
        t_end_a =[]
        
        Br_c_0 =[]
        Br_c_2 =[]
        Br_c_3 =[]
        Btheta_c_0 =[]
        Btheta_c_2 =[]
        Btheta_c_3 =[]
        Bphi_c_0 =[]
        Bphi_c_2 =[]
        Bphi_c_3 =[]
        
        while(1):
            
            start_d = move_d
            
            start_d = nearest(date1_a, start_d)
            
            if start_d == start_previous:
                start_d = date1_a[idx_s+1]
            
            end_d = start_d+relativedelta(months=+18) 
            
            last_date = nearest(date1_a, end_d)
            
            # Find index of first date & last date
            idx_s = date1_a.index(start_d)
            idx_e = date1_a.index(last_date)
            
            Br_fit = Br_plot[idx_s:idx_e+1]
            Btheta_fit = Btheta_plot[idx_s:idx_e+1]
            Bphi_fit = Bphi_plot[idx_s:idx_e+1]
            date_fit = t_s_plot_norm[idx_s:idx_e+1]     # Based on normalized array
            date_fit_a = np.asarray(date_fit)
            date_fit_norm = (date_fit_a - min(date_fit_a))/(max(date_fit_a) - min(date_fit_a))
            
            t_start = t_s_plot_norm[idx_s]
            t_start_a.append(t_start)
            t_end = t_s_plot_norm[idx_e]
            t_end_a.append(t_end)
            
            start_previous = start_d
            move_d = start_d+relativedelta(months=+1)
            
            # Fitting with 3rd degree polynomial
            
            if (len(date_fit)<10): # Ensure that there is ample data to do a fit (tried 5 and it works just as well)
                continue
            
            n_p = 3
            
            Br_polyfit_c = polyfit(date_fit_norm,Br_fit,n_p)
            Br_polyfit = Br_polyfit_c[0]+Br_polyfit_c[1]*(date_fit_norm)+Br_polyfit_c[2]*(date_fit_norm**2)+Br_polyfit_c[3]*(date_fit_norm**3)
            Br_fit_slope.append(Br_polyfit_c[1])
            Br_c_0.append(Br_polyfit_c[0])
            Br_c_2.append(Br_polyfit_c[2])
            Br_c_3.append(Br_polyfit_c[3])
            
            Btheta_polyfit_c = polyfit(date_fit_norm,Btheta_fit,n_p)
            Btheta_polyfit = Btheta_polyfit_c[0]+Btheta_polyfit_c[1]*(date_fit_norm)+Btheta_polyfit_c[2]*(date_fit_norm**2)+Btheta_polyfit_c[3]*(date_fit_norm**3)
            Btheta_fit_slope.append(Btheta_polyfit_c[1])
            Btheta_c_0.append(Btheta_polyfit_c[0])
            Btheta_c_2.append(Btheta_polyfit_c[2])
            Btheta_c_3.append(Btheta_polyfit_c[3])
            
            Bphi_polyfit_c = polyfit(date_fit_norm,Bphi_fit,n_p)
            Bphi_polyfit = Bphi_polyfit_c[0]+Bphi_polyfit_c[1]*(date_fit_norm)+Bphi_polyfit_c[2]*(date_fit_norm**2)+Bphi_polyfit_c[3]*(date_fit_norm**3)
            Bphi_fit_slope.append(Bphi_polyfit_c[1])
            Bphi_c_0.append(Bphi_polyfit_c[0])
            Bphi_c_2.append(Bphi_polyfit_c[2])
            Bphi_c_3.append(Bphi_polyfit_c[3])
            
            # Mid point fo the fit (Check this please)
            t_slope = date_fit[int(len(date_fit)/2)]
            t_slope_a.append(t_slope)
            
            # Plot the fits
            figure(1)
            plt.plot(date_fit, Br_polyfit, '-')
            plt.savefig('TimeSeriesPlots/FilteredBrDataTimeSeries_fit_lat_%d_long_%d.png'%(lat_oc,long_oc), dpi=400)
            
            figure(2)
            plt.plot(date_fit, Btheta_polyfit, '-')
            plt.savefig('TimeSeriesPlots/FilteredBthetaDataTimeSeries_fit_lat_%d_long_%d.png'%(lat_oc,long_oc), dpi=400)
           
            figure(3)
            plt.plot(date_fit, Bphi_polyfit, '-')
            plt.savefig('TimeSeriesPlots/FilteredBphiDataTimeSeries_fit_lat_%d_long_%d.png'%(lat_oc,long_oc), dpi=400)
            
            if (move_d > date_end_f):
                break
            
        # Plot the vatriation of the linear term
        
        Br_fit_slope = np.asarray(Br_fit_slope)
        Btheta_fit_slope = np.asarray(Btheta_fit_slope)
        Bphi_fit_slope = np.asarray(Bphi_fit_slope)
        t_slope_a = np.asarray(t_slope_a)
        
        Br_c_0 = np.asarray(Br_c_0)
        Br_c_2 = np.asarray(Br_c_2)
        Br_c_3 = np.asarray(Br_c_3)
        Btheta_c_0 = np.asarray(Btheta_c_0)
        Btheta_c_2 = np.asarray(Btheta_c_2)
        Btheta_c_3 = np.asarray(Btheta_c_3)
        Bphi_c_0 = np.asarray(Bphi_c_0)
        Bphi_c_2 = np.asarray(Bphi_c_2)
        Bphi_c_3 = np.asarray(Bphi_c_3)
        
        # Reduce noise by finding the average of the quartiles       
        
        Br_mean_seg_a = []
        Br_slope_mean_seg_a = []
        Btheta_mean_seg_a = []
        Btheta_slope_mean_seg_a = []
        Bphi_mean_seg_a = []
        Bphi_slope_mean_seg_a = []
        
        for i in range(len(t_slope_a)):
            
            t_inq = t_slope_a[i]
            
            idx = where( (t_start_a < t_inq) & (t_end_a > t_inq) )[0]
            
            Br_fit_t = []
            Br_slope_t = []
            Btheta_fit_t = []
            Btheta_slope_t = []
            Bphi_fit_t = []
            Bphi_slope_t = []
            
            for k in idx:#t_up_lim, t_low_lim):
                
                Br_c_0_t = Br_c_0[k]
                Br_c_1_t = Br_fit_slope[k]
                Br_c_2_t = Br_c_2[k]
                Br_c_3_t = Br_c_3[k]
                Br_t = Br_c_0_t+Br_c_1_t*(t_inq)+Br_c_2_t*(t_inq**2)+Br_c_3_t*(t_inq**3)
                
                Br_slope_t.append(Br_c_1_t)
                Br_fit_t.append(Br_t)
                
                Btheta_c_0_t = Btheta_c_0[k]
                Btheta_c_1_t = Btheta_fit_slope[k]
                Btheta_c_2_t = Btheta_c_2[k]
                Btheta_c_3_t = Btheta_c_3[k]
                Btheta_t = Btheta_c_0_t+Btheta_c_1_t*(t_inq)+Btheta_c_2_t*(t_inq**2)+Btheta_c_3_t*(t_inq**3)
                
                Btheta_slope_t.append(Btheta_c_1_t)
                Btheta_fit_t.append(Btheta_t)
                
                Bphi_c_0_t = Bphi_c_0[k]
                Bphi_c_1_t = Bphi_fit_slope[k]
                Bphi_c_2_t = Bphi_c_2[k]
                Bphi_c_3_t = Bphi_c_3[k]
                Bphi_t = Bphi_c_0_t+Bphi_c_1_t*(t_inq)+Bphi_c_2_t*(t_inq**2)+Bphi_c_3_t*(t_inq**3)
                
                Bphi_slope_t.append(Bphi_c_1_t)
                Bphi_fit_t.append(Bphi_t)
                                                                      
            try:
                Br_q1 = np.percentile(Br_fit_t, 25)
                Br_q3 = np.percentile(Br_fit_t, 75)
                Br_fit_t = np.sort(Br_fit_t).tolist()
                Br_q1 = nearest(Br_fit_t, Br_q1)
                Br_q3 = nearest(Br_fit_t, Br_q3)
            except:
               break
            
            idx_q1_r = Br_fit_t.index(Br_q1)
            idx_q3_r = Br_fit_t.index(Br_q3)
            
            Br_mean_seg = mean(Br_fit_t[idx_q1_r:idx_q3_r])
            Br_mean_seg_a.append(Br_mean_seg)
            
            Br_slope_q1 = np.percentile(Br_slope_t, 25)
            Br_slope_q3 = np.percentile(Br_slope_t, 75)
            Br_slope_t = np.sort(Br_slope_t).tolist()
            Br_slope_q1 = nearest(Br_slope_t, Br_slope_q1)
            Br_slope_q3 = nearest(Br_slope_t, Br_slope_q3)
            
            idx_q1_r_slope = Br_slope_t.index(Br_slope_q1)
            idx_q3_r_slope = Br_slope_t.index(Br_slope_q3)
            
            Br_slope_mean_seg = mean(Br_slope_t[idx_q1_r_slope:idx_q3_r_slope])
            Br_slope_mean_seg_a.append(Br_slope_mean_seg)
            
            Btheta_q1 = np.percentile(Btheta_fit_t, 25)
            Btheta_q3 = np.percentile(Btheta_fit_t, 75)
            Btheta_fit_t = np.sort(Btheta_fit_t).tolist()
            Btheta_q1 = nearest(Btheta_fit_t, Btheta_q1)
            Btheta_q3 = nearest(Btheta_fit_t, Btheta_q3)
            
            idx_q1_theta = Btheta_fit_t.index(Btheta_q1)
            idx_q3_theta = Btheta_fit_t.index(Btheta_q3)
            
            Btheta_mean_seg = mean(Btheta_fit_t[idx_q1_theta:idx_q3_theta])
            Btheta_mean_seg_a.append(Btheta_mean_seg)
            
            Btheta_slope_q1 = np.percentile(Btheta_slope_t, 25)
            Btheta_slope_q3 = np.percentile(Btheta_slope_t, 75)
            Btheta_slope_t = np.sort(Btheta_slope_t).tolist()
            Btheta_slope_q1 = nearest(Btheta_slope_t, Btheta_slope_q1)
            Btheta_slope_q3 = nearest(Btheta_slope_t, Btheta_slope_q3)
            
            idx_q1_theta_slope = Btheta_slope_t.index(Btheta_slope_q1)
            idx_q3_theta_slope = Btheta_slope_t.index(Btheta_slope_q3)
            
            Btheta_slope_mean_seg = mean(Btheta_slope_t[idx_q1_theta_slope:idx_q3_theta_slope])
            Btheta_slope_mean_seg_a.append(Btheta_slope_mean_seg)
            
            Bphi_q1 = np.percentile(Bphi_fit_t, 25)
            Bphi_q3 = np.percentile(Bphi_fit_t, 75)
            Bphi_fit_t = np.array(Bphi_fit_t)
            Bphi_mean_seg = mean(Bphi_fit_t[(Bphi_fit_t > Bphi_q1) & (Bphi_fit_t < Bphi_q3)])
            Bphi_mean_seg_a.append(Bphi_mean_seg)
            
            Bphi_slope_q1 = np.percentile(Bphi_slope_t, 25)
            Bphi_slope_q3 = np.percentile(Bphi_slope_t, 75)
            Bphi_slope_t = np.sort(Bphi_slope_t).tolist()
            Bphi_slope_q1 = nearest(Bphi_slope_t, Bphi_slope_q1)
            Bphi_slope_q3 = nearest(Bphi_slope_t, Bphi_slope_q3)
            
            idx_q1_phi_slope = Bphi_slope_t.index(Bphi_slope_q1)
            idx_q3_phi_slope = Bphi_slope_t.index(Bphi_slope_q3)
            
            Bphi_slope_mean_seg = mean(Bphi_slope_t[idx_q1_phi_slope:idx_q3_phi_slope])
            Bphi_slope_mean_seg_a.append(Bphi_slope_mean_seg)
        
        Br_mean_seg_a = np.asarray(Br_mean_seg_a)
        Btheta_mean_seg_a = np.asarray(Btheta_mean_seg_a)
        Bphi_mean_seg_a = np.asarray(Bphi_mean_seg_a)
        Br_slope_mean_seg_a = np.asarray(Br_slope_mean_seg_a)
        Btheta_slope_mean_seg_a = np.asarray(Btheta_slope_mean_seg_a)
        Bphi_slope_mean_seg_a = np.asarray(Bphi_slope_mean_seg_a)
        
        figure(4)
        plt.plot(t_slope_a, Br_mean_seg_a,'o-')
        xlabel('Time')
        ylabel(r"$B_{r}$, nT")
        xticks(year_norm, ['2010','2011','2012','2013','2014','2015','2016'])
        plt.title(r'$B_{r}$ From Polynomial Fit Segments for latitude %d & longitude %d'%(lat_oc,long_oc))
        plt.savefig('TimeSeriesPlots/Br_fit_seg_lat_%d_long_%d.png'%(lat_oc,long_oc), dpi=400)
        
        figure(5)
        plt.plot(t_slope_a, Btheta_mean_seg_a,'o-')
        xlabel('Time')
        ylabel(r"$B_{\theta}$, nT")
        xticks(year_norm, ['2010','2011','2012','2013','2014','2015','2016'])
        plt.title(r'$B_{\theta}$ From Polynomial Fit Segments for latitude %d & longitude %d'%(lat_oc,long_oc))
        plt.savefig('TimeSeriesPlots/Btheta_fit_seg_lat_%d_long_%d.png'%(lat_oc,long_oc), dpi=400)
        
        figure(6)
        plt.plot(t_slope_a, Bphi_mean_seg_a,'o-')
        xlabel('Time')
        ylabel(r"$B_{\phi}$, nT")
        xticks(year_norm, ['2010','2011','2012','2013','2014','2015','2016'])
        plt.title(r'$B_{\phi}$ From Polynomial Fit Segments for latitude %d & longitude %d'%(lat_oc,long_oc))
        plt.savefig('TimeSeriesPlots/Bphi_fit_seg_lat_%d_long_%d.png'%(lat_oc,long_oc), dpi=400)
        
        t_slope_a_6 = t_slope_a*6
        year_norm_6 = [-0.00023242300987797793*6, 0.1694363742010459*6, 0.3391051714119698*6, 0.5092388146426496*6, 0.6789076118535735*6, 0.8485764090644974*6, 1.0182452062754213*6]	    
	
        figure(7)
        plt.plot(t_slope_a_6, gradient(Br_mean_seg_a,t_slope_a_6),'o-')
        xlabel('Time')
        ylabel(r"$B_{r}$, nT/yr")
        xticks(year_norm_6, ['2010','2011','2012','2013','2014','2015','2016'])
        plt.ylim(-600,350)
        plt.title(r'$B_{r}$ slope for latitude %d & longitude %d'%(lat_oc,long_oc))
        plt.savefig('TimeSeriesPlots/Br_grad_lat_%d_long_%d_norm.png'%(lat_oc,long_oc), dpi=400)
        
        
        figure(8)
        plt.plot(t_slope_a_6, gradient(Btheta_mean_seg_a,t_slope_a_6),'o-')
        xlabel('Time')
        ylabel(r"$B_{\theta}$, nT/yr")
        xticks(year_norm_6, ['2010','2011','2012','2013','2014','2015','2016'])
        plt.ylim(-250,340)
        plt.title(r'$B_{\theta}$ slope for latitude %d & longitude %d'%(lat_oc,long_oc))
        plt.savefig('TimeSeriesPlots/Btheta_grad_lat_%d_long_%d_norm.png'%(lat_oc,long_oc), dpi=400)
        
        
        figure(9)
        plt.plot(t_slope_a_6, gradient(Bphi_mean_seg_a,t_slope_a_6),'o-')
        xlabel('Time')
        ylabel(r"$B_{\phi}$, nT/yr")
        xticks(year_norm_6, ['2010','2011','2012','2013','2014','2015','2016'])
        plt.ylim(-200,275)
        plt.title(r'$B_{\phi}$ slope for latitude %d & longitude %d'%(lat_oc,long_oc))
        plt.savefig('TimeSeriesPlots/Bphi_grad_lat_%d_long_%d_norm.png'%(lat_oc,long_oc), dpi=400)
        
        
        Br_grad = gradient(Br_mean_seg_a,t_slope_a_6)
        Btheta_grad = gradient(Btheta_mean_seg_a,t_slope_a_6)
        Bphi_grad = gradient(Bphi_mean_seg_a,t_slope_a_6)

        savemat('grad_timeseries_b_%d_%d.mat' %(lat_oc, long_oc), {'Br_grad':Br_grad, 'Btheta_grad':Btheta_grad, 'Bphi_grad':Bphi_grad,'time':t_slope_a})        

        plt.close('all')
    plt.close('all')

