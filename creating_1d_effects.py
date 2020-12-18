# -*- coding: utf-8 -*-
"""
Created on Mon Jul  6 21:07:37 2020

Creating 1D effects for 6 biomechanical datasets plus noise models.

@author: Mark Robinson, m.a.robinson@ljmu.ac.uk
"""
import os
import numpy as np
import power1d
from matplotlib import pyplot as plt
from scipy import interpolate, stats

os.chdir('C:/Users/...') # data file directory

#############################################################################
# Functions used

def resample101(x,y):
    "resample an x-y signal. Requires x to start at 0, end at 100"
    f = interpolate.interp1d(x,y)
    xnew = np.arange(0,101,1)
    ynew = f(xnew)
    
    plt.plot(x,y, 'bo')
    plt.plot(xnew,ynew, 'rx')
    plt.show()
    
    return ynew

def plot_dataset(signal, baseline, noise, null, alt, ymin, ymax):
    "Plot the model and noise inputs required for power analysis"
    labels = ['signal', 'baseline', 'noise', 'null model', 'alternative', 'null + alt']

    ax0      = plt.subplot(231)
    ax1      = plt.subplot(232)
    ax2      = plt.subplot(233)
    ax3      = plt.subplot(234)
    ax4      = plt.subplot(235)
    ax5      = plt.subplot(236) 
    AX       = [ax0, ax1, ax2, ax3, ax4, ax5]
    
    signal.plot(ax=ax0, color='g')
    baseline.plot(ax=ax1, color='y')
    noise.plot(ax=ax2)    
    null.plot(ax=ax3)
    alt.plot(ax=ax4, color='m')
    null.plot(ax=ax5)
    alt.plot(ax=ax5, color='m')
    
    for ax,s in zip( AX , labels ):
            ax.text(0.05, 0.9, s, transform=ax.transAxes, bbox=dict(facecolor='w'))
            ax.set_ylim([ymin, ymax])
            
    ax2.set_ylim([-2,2])
    
    plt.tight_layout()
    plt.show()


#############################################################################
# Effect Datasets
J  = 8
Q  = 101
JJ = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]

# Barrios 2017
# Barrios report SEM not SD
# SEM of peak medial TF JCF = 0.089, n = 25 participants
data0 = np.loadtxt('barrios2017data.txt') 
lat_wedge = resample101(data0[:,0],data0[:,1])
no_wedge = resample101(data0[:144,2],data0[:144,3])
data0sd =  0.089 * np.sqrt(25) # calc SD from SEM = 0.445

data0_bl = power1d.geom.Continuum1D( lat_wedge ) # baseline
data0_sig0  = power1d.geom.Null( Q ) # null signal - needed for null power model
data0_sig1  = power1d.geom.GaussianPulse( Q , q=20  , fwhm=8 , amp= 0.25 ) #signal
data0_nse    = power1d.noise.SmoothGaussian( J , Q , mu = 0 , sigma = 0.445 , fwhm = 20 )
data0_mod0   = power1d.models.DataSample( data0_bl, data0_sig0, data0_nse, J ) # null for power analysis
data0_mod1   = power1d.models.DataSample( data0_bl, data0_sig1, data0_nse, J ) # alternative

plot_dataset(data0_sig1, data0_bl, data0_nse, data0_mod0, data0_mod1, 0, 6)


# Phan et al. 2017 - 26 participants
# SD reported in the abstract = 0.05 BW
data1 = np.loadtxt('phan2017data.txt') 
quiet = resample101(data1[:,0],data1[:,1])
normal = resample101(data1[:,2],data1[:,3])
data1sd =  0.6 # estimate based on impact peak (GRF peak SD was 0.38)

data1_bl = power1d.geom.Continuum1D( quiet ) # baseline
data1_sig0  = power1d.geom.Null( Q ) # null signal - needed for null power model
data1_sig1  = power1d.geom.GaussianPulse( Q , q=8  , fwhm=10 , amp= 0.8 ) #signal
data1_nse    = power1d.noise.SmoothGaussian( J , Q , mu = 0 , sigma = 0.6 , fwhm = 20 )
data1_mod0   = power1d.models.DataSample( data1_bl, data1_sig0, data1_nse, J ) # null for power analysis
data1_mod1   = power1d.models.DataSample( data1_bl, data1_sig1, data1_nse, J ) # alternative

plot_dataset(data1_sig1, data1_bl, data1_nse, data1_mod0, data1_mod1, 0, 6)


# Bovi 2011 - 40 participants
# Average SD over time of the Young group, Natural condition = 0.14
data2 = np.loadtxt('bovi2011data.txt') 
adult = data2[:,0]
young = np.array(data2[:,1])
data2sd =  0.14

data2_bl = power1d.geom.Continuum1D( young ) # baseline
data2_sig0  = power1d.geom.Null( Q ) # null signal - needed for null power model
data2_sig1  = power1d.geom.GaussianPulse( Q , q=38  , fwhm=12 , amp= 0.07 ) #signal
data2_nse    = power1d.noise.SmoothGaussian( J , Q , mu = 0 , sigma = 0.14 , fwhm = 20 )
data2_mod0   = power1d.models.DataSample( data2_bl, data2_sig0, data2_nse, J ) # null for power analysis
data2_mod1   = power1d.models.DataSample( data2_bl, data2_sig1, data2_nse, J ) # alternative

plot_dataset(data2_sig1, data2_bl, data2_nse, data2_mod0, data2_mod1, 0, 0.5)



# Bakke - n=5
# Not clear from the data but very small - suggest to use '2' as a best guess 
data3 = np.loadtxt('bakke2020data.txt') 
con1 = resample101(data3[:64,0],data3[:64,1])
con2 = resample101(data3[:,2],data3[:,3])
data3sd =  2

data3_bl = power1d.geom.Continuum1D( con1 ) # baseline
data3_sig0  = power1d.geom.Null( Q ) # null signal - needed for null power model
data3_sig1  = power1d.geom.Constant( Q , amp= 6.5 ) #signal
data3_nse    = power1d.noise.SmoothGaussian( J , Q , mu = 0 , sigma = 2 , fwhm = 20 )
data3_mod0   = power1d.models.DataSample( data3_bl, data3_sig0, data3_nse, J ) # null for power analysis
data3_mod1   = power1d.models.DataSample( data3_bl, data3_sig1, data3_nse, J ) # alternative

plot_dataset(data3_sig1, data3_bl, data3_nse, data3_mod0, data3_mod1, -50, 10)



# Robinson - n=34
# Average SD of the moment Y DK signal was 0.5
data4 = np.loadtxt('robinson2014data.txt')
dk = data4[:,0]
ik = data4[:,1]
data4sd =  0.5

data4_bl = power1d.geom.Continuum1D( ik ) # baseline
data4_sig0  = power1d.geom.Null( Q ) # null signal - needed for null power model
data4_sig1  = power1d.geom.GaussianPulse( Q , q=20  , fwhm=12 , amp= 0.7 ) #signal
data4_sig2  = power1d.geom.GaussianPulse( Q , q=42  , fwhm=40 , amp= 0.9 ) #signal
dbl_sig = data4_sig1.toarray() + data4_sig2.toarray() # convert signals to arrays
data4_sig1.value = dbl_sig # re-store into first signal
data4_nse    = power1d.noise.SmoothGaussian( J , Q , mu = 0 , sigma = 0.5 , fwhm = 20 )
data4_mod0   = power1d.models.DataSample( data4_bl, data4_sig0, data4_nse, J ) # null for power analysis
data4_mod1   = power1d.models.DataSample( data4_bl, data4_sig1, data4_nse, J ) # alternative

plot_dataset(data4_sig1, data4_bl, data4_nse, data4_mod0, data4_mod1, -2, 2)



# Gomez 2017 - n=30
# In table 4 the soleus SD is 40 for the control group (at 40-60%)
data5 = np.loadtxt('gomes2017data.txt') 
con = resample101(data5[:61,0],data5[:61,1])
diab = resample101(data5[:,2],data5[:,3])
data5sd =  40

data5_bl = power1d.geom.Continuum1D( con ) # baseline
data5_sig0  = power1d.geom.Null( Q ) # null signal - needed for null power model
data5_sig1  = power1d.geom.GaussianPulse( Q , q=46  , fwhm=25 , amp= 30 ) #signal
data5_nse    = power1d.noise.SmoothGaussian( J , Q , mu = 0 , sigma = 40 , fwhm = 20 )
data5_mod0   = power1d.models.DataSample( data5_bl, data5_sig0, data5_nse, J ) # null for power analysis
data5_mod1   = power1d.models.DataSample( data5_bl, data5_sig1, data5_nse, J ) # alternative

plot_dataset(data5_sig1, data5_bl, data5_nse, data5_mod0, data5_mod1, 0, 250)
