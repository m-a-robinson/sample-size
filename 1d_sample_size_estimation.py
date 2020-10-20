# -*- coding: utf-8 -*-
"""
Sample size calculation for 6 biomechanical datasets with noise.
Note: Multiple simulations are generated so this script may take a few minutes to run

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

def power_analysis(null, alt):
    "Power calculation - input models.DataSameple objects with plot"
    # create experiment models:
    teststat = power1d.stats.t_1sample
    emodel0  = power1d.models.Experiment( null , teststat )    # null hypothesis
    emodel1  = power1d.models.Experiment( alt , teststat )    # alternative hypothesis   

    # simulate the experiments:
    sim      = power1d.ExperimentSimulator( emodel0 , emodel1 )
    results  = sim.simulate( 5000, progress_bar=True )

    # visualize the power results:
    results.plot()     
    plt.show()
    
    return emodel0, emodel1
        
def sample_size_calc(model0, model1, emodel0, emodel1, q):
    "Simulate sample size calculation with plot"
    
    np.random.seed(0)    #seed the random number generator
    JJ         = [5, 10, 15, 20, 25, 30, 35, 40, 45, 50]  #sample sizes
    tstat      = power1d.stats.t_1sample  #test statistic function
    emodel0    = power1d.models.Experiment(model0, tstat) # null
    emodel1    = power1d.models.Experiment(model1, tstat) # alternative
    sim        = power1d.ExperimentSimulator(emodel0, emodel1)
    
    ### loop through the different sample sizes:
    power_omni = []
    power_coi  = []
    coir       = 3
    for J in JJ:
            emodel0.set_sample_size( J )
            emodel1.set_sample_size( J )
            results = sim.simulate( 5000 )
            results.set_coi( ( q , coir ) )  #create a COI at the signal location
            power_omni.append( results.p_reject1 )  #omnibus power
            power_coi.append( results.p_coi1[0] )   #coi power
    
    #(3) Plot the results:
    ax = plt.axes()
    ax.plot(JJ, power_omni, 'o-', label='Omnibus')
    ax.plot(JJ, power_coi,  'o-', label='COI (radius=%d)' %coir)
    ax.axhline(0.8, color='k', linestyle='--')
    ax.set_xlabel('Sample size', size=14)
    ax.set_ylim([0,1])
    ax.set_ylabel('Power', size=14)
    ax.legend()
    plt.show()
    
    return power_omni, power_coi


def effect_size(model1,baseline,noise):
    "1D effect size based on power1d models"
    es = ( model1.value0 - baseline.toarray() ) / noise.sigma
    plt.plot(es)
    return(es)

def power_1sample(n, effect, alpha=0.05):
	delta  = effect * n ** 0.5   # non-centrality parameter
	u      = stats.t.isf( alpha , n-1 )
	return stats.nct.sf( u , n-1 , delta )

def sample_size_1sample(effect, alpha=0.05, target_power=0.8, n_range=(3,50)):
	'''
	Adjust n_range to a broader sample size range if necessary
	'''
	n   = np.arange( *n_range )
	p   = np.array([power_1sample(nn, effect, alpha) for nn in n])
	ind = np.argwhere(p > target_power).flatten()[0]
	return n[ind]


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
data0_es = effect_size(data0_mod1, data0_bl, data0_nse) # 1D effect size
data0_ss_0d = sample_size_1sample(data0_es[20],alpha=0.05, target_power=0.8, n_range=(3,50)) # 0D sample size


plot_dataset(data0_sig1, data0_bl, data0_nse, data0_mod0, data0_mod1, 0, 6)
(data0_emod0, data0_emod1) = power_analysis(data0_mod0, data0_mod1)
(data0_pwr_omni, data0_pwr_coi) = sample_size_calc(data0_mod0, data0_mod1, data0_emod0, data0_emod1,20)


# Phan et al. 2017 - 26 participants
# SD reported in the abstract = 0.05 BW
data1 = np.loadtxt('phan2017data.txt') 
quiet = resample101(data1[:,0],data1[:,1])
normal = resample101(data1[:,2],data1[:,3])
data1sd =  0.6 # estimate based on impact peak (GRF peak SD was 0.38)

data1_bl = power1d.geom.Continuum1D( quiet ) # baseline
data1_sig0  = power1d.geom.Null( Q ) # null signal - needed for null power model
data1_sig1  = power1d.geom.GaussianPulse( Q , q=8  , fwhm=10 , amp= 0.8 ) #signal
data1_nse    = power1d.noise.SmoothGaussian( J , Q , mu = 0 , sigma = 0.05 , fwhm = 20 )
data1_mod0   = power1d.models.DataSample( data1_bl, data1_sig0, data1_nse, J ) # null for power analysis
data1_mod1   = power1d.models.DataSample( data1_bl, data1_sig1, data1_nse, J ) # alternative
data1_es = effect_size(data1_mod1, data1_bl, data1_nse) # 1D effect size
data1_ss_0d = sample_size_1sample(data1_es[8],alpha=0.05, target_power=0.8, n_range=(3,50)) # 0D sample size

plot_dataset(data1_sig1, data1_bl, data1_nse, data1_mod0, data1_mod1, 0, 6)
(data1_emod0, data1_emod1) = power_analysis(data1_mod0, data1_mod1)
(data1_pwr_omni, data1_pwr_coi) = sample_size_calc(data1_mod0, data1_mod1, data1_emod0, data1_emod1,20)


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
data2_es = effect_size(data2_mod1, data2_bl, data2_nse) # 1D effect size
data2_ss_0d = sample_size_1sample(data2_es[38],alpha=0.05, target_power=0.8, n_range=(3,50)) # 0D sample size

plot_dataset(data2_sig1, data2_bl, data2_nse, data2_mod0, data2_mod1, 0, 0.5)
(data2_emod0, data2_emod1) = power_analysis(data2_mod0, data2_mod1)
(data2_pwr_omni, data2_pwr_coi) = sample_size_calc(data2_mod0, data2_mod1, data2_emod0, data2_emod1,38)


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
data3_es = effect_size(data3_mod1, data3_bl, data3_nse) # 1D effect size
data3_ss_0d = sample_size_1sample(data3_es[52],alpha=0.05, target_power=0.8, n_range=(3,50)) # 0D sample size

plot_dataset(data3_sig1, data3_bl, data3_nse, data3_mod0, data3_mod1, -50, 10)
(data3_emod0, data3_emod1) = power_analysis(data3_mod0, data3_mod1)
(data3_pwr_omni, data3_pwr_coi) = sample_size_calc(data3_mod0, data3_mod1, data3_emod0, data3_emod1,50) # What should q be for constant effect?


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
data4_es = effect_size(data4_mod1, data4_bl, data4_nse) # 1D effect size
data4_ss_0d = sample_size_1sample(data4_es[20],alpha=0.05, target_power=0.8, n_range=(3,50)) # 0D sample size

plot_dataset(data4_sig1, data4_bl, data4_nse, data4_mod0, data4_mod1, -2, 2)
(data4_emod0, data4_emod1) = power_analysis(data4_mod0, data4_mod1)
(data4_pwr_omni, data4_pwr_coi) = sample_size_calc(data4_mod0, data4_mod1, data4_emod0, data4_emod1,25) # What should q be for this effect?


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
data5_es = effect_size(data5_mod1, data5_bl, data5_nse) # 1D effect size
data5_ss_0d = sample_size_1sample(data5_es[46],alpha=0.05, target_power=0.8, n_range=(3,50)) # 0D sample size

plot_dataset(data5_sig1, data5_bl, data5_nse, data5_mod0, data5_mod1, 0, 250)
(data5_emod0, data5_emod1) = power_analysis(data5_mod0, data5_mod1)
(data5_pwr_omni, data5_pwr_coi) = sample_size_calc(data5_mod0, data5_mod1, data5_emod0, data5_emod1,46)



