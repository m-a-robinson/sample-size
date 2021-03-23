# Figure 2 - Creating example vGRF signals

import os
from matplotlib import pyplot as plt
import power1d
import numpy as np


# Vertical GRF data:
os.chdir('ADD YOUR WORKING DIRECTORY')
grf = np.loadtxt('vgrf_fig2.txt') # grf data 13 trials
m   = grf.mean( axis=1 )        # mean continuum
baseline = power1d.geom.Continuum1D( m )

J        = 8    # sample size
Q        = 101  # continuum size
q        = 8    # location of effect in time

# Constant addition
signal0  = power1d.geom.SquarePulse( Q , q0=5, q1=101, x0=0, x1=200)#signal
noise    = power1d.noise.SmoothGaussian( J , Q , mu = 0 , sigma = 200 , fwhm = 40 )
model0   = power1d.models.DataSample( baseline, signal0, noise, J=J ) # alternative

# Constant subtraction
signal00  = power1d.geom.SquarePulse( Q , q0=3, q1=101, x0=0, x1=-200)#signal
noise    = power1d.noise.SmoothGaussian( J , Q , mu = 0 , sigma = 200 , fwhm = 40 )
model00   = power1d.models.DataSample( baseline, signal00, noise, J=J ) # alternative

# Pulse addition
signal1  = power1d.geom.GaussianPulse( Q , q=q  , fwhm=8 , amp=250 ) #signal
noise    = power1d.noise.SmoothGaussian( J , Q , mu = 0 , sigma = 200 , fwhm = 40 )
model1   = power1d.models.DataSample( baseline, signal1, noise, J=J ) # alternative

# Pulse subtraction
signal2  = power1d.geom.GaussianPulse( Q , q=q  , fwhm=8 , amp=-250 ) #signal
noise    = power1d.noise.SmoothGaussian( J , Q , mu = 0 , sigma = 200 , fwhm = 40 )
model2   = power1d.models.DataSample( baseline, signal2, noise, J=J ) # alternative

# Double pulse
signal3  = power1d.geom.GaussianPulse( Q , q=q  , fwhm=8 , amp=300 ) #signal
signal4  = power1d.geom.GaussianPulse( Q , q=45  , fwhm=25 , amp=200 ) #signal
dbl_sig = signal3.toarray() + signal4.toarray() # convert signals to arrays
signal3.value = dbl_sig # re-store into first signal
noise    = power1d.noise.SmoothGaussian( J , Q , mu = 0 , sigma = 200 , fwhm = 40 )
model3   = power1d.models.DataSample( baseline, signal3, noise, J=J ) # alternative

# Double pulse together
signal5  = power1d.geom.GaussianPulse( Q , q=20  , fwhm=8 , amp=-300 ) #signal
signal6  = power1d.geom.GaussianPulse( Q , q=45  , fwhm=25 , amp=200 ) #signal
dbl_sig = signal5.toarray() + signal6.toarray() # convert signals to arrays
signal5.value = dbl_sig # re-store into first signal
noise    = power1d.noise.SmoothGaussian( J , Q , mu = 0 , sigma = 200 , fwhm = 40 )
model4   = power1d.models.DataSample( baseline, signal5, noise, J=J ) # alternative
plt.plot(dbl_sig)


# Plot
plt.close('all')
f, axs = plt.subplots(3,2, subplot_kw={'xticks': np.arange(0,101,20), 'ylim': [0,2000], 'sharey':True, 'sharex':True})

axs[0,0].plot(m,'k', label = 'baseline',lw=2)
model0.plot(ax=axs[0,0], color='m', with_noise=False,lw=2)
axin0 = axs[0,0].inset_axes([0.65,0.65,0.3,0.3],frameon=False, xticks=[],yticks=[],ylim=[-10, 400])
axin0.set_xticklabels('')
axin0.set_yticklabels('')
signal0.plot(ax=axin0,color='m', lw=1)
axs[0,0].set_ylabel('Force (N)')

axs[0,1].plot(m,'k', label = 'baseline',lw=2)
model00.plot(ax=axs[0,1], color='b', with_noise=False,lw=2)
axin1 = axs[0,1].inset_axes([0.65,0.65,0.3,0.3],frameon=False, xticks=[],yticks=[],ylim=[-400, 10])
axin1.set_xticklabels('')
axin1.set_yticklabels('')
signal00.plot(ax=axin1,color='b', lw=1)

axs[1,0].plot(m,'k', label = 'baseline',lw=2)
model1.plot(ax=axs[1,0], color='y', with_noise=False,lw=2)
axin2 = axs[1,0].inset_axes([0.65,0.65,0.3,0.3],frameon=False, xticks=[],yticks=[],ylim=[-10, 400])
axin2.set_xticklabels('')
axin2.set_yticklabels('')
signal1.plot(ax=axin2,color='y', lw=1)
axs[1,0].set_ylabel('Force (N)')


axs[1,1].plot(m,'k', label = 'baseline',lw=2)
model2.plot(ax=axs[1,1], color='c', with_noise=False,lw=2)
axin3 = axs[1,1].inset_axes([0.65,0.65,0.3,0.3],frameon=False, xticks=[],yticks=[],ylim=[-300, 10])
axin3.set_xticklabels('')
axin3.set_yticklabels('')
signal2.plot(ax=axin3,color='c', lw=1)


axs[2,0].plot(m,'k', label = 'baseline',lw=2)
model3.plot(ax=axs[2,0], color='r', with_noise=False,lw=2)
axin4 = axs[2,0].inset_axes([0.65,0.65,0.3,0.3],frameon=False, xticks=[],yticks=[],ylim=[-10, 400])
axin4.set_xticklabels('')
axin4.set_yticklabels('')
signal3.plot(ax=axin4,color='r', lw=1)
axs[2,0].set_xlabel('Stance (%)')
axs[2,0].set_ylabel('Force (N)')


axs[2,1].plot(m,'k', label = 'baseline',lw=2)
model4.plot(ax=axs[2,1], color='C1', with_noise=False,lw=2)
axin5 = axs[2,1].inset_axes([0.65,0.65,0.3,0.3],frameon=False, xticks=[],yticks=[],ylim=[-300, 300])
axin5.set_xticklabels('')
axin5.set_yticklabels('')
signal5.plot(ax=axin5,color='C1', lw=1)
axs[2,1].set_xlabel('Stance (%)')

plt.tight_layout()
plt.show()
