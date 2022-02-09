import numpy as np
from os.path import dirname, abspath, join as pjoin
from Corrfunc.theory.DD import DD
from Corrfunc.io import read_catalog
from Corrfunc.utils import convert_3d_counts_to_cf
from Corrfunc.mocks import *
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import glob, time 
from astropy.io import fits as fits
from astropy.table import Table, vstack
from astropy import constants as const
from astropy import units as u
from astropy.table import QTable
import math
import sys
sys.path.append('/global/homes/l/lbigwood/LSS/py')
import LSS
import LSS.SV3
import LSS.SV3.cattools as cattools
from desitarget.targets import desi_mask, bgs_mask, mws_mask 
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from Corrfunc.utils import convert_3d_counts_to_cf


#make big file of bright targets

def load_cat():

    to_grab=glob.glob('/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/bright/targets-bright-hp-*.fits') 

    # very good practice to apply sorted, otherwise the file ordering will be random 	and non-repeatable.  
    to_grab = sorted(to_grab) 

    hp_stack = []

    #do timer as takes a while
    start = time.time() 

    #total number of pixels, not quite sure where this has come from as npix is less than this above 


    mmask = 'BGS_TARGET'
    ttype = 'BGS_BRIGHT'

    #loop through pixels
    for i, x in enumerate(to_grab):
        x = fits.open(x)
        f = np.array(x[1].data)[['RA','DEC','TARGETID', 'BGS_TARGET', 'MWS_TARGET','PHOTSYS']]
        #mask for bgs objects

        is_bgs = (f[mmask] & bgs_mask[ttype]) != 0
        #idx = np.arange(len(x))[is_bgs]
        #x = x.iloc[idx] 
        hp_stack.append(f[is_bgs])

        #more timing stuff
        if (i % 20) == 0:
            runtime = (time.time() - start)

            print('Runtime of {:.6f} seconds after {:d} pixels'.format(runtime, i))

        

    data_stack = np.concatenate(hp_stack)      
    data_stack = Table(data_stack)

    #only select decals area
    data_stack = data_stack[(data_stack['PHOTSYS']=='S')]
    
    return data_stack

def angular_corrfunc(data_stack, ra_min, ra_max, dec_min, dec_max, num_threads):
    data_stack_small = data_stack[((data_stack['RA']>ra_min) & (data_stack['RA']<ra_max) & (data_stack['DEC']>dec_min)&(data_stack['DEC']<dec_max))]

    RA = data_stack_small['RA']
    DEC = data_stack_small['DEC']
    N = len(RA)

    # Read the supplied randoms catalog
    f = fits.open('/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-1-0.fits')
    random1=f[1].data

    random1_small = random1[((random1['RA']>ra_min) & (random1['RA']<ra_max) & (random1['DEC']>dec_min) & (random1['DEC']<dec_max))]

    rand_RA = random1_small['RA']
    rand_DEC = random1_small['DEC']
    rand_N = len(rand_RA)

    # Setup the bins
    nbins = 30
    bins = np.logspace(-3, 1, nbins + 1, base=10)
    #bins = np.linspace(0.001, 10.0, nbins + 1)

    # Number of threads to use
    nthreads = 2

    # Auto pair counts in DD
    autocorr=1
    DD_counts = DDtheta_mocks(autocorr, nthreads, bins,RA, DEC)

    # Cross pair counts in DR
    autocorr=0
    DR_counts = DDtheta_mocks(autocorr, nthreads, bins,RA, DEC,RA2=rand_RA, DEC2=rand_DEC)

    # Auto pairs counts in RR
    autocorr=1
    RR_counts = DDtheta_mocks(autocorr, nthreads, bins, rand_RA, rand_DEC)

    # All the pair counts are done, get the angular correlation function
    wtheta = convert_3d_counts_to_cf(N, N, rand_N, rand_N, DD_counts, DR_counts,DR_counts, RR_counts)
    
    plt.rc('xtick',direction='in',labelsize=22,top=True)
    plt.rc('ytick',direction='in',labelsize=22, right = True)
    plt.rc('xtick.major',size = 22)
    plt.rc('xtick.minor',size = 6)
    plt.rc('ytick.major',size = 22)
    plt.rc('ytick.minor',size = 6)
    plt.rc('axes', labelsize = 22)
    plt.rc('legend',fontsize=22)
    plt.rc('font', family='serif',size=20)


    bins = np.logspace(-3, 1, nbins + 1, base=10)
    x_axis = []
    for i in range(len(bins)-1):
        x_axis.append((bins[i]+bins[i+1])/2)

    y_axis = []
    for i in range(len(x_axis)):
        y_axis.append(wtheta[i]*x_axis[i]**(1.8-1))

    plt.figure(figsize=(12,10))
    plt.plot(x_axis, y_axis,color='#1e3799')
    plt.xlabel(r'$\theta$ / deg')
    plt.ylabel(r'$\omega (\theta)\times \theta^{-(1-\gamma)}$')
    plt.xscale('log')
    plt.yscale('log')
    plt.ylim(1e-4,5e-2)
    plt.savefig('corr_func_new.png')
    plt.show()
    
    
    
    