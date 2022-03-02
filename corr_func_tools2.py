import numpy as np
import Corrfunc
from Corrfunc.mocks.DDtheta_mocks import DDtheta_mocks
from Corrfunc.utils import convert_3d_counts_to_cf

def calc_wtheta(dat,ran,nbins):

    RA = dat['RA']
    DEC = dat['DEC']
    N = len(RA)


    rand_RA = ran['RA']
    rand_DEC = ran['DEC']
    rand_N = len(rand_RA)

    # Setup the bins
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
    return convert_3d_counts_to_cf(N, N, rand_N, rand_N, DD_counts, DR_counts,DR_counts, RR_counts)

def create_axes(wtheta,nbins):
    
    bins = np.logspace(-3, 1, nbins + 1, base=10)
    x_axis = []
    for i in range(len(bins)-1):
        x_axis.append((bins[i]+bins[i+1])/2)

    y_axis = []
    for i in range(len(x_axis)):
        y_axis.append(wtheta[i]*x_axis[i]**(1-1.8))
        
    return x_axis,y_axis
