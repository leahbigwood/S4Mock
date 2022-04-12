import sys
sys.path.append('/global/homes/l/lbigwood/S4Mock/')
sys.path.append('/global/project/projectdirs/desi/mocks/bgs/MXXL/one_percent')
import footprint 
import h5py
import S4Mock_io
from corr_func_tools import calc_wtheta,create_axes
import numpy as np
from astropy.table import Table

def read_mxxl_real():
    mxxl = S4Mock_io.read_mxxl(small=False)

    root  = "/global/cfs/cdirs/desi/mocks/bgs/MXXL/one_percent/"
    fpath = root + "one_percent_v2.hdf5"
    f     = h5py.File(fpath, mode='r')

    # assert 
    
    mxxl['NMOCK'] = f['nmock'][:]

    return mxxl

def read_mxxl_real_rand():
    version =2 
    Ngal=3000000

    if version==2:
        tiles_file = "/project/projectdirs/desi/mocks/bgs/MXXL/one_percent/tiles/v2/tiles-sv3.ecsv"
        Nmock = 36
    elif version==1:
        tiles_file = "/project/projectdirs/desi/mocks/bgs/MXXL/one_percent/tiles/v1/onepercent.fits"
        Nmock = 96

    # generate random galaxy ra, dec, then find galaxies in footprint                                                                             
    #Ngal = 1000000                                                                                                                               
    ra = np.random.rand(Ngal) * 360
    sin_dec = np.random.rand(Ngal) * 2 - 1
    dec = np.arcsin(sin_dec) * 180/np.pi

    ra_array = []
    dec_array = []
    nmock_array = []

    count = int(0)
    for i in range(Nmock):
        in_footprint = footprint.in_onepercent_footprint(ra, dec, tiles_file, nmock=i, version=version)
        ra_array.append(ra[in_footprint])
        dec_array.append(dec[in_footprint])
        nmock_array.append([count]*len(ra[in_footprint]))
        count +=1

    ra_array = np.array(ra_array)
    dec_array = np.array(dec_array)
    nmock_array = np.array(nmock_array)


    ra_array = np.concatenate(ra_array)
    dec_array = np.concatenate(dec_array)
    nmock_array = np.concatenate(nmock_array)

    return Table(np.c_[ra_array, dec_array,nmock_array], names=['RA', 'DEC','NMOCK'],dtype=[np.float32,np.float32,np.int32])


def calc_mxxl_errors():
    mxxl = read_mxxl_real()
    rand = read_mxxl_real_rand()
    for i in range(36):
        mxxl_n = mxxl[(mxxl['NMOCK']==i)]
        rand_n = rand[(rand['NMOCK']==i)]
        calc_wtheta(mxxl_n,rand_n)