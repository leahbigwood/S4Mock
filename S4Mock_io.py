import time
import glob
import desitarget
import numpy as np
from astropy.io import fits as fits
from desitarget.targets import desi_mask, bgs_mask, mws_mask #sv3!
import pandas as pd
import fitsio
import h5py
import healpy as hp
from desimodel.footprint import is_point_in_desi, tiles2pix
from desitarget.geomask import pixarea2nside
import sys
import healpy as hp
sys.path.append('/global/homes/l/lbigwood/S4Mock/')
import geometry
from   astropy.table import Table, vstack, unique
from ros_tools import tile2rosette


def read_mainsurvey_targets_bright(pixlist=None):
    
    
    # What's the original path.
    # select bright targets.
    # calculate mags and ebv correction.
    root = '/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/bright/'

    if pixlist == None:
        to_grab=glob.glob(root + 'targets-bright-hp-*.fits') 
    else:
        to_grab=[root + f'targets-bright-hp-{pix}.fits' for pix in pixlist]
        
    # very good practice to apply sorted, otherwise the file ordering will be random 	and non-repeatable.  
    to_grab = sorted(to_grab) 

    hp_stack = []

    #do timer as takes a while
    start = time.time() 

    #total number of pixels, not quite sure where this has come from as npix is less than this above 
    mmask = 'BGS_TARGET' # NOTE: SV3 targets.  
    ttype = 'BGS_BRIGHT' 

    min_cols = ['RA','DEC','TARGETID', 'BGS_TARGET', 'MWS_TARGET','PHOTSYS']
    
    #loop through pixels
    for i, x in enumerate(to_grab):
        # f = np.array(x[1].data)[]
        f = fitsio.read(x,columns=min_cols)
        
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
    mask,idx = np.unique(data_stack['TARGETID'],return_index=True)
    data_stack = data_stack[idx]
    
    # MW_TRANSMISSION_GRZ, EBV, 
    # data_stack['RMAG_DRED'] = 22.5 - 2.5 * np.log10(flux / mwtrans)
    # RMAG, GMAG, W1, ZMAG, RFIBMAG, 
    return data_stack

def read_mainsurvey_targets_faint(pixlist=None):
    
    
    # What's the original path.
    # select bright targets.
    # calculate mags and ebv correction.
    root = '/global/cfs/cdirs/desi/target/catalogs/dr9/1.1.1/targets/main/resolve/bright/'

    if pixlist == None:
        to_grab=glob.glob(root + 'targets-bright-hp-*.fits') 
    else:
        to_grab=[root + f'targets-bright-hp-{pix}.fits' for pix in pixlist]
        
    # very good practice to apply sorted, otherwise the file ordering will be random 	and non-repeatable.  
    to_grab = sorted(to_grab) 

    hp_stack = []

    #do timer as takes a while
    start = time.time() 

    #total number of pixels, not quite sure where this has come from as npix is less than this above 
    mmask = 'BGS_TARGET' # NOTE: SV3 targets.  
    ttype = 'BGS_FAINT' 

    min_cols = ['RA','DEC','TARGETID', 'BGS_TARGET', 'MWS_TARGET','PHOTSYS']
    
    #loop through pixels
    for i, x in enumerate(to_grab):
        # f = np.array(x[1].data)[]
        f = fitsio.read(x,columns=min_cols)
        
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
    mask,idx = np.unique(data_stack['TARGETID'],return_index=True)
    data_stack = data_stack[idx]
    
    # MW_TRANSMISSION_GRZ, EBV, 
    # data_stack['RMAG_DRED'] = 22.5 - 2.5 * np.log10(flux / mwtrans)
    # RMAG, GMAG, W1, ZMAG, RFIBMAG, 
    return data_stack


def read_mainsurvey_ledgers(data=True, uniquetargs=True):
    # desitarget.io.function
    # altmtls/ledger/
    if data:
        to_grab=glob.glob('/global/cfs/cdirs/desi/survey/ops/surveyops/trunk/mtl/main/bright/mtl-bright-hp-*.ecsv') 

    # very good practice to apply sorted, otherwise the file ordering will be random 	and non-repeatable.  
        to_grab = sorted(to_grab) 

        hp_stack = []

        #do timer as takes a while
        start = time.time() 

        #total number of pixels, not quite sure where this has come from as npix is less than this above 
        npix_todo = 200000

        mmask = 'BGS_TARGET'
        ttype = 'BGS_BRIGHT'

        #loop through pixels
        for i, x in enumerate(to_grab):
            x = pd.read_csv(x, comment='#', delimiter='\s+', usecols=['RA', 'DEC', 'TARGETID', 'BGS_TARGET', 'MWS_TARGET'])

            #mask for bgs objects
            is_bgs = (x[mmask] & bgs_mask[ttype]) != 0
            idx = np.arange(len(x))[is_bgs]
            x = x.iloc[idx] 
            hp_stack.append(x)

            #more timing stuff
            if (i % 100) == 0:
                runtime = (time.time() - start)

                print('Runtime of {:.6f} seconds after {:d} pixels'.format(runtime, i))

            if i > npix_todo:
                break

        # Create a big table from the list of tables.  
        data_stack = pd.concat(hp_stack, ignore_index=True)

        #unique targets only and put it in right table format
        
        if uniquetargs==True:
            mask,idx = np.unique(data_stack['TARGETID'],return_index=True)
            data_stack = data_stack.iloc[idx]
        
        data_stack = Table.from_pandas(data_stack)

        #more timing stuff
        runtime = (time.time() - start)
        print('\n\nTotal runtime of {:.6f} seconds after {:d} pixels'.format(runtime, npix_todo))

    return data_stack


def read_sv3_ledgers(mock=True, uniquetargs=True):
    # desitarget.io.function
    # altmtls/ledger/
    if mock==False:
        to_grab=glob.glob('/global/cscratch1/sd/mjwilson/S4MOCK/SV3REAL/SV3REALLEDGER/bright/sv3mtl-bright-hp-*.ecsv') 

    # very good practice to apply sorted, otherwise the file ordering will be random 	and non-repeatable.  
        to_grab = sorted(to_grab) 

        hp_stack = []

        #do timer as takes a while
        start = time.time() 

        #total number of pixels, not quite sure where this has come from as npix is less than this above 
        npix_todo = 200000

        #loop through pixels
        for i, x in enumerate(to_grab):
            x = pd.read_csv(x, comment='#', delimiter='\s+') #usecols=['RA', 'DEC', 'TARGETID', 'BGS_TARGET', 'MWS_TARGET'])

            hp_stack.append(x)

            #more timing stuff
            if (i % 100) == 0:
                runtime = (time.time() - start)

                print('Runtime of {:.6f} seconds after {:d} pixels'.format(runtime, i))

            if i > npix_todo:
                break

        # Create a big table from the list of tables.  
        data_stack = pd.concat(hp_stack, ignore_index=True)

        #unique targets only and put it in right table format
        
        data_stack = Table.from_pandas(data_stack)

        #more timing stuff
        runtime = (time.time() - start)
        print('\n\nTotal runtime of {:.6f} seconds after {:d} pixels'.format(runtime, npix_todo))

    if mock==True:
        to_grab=glob.glob('/global/cscratch1/sd/lbigwood/S4MOCK/mockledger/sv3/bright/sv3mtl-bright-hp-*.ecsv') 

    # very good practice to apply sorted, otherwise the file ordering will be random 	and non-repeatable.  
        to_grab = sorted(to_grab) 

        hp_stack = []

        #do timer as takes a while
        start = time.time() 

        #total number of pixels, not quite sure where this has come from as npix is less than this above 
        npix_todo = 200000

        #loop through pixels
        for i, x in enumerate(to_grab):
            x = pd.read_csv(x, comment='#', delimiter='\s+') #usecols=['RA', 'DEC', 'TARGETID', 'BGS_TARGET', 'MWS_TARGET'])

            hp_stack.append(x)

            #more timing stuff
            if (i % 100) == 0:
                runtime = (time.time() - start)

                print('Runtime of {:.6f} seconds after {:d} pixels'.format(runtime, i))

            if i > npix_todo:
                break

        # Create a big table from the list of tables.  
        data_stack = pd.concat(hp_stack, ignore_index=True)

        #unique targets only and put it in right table format
        
        data_stack = Table.from_pandas(data_stack)

        #more timing stuff
        runtime = (time.time() - start)
        print('\n\nTotal runtime of {:.6f} seconds after {:d} pixels'.format(runtime, npix_todo))

    
    
    return data_stack





def read_desitargetrandoms(number=1):
    ns = np.arange(1, number+1, 1)
    
    min_cols = ['RA' , 'DEC', 'MASKBITS', 'NOBS_G', 'NOBS_R', 'NOBS_Z']
    
    rand = np.hstack([fitsio.read(f'/global/cfs/cdirs/desi/target/catalogs/dr9/0.49.0/randoms/resolve/randoms-{nn}-0.fits', columns=min_cols) for nn in ns])
    return rand 

def tile2rosette(tile):
    if tile < 433:
        return (tile-1)//27
    else:
        if tile >= 433 and tile < 436:
            return 13
        if tile >= 436 and tile < 439:
            return 14
        if tile >= 439 and tile < 442:
            return 15
        if tile >= 442 and tile <=480:
            return (tile-442)//3
            
        if tile > 480:
            return tile//30    
    return 999999 #shouldn't be any more?


def read_sv3tiles():
    tiles = Table.read('/global/cscratch1/sd/mjwilson/S4MOCK/tiles-sv3.ecsv')
    tiles = tiles[(tiles['STATUS'] == 'done') & (tiles['PROGRAM']=='BRIGHT')]
    tiles['ROSETTE'] = np.array([tile2rosette(x) for x in tiles['TILEID']])
    return tiles

def read_clustering_dat():
    # TODO: cut to dec. 32.375 
    x = fits.open('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/LSScats/2.1/BGS_BRIGHT_clustering.dat.fits')
    dat = x[1].data
    dat = dat[(dat['DEC']<32.375)]
    
    return dat
        


def read_clustering_ran(): 
    x = fits.open('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/everest/LSScats/2.1/BGS_BRIGHT_S_9_clustering.ran.fits')
    return x[1].data


def read_mxxl(small=True,nside=32):

    if small == True:
        root  = "/global/project/projectdirs/desi/mocks/bgs/MXXL/small/"

        fpath = root + "galaxy_catalogue_small.hdf5"
    
    else:
        fpath = '/project/projectdirs/desi/mocks/bgs/MXXL/full_sky/v0.0.4/BGS_r20.6.hdf5'
    
    
    f   = h5py.File(fpath, mode='r')

    ra  = f["Data/ra"][...]
    dec = f["Data/dec"][...]
    z   = f["Data/z_obs"][...]
    r   = f["Data/app_mag"][...]

    f.close()
    
    temp = np.c_[ra, dec, z, r]

    mxxl = Table(temp, names=['RA', 'DEC','Z_OBS','APP_MAG'])

    mxxl  = mxxl
    
    theta = np.pi / 2. - np.radians(mxxl['DEC'].data)
    phi   = np.radians(mxxl['RA'].data)

    mxxl['HPX'] = hp.ang2pix(nside, theta, phi,nest=True, lonlat=False)

    # single_pixel_mxxl['BGS_BRIGHT'] = single_pixel_mxxl['RMAG_DRED'] <= 19.5
    # single_pixel_mxxl['BGS_FAINT']  = (single_pixel_mxxl['RMAG_DRED'] > 19.5) & (single_pixel_mxxl['RMAG_DRED'] <= 20.175)

    return  mxxl

def read_mxxl_v2():

    #read targets in 
    #get tiles
    tiles = read_sv3tiles()
    # closest nside to DESI tile area of ~7 deg
    nside = pixarea2nside(7.)
    # ADM determine the pixels that touch the tiles.
    pixlist = tiles2pix(nside, tiles=tiles)
    #read in mxxl
    mxxl =read_mxxl(small=False,nside=nside)
    #read in our mxxl targets but having this nside and this pixlist 
    targets = mxxl[np.in1d(mxxl['HPX'],pixlist)]
    #restrict only to targets in the requested tiles...
    ii = is_point_in_desi(tiles, targets["RA"], targets["DEC"])
    targets = targets[ii]
    #now get pixlist in nside=32
    pix32 = geometry.radec2pix(targets,nside=32)
    targets['HPX']=pix32
    
    return targets

nights    = [x.split('/')[-1] for x in sorted(glob.glob('/global/cscratch1/sd/mjwilson/altmtls/ledger/initial/Univ000/fa/SV3' + '/*'))]

def read_favail_mock():
    return vstack([Table(fits.open(x)['FAVAIL'].data) for night in nights for x in glob.glob(f'/global/cscratch1/sd/mjwilson/altmtls/ledger/initial/Univ000/fa/SV3/{night}/fba-*.fits')])

def read_fassign_mock():
    return vstack([Table(fits.open(x)['FASSIGN'].data) for night in nights for x in glob.glob(f'/global/cscratch1/sd/mjwilson/altmtls/ledger/initial/Univ000/fa/SV3/{night}/fba-*.fits')])

def read_targ_mock():
    return vstack([Table.read(x) for night in nights for x in glob.glob(f'/global/cscratch1/sd/mjwilson/altmtls/ledger/initial/Univ000/fa/SV3/{night}//*-targ.fits')])

def read_sv3_randoms(number=1):
    if number ==1:
        rand = Table.read('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random0/rancomb_brightwdup_Alltiles.fits')
        rand             = unique(rand, keys='TARGETID')
        return rand
    else:
        return vstack([Table.read(x) for x in glob.glob('/global/cfs/cdirs/desi/survey/catalogs/SV3/LSS/random*/rancomb_brightwdup_Alltiles.fits')])

def read_mock_ledger():
    return vstack([Table.read(x) for x in glob.glob('/global/cscratch1/sd/mjwilson/altmtls/ledger/initial/Univ000/sv3/bright/sv3mtl-bright-hp-*.ecsv')])

def read_init_ledger():
    return vstack([Table.read(x) for x in glob.glob('/global/cscratch1/sd/mjwilson/altmtls/ledger/initial/Univ000/sv3/bright/initial/sv3mtl-bright-hp-*.ecsv')])