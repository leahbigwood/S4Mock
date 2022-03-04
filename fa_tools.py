from astropy.table import Table
from astropy.io import fits as fits

import glob
import pylab as pl
import numpy as np
import matplotlib
import astropy.io.fits as fits
import matplotlib.pyplot as plt
import desimodel
import desimeter

from   matplotlib.patches import Ellipse
from   astropy.time import Time
from   astropy.table import Table, vstack, join
from   desimodel.focalplane.geometry import xy2radec 
from   desimodel.io import load_fiberpos 
from   desitarget.geomask import circles
from   desitarget.sv3.sv3_targetmask import desi_mask, bgs_mask, mws_mask, scnd_mask
from   desimeter.fiberassign import fiberassign_flat_xy2radec, radec2tan
from   fiberassign.hardware import xy2radec
from   fiberassign.hardware import load_hardware
from   desispec.maskbits import fibermask as fmsk
from   desispec.fiberbitmasking import get_all_fiberbitmask_with_amp

def fa_plot(fba_path,ledger):
    colors = plt.rcParams['axes.prop_cycle'].by_key()['color']
    
    fp            = load_fiberpos() 
    fp.sort('LOCATION')
    
    patrol_radii  = 1.48/60. # degrees 
    
    tile_radii    = desimodel.focalplane.geometry.get_tile_radius_deg()
    
    def radec2standard(ra, dec, tile_ra, tile_dec):

        _ra       = np.radians(ra)
        _dec      = np.radians(dec)
        _tra      = np.radians(tile_ra)
        _tdec     = np.radians(tile_dec)

        xi        = np.sin(_ra - _tra) / (np.sin(_tdec) * np.tan(_dec) + np.cos(_tdec) * np.cos(_ra - _tra))
        eta       = (np.tan(_dec) - np.tan(_tdec) * np.cos(_ra - _tra)) / (np.tan(_tdec) * np.tan(_dec) + np.cos(_ra - _tra)) 

        return xi, eta

    rundate        = '2019-03-17T23:20:01' 
    # rundate      = '2021-03-17T23:20:01'

    hw             = load_hardware(rundate=rundate)
    
    fig = plt.figure(figsize=(15,15))
    axes = fig.add_subplot(111)

    fpath = fba_path

    hdr       = fits.open(fpath)[1].header

    tra       = hdr['TILERA']                                                  
    tdec      = hdr['TILEDEC']                                           
    fieldrot  = hdr['FIELDROT']                                                
    fa_plan   = hdr['FA_PLAN']
    fa_run    = hdr['FA_RUN']
    fa_ha     = hdr['FA_HA']                                                           

    hw        = load_hardware(rundate=fa_run)

    tt        = Time(fa_run, format='isot', scale='utc')
    fba_mjd   = tt.mjd 

    tra       = fits.open(fpath)[1].header['TILERA']
    tdec      = fits.open(fpath)[1].header['TILEDEC']

    # targ    = Table.read(tpath)
    fba       = Table.read(fpath, hdu='FASSIGN')
    ftarg     = Table.read(fpath, hdu='FTARGETS')
    favl      = Table.read(fpath, hdu='FAVAIL')

    fba       = fba[fba['DEVICE_TYPE'] != 'ETC']

    # Unbelieveably, join alters the row order.  Not sure I knew this ...
    fba       = join(fba,  ledger, keys='TARGETID', join_type='left')
    favl      = join(favl, ledger, keys='TARGETID', join_type='left')

    in_ledger = ~fba['RA'].mask

    fba.sort('LOCATION')

    assert np.all(fp['LOCATION'].data == fba['LOCATION'].data)

    # pretty sure above does same thing so why go to the effort of doing this 
    targ_ras, targ_decs = xy2radec(hw, tra, tdec, tt, fieldrot, fa_ha, fba['FIBERASSIGN_X'], fba['FIBERASSIGN_Y'], use_cs5=False)

    # Note broken fibers can be used as skies, so we don't assume fiberstatus == 0.
    standard   = (fba['FA_TYPE'].data & 2)  != 0
    sky        = (fba['FA_TYPE'].data & 4)  != 0  
    safe       = (fba['FA_TYPE'].data & 8)  != 0  
    supp       = (fba['FA_TYPE'].data & 16) != 0

    #secondary  = (fba['SV3_SCND_TARGET'].data > 0)
    bgs        = (fba['SV3_DESI_TARGET'] & desi_mask['BGS_ANY']) != 0

    bgs_bright = (fba['SV3_BGS_TARGET'] & bgs_mask['BGS_BRIGHT']) != 0
    bgs_faint  = (fba['SV3_BGS_TARGET'] & bgs_mask['BGS_FAINT']) != 0
    bgs_hip    = (fba['SV3_BGS_TARGET'] & bgs_mask['BGS_FAINT_HIP']) != 0

    mws        = (fba['SV3_DESI_TARGET'] & desi_mask['MWS_ANY']) != 0

    unassigned = (fba['FIBERSTATUS'].data & 2**fmsk.bitnum('UNASSIGNED')) != 0
    stuck      = (fba['FIBERSTATUS'].data & 2**fmsk.bitnum('STUCKPOSITIONER')) != 0
    broken     = (fba['FIBERSTATUS'].data & 2**fmsk.bitnum('BROKENFIBER')) != 0
    restricted = (fba['FIBERSTATUS'].data & 2**fmsk.bitnum('RESTRICTED')) != 0

    good_fba   = (fba['FIBERSTATUS'] == 0) | (restricted & ~stuck & ~broken)

    unassigned = unassigned & ~stuck & ~broken
    restricted = restricted & ~stuck & ~broken & ~unassigned

    # ra, dec  = xy2radec(tra, tdec, fp["X"], fp["Y"])

    # assumes x,y are flat focal plane coordinates
    # ras, decs = fiberassign_flat_xy2radec(fp['X'], fp['Y'], tra, tdec, fba_mjd, fa_ha, fieldrot, from_platemaker=False)

    # https://github.com/desihub/fiberassign/blob/a3abb8758ddff19f7d256885ad4dfa83f787bd7f/py/fiberassign/hardware.py#L408
    # xy2radec(hw, tile_ra, tile_dec, tile_obstime, tile_obstheta, tile_obsha, x, y, use_cs5, threads=0)

    ells = []

    locs = sorted(hw.loc_pos_cs5_mm.keys())
    locs = [loc for loc in locs if hw.loc_device_type[loc] == 'POS']

    xys  = np.array([hw.loc_pos_cs5_mm[loc] for loc in locs])

    # Here we take the model hardware state for an early run date 2019-03-17T23:20:01 (to get home positions).
    ras, decs       = xy2radec(hw, tra, tdec, tt, fieldrot, fa_ha, xys[:,0], xys[:,1], use_cs5=True)

    # Here we loop over fiber home locations, and plot according to their status. 
    for loccount, loc in enumerate(locs):
        ra          = ras[loccount]
        dec         = decs[loccount]

        # https://en.wikipedia.org/wiki/Angular_distance
        #
        # At fixed dec, change in ra for patrol radii is:
        patrol_dra  = patrol_radii / np.cos(np.radians(dec))

        # At fixed ra, change in dec for patrol radii is:
        patrol_ddec = patrol_radii

        e           = Ellipse(xy=(ra, dec),
                              width=2. * patrol_dra, height=2. * patrol_ddec,
                              angle=0.0)

        ells.append(e)

        axes.add_artist(e)

        e.set_clip_box(axes.bbox)
        e.set_alpha(0.3)
        e.set_facecolor('None')

        # working fiber, but restricted.  successfully assigned.
        if restricted[loccount]:
            e.set_edgecolor('g')
            e.set_linestyle('--')

        # working fiber.  not successfully assigned.
        elif unassigned[loccount]:
            e.set_edgecolor('r')
            e.set_linestyle('--')

        # working fiber.  successfully assigned and not restricted.
        elif good_fba[loccount]:
            e.set_edgecolor('g')

        else:
            e.set_edgecolor('red')

        e.set_linewidth(3)

        # If sky, we override with blue.  Can happen for both good and bad fibers.
        if sky[loccount] | supp[loccount]:
            e.set_facecolor('blue')

        if standard[loccount]:
            e.set_facecolor('gold')

        if safe[loccount]:
            e.set_facecolor('cyan')

        tids = []

        if in_ledger[loccount]:
            # --- Targets ---
            if bgs_bright[loccount]:
                tids.append(fba['TARGETID'].data[loccount])

                axes.plot(targ_ras[loccount], targ_decs[loccount], marker='.', lw=0.0, markersize=10, alpha=1., c=colors[0], label='BGS BRIGHT')

            elif bgs_faint[loccount]:
                tids.append(fba['TARGETID'].data[loccount])

                axes.plot(targ_ras[loccount], targ_decs[loccount], marker='.', lw=0.0, markersize=10, alpha=1., c=colors[1], label='BGS FAINT')

            elif mws[loccount]:
                tids.append(fba['TARGETID'].data[loccount])

                axes.plot(targ_ras[loccount], targ_decs[loccount], marker='.', lw=0.0, markersize=10, alpha=1., c=colors[2], label='MWS')

            # To be secondaries here, they were found in the ledger to get scnd. bit mask, so also main survey style targets. 
            #elif secondary[loccount]:
                #tids.append(fba['TARGETID'].data[loccount])

                #axes.plot(targ_ras[loccount], targ_decs[loccount], marker='.', lw=0.0, markersize=10, alpha=1., c=colors[3], label='SCND')

            else:
                pass

        else:
            pass

    tids = np.array(tids)

    for itype, ttype in enumerate(['BRIGHT', 'FAINT']):
        is_type = (favl['SV3_BGS_TARGET'].data & bgs_mask['BGS_{}'.format(ttype)]) != 0
        is_type = is_type & ~favl['RA'].mask
        is_type = is_type & ~np.isin(favl['TARGETID'].data, tids)

        axes.plot(favl['RA'].data[is_type].data, favl['DEC'].data[is_type].data, marker='x', lw=0.0, markersize=5, alpha=1., c=colors[itype], label='')

    left=tra-tile_radii / np.cos(np.radians(tdec))
    bottom=tdec-tile_radii

    #axes.set_xlim(right=tra, left=left)
    #axes.set_ylim(top=tdec, bottom=bottom)        

    handles, labels = axes.get_legend_handles_labels()
    unique = [(h, l) for i, (h, l) in enumerate(zip(handles, labels)) if l not in labels[:i]]
    axes.legend(*zip(*unique), frameon=False)

    pl.xlabel('Right ascension [deg.]')
    pl.ylabel('Declination [deg.]')        
    #pl.title('Tiles {} on night {}'.format(ts, night))
    pl.show()
    pl.clf()