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

