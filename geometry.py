import healpy as hp
import numpy as np
from desitarget.geomask import get_imaging_maskbits


def radec2pix(ra, dec, nside=32):
    theta = np.pi / 2. - np.radians(dec)
    phi   = np.radians(ra)
    
    all_pixel_indices = hp.ang2pix(nside, theta, phi,nest=True, lonlat=False)

    return all_pixel_indices

def hp_props(nside):
    npix = hp.nside2npix(nside)
    pixel_area = hp.nside2pixarea(nside,degrees=True)

    return npix, pixel_area

def targ_hpmap(targs, norm=None, nside=32, return_pix=False):
    pix = radec2pix(targs['RA'], targs['DEC'], nside=nside)
    
    #indice of filled pixels and corrosponding targets in pixel
    filled_pixel_index, filled_targets_per_pixel = np.unique(pix, return_counts=True) 

    #no. targets per pixel, initially 0 
    target_pixel_density = np.zeros(hp.nside2npix(nside))

    #update no. targets per pixel 
    target_pixel_density[filled_pixel_index] = filled_targets_per_pixel 
    target_pixel_density[target_pixel_density == 0] = np.NaN

    if norm is not None:
        
        for i in range(len(norm)):
            target_pixel_density[i] /= norm[i]
    
    if return_pix:
        return target_pixel_density,np.unique(pix)
    
    else:
        return  target_pixel_density

def rand_inrect(ra_lower, ra_upper, dec_lower,dec_upper,nside=32):
    ra_rand = np.random.uniform(ra_lower,ra_upper,100)
    dec_rand = np.random.uniform(dec_lower,dec_upper,100)
    
    theta = np.pi / 2. - np.radians(dec_rand.data)
    phi = np.radians(ra_rand.data)

    #indices of pixels with non-zero density, unorganised list
    all_pixel_indices = hp.ang2pix(nside, theta, phi,nest=True, lonlat=False)
    
    return np.unique(all_pixel_indices)

def bgs_mask_randoms(random):    

    # Apply custom imaging mask around bright stars etc.   
    bitnamelist = ["BRIGHT", "CLUSTER"] 

    bits = get_imaging_maskbits(bitnamelist) 

    print(bits)
    
    retain_random = np.ones(len(random['MASKBITS']), dtype=bool) #got rid of .data as didnt work below 

    for bit, ttype in zip(bits, bitnamelist): 
        # Keep random if bit not set for bits corresponding to BRIGHT and CLUSTER. 
        retain_random &= ((random['MASKBITS'] & 2**bit) == 0)  

        print(ttype, bit, np.mean(retain_random))

    #other cuts
    NOBS_mask = ((random['NOBS_G'] > 0) | (random['NOBS_R'] > 0) | (random['NOBS_Z'] > 0))
    
    retain_random = retain_random & NOBS_mask
 
    print('NOBS', np.mean(retain_random))

    return random[retain_random]
    