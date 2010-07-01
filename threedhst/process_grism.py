"""
process_grism.py - Pipeline to process a set of grism+direct exposures

The overall outline of the process is as follows:

o Copy a fresh set of raw flt files that haven't been touched by 
  multidrizzle
  
o Get shifts based on direct image.
  
    Start with tweakshifts, will eventually need to refine shifts to align
    to WCS reference frame, such as ACS/GOODS

o Run Multidrizzle on direct images to make a detection image

o Run sextractor on the combined direct image to make a catalog

    Need to process SExtractor output catalog to have MAG_F[LAMBDA]W column
    Make a region file of the detected sources
    
    Also make a region file of the current mosaic at this point
    
o Run aXe/tiolprep to process the mosaic catalog for the individual exposures

    Generate a "prep.lis" file for the image/catalog combinations 
    make_aXe_lis(asn_grism, asn_direct)
    
o Run aXe/taxeprep to subtract the sky background from the grism images

o Run Multidrizzle on grism images to flag cosmic rays in the DQ header

o CHECK if newer calwfc3 pipeline flags all of the bad pixels.  If not, 
  could need to apply a static bad pixel mask by hand.

o Prepare a fluxcube for either a single band, or for multiple bands as
  they become available from CANDELS
  
o Run aXe/taxecore for most of the aXe processing (contamination, PET, etc.)

o Run aXe/tdrzprep to prepare the Drizzle stage

o Run aXe/taxedrizzle for the final drizzle / spectra generation

o Clean up intermediate files

o Additional catalog/high level steps to prepare spectra for actual use.

"""

#   $URL$
#   $Rev$
#   $Author$
#   $Date$

import os
import pyfits
from pyraf import iraf
no = iraf.no
yes = iraf.yes
import threedhst
   
def set_aXe_environment(grating='G141'):
    """
set_aXe_environment(grating='G141')
    
Setup aXe environment variables:
    AXE_IMAGE_PATH   = ./DATA
    AXE_CONFIG_PATH  = ./CONF
    AXE_OUTPUT_PATH  = ./OUTPUT_[GRATING]
    AXE_DRIZZLE_PATH = ./DRIZZLE_[GRATING]
    
CONF should by symlinked from /research/HST/GRISM/CONF
    
    """
    os.environ['AXE_IMAGE_PATH'] = './DATA/'
    print '--> variable AXE_IMAGE_PATH   set to "./DATA"'
    
    os.environ['AXE_CONFIG_PATH'] = './CONF/'
    print '--> variable AXE_CONFIG_PATH  set to "./CONF/"'
     
    os.environ['AXE_OUTPUT_PATH'] = './OUTPUT_'+grating.upper()+'/'
    print '--> variable AXE_OUTPUT_PATH  set to "./OUTPUT_'+grating.upper()+'/"'
    
    os.environ['AXE_DRIZZLE_PATH'] = './DRIZZLE_'+grating.upper()+'/'
    print '--> variable AXE_DRIZZLE_PATH set to "./DRIZZLE_'+grating.upper()+'/"'


def check_3dhst_environment(makeDirs=False):
    """
check_3dhst_environment(makeDirs=False)
    
Check that all of the expected directories exist for 
3D-HST data reduction.
    
If makeDirs is True, then mkdir any that isn't found in ./
    
    """    
    directories = ['DATA','CAT','RAW','OUTPUT_G141','DRIZZLE_G141']
    for dir in directories:
        if not os.path.exists(dir):
            if makeDirs:
                os.mkdir(dir)
            else:
                raise IOError('Directory %s doesn\'t exist in %s.' %(dir,os.getcwd()))


def process_grism(asn_grism, asn_direct):
    """
process_grism(asn_grism, asn_direct)
    
Pipeline to process a set of grism/direct exposures.
    
    """
    root_grism = asn_grism.split('_asn.fits')[0]
    root_direct = asn_direct.split('_asn.fits')[0]
    
    os.chdir('./DATA')
    
    #### Compute shifts
    threedhst.compute_shifts(asn_direct)
    
    #### Make a shiftfile for the GRISM ASN, with 
    #### same shifts taken from the direct images
    threedhst.make_grism_shiftfile(asn_direct,asn_grism)
    #### First Multidrizzle run on DIRECT image for detection image
    iraf.unlearn('multidrizzle')
    iraf.multidrizzle ( input = asn_direct, shiftfile = root_direct + '_shifts.txt', \
       output = '', final_scale = INDEF, final_pixfrac = 1.0)
       
    direct_mosaic = root_direct.upper()+'_drz.fits'


                    
def multidrizzle_defaults():
    multidrizzle ( input = 'ib3714060_asn.fits', output = '', mdriztab = no, \
       refimage = '', runfile = '', workinplace = no, updatewcs = yes, \
       proc_unit = 'native', coeffs = 'header', context = no, clean = no, \
       group = '', ra = , dec = , build = yes, shiftfile = '', staticfile = '', \
       static = yes, static_sig = 4.0, skysub = yes, skywidth = 0.1, \
       skystat = 'median', skylower = INDEF, skyupper = INDEF, skyclip = 5, \
       skylsigma = 4.0, skyusigma = 4.0, skyuser = '', driz_separate = yes, \
       driz_sep_outnx = , driz_sep_outny = , driz_sep_kernel = 'turbo', \
       driz_sep_wt_scl = 'exptime', driz_sep_scale = INDEF, driz_sep_pixfrac = 1.0, \
       driz_sep_rot = INDEF, driz_sep_fillval = 'INDEF', driz_sep_bits = 0, \
       median = yes, median_newmasks = yes, combine_maskpt = 0.7, \
       combine_type = 'minmed', combine_nsigma = '4 3', combine_nlow = 0, \
       combine_nhigh = 1, combine_lthresh = 'INDEF', combine_hthresh = 'INDEF', \
       combine_grow = 1, blot = yes, blot_interp = 'poly5', blot_sinscl = 1.0, \
       driz_cr = yes, driz_cr_corr = no, driz_cr_snr = '3.5 3.0', \
       driz_cr_grow = 1, driz_cr_ctegrow = 0, driz_cr_scale = '1.2 0.7', \
       driz_combine = yes, final_wht_type = 'EXP', final_outnx = , \
       final_outny = , final_kernel = 'square', final_wt_scl = 'exptime', \
       final_scale = INDEF, final_pixfrac = 1.0, final_rot = 0.0, \
       final_fillval = 'INDEF', final_bits = 0, final_units = 'cps', \
       gain = , gnkeyword = , rdnoise = , rnkeyword = , exptime = , \
    )
    

    

    