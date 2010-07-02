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
import pyraf
from pyraf import iraf
from iraf import stsdas,dither
no = iraf.no
yes = iraf.yes
INDEF = iraf.INDEF

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


def process_grism(asn_grism_file, asn_direct_file):
    """
process_grism(asn_grism, asn_direct)
    
Pipeline to process a set of grism/direct exposures.
    
    """
    import shutil
    import aXe2html.sexcat.sextractcat
	
    asn_grism_file  = 'ib3701060_asn.fits'
    asn_direct_file = 'ib3701050_asn.fits'
    
    #### Check that we're in the home directory of a 3D-HST field
    check_3dhst_environment(makeDirs=False)
    
    #### ASN root names
    root_grism = asn_grism_file.split('_asn.fits')[0]
    root_direct = asn_direct_file.split('_asn.fits')[0]
    
    #### Read ASN files
    asn_grism  = threedhst.utils.ASNFile(file='RAW/'+asn_grism_file)
    asn_direct = threedhst.utils.ASNFile(file='RAW/'+asn_direct_file)
        
    os.chdir('./DATA')
    
    #### Copy ASN files from RAW
    shutil.copy('../RAW/'+asn_grism_file,'./')
    shutil.copy('../RAW/'+asn_direct_file,'./')
    #### Copy FLT files
    explist = asn_grism.exposures
    explist.extend(asn_direct.exposures)
    for exp in explist:
        fits_file = threedhst.utils.find_fits_gz('../RAW/'+exp+'_flt.fits')
        print fits_file
        fi = pyfits.open(fits_file)
        try:
            os.remove('./'+exp+'_flt.fits')
        except:
            pass
            
        fi.writeto('./'+exp+'_flt.fits', clobber=True)
        
    #### Compute shifts
    try:
        os.remove(root_direct+'_tweak.fits')
    except:
        pass
        
    threedhst.shifts.compute_shifts(asn_direct_file)
    #### Check to make sure that every exposure in the ASN file
    #### has a line in the shiftfile.  If something goes wrong, 
    #### tweakshifts will usually omit an exposure from the shiftfile
    threedhst.shifts.checkShiftfile(asn_direct_file)
    #### Make a shiftfile for the GRISM ASN, with 
    #### same shifts taken from the direct images
    threedhst.shifts.make_grism_shiftfile(asn_direct_file,asn_grism_file)
    
    #### First Multidrizzle run on DIRECT images to create a detection image
    iraf.unlearn('multidrizzle')
    dither.multidrizzle ( input = asn_direct_file, shiftfile = root_direct + '_shifts.txt', \
       output = '', final_scale = INDEF, final_pixfrac = 1.0)
    cleanMultidrizzleOutput()
    
    direct_mosaic = root_direct.upper()+'_drz.fits'
    
    #### Get out SCI and WHT extensions of the Multidrizzle mosaic
    im = pyfits.open(direct_mosaic)
    SCI = im[1]
    SCI.writeto(root_direct.upper()+'_SCI.fits')
    WHT = im[2]
    WHT.writeto(root_direct.upper()+'_WHT.fits')
    
    #### Run SExtractor on the direct image, with the WHT extension as a weight image
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.options['CATALOG_NAME']    = root_direct.upper()+'_SCI.cat'
    se.options['CHECKIMAGE_NAME'] = root_direct.upper()+'_seg.fits'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = root_direct.upper()+'_WHT.fits'
    se.options['FILTER']    = 'Y'
    
    se.options['DETECT_THRESH']    = '3'   ## Default 1.5
    se.options['ANALYSIS_THRESH']    = '3' ## Default 1.5
    
    se.copyConvFile()
    se.overwrite = True
    status = se.sextractImage(root_direct.upper()+'_SCI.fits')
    #### Make region file for SExtracted catalog
    threedhst.sex.sexcatRegions(root_direct.upper()+'_SCI.cat', root_direct.upper()+'_SCI.reg', format=2)
    #### Read catalog to keep around
    sexCat = aXe2html.sexcat.sextractcat.SexCat(root_direct.upper()+'_SCI.cat')
	
def prep_name(input_asn):
    """
    make_prep_name(input_asn)
    
    Example: 
        >>> prepFile = make_prep_name("ib3714060_asn.fits")
        >>> print prepFile
        ib3714060_prep.lis
    """
    return input_asn.split('_asn.fits')[0] + '_prep.lis'


def make_aXe_lis(asn_grism_file, asn_direct_file):
    """
    status = make_aXe_lis(asn_grism_file, asn_direct_file)
    
    Make "inlist" file for aXe routines, with format
    
        grismA_flt.fits directA_flt.1.cat directA_flt.fits 0.0
        grismB_flt.fits directB_flt.1.cat directB_flt.fits 0.0
        ...
        
    Returns "True" if executes correctly, "False" on an error
    
    """
    from threedhst.utils import ASNFile
    
    asn_grism = ASNFile(asn_grism_file)
    asn_direct = ASNFile(asn_direct_file)
    if len(asn_grism.exposures) != len(asn_direct.exposures):
        print """
3D-HST / make_aXe_lis: Number of grism exposures (%d) in %s is different from the
                     : number of direct images (%d) in %s.
              """ %(len(grism_files), asn_grism, len(direct_files), asn_direct)
        return False
    
    NFILE = len(asn_grism.exposures)
    outfile = prep_name(asn_grism_file)
    fp = open(outfile,'w')
    for i in range(NFILE):
        fp.write("%s_flt.fits %s_flt_1.cat %s_flt.fits 0.0\n" 
              %(asn_grism.exposures[i], asn_direct.exposures[i], asn_direct.exposures[i]))
    fp.close()
    print "3D-HST / make_aXe_lis: Created %s\n" %outfile
    return True

def cleanMultidrizzleOutput():
    """
    cleanMultidrizzleOutput()
    
    Remove *single_[sci/wht].fits, *sci1_blt.fits, *flt*mask1.fits, *coeffs1.dat
    """
    import os,glob
    files = glob.glob('*single_???.fits')
    files.extend(glob.glob('*sci1_blt.fits'))
    files.extend(glob.glob('*flt*mask1.fits'))
    files.extend(glob.glob('*coeffs1.dat'))
    files.extend(glob.glob('IB*.run'))
    files.extend(glob.glob('IB*_med.fits'))
    for file in files:
        os.remove(file)
        
def multidrizzle_defaults():
    multidrizzle ( input = 'ib3714060_asn.fits', output = '', mdriztab = no, \
       refimage = '', runfile = '', workinplace = no, updatewcs = yes, \
       proc_unit = 'native', coeffs = 'header', context = no, clean = no, \
       group = '', ra = INDEF, dec = INDEF, build = yes, shiftfile = '', staticfile = '', \
       static = yes, static_sig = 4.0, skysub = yes, skywidth = 0.1, \
       skystat = 'median', skylower = INDEF, skyupper = INDEF, skyclip = 5, \
       skylsigma = 4.0, skyusigma = 4.0, skyuser = '', driz_separate = yes, \
       driz_sep_outnx = INDEF, driz_sep_outny = INDEF, driz_sep_kernel = 'turbo', \
       driz_sep_wt_scl = 'exptime', driz_sep_scale = INDEF, driz_sep_pixfrac = 1.0, \
       driz_sep_rot = INDEF, driz_sep_fillval = 'INDEF', driz_sep_bits = 0, \
       median = yes, median_newmasks = yes, combine_maskpt = 0.7, \
       combine_type = 'minmed', combine_nsigma = '4 3', combine_nlow = 0, \
       combine_nhigh = 1, combine_lthresh = 'INDEF', combine_hthresh = 'INDEF', \
       combine_grow = 1, blot = yes, blot_interp = 'poly5', blot_sinscl = 1.0, \
       driz_cr = yes, driz_cr_corr = no, driz_cr_snr = '3.5 3.0', \
       driz_cr_grow = 1, driz_cr_ctegrow = 0, driz_cr_scale = '1.2 0.7', \
       driz_combine = yes, final_wht_type = 'EXP', final_outnx = INDEF, \
       final_outny = INDEF, final_kernel = 'square', final_wt_scl = 'exptime', \
       final_scale = INDEF, final_pixfrac = 1.0, final_rot = 0.0, \
       final_fillval = 'INDEF', final_bits = 0, final_units = 'cps', \
       gain = INDEF, gnkeyword = INDEF, rdnoise = INDEF, rnkeyword = INDEF, exptime = INDEF, \
    )
    

    

    