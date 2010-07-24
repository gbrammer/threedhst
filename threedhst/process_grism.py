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

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import os
import pyfits
import pyraf
from pyraf import iraf
from iraf import stsdas,dither,slitless,axe 

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
    print '--> variable AXE_DRIZZLE_PATH set to' + \
               '"./DRIZZLE_'+grating.upper()+'/"'


def check_3dhst_environment(makeDirs=False):
    """
check_3dhst_environment(makeDirs=False)
    
Check that all of the expected directories exist for 
3D-HST data reduction.
    
If makeDirs is True, then mkdir any that isn't found in ./
    
    """    
    directories = ['DATA','RAW','OUTPUT_'+threedhst.options['GRISM_NAME'],
                   'DRIZZLE_'+threedhst.options['GRISM_NAME'],'CONF']
    for dir in directories:
        if not os.path.exists(dir):
            if makeDirs:
                os.mkdir(dir)
            else:
                raise IOError('Directory %s doesn\'t exist in %s.'
                               %(dir,os.getcwd()))
    
    if not os.path.exists('CONF/'+threedhst.options['CONFIG_FILE']):
        raise IOError("options['CONFIG_FILE'] doesn't exist: CONF/" +
                      threedhst.options['CONFIG_FILE'])
    
    if not os.path.exists('CONF/'+threedhst.options['SKY_BACKGROUND']):
        raise IOError("options['SKY_BACKGROUND'] doesn't exist:" +   
                      "CONF/"+threedhst.options['SKY_BACKGROUND'])
                      

def reduction_script(asn_grism_file=None, asn_direct_file=None):
    """
process_grism(asn_grism, asn_direct)
    
Pipeline to process a set of grism/direct exposures.
    
    """
    import shutil
    import aXe2html.sexcat.sextractcat
    import glob
    import pywcs
    import numpy as np
    
    ##########################################
    ####   Bookkeeping, set up DATA directory
    ##########################################
    
    #### Check that we're in the home directory of a 3D-HST field
    check_3dhst_environment(makeDirs=False)
    
    os.chdir('./DATA')
    
    if not asn_grism_file:
        asn_grism_file  = 'ib3721060_asn.fits'
    if not asn_direct_file:
        asn_direct_file = 'ib3721050_asn.fits'
    
    #### ASN root names
    root_grism = asn_grism_file.split('_asn.fits')[0]
    root_direct = asn_direct_file.split('_asn.fits')[0]
    
    threedhst.currentRun['asn_grism_file'] = asn_grism_file
    threedhst.currentRun['asn_direct_file'] = asn_direct_file
    threedhst.currentRun['root_grism'] = root_grism
    threedhst.currentRun['root_direct'] = root_direct
        
    #### Copy ASN files from RAW, if they don't already exist
    if not os.path.exists(asn_grism_file):
        shutil.copy('../RAW/'+asn_grism_file,'./')
    if not os.path.exists(asn_direct_file):
        shutil.copy('../RAW/'+asn_direct_file,'./')
    
    #### Read ASN files
    asn_grism  = threedhst.utils.ASNFile(file=asn_grism_file)
    asn_direct = threedhst.utils.ASNFile(file=asn_direct_file)
    ### Force PROD-DTH to be the root of the input ASN file
    asn_grism.product = root_grism
    asn_direct.product = root_direct
    asn_grism.writeToFile(asn_grism_file)
    asn_direct.writeToFile(asn_direct_file)
    
    #### Copy FLT files
    explist = []
    explist.extend(asn_grism.exposures)
    explist.extend(asn_direct.exposures)
    for exp in explist:
        fits_file = threedhst.utils.find_fits_gz('../RAW/'+exp+'_flt.fits')
        fi = pyfits.open(fits_file)
        try:
            os.remove('./'+exp+'_flt.fits')
        except:
            dummy = 1
        print exp
        fi.writeto('./'+exp+'_flt.fits', clobber=True)
        #### Apply DQ mask (.mask.reg), if it exists
        threedhst.dq.apply_dq_mask(fits_file, addval=2048)
        
    threedhst.currentRun['step'] = 'COPY_FROM_RAW'
    
    #########################################
    ####   Compute shifts
    #########################################
    
    try:
        os.remove(root_direct+'_tweak.fits')
    except:
        dummy = 1
        
    threedhst.shifts.compute_shifts(asn_direct_file)
    #### Check to make sure that every exposure in the ASN file
    #### has a line in the shiftfile.  If something goes wrong, 
    #### tweakshifts will usually omit an exposure from the shiftfile
    threedhst.shifts.checkShiftfile(asn_direct_file)
    #### Make a shiftfile for the GRISM ASN, with 
    #### same shifts taken from the direct images
    threedhst.shifts.make_grism_shiftfile(asn_direct_file,asn_grism_file)
    
    threedhst.currentRun['step'] = 'SHIFTS'
    
    #########################################
    ####   Make detection image and generate a catalog
    #########################################
    #dither
    
    #### First Multidrizzle run on DIRECT images to create a detection image
    iraf.unlearn('multidrizzle')
    iraf.multidrizzle(input=asn_direct_file, \
       shiftfile=root_direct + '_shifts.txt', \
       output = '', final_scale = INDEF, final_pixfrac = 1.0)
    cleanMultidrizzleOutput()
    threedhst.currentRun['step'] = 'DETECTION_IMAGE'
    
    direct_mosaic = root_direct+'_drz.fits'
    
    #### Get out SCI and WHT extensions of the Multidrizzle mosaic
    iraf.imcopy(input=direct_mosaic+'[1]',output=root_direct+'_SCI.fits')
    iraf.imcopy(input=direct_mosaic+'[2]',output=root_direct+'_WHT.fits')
    #### Run SExtractor on the direct image, with the WHT 
    #### extension as a weight image
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CATALOG_NAME']    = root_direct+'_drz.cat'
    se.options['CHECKIMAGE_NAME'] = root_direct+'_seg.fits'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = root_direct+'_WHT.fits'
    se.options['FILTER']    = 'Y'
    ## Detect thresholds (default = 1.5)
    se.options['DETECT_THRESH']    = str(threedhst.options['DETECT_THRESH']) 
    se.options['ANALYSIS_THRESH']  = str(threedhst.options['ANALYSIS_THRESH']) 
    se.options['MAG_ZEROPOINT'] = str(threedhst.options['MAG_ZEROPOINT'])
    
    ## Run SExtractor
    status = se.sextractImage(root_direct+'_SCI.fits')
    
    threedhst.currentRun['step'] = 'SEXTRACT_CATALOG'
    
    #### Read catalog to keep around
    sexCat = threedhst.sex.mySexCat(root_direct+'_drz.cat')
    threedhst.currentRun['sexCat'] = sexCat
    
    #### Replace MAG_AUTO column name to MAG_F1392W
    sexCat.change_MAG_AUTO_for_aXe(filter='F1392W')
    sexCat.writeToFile()
    #### Make region file for SExtracted catalog and the pointing itself
    threedhst.sex.sexcatRegions(root_direct+'_drz.cat', 
                                root_direct+'_drz.reg', format=2)
    threedhst.utils.asn_region(asn_direct_file)
    
    #### Make zeroth order region file
    x_world = sexCat.columns[sexCat.searchcol('X_WORLD')].entry
    y_world = sexCat.columns[sexCat.searchcol('Y_WORLD')].entry
    a_col = sexCat.columns[sexCat.searchcol('A_WORLD')].entry
    b_col = sexCat.columns[sexCat.searchcol('B_WORLD')].entry
    theta_col = sexCat.columns[sexCat.searchcol('THETA_WORLD')].entry
    asec = 3600.
    pp = '"'
    theta_sign = -1
    
    flt = pyfits.open(asn_direct.exposures[0]+'_flt.fits')
    wcs = pywcs.WCS(flt[1].header)
    world_coord = []
    for i in range(len(x_world)):
        world_coord.append([np.float(x_world[i]),np.float(y_world[i])])
    xy_coord = wcs.wcs_sky2pix(world_coord,0)
    ## conf: XOFF_B -192.2400520   -0.0023144    0.0111089
    for i in range(len(x_world)):
        xy_coord[i][0] += (-192.2400520 - 0.0023144*xy_coord[i][0] +
                            0.0111089*xy_coord[i][1])
    world_coord = wcs.wcs_pix2sky(xy_coord,0)
    fp = open(root_grism+'_zeroth.reg','w')
    fp.write('fk5\n')
    for i in range(len(x_world)):
        line = "ellipse(%s, %s, %6.2f%s, %6.2f%s, %6.2f)\n" %(world_coord[i][0],
              world_coord[i][1], 
              float(a_col[i])*asec, pp,
              float(b_col[i])*asec, pp, float(theta_col[i])*theta_sign)
        fp.write(line)
    fp.close()
    
    threedhst.currentRun['step'] = 'PROCESS_CATALOG'
    
    #########################################
    ####   Run aXe scripts
    #########################################
    #taxe20
    #### initialize parameters
    conf = Conf(threedhst.options['CONFIG_FILE'])
    threedhst.currentRun['conf'] = conf
    
    conf.params['DRZROOT'] = root_direct
    ## Workaround to get 0th order contam. from fluxcube
    conf.params['BEAMB'] = '-220 220'    
    conf.writeto(root_direct+'.conf')
    CONFIG = root_direct+'.conf'
    #CONFIG = 'WFC3.IR.G141.V1.0.conf'
    SKY = threedhst.options['SKY_BACKGROUND']
    
    set_aXe_environment(grating=threedhst.options['GRISM_NAME'])
    
    ## Make 'lis' file for input into aXe
    status = make_aXe_lis(asn_grism_file, asn_direct_file)
    shutil.move(prep_name(asn_grism_file),'..')
    ## taxe20.tiolprep
    flprMulti()
    iraf.iolprep(mdrizzle_ima=root_direct+'_drz.fits',
                 input_cat=root_direct+'_drz.cat')
    ## taxe20.taxeprep
    os.chdir('../')
    flprMulti()
    iraf.axeprep(inlist=prep_name(asn_grism_file), configs=CONFIG,
        backgr="YES", backims=SKY, mfwhm="2.0",norm="NO")
    threedhst.currentRun['step'] = 'AXEPREP'
    
    #### Multidrizzle the grism images to flag cosmic rays
    os.chdir('./DATA')
    flprMulti()
    iraf.multidrizzle(input=asn_grism_file, 
                      shiftfile=root_grism +  '_shifts.txt', 
                      output = '', final_scale = INDEF, 
                      final_pixfrac = 1.0, skysub = no)
    cleanMultidrizzleOutput()
    threedhst.currentRun['step'] = 'MULTIDRIZZLE_GRISM'
    
    #### Prepare fluxcube image
    #### (Force flat spectrum in f_nu)
    fp = open('zeropoints.lis','w')
    fp.writelines([root_direct+'_drz.fits, 1792.3, 26.46\n',
                   root_direct+'_drz.fits, 1392.3, 26.46\n',
                   root_direct+'_drz.fits, 792.3, 26.46\n'])
    fp.close()
    flprMulti()
    iraf.fcubeprep(grism_image = root_grism+'_drz.fits',
       segm_image = root_direct+'_seg.fits',
       filter_info = 'zeropoints.lis', AB_zero = yes, 
       dimension_info = '0,0,0,0')
    
    threedhst.currentRun['step'] = 'FCUBEPREP'
    
    ####################       Main aXe run
    os.chdir('../')
    flprMulti()
    #### taxecore
    iraf.axecore(inlist=prep_name(asn_grism_file), configs=CONFIG,
        back="NO",extrfwhm=4.1, drzfwhm=4, backfwhm=4.1,
        slitless_geom="YES", orient="YES", exclude="NO",lambda_mark=1392.0, 
        cont_model="fluxcube", model_scale=3, lambda_psf=1392.0,
        inter_type="linear", np=10,interp=-1,smooth_lengt=1, smooth_fwhm=1.0,
        spectr="NO", adj_sens="NO", weights="NO", sampling="drizzle")
    
    threedhst.currentRun['step'] = 'AXECORE'
    
    #### tdrzprep
    flprMulti()
    iraf.drzprep(inlist=prep_name(asn_grism_file), configs=CONFIG,
        opt_extr="YES", back="NO")
    
    threedhst.currentRun['step'] = 'DRZPREP'
        
    #### taxedrizzle
    flprMulti()
    iraf.axedrizzle(inlist=prep_name(asn_grism_file), configs=CONFIG,
                    infwhm=4.1,outfwhm=4, back="NO", makespc="YES",
                    opt_extr="YES",adj_sens="YES")
    
    threedhst.currentRun['step'] = 'AXEDRIZZLE'
    
    #### Make a multidrizzled contamination image
    os.chdir('./DATA')
    for expi in asn_grism.exposures:
        ### copy ../OUTPUT_G141/CONT.fits into existing grism FLT files
        flt = pyfits.open(expi+'_flt.fits','update')
        cont = pyfits.open('../OUTPUT_G141/'+expi+'_flt_2.CONT.fits')
        flt[1].data = cont[1].data.copy()
        flt.flush()
    asn_grism.product = None #root_grism+'CONT'
    asn_grism.writeToFile(root_grism+'CONT_asn.fits')
    flprMulti()
    iraf.multidrizzle(input = root_grism+'CONT_asn.fits', output = '', \
       shiftfile = root_grism + '_shifts.txt', \
       final_scale = INDEF, final_pixfrac = 1.0, skysub = no,
       static=no, driz_separate=no,median=no, blot=no, driz_cr=no)
    cleanMultidrizzleOutput()
    
    threedhst.currentRun['step'] = 'MULTIDRIZZLE_CONTAMINATION'
    
    #### Clean up
    rmfiles = glob.glob('*FLX.fits')
    rmfiles.extend(glob.glob('*_flt.fits'))
    rmfiles.extend(glob.glob('*flt_1.cat'))
    rmfiles.extend(glob.glob(root_direct+'*[SW][CH]?.fits'))
    rmfiles.extend(glob.glob('*coeffs1.dat'))
    rmfiles.extend(glob.glob('threedhst_auto.*'))
    rmfiles.extend(glob.glob('zeropoints.lis'))
    rmfiles.extend(glob.glob('default.conv'))
    for expi in asn_grism.exposures:
        rmfiles.extend(glob.glob('../OUTPUT_G141/'+expi+'*'))
         
    if len(rmfiles) > 0:
        for rmfile in rmfiles:
            os.remove(rmfile)

    threedhst.currentRun['step'] = 'CLEANUP_DATA'
    
    #### Make output webpages with spectra thumbnails    
    try:
        os.mkdir('../HTML/scripts')
    except:
        print("""
        WARNING: ../HTML/scripts/ not found.  Download from
                 http://code.google.com/p/threedhst/
              """)
              
    try:
        os.mkdir('../HTML/images')
    except:
        pass
        
    SPC = threedhst.plotting.SPCFile(root_direct+'_2_opt.SPC.fits')
    threedhst.currentRun['SPC'] = SPC
    
    print '\nthreedhst.plotting.makeThumbs: Creating direct image ' + \
          'thumbnails...\n\n'
    threedhst.plotting.makeThumbs(SPC, sexCat, path='../HTML/images/')
    
    print '\nthreedhst.plotting.makeSpec1dImages: Creating 1D spectra '+ \
          'thumbnails...\n\n'
    threedhst.plotting.makeSpec1dImages(SPC, path='../HTML/images/')

    print '\nthreedhst.plotting.makeSpec1dImages: Creating 2D spectra '+ \
          'thumbnails...\n\n'
    threedhst.plotting.makeSpec2dImages(SPC, path='../HTML/images/')
    
    threedhst.currentRun['step'] = 'MAKE_THUMBNAILS'
    
    #### Make tiles for Google map
    try:
        os.mkdir('../HTML/tiles')
    except:
        pass
    
    threedhst.gmap.makeCatXML(catFile=root_direct+'_drz.cat',
                              xmlFile='../HTML/'+root_direct+'.xml')
    
    threedhst.gmap.makeCirclePNG(outfile='../HTML/circle.php')            
                  
    mapParamsD = threedhst.gmap.makeGMapTiles(fitsfile=root_direct+'_drz.fits',
                                             outPath='../HTML/tiles/',
                                             tileroot='direct')
    
    mapParamsG = threedhst.gmap.makeGMapTiles(fitsfile=root_grism+'_drz.fits',
                                             outPath='../HTML/tiles/',
                                             tileroot='grism')
    
    threedhst.currentRun['step'] = 'MAKE_GMAP_TILES'
    
    out_web = '../HTML/'+root_direct+'_index.html'
    print '\nthreedhst.plotting.makeHTML: making webpage: %s\n' %out_web
    threedhst.plotting.makeHTML(SPC, sexCat, mapParamsD, output=out_web)
    
    threedhst.currentRun['step'] = 'MAKE_HTML'
    
    #### Done!
    print 'threedhst: cleaned up and Done!\n'
    
    threedhst.currentRun['step'] = 'DONE'
    

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
3D-HST / make_aXe_lis: Number of grism exposures (%d) in %s is different from
                       the number of direct images (%d) in %s.
              """ %(len(grism_files), asn_grism, len(direct_files), asn_direct)
        return False
    
    NFILE = len(asn_grism.exposures)
    outfile = prep_name(asn_grism_file)
    fp = open(outfile,'w')
    for i in range(NFILE):
        fp.write("%s_flt.fits %s_flt_1.cat %s_flt.fits 0.0\n" 
              %(asn_grism.exposures[i], asn_direct.exposures[i], 
                asn_direct.exposures[i]))
    fp.close()
    print "3D-HST / make_aXe_lis: Created %s\n" %outfile
    return True

def cleanMultidrizzleOutput():
    """
    cleanMultidrizzleOutput()
    
    Remove *single_[sci/wht].fits, *sci1_blt.fits, *flt*mask1.fits, *coeffs1.dat
    """
    import os,glob
    rmfiles = glob.glob('*single_???.fits')
    rmfiles.extend(glob.glob('*sci1_blt.fits'))
    rmfiles.extend(glob.glob('*flt*mask1.fits'))
    #files.extend(glob.glob('*coeffs1.dat'))
    rmfiles.extend(glob.glob('ib*.run'))
    rmfiles.extend(glob.glob('ib*_med.fits'))
    if len(rmfiles) > 0:
        for rmfile in rmfiles:
            os.remove(rmfile)

def flprMulti(n=3):
    """
    flprMulti(n=3)
    
    Run iraf.flpr() ``n`` times.
    """
    if n < 1:
        n=1
    for i in range(n):
        iraf.flpr()

def multidrizzle_defaults():
    multidrizzle(input = 'ib3714060_asn.fits', output = '', mdriztab = no, \
       refimage = '', runfile = '', workinplace = no, updatewcs = yes, \
       proc_unit = 'native', coeffs = 'header', context = no, clean = no, \
       group = '', ra = INDEF, dec = INDEF, build = yes, shiftfile = '', 
       staticfile = '', \
       static = yes, static_sig = 4.0, skysub = yes, skywidth = 0.1, \
       skystat = 'median', skylower = INDEF, skyupper = INDEF, skyclip = 5, \
       skylsigma = 4.0, skyusigma = 4.0, skyuser = '', driz_separate = yes, \
       driz_sep_outnx = INDEF, driz_sep_outny = INDEF, 
       driz_sep_kernel = 'turbo', \
       driz_sep_wt_scl = 'exptime', driz_sep_scale = INDEF, 
       driz_sep_pixfrac = 1.0, \
       driz_sep_rot = INDEF, driz_sep_fillval = 'INDEF', driz_sep_bits = 0, \
       median = yes, median_newmasks = yes, combine_maskpt = 0.7, \
       combine_type = 'minmed', combine_nsigma = '4 3', combine_nlow = 0, \
       combine_nhigh = 1, combine_lthresh = 'INDEF', combine_hthresh = 'INDEF',\
       combine_grow = 1, blot = yes, blot_interp = 'poly5', blot_sinscl = 1.0, \
       driz_cr = yes, driz_cr_corr = no, driz_cr_snr = '3.5 3.0', \
       driz_cr_grow = 1, driz_cr_ctegrow = 0, driz_cr_scale = '1.2 0.7', \
       driz_combine = yes, final_wht_type = 'EXP', final_outnx = INDEF, \
       final_outny = INDEF, final_kernel = 'square', final_wt_scl = 'exptime', \
       final_scale = INDEF, final_pixfrac = 1.0, final_rot = 0.0, \
       final_fillval = 'INDEF', final_bits = 0, final_units = 'cps', \
       gain = INDEF, gnkeyword = INDEF, rdnoise = INDEF, rnkeyword = INDEF,
       exptime = INDEF)
    
class Conf(object):
    """
    Conf(infile='WFC3.IR.G141.V1.0.conf')
    
    Read an aXe configuration file and easily change parameters.
    
    """            
    def _getPath(self):
        """
        _getPath()
        
        Figure out if we're in the root directory or in DATA
        """
        import os
        if os.getcwd().split('/')[-1] == 'DATA':
            self.path = '../CONF/'
        else:
            self.path = './CONF/'
        
    def _processLines(self):
        self.nlines = len(self.lines)
        self.params = {}
        self._pline = {}
        for i,line in enumerate(self.lines):
            if (line[0] is not '#') & (line.strip() is not  ''):
                spl = line.split()
                self.params[spl[0]] = ' '.join(spl[1:])
                self._pline[spl[0]] = i
        self.nkeys = self.params.keys().__len__()
        
    def _readLines(self):
        """
        _readlines()
        """
        self._getPath()
        fp = open(self.path+self.infile,'r')
        self.lines = fp.readlines()
        fp.close()
        
    def __init__(self, infile='WFC3.IR.G141.V1.0.conf'):
        self.infile = infile
        self._getPath()
        self._readLines()
        self._processLines()
        
    def _assignPars(self):
        """
        _assignPars()
        
        Apply parameter values to self.lines
        """
        for key in self.params.keys():
            self.lines[self._pline[key]] = key + ' ' + self.params[key]+'\n'
    
    def writeto(self, output='tmp.conf'):
        """
        writeto(output='tmp.conf')
        """ 
        self._getPath()
        self._assignPars()
        fp = open(self.path+output,'w')
        fp.writelines(self.lines)
        fp.close()
    
    

    