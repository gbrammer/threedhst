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

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import numpy as np

import threedhst

try:
    import drizzlepac
    drizzle_function = drizzlepac.astrodrizzle.AstroDrizzle
except:
    import pyraf
    from pyraf import iraf
    drizzle_function = iraf.multidrizzle
    
def set_aXe_environment(grating='G141'):
    """
set_aXe_environment(grating='G141')
    
Setup aXe environment variables:

    AXE_IMAGE_PATH   = ./DATA
    AXE_CONFIG_PATH  = ./CONF
    AXE_OUTPUT_PATH  = ./OUTPUT_[GRATING]
    AXE_DRIZZLE_PATH = ./DRIZZLE_[GRATING]
    
CONF can be symlinked from e.g. /research/HST/GRISM/CONF
    
    """
    os.environ['AXE_IMAGE_PATH'] = './DATA/'
    print('--> variable AXE_IMAGE_PATH   set to "./DATA"')
    
    os.environ['AXE_CONFIG_PATH'] = './CONF/'
    print('--> variable AXE_CONFIG_PATH  set to "./CONF/"')
     
    os.environ['AXE_OUTPUT_PATH'] = './OUTPUT_'+grating.upper()+'/'
    print('--> variable AXE_OUTPUT_PATH  set to "./OUTPUT_'+grating.upper()+'/"')
    
    os.environ['AXE_DRIZZLE_PATH'] = './DRIZZLE_'+grating.upper()+'/'
    print('--> variable AXE_DRIZZLE_PATH set to' + \
               '"./DRIZZLE_'+grating.upper()+'/"')


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
    
    if not threedhst.options['ACS_G800L']:
        if not os.path.exists('CONF/'+threedhst.options['CONFIG_FILE']):
            raise IOError("options['CONFIG_FILE'] doesn't exist: CONF/" +
                          threedhst.options['CONFIG_FILE'])

        if threedhst.options['SKY_BACKGROUND'] is not None:
            if not os.path.exists('CONF/'+threedhst.options['SKY_BACKGROUND']):
                raise IOError("options['SKY_BACKGROUND'] doesn't exist:" +   
                          "CONF/"+threedhst.options['SKY_BACKGROUND'])

def write_to_log(lines, init=False):
    """
write_to_log(lines, init=False)
    
    Write output of IRAF tasks to options['ROOT_GRISM']+'.iraf.log'
    """
    if init:
        fp = open(threedhst.options['ROOT_GRISM']+'.iraf.log','w')
    else:
        fp = open(threedhst.options['ROOT_GRISM']+'.iraf.log','a')
        
    fp.writelines(lines)
    fp.close()
        
################################################################################
####    This is the main routine for running the grism data reduction
################################################################################
def reduction_script(asn_grism_file=None, asn_direct_file=None):
    """
process_grism(asn_grism, asn_direct)
    
Pipeline to process a set of grism/direct exposures.
    
    """
    import shutil
    import glob
    
    import pyraf
    import aXe2html.sexcat.sextractcat
    
    from pyraf import iraf
    from iraf import stsdas,dither,slitless,axe 
    
    no = iraf.no
    yes = iraf.yes
    INDEF = iraf.INDEF
    
    #import threedhst.dq
    if threedhst.options['USE_TAXE']:
        try:
            from iraf import taxe21
        except:
            threedhst.showMessage("USE_TAXE is True, but no iraf.taxe21 found.", warn=True)
            return False
            
    ############################################################################
    ####   Bookkeeping, set up DATA directory
    ############################################################################
    
    #### Check that we're in the home directory of a 3D-HST field
    check_3dhst_environment(makeDirs=False)
    
    os.chdir('./DATA')
    
    if not asn_grism_file:
        asn_grism_file  = 'ib3721060_asn.fits'
    if not asn_direct_file:
        asn_direct_file = 'ib3721050_asn.fits'
            
    #### ASN root names
    ROOT_GRISM = asn_grism_file.split('_asn.fits')[0]
    threedhst.currentRun['asn_grism_file'] = asn_grism_file
    threedhst.options['ROOT_GRISM'] = ROOT_GRISM

    #### Copy ASN files from RAW, if they don't already exist
    if not os.path.exists(asn_grism_file):
        shutil.copy(threedhst.options['PATH_TO_RAW']+asn_grism_file,'./')
    
    if threedhst.options['PREFAB_DIRECT_IMAGE'] is not None:
        ROOT_DIRECT = os.path.basename(
             threedhst.options['PREFAB_DIRECT_IMAGE'].split('_drz.fits')[0])
        threedhst.options['ROOT_DIRECT'] = ROOT_DIRECT     
    else:
        ROOT_DIRECT = asn_direct_file.split('_asn.fits')[0]
        threedhst.currentRun['asn_direct_file'] = asn_direct_file
        threedhst.options['ROOT_DIRECT'] = ROOT_DIRECT
        if not os.path.exists(asn_direct_file):
            shutil.copy(threedhst.options['PATH_TO_RAW']+asn_direct_file,'./')
    
    #print "\n3DHST.process_grism: Fresh ASN and FLT files...\n"
    threedhst.showMessage("Fresh ASN and FLT files...")
    
    #### Read ASN files
    asn_grism  = threedhst.utils.ASNFile(file=asn_grism_file)
    ### Force PROD-DTH to be the root of the input ASN file
    asn_grism.product = ROOT_GRISM
    asn_grism.write(asn_grism_file)
    #### Copy fresh FLT files from RAW
    threedhst.process_grism.fresh_flt_files(asn_grism_file, 
                      from_path=threedhst.options['PATH_TO_RAW'])
    
    #### Do it now for the direct images
    if threedhst.options['PREFAB_DIRECT_IMAGE'] is None:
        asn_direct = threedhst.utils.ASNFile(file=asn_direct_file)
        asn_direct.product = ROOT_DIRECT
        asn_direct.write(asn_direct_file)                            
        threedhst.process_grism.fresh_flt_files(asn_direct_file, 
                      from_path=threedhst.options['PATH_TO_RAW'])
            
    threedhst.currentRun['step'] = 'COPY_FROM_RAW'
    
    #########################################################################
    ####   Compute shifts
    #########################################################################
    if os.path.exists(ROOT_DIRECT + '_shifts.txt') is False:
        threedhst.options['MAKE_SHIFTFILES'] = True
    
    #### If a direct image is provided, assume shiftfiles exist
    if threedhst.options['PREFAB_DIRECT_IMAGE'] is not None:
        threedhst.options['MAKE_SHIFTFILES'] = False
        if not os.path.exists(asn_grism_file.split('_asn')[0]+'_shifts.txt'):
            threedhst.showMessage('No grism shiftfile found, making an empty one...', warn=True)
            threedhst.shifts.make_blank_shiftfile(asn_file=asn_grism_file)
            
    if threedhst.options['MAKE_SHIFTFILES']:
        #print "\n3DHST.process_grism: Getting initial shifts...\n"
        threedhst.showMessage("Getting initial shifts...")
        
        #### Run iraf.tweakshifts on the diret images
        threedhst.shifts.run_tweakshifts(asn_direct_file)
        
        #### Check to make sure that every exposure in the ASN file
        #### has a line in the shiftfile.  If something goes wrong, 
        #### tweakshifts will usually omit an exposure from the shiftfile
        threedhst.shifts.checkShiftfile(asn_direct_file)
        
        #### Make a shiftfile for the GRISM ASN, with 
        #### same shifts taken from the direct images
        threedhst.shifts.make_grism_shiftfile(asn_direct_file, asn_grism_file)
        
        threedhst.currentRun['step'] = 'SHIFTS'
        
    #########################################################################
    ####   Make detection image and generate a catalog
    #########################################################################
    
    #### First Multidrizzle run on DIRECT images to create a detection image
    #print "\n3DHST.process_grism: Creating MultiDrizzled detection image\n"
    threedhst.options['USED_PIXFRAC'] = threedhst.options['PIXFRAC']
    
    if threedhst.options['PREFAB_DIRECT_IMAGE'] is None:
    
        threedhst.showMessage('Creating MultiDrizzled detection image.')

        iraf.unlearn('multidrizzle')

        if len(asn_direct.exposures) < threedhst.options['NFRAME_FOR_PIXFRAC']:
            threedhst.options['USED_PIXFRAC']=1.0
        else:
            threedhst.options['USED_PIXFRAC'] = threedhst.options['PIXFRAC']

        if threedhst.options['SKY_BACKGROUND'] is None:
            MDRIZ_SKYSUB=no
        else:
            MDRIZ_SKYSUB=yes

        #threedhst.process_grism.multidrizzle_defaults(asn_direct_file)
        
        try:
            status = drizzle_function(input=asn_direct_file, \
            shiftfile=ROOT_DIRECT + '_shifts.txt', skysub=MDRIZ_SKYSUB, \
            output = '', final_scale = threedhst.options['DRZSCALE'], 
            final_pixfrac = threedhst.options['USED_PIXFRAC'], Stdout=1)
        except:
            status = drizzle_function(input=asn_direct_file, \
            skysub=MDRIZ_SKYSUB, \
            output = '', final_scale = threedhst.options['DRZSCALE'], 
            final_pixfrac = threedhst.options['USED_PIXFRAC'], build=False)
                
        ### Read the '.run' file output by MultiDrizzle
        #mdrzRun = threedhst.prep_flt_files.MultidrizzleRun(root=ROOT_DIRECT)

        cleanMultidrizzleOutput()
        threedhst.currentRun['step'] = 'DETECTION_IMAGE'
    
        #############################################
        #### Align to reference image
        #### e.g. threedhst.options['ALIGN_IMAGE'] = '../ACS/h_sz*drz_img.fits'
        #############################################
        if threedhst.options['ALIGN_IMAGE']:
          #### Have to run this twice because I'm note sure 
          #### how to translate the reference point for rotation
          #### into the FLT frame
          for geom in ['rxyscale','shift']:  

            threedhst.shifts.refine_shifts(ROOT_DIRECT=ROOT_DIRECT, 
                      ALIGN_IMAGE=threedhst.options['ALIGN_IMAGE'], 
                      fitgeometry=geom, clean=True)
            
            threedhst.showMessage('MultiDrizzle detection image, refined shifts')
            
            threedhst.process_grism.multidrizzle_defaults(asn_direct_file)
            
            status = drizzle_function(input=asn_direct_file, 
               shiftfile=ROOT_DIRECT + '_shifts.txt', skysub=MDRIZ_SKYSUB, 
               output = '', final_scale = threedhst.options['DRZSCALE'],
               final_pixfrac = threedhst.options['USED_PIXFRAC'],
               Stdout=1, updatewcs=iraf.no)
            
            cleanMultidrizzleOutput()
          
          #### Make the grism shiftfile with the same shifts as 
          #### for the direct images
          threedhst.shifts.make_grism_shiftfile(asn_direct_file, asn_grism_file)
          threedhst.currentRun['step'] = 'ALIGN_TO_REFERENCE'
        
        direct_mosaic = ROOT_DIRECT+'_drz.fits'
    else:
        direct_mosaic = threedhst.options['PREFAB_DIRECT_IMAGE']
    
    threedhst.options['DIRECT_MOSAIC'] = direct_mosaic
    
    #############################################
    #### Make a catalog with SExtractor
    #############################################
    
    #### If a catalog is not specified in 'FORCE_CATALOG', 
    #### make one with SExtractor
    if threedhst.options['FORCE_CATALOG'] is None:
    
        threedhst.showMessage('Making the catalog')

        #### Get out SCI and WHT extensions of the Multidrizzle mosaic, 
        #### SExtractor was choking on multi-extension FITS files
        for ext in ['SCI','WHT']:
            try:
                os.remove(ROOT_GRISM+'_'+ext+'.fits')
            except:
                pass
        
        #print 'IMCOPY: '+ROOT_GRISM, os.getcwd()
        os.chdir('../DATA')
        iraf.imcopy(direct_mosaic+'[1]','../DATA/'+ROOT_GRISM+'_SCI.fits')
        iraf.imcopy(direct_mosaic+'[2]','../DATA/'+ROOT_GRISM+'_WHT.fits')
            
        #### Run SExtractor on the direct image, with the WHT 
        #### extension as a weight image
        se = threedhst.sex.SExtractor()

        ## Set the output parameters required for aXe 
        ## (stored in [threedhst source]/data/aXe.param) 
        se.aXeParams()

        ## XXX add test for user-defined .conv file
        se.copyConvFile()

        se.overwrite = True
        se.options['CATALOG_NAME']    = ROOT_GRISM+'_drz.cat'
        se.options['CHECKIMAGE_NAME'] = ROOT_GRISM+'_seg.fits'
        se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
        se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
        se.options['WEIGHT_IMAGE']    = ROOT_GRISM+'_WHT.fits'
        se.options['FILTER']    = 'Y'

        #### Detect thresholds (default = 1.5)
        se.options['DETECT_THRESH']    = str(threedhst.options['DETECT_THRESH']) 
        se.options['ANALYSIS_THRESH']  = str(threedhst.options['ANALYSIS_THRESH']) 
        se.options['MAG_ZEROPOINT'] = str(threedhst.options['MAG_ZEROPOINT'])
    
        #### Run SExtractor
        status = se.sextractImage(ROOT_GRISM+'_SCI.fits')

        threedhst.currentRun['step'] = 'SEXTRACT_CATALOG'

        #### Read catalog to keep around
        sexCat = threedhst.sex.mySexCat(ROOT_GRISM+'_drz.cat')

        #### Trim faint sources from the catalog (will still be in the seg. image)
        mag = np.cast[float](sexCat.MAG_AUTO)
        q = np.where(mag > threedhst.options['LIMITING_MAGNITUDE'])[0]
        if len(q) > 0:
            threedhst.showMessage('Trimming objects with direct M < %6.2f.' 
                                  %(threedhst.options['LIMITING_MAGNITUDE']))
            numbers = np.cast[int](sexCat.NUMBER)[q]
            sexCat.popItem(numbers)

        #### Replace MAG_AUTO column name to MAG_F1392W
        threedhst.options['FILTWAVE'] = np.float(threedhst.options['FILTWAVE'])
        sexCat.change_MAG_AUTO_for_aXe(filter='F%0dW'
                    %(np.round(threedhst.options['FILTWAVE'])))
    
    else:
        #### Assumes that the MAG_AUTO column has already been changed 
        #### appropriately to something like MAG_F1392W
        sexCat = threedhst.sex.mySexCat(threedhst.options['FORCE_CATALOG'])
        col_name = 'MAG_F%0dW' %(np.round(threedhst.options['FILTWAVE']))
        mag = np.cast[float](sexCat[col_name])
        q = np.where(mag > threedhst.options['LIMITING_MAGNITUDE'])[0]
        if len(q) > 0:
            threedhst.showMessage('Trimming objects with direct M < %6.2f.' 
                                  %(threedhst.options['LIMITING_MAGNITUDE']))
            numbers = np.cast[int](sexCat.NUMBER)[q]
            sexCat.popItem(numbers)
        ##### Note: segmentation image has to accompany the force_catalog!
        
    sexCat.write(ROOT_GRISM+'_drz.cat')
    
    ##### Subset of galaxies for z8.6 test
    #sexCat = threedhst.sex.mySexCat('z86.cat')
    #sexCat.write(outfile=ROOT_DIRECT+'_drz.cat')
    
    #### Trim objects on the edge of the detection image whose
    #### spectra will fall off of the grism image [turned off]
    #threedhst.regions.trim_edge_objects(sexCat)
    
    threedhst.currentRun['sexCat'] = sexCat
    
    #### Make region file for SExtracted catalog
    threedhst.sex.sexcatRegions(ROOT_GRISM+'_drz.cat', 
                                ROOT_GRISM+'_drz.reg', format=2)
    
    #### Make a region file for the pointing itself
    #if threedhst.options['PREFAB_DIRECT_IMAGE'] is None: 
    threedhst.regions.asn_region(asn_grism_file)
    
    #### Make zeroth order region file [turned off]
    #threedhst.regions.make_zeroth(sexCat, outfile=ROOT_GRISM+'_zeroth.reg')
        
    threedhst.currentRun['step'] = 'PROCESS_CATALOG'
    
    ############################################################################
    ####   Run aXe scripts
    ############################################################################

    #### Initialize parameters, update the config file in CONF
    chip = 1
    for conf_file in threedhst.options['CONFIG_FILE'].split(','):
        conf = Conf(conf_file)
        threedhst.currentRun['conf'] = conf

        #### Need to scale the 0th order sensitivity curve
        # if conf.params['SENSITIVITY_B'] == 'wfc3_abscal_IRg141_0th_sens.fits':
        zeroth_list = ['wfc3_abscal_IRg141_0th_sens.fits', 'WFC3.IR.G141.0th.sens.1.fits']
        if conf.params['SENSITIVITY_B'] in zeroth_list:
            zeroth_file = pyfits.open(conf.path+'/'+conf.params['SENSITIVITY_B'])
            zeroth_data = zeroth_file[1].data
            sens = zeroth_data.field('SENSITIVITY')
            err = zeroth_data.field('ERROR')
            scale_factor = 3.6
            sens *= scale_factor
            err *= scale_factor
            zeroth_file.writeto(conf.path+'/'+'WFC3_G141_0th_SCALED.fits',
                                clobber=True)
            conf.params['SENSITIVITY_B'] = 'WFC3_G141_0th_SCALED.fits'

        ##### Parameters for aXe
        conf.params['DRZROOT'] = ROOT_GRISM
        conf.params['DRZRESOLA'] = threedhst.options['DRZRESOLA']
        conf.params['DRZSCALE'] = threedhst.options['DRZSCALE']
        conf.params['DRZPFRAC'] = threedhst.options['PIXFRAC']
        
        ## Try expanding the SMFACTOR to account for different pixel scales
        ## in the sensitivity smoothing.  Bug in aXe???
        if threedhst.options['GRISM_NAME'] == 'G141':
            conf.params['SMFACTOR'] = '%.3f' %(0.128254/np.float(threedhst.options['DRZSCALE']))
        
        #### Make sure 4096 is set in CONF.DQMASK, which is the value 
        #### DQ uses for flagging problem regions
        #print conf.params['DQMASK']
        
        conf.params['DQMASK'] = np.str(np.int(conf.params['DQMASK'].split()[0]) | 4096 | 2048)
                
        #### Workaround to get 0th order contam. in the right place for the fluxcube
        if threedhst.options['CONFIG_FILE'] == 'WFC3.IR.G141.V1.0.conf':
            conf.params['BEAMB'] = '-220 220'    
        
        #
        if threedhst.options['CONFIG_FILE'].startswith('WFC3.IR'):
            #conf.params['YOFF_A'] = '2.0'
            pass
            
        conf.writeto(ROOT_GRISM+'_%0d_full.conf' %chip)
        
        chip+=1
        
    if threedhst.options['ACS_G800L']:
        #### The order here is OK because the output order is the same as the 
        #### input.  Though note that now "_1_full.conf" corresponds to 
        #### WFC/Chip2 because CONFIG is supplied as "Chip2,Chip1"
        CONFIG = ROOT_GRISM+'_1_full.conf,'+ROOT_GRISM+'_2_full.conf'
    else:
        CONFIG = ROOT_GRISM+'_1_full.conf'
    
    #CONFIG = 'WFC3.IR.G141.V1.0.conf'
    if threedhst.options['SKY_BACKGROUND'] is None:
        BACKGR = "NO"
        SKY=""
    else:
        BACKGR = "YES"
        SKY = threedhst.options['SKY_BACKGROUND']
        if SKY == 'WFC3.IR.G141.sky.V1.0.fits':
            threedhst.process_grism.fix_g141_sky_image(conf_path=conf.path,
                verbose=True)
            SKY = 'g141_fixed_sky.fits'
            
    #### Multidrizzle the grism images to flag cosmic rays
    if threedhst.options['PREFAB_DIRECT_IMAGE'] is not None:
        threedhst.options['CATALOG_MODE'] = 'GRISM'
    
    #### If you want to use the grism images for making the object catalogs
    #### (CATALOG_MODE != 'DIRECT'), then might need to prepare a grism mosaic 
    #### here.
    redo_second = iraf.yes
    if ((threedhst.options['PREFAB_GRISM_IMAGE'] is None) & 
       (threedhst.options['CATALOG_MODE'] != 'DIRECT')):
        
        redo_second = iraf.no
        flprMulti()
        threedhst.showMessage('Running MultiDrizzle on the grism exposures to flag cosmic rays\nand make grism mosaic.')
        
        if threedhst.options['PREFAB_DIRECT_IMAGE'] is None:
            drizzle_function(input=asn_grism_file, 
                          shiftfile=ROOT_GRISM + '_shifts.txt', 
                          output = '', 
                          final_scale = threedhst.options['DRZSCALE'], 
                          final_pixfrac = threedhst.options['USED_PIXFRAC'], 
                          skysub = BACKGR, static=BACKGR, updatewcs=BACKGR,
                          driz_separate=BACKGR, median=BACKGR,
                          blot=BACKGR, driz_cr=BACKGR)
        else:
            #### Match the grism mosaic WCS to the direct image
            # header = pyfits.getheader(threedhst.options['PREFAB_DIRECT_IMAGE'],1)
            # NX = header.get('NAXIS1')
            # NY = header.get('NAXIS2')
            # RA_CENTER = header.get('CRVAL1')
            # DEC_CENTER = header.get('CRVAL2')
            # #SCALE = np.abs(header.get('CD1_1'))*3600.
            # SCALE = np.sqrt((header.get('CD1_1')*3600.)**2+
            #                 (header.get('CD1_2')*3600.)**2)
            # 
            # ANGLE = np.arctan2(header.get('CD2_1'),
            #                    header.get('CD2_2'))/2/np.pi*360.
            # 
            # print 'NX NY RA DEC SCALE'
            # print '%d %d %13.6f %13.6f %6.3f' %(NX, NY, RA_CENTER, DEC_CENTER, SCALE)
            # iraf.multidrizzle(input=asn_grism_file, 
            #               shiftfile=ROOT_GRISM + '_shifts.txt', 
            #               output = '', final_scale = SCALE, 
            #               final_pixfrac = threedhst.options['USED_PIXFRAC'], 
            #               skysub = BACKGR, final_outnx=NX, final_outny=NY, 
            #               final_rot=0., ra=RA_CENTER, dec=DEC_CENTER,
            #               static=BACKGR, updatewcs=BACKGR,
            #               driz_separate=BACKGR, median=BACKGR,
            #               blot=BACKGR, driz_cr=BACKGR)
            drizzle_function(input=asn_grism_file, 
                          shiftfile=ROOT_GRISM + '_shifts.txt', 
                          output = '', 
                          final_pixfrac = threedhst.options['USED_PIXFRAC'], 
                          skysub = BACKGR, 
                    refimage=threedhst.options['PREFAB_DIRECT_IMAGE']+'[1]',
                          static=BACKGR, updatewcs=BACKGR,
                          driz_separate=BACKGR, median=BACKGR,
                          blot=BACKGR, driz_cr=BACKGR)
        
        
        cleanMultidrizzleOutput()
    threedhst.currentRun['step'] = 'MULTIDRIZZLE_GRISM'
    
    #### Set the aXe environment variables
    set_aXe_environment(grating=threedhst.options['GRISM_NAME'])
    
    #### Make 'lis' file for input into aXe
    status = make_aXe_lis(asn_grism_file, asn_direct_file,
                          mode=threedhst.options['CATALOG_MODE'])
    try:
        os.remove('../'+prep_name(asn_grism_file))
    except:
        pass
    shutil.move(prep_name(asn_grism_file),'..')
    
    #############################################
    #### Run aXe.iolprep to process the catalog as needed [[ in DATA directory]]
    #############################################
    
    threedhst.showMessage('Running aXe.iolprep')    
    IOL_IMAGE = direct_mosaic
    if threedhst.options['CATALOG_MODE'] != 'DIRECT':
        if threedhst.options['PREFAB_GRISM_IMAGE'] is None:
            IOL_IMAGE = ROOT_GRISM+'_drz.fits'
        else:
            IOL_IMAGE = threedhst.options['PREFAB_GRISM_IMAGE']
    
    # threedhst.showMessage('IOL_IMAGE: %s' %(IOL_IMAGE), warn=True)
    flprMulti()
    if threedhst.options['USE_TAXE']:
        taxe21.tiolprep(mdrizzle_ima=IOL_IMAGE,
            input_cat=ROOT_GRISM+'_drz.cat', 
            dimension_in=threedhst.options["AXE_EDGES"])
    else:
        iraf.iolprep(mdrizzle_ima=IOL_IMAGE,
            input_cat=ROOT_GRISM+'_drz.cat', 
            dimension_in=threedhst.options["AXE_EDGES"])
                 
    threedhst.currentRun['step'] = 'IOLPREP'
    
    #############################################
    #### In root directory, run aXe.axeprep to subtract the sky background
    #### from the grism images
    #############################################
    os.chdir('../')
    flprMulti()
    threedhst.showMessage('Running aXe.axeprep')
    if threedhst.options['USE_TAXE']:
        status = taxe21.taxeprep(inlist=prep_name(asn_grism_file),
            configs=CONFIG, backgr=BACKGR, backims=SKY, mfwhm=3.0,norm="NO",
            Stdout=1)
    else:
        status = iraf.axeprep(inlist=prep_name(asn_grism_file), configs=CONFIG,
            backgr=BACKGR, backims=SKY, mfwhm=3.0,norm="NO", Stdout=1)
        
    threedhst.currentRun['step'] = 'AXEPREP'
    
    os.chdir('./DATA')
    
    #### If you used aXe to subtract the background from the FLT files, need
    #### to remake the grism mosaic
    if (threedhst.options['PREFAB_GRISM_IMAGE'] is None) & (BACKGR == "YES"):
        flprMulti()
        threedhst.showMessage('Second pass Multidrizzle on sky-subtracted grism exposures.')
        
        if threedhst.options['PREFAB_DIRECT_IMAGE'] is None:
            status = drizzle_function(input=asn_grism_file, 
                          shiftfile=ROOT_GRISM +  '_shifts.txt', 
                          output = '', 
                          final_scale = threedhst.options['DRZSCALE'], 
                          final_pixfrac = threedhst.options['USED_PIXFRAC'], 
                          skysub = no, static=redo_second,
                          driz_separate=redo_second, median=redo_second,
                          blot=redo_second, driz_cr=redo_second, Stdout=1)
        else:
            #### Match the grism mosaic WCS to the direct image
            # header = pyfits.getheader(threedhst.options['PREFAB_DIRECT_IMAGE'],1)
            # NX = header.get('NAXIS1')
            # NY = header.get('NAXIS2')
            # #SCALE = np.abs(header.get('CD1_1'))*3600.
            # SCALE = np.sqrt((header.get('CD1_1')*3600.)**2+
            #                 (header.get('CD1_2')*3600.)**2)
            # 
            # ANGLE = np.arctan2(header.get('CD2_1'),
            #                    header.get('CD2_2'))/2/np.pi*360.
            # 
            # print 'NX NY RA DEC SCALE'
            # print '%d %d %13.6f %13.6f %6.3f' %(NX, NY, RA_CENTER, DEC_CENTER, SCALE)
            
            # status = drizzle_function(input=asn_grism_file, 
            #               shiftfile=ROOT_GRISM +  '_shifts.txt', 
            #               output = '', final_scale = SCALE, 
            #               final_pixfrac = threedhst.options['USED_PIXFRAC'], 
            #               skysub = no, final_outnx=NX, final_outny=NY, 
            #               ra=RA_CENTER, dec=DEC_CENTER, final_rot=0., 
            #               static=redo_second, 
            #               driz_separate=redo_second,median=redo_second, 
            #               blot=redo_second, driz_cr=redo_second,
            #               Stdout=1)
            status = drizzle_function(input=asn_grism_file, 
                          shiftfile=ROOT_GRISM +  '_shifts.txt', 
                          output = '',
                          final_pixfrac = threedhst.options['USED_PIXFRAC'], 
                          skysub = no, static=redo_second, 
                          driz_separate=redo_second,median=redo_second, 
                          blot=redo_second, driz_cr=redo_second,
                     refimage=threedhst.options['PREFAB_DIRECT_IMAGE']+'[1]',
                          Stdout=1)
        
        cleanMultidrizzleOutput()
        
    #############################################
    #### Prepare the fluxcube for the contamination model
    #############################################
    swarpOtherBands()
    mag_zeropoint='%5.2f' %(np.cast[float](threedhst.options['MAG_ZEROPOINT']))
    
    if threedhst.options['OTHER_BANDS']:
        #### If OTHER_BANDS are available, use them for the fluxcube, start
        #### with the F140W image.
        lines = [direct_mosaic+', %6.1f, %s\n' %(threedhst.options['FILTWAVE'], mag_zeropoint)]
        for band in threedhst.options['OTHER_BANDS']:
            ## 'band' like ['../ACS/h_sz*drz_img.fits','F850LP',903.,24.84]
            if len(band) == 4:
                lines.append(ROOT_GRISM+'_%s_drz.fits,' %band[1] +
                             ' %6.1f, %5.2f\n' %(band[2],band[3]))
    else:
        #### If no other bands present for the fluxcube, force flat spectrum
        #### in f_nu (constant AB mag)

        ## Make dummy copies of the direct image for the fluxcube
        shutil.copy(direct_mosaic,'flux1.fits')
        shutil.copy(direct_mosaic,'flux2.fits')
        lines = ['flux1.fits, %6.1f, %s\n' %(threedhst.options['FILTWAVE']+50,mag_zeropoint),
                   direct_mosaic+', %6.1f, %s\n' %(threedhst.options['FILTWAVE'], mag_zeropoint)] #,
#                   'flux2.fits, %6.1f, %s\n' %(threedhst.options['FILTWAVE']-50,mag_zeropoint)]
    
    #### Make a 'zeropoints.lis' file needed for fcubeprep   
    fp = open('zeropoints.lis','w')
    fp.writelines(lines)
    fp.close()
    
    #### Run aXe.fcubeprep
    # print "\n3DHST.proces_grism: iraf.FCUBEPREP\n"
    threedhst.showMessage('Prepare FLUXCUBE images with iraf.fcubeprep')
    
    FLX_files = glob.glob('*FLX.fits')
    for file in FLX_files:
        os.remove(file)
        
    flprMulti()
    # iraf.fcubeprep(grism_image = ROOT_GRISM+'_drz.fits',
    #    segm_image = ROOT_DIRECT+'_seg.fits',
    #    filter_info = 'zeropoints.lis', AB_zero = yes, 
    #    dimension_info = '400,400,0,0', interpol="poly5")

    if threedhst.options['USE_TAXE']:
        taxe21.tfcubeprep(grism_image = ROOT_GRISM+'_drz.fits',
            segm_image = ROOT_GRISM+'_seg.fits',
            filter_info = 'zeropoints.lis', AB_zero = yes, 
            dimension_info =threedhst.options["AXE_EDGES"], interpol="poly5")
    else:
        iraf.fcubeprep(grism_image = ROOT_GRISM+'_drz.fits',
            segm_image = ROOT_GRISM+'_seg.fits',
            filter_info = 'zeropoints.lis', AB_zero = yes, 
            dimension_info =threedhst.options["AXE_EDGES"], interpol="poly5")
    
    #### Try matching the FLX WCS to the grism images, could need to add shifts
    #threedhst.process_grism.update_FLX_WCS(path_to_FLT='./')
    
    #### Need to get this file out to avoid permissions conflicts
    try:
        sys_result = os.system('rm tmp_coeffs1.dat')
    except:
        pass
        
    threedhst.currentRun['step'] = 'FCUBEPREP'
    
    #############################################
    ####  aXe.axecore to produce the PETs and spectral models
    #############################################
    os.chdir('../')
    # print "\n3DHST.proces_grism: iraf.AXECORE\n"
    threedhst.showMessage('Running iraf.axecore')
    flprMulti()
    
    if threedhst.options['FULL_EXTRACTION_GEOMETRY'] is True:
        geom_yn = "YES"
    else:
        geom_yn = "NO"
        
    # iraf.axecore(inlist=prep_name(asn_grism_file), configs=CONFIG,
    #     back="NO",extrfwhm=5.0, drzfwhm=4.0, backfwhm=0.0,
    #     slitless_geom=geom_yn, orient=geom_yn, exclude="NO", 
    #     lambda_mark=threedhst.options['FILTWAVE'], 
    #     cont_model="fluxcube", model_scale=4.0, 
    #     lambda_psf=threedhst.options['FILTWAVE'],
    #     inter_type="linear", np=10, interp=-1, smooth_lengt=0, smooth_fwhm=0.0,
    #     spectr="NO", adj_sens=threedhst.options['AXE_ADJ_SENS'], weights="NO",
    #     sampling="drizzle")
        
    #### Local background
    #LOCAL_BACKGROUND="YES"
    LOCAL_BACKGROUND="NO"
    
    if threedhst.options['USE_TAXE']:
        taxe21.taxecore(inlist=prep_name(asn_grism_file), configs=CONFIG,
            back=LOCAL_BACKGROUND,extrfwhm=4.0, drzfwhm=3.0, backfwhm=4.0,
            slitless_geom=geom_yn, orient=geom_yn, exclude="NO", 
            lambda_mark=threedhst.options['FILTWAVE'], 
            cont_model="fluxcube", model_scale=4.0, 
            lambda_psf=threedhst.options['FILTWAVE'],
            inter_type="linear", np=10, interp=0, smooth_lengt=0,
            smooth_fwhm=0.0,
            spectr="NO", adj_sens=threedhst.options['AXE_ADJ_SENS'],
            weights="NO",
            sampling="drizzle")
    else:
        iraf.axecore(inlist=prep_name(asn_grism_file), configs=CONFIG,
            back=LOCAL_BACKGROUND,extrfwhm=4.0, drzfwhm=3.0, backfwhm=4.0,
            slitless_geom=geom_yn, orient=geom_yn, exclude="NO", 
            lambda_mark=threedhst.options['FILTWAVE'], 
            cont_model="fluxcube", model_scale=4.0, 
            lambda_psf=threedhst.options['FILTWAVE'],
            inter_type="linear", np=10, interp=0, smooth_lengt=0,
            smooth_fwhm=0.0,
            spectr="NO", adj_sens=threedhst.options['AXE_ADJ_SENS'],
            weights="NO",
            sampling="drizzle")
    
    threedhst.currentRun['step'] = 'AXECORE'
        
    #############################################
    #### aXe.tdrzprep - prepare for aXe drizzling
    #############################################
    
    flprMulti()
    #print "\n3DHST.proces_grism: iraf.DRZPREP\n"
    threedhst.showMessage('Running iraf.drzprep')
    
    if threedhst.options['USE_TAXE']:
        status = taxe21.tdrzprep(inlist=prep_name(asn_grism_file),
            configs=CONFIG,
            opt_extr="YES", back=LOCAL_BACKGROUND, Stdout=1)
    else:
        status = iraf.drzprep(inlist=prep_name(asn_grism_file), 
            configs=CONFIG,
            opt_extr="YES", back=LOCAL_BACKGROUND, Stdout=1)
    
    threedhst.currentRun['step'] = 'DRZPREP'
        
    #############################################
    #### aXe.taxedrizzle - drizzle combine the dithers.  This step also 
    #### produces the 1D extraction
    #############################################
    threedhst.options['DRIZZLE_PATH'] = os.environ['AXE_DRIZZLE_PATH']
    #print "\n3DHST.proces_grism: iraf.AXEDRIZZLE (#1)\n"
    threedhst.showMessage('Running iraf.axedrizzle')
    
    flprMulti()
    if threedhst.options['USE_TAXE']:
        status = taxe21.taxedrizzle(inlist=prep_name(asn_grism_file),
                        configs=CONFIG,
                        infwhm=4.0,outfwhm=3.0, back=LOCAL_BACKGROUND,
                        makespc="YES",
                        opt_extr=threedhst.options['AXE_OPT_EXTR'],
                        adj_sens=threedhst.options['AXE_ADJ_SENS'], 
                        driz_separate='NO', Stdout=1)
    else:
        status = iraf.axedrizzle(inlist=prep_name(asn_grism_file),
                        configs=CONFIG,
                        infwhm=4.0,outfwhm=3.0, back=LOCAL_BACKGROUND,
                        makespc="YES",
                        opt_extr=threedhst.options['AXE_OPT_EXTR'],
                        adj_sens=threedhst.options['AXE_ADJ_SENS'], 
                        driz_separate='NO', Stdout=1)
    
    #### Run drizzle again for internal CR rejection.  Perhaps all CRs are
    #### already flagged from the earlier MultiDrizzle run
    if threedhst.options['RUN_DRZREJ']:
        #### Check if DRZREJ directory exists
        try:
            os.mkdir('./DRZREJ_'+
                        threedhst.options['GRISM_NAME']+'/')
        except:
            pass

        os.environ['AXE_DRIZZLE_PATH'] = ('./DRZREJ_'+
                    threedhst.options['GRISM_NAME']+'/')
        threedhst.options['DRIZZLE_PATH'] = os.environ['AXE_DRIZZLE_PATH']
        
        #### Run axedrizzle again
        #print "\n3DHST.proces_grism: iraf.AXEDRIZZLE (#2)\n"
        threedhst.showMessage('Running iraf.axedrizzle, second pass.')
        
        flprMulti()
        if threedhst.options['USE_TAXE']:
            taxe21.taxedrizzle(inlist=prep_name(asn_grism_file), configs=CONFIG,
                            infwhm=4.0,outfwhm=3.0, back=LOCAL_BACKGROUND,
                            makespc="YES",
                            opt_extr=threedhst.options['AXE_OPT_EXTR'], 
                            adj_sens=threedhst.options['AXE_ADJ_SENS'], 
                            driz_separate='YES',
                            combine_type='median', combine_maskpt=0.7,
                            combine_nsigmas="4.0 3.0", combine_nlow=0, 
                            combine_nhigh=0, combine_lthresh="INDEF",
                            combine_hthresh="INDEF", combine_grow=1.0, 
                            blot_interp='poly5', blot_sinscl=1.0, 
                            driz_cr_snr="5.0 4.0", driz_cr_grow=1,
                            driz_cr_scale="1.2 0.7")
        else:
            iraf.axedrizzle(inlist=prep_name(asn_grism_file), configs=CONFIG,
                        infwhm=4.0,outfwhm=3.0, back=LOCAL_BACKGROUND,
                        makespc="YES",
                        opt_extr=threedhst.options['AXE_OPT_EXTR'], 
                        adj_sens=threedhst.options['AXE_ADJ_SENS'], 
                        driz_separate='YES',
                        combine_type='median', combine_maskpt=0.7,
                        combine_nsigmas="4.0 3.0", combine_nlow=0, 
                        combine_nhigh=0, combine_lthresh="INDEF",
                        combine_hthresh="INDEF", combine_grow=1.0, 
                        blot_interp='poly5', blot_sinscl=1.0, 
                        driz_cr_snr="5.0 4.0", driz_cr_grow=1,
                        driz_cr_scale="1.2 0.7")
    
    threedhst.currentRun['step'] = 'AXEDRIZZLE'
    
    #return
    
    #############################################
    #### Make a multidrizzled contamination image, which can be 
    #### compared to the grism multidrizzle mosaic
    #############################################
    #print "\n3DHST.proces_grism: iraf.Multidrizzle GRISM mosaic\n"
    threedhst.showMessage('Run multidrizzle to create a mosaic of the grism exposures')
    
    os.chdir('./DATA')
    for expi in asn_grism.exposures:
        ### copy ../OUTPUT_G141/CONT.fits into existing grism FLT files
        flt = pyfits.open(expi+'_flt.fits','update')
        cont = pyfits.open('../OUTPUT_'+threedhst.options['GRISM_NAME']+
                           '/'+expi+'_flt_2.CONT.fits')
        flt[1].data  = cont[1].data
        
        if threedhst.options['GRISM_NAME'] == 'G800L':
            cont = pyfits.open('../OUTPUT_'+threedhst.options['GRISM_NAME']+
                               '/'+expi+'_flt_5.CONT.fits')
            ##### Contamination is in cts/s, normal ACS FLT is just cts
            flt[1].data *= flt[0].header['EXPTIME']
            flt[4].data = cont[1].data*flt[0].header['EXPTIME']
            # flt[4].data = cont[1].data*flt[0].header['EXPTIME']
            # flt[1].header.update('BUNIT','ELECTRONS/S')
            # flt[4].header.update('BUNIT','ELECTRONS/S')

        flt.flush()

    asn_grism.product = ROOT_GRISM+'CONT' #ROOT_GRISM+'CONT'
    asn_grism.write(ROOT_GRISM+'CONT_asn.fits')
    flprMulti()
    
    if threedhst.options['PREFAB_DIRECT_IMAGE'] is None:
        status = drizzle_function(input=ROOT_GRISM+'CONT_asn.fits', 
                      shiftfile=ROOT_GRISM +  '_shifts.txt', 
                      output = '', 
                      final_scale = threedhst.options['DRZSCALE'], 
                      final_pixfrac = threedhst.options['USED_PIXFRAC'], 
                      skysub = no, static=no, driz_separate=no, median=no, 
                      blot=no, driz_cr=no, Stdout=1)
    else:
        #### Match the grism mosaic WCS to the direct image
        # header = pyfits.getheader(threedhst.options['PREFAB_DIRECT_IMAGE'],1)
        # NX = header.get('NAXIS1')
        # NY = header.get('NAXIS2')
        # RA_CENTER = header.get('CRVAL1')
        # DEC_CENTER = header.get('CRVAL2')
        # #SCALE = np.abs(header.get('CD1_1'))*3600.
        # SCALE = np.sqrt((header.get('CD1_1')*3600.)**2+
        #                 (header.get('CD1_2')*3600.)**2)
        # 
        # ANGLE = np.arctan2(header.get('CD2_1'),
        #                    header.get('CD2_2'))/2/np.pi*360.
        # 
        # print 'NX NY RA DEC SCALE'
        # print '%d %d %13.6f %13.6f %6.3f' %(NX, NY, RA_CENTER, DEC_CENTER, SCALE)
        # 
        # iraf.multidrizzle(input=ROOT_GRISM+'CONT_asn.fits', 
        #               shiftfile=ROOT_GRISM + '_shifts.txt', 
        #               output = '', final_scale = SCALE, 
        #               final_pixfrac = threedhst.options['USED_PIXFRAC'], 
        #               skysub = no, final_outnx=NX, final_outny=NY, 
        #               ra=RA_CENTER, dec=DEC_CENTER, final_rot=0.,
        #               static=no, 
        #               driz_separate=no, median=no, blot=no, driz_cr=no)

        drizzle_function(input=ROOT_GRISM+'CONT_asn.fits', 
                      shiftfile=ROOT_GRISM + '_shifts.txt', 
                      output = '', 
                      final_pixfrac = threedhst.options['USED_PIXFRAC'], 
                      refimage=threedhst.options['PREFAB_DIRECT_IMAGE']+'[1]',
                      skysub = no, static=no, 
                      driz_separate=no, median=no, blot=no, driz_cr=no)
    
    # status = iraf.multidrizzle(input = ROOT_GRISM+'CONT_asn.fits', \
    #    output = '', shiftfile = ROOT_GRISM + '_shifts.txt', \
    #    final_scale = threedhst.options['DRZSCALE'], 
    #    final_pixfrac = threedhst.options['USED_PIXFRAC'], skysub = no,
    #    static=no, driz_separate=no,median=no, blot=no, driz_cr=no, Stdout=1)
    # status = iraf.multidrizzle(input = ROOT_GRISM+'CONT_asn.fits', \
    #    output = '', shiftfile = ROOT_GRISM + '_shifts.txt', \
    #    final_scale = iraf.INDEF, 
    #    final_pixfrac = 1.0, skysub = no,
    #    static=no, driz_separate=no,median=no, blot=no, driz_cr=no, Stdout=1)

    cleanMultidrizzleOutput()
    
    threedhst.currentRun['step'] = 'MULTIDRIZZLE_CONTAMINATION'
    
    #############################################
    #### Clean up after multidrizzle and aXe routines
    #############################################
    
    if threedhst.options['CLEAN_UP']:
        rmfiles = []
        #rmfiles.extend(glob.glob('*FLX.fits'))
        rmfiles.extend(glob.glob('flux?.fits'))
        rmfiles.extend(glob.glob('*_flt.fits'))
        rmfiles.extend(glob.glob('*flt_1.cat'))
        rmfiles.extend(glob.glob(ROOT_GRISM+'*[SW][CH]?.fits'))
        #rmfiles.extend(glob.glob('*coeffs1.dat'))
        #rmfiles.extend(glob.glob('threedhst_auto.*'))
        #rmfiles.extend(glob.glob('zeropoints.lis'))
        rmfiles.extend(glob.glob('default.conv'))
        #for expi in asn_grism.exposures:
        #    rmfiles.extend(glob.glob('../OUTPUT_G141/'+expi+'*'))
        
        if len(rmfiles) > 0:
            for rmfile in rmfiles:
                os.remove(rmfile)
        
        threedhst.currentRun['step'] = 'CLEANUP_DATA'
    
    #############################################
    #### Make output webpages with spectra thumbnails    
    #############################################
    if threedhst.options['MAKE_WEBPAGE']:
        threedhst.plotting.make_data_products(ROOT_DIRECT, ROOT_GRISM)
            
    #############################################
    #### Write parameter log
    #############################################
    os.chdir('../')
    logfile = ROOT_GRISM+'.threedhst.param'
    print('\nthreedhst: Parameter log in <%s>.\n' %logfile)
    
    threedhst.showOptions(to_file = logfile)
    if threedhst.options['MAKE_WEBPAGE']:
        threedhst.showOptions(to_file = "./HTML/"+logfile)

    #############################################
    #### Add HAS_SPEC and FCONTAM columns to the catalog
    #############################################
    if threedhst.options['GRISM_NAME'] == 'G141':
        CONT_LAM = 1.4e4
    if threedhst.options['GRISM_NAME'] == 'G102':
        CONT_LAM = 0.95e4
    if threedhst.options['GRISM_NAME'] == 'G800L':
        CONT_LAM = 0.7e4
        
    update_catalogs(root=ROOT_GRISM, HTML_DIR='./HTML/', DATA_DIR='./DATA/', CONT_LAM = CONT_LAM)
    
    #############################################
    #### Done!
    #############################################
    
    print('\nthreedhst: cleaned up and Done!\n')
    threedhst.currentRun['step'] = 'DONE'

def safe_chdir(dir):
    """
safe_chdir(dir)
    
    Wrap os.chdir within a try/except clause to avoid dying on a
    directory not found.
    """
    import os
    try:
        os.chdir(dir)
    except:
        threedhst.showMessage('Directory %s not found!' %(dir), warn=True)
        pass

def clean_flt_files(asn_filename):
    """
    Delete FLT files defined in an ASN table
    """
    import threedhst    
    ASN  = threedhst.utils.ASNFile(file=asn_filename)
    for exp in ASN.exposures:
        os.remove(exp+'_flt.fits')
        
def fresh_flt_files(asn_filename, from_path="../RAW", preserve_dq = False):
    """
    """ 
    #import threedhst.dq
    import threedhst
    import threedhst.prep_flt_files
    
    ASN  = threedhst.utils.ASNFile(file=asn_filename)
    
    #### Copy fresh FLT files from RAW
    explist = []
    explist.extend(ASN.exposures)
    for exp in explist:
        fits_file = threedhst.utils.find_fits_gz(from_path+'/'+exp+'_flt.fits')
        fi = pyfits.open(fits_file)
        
        try:
            os.remove('./'+exp+'_flt.fits')
        except:
            pass
        
        print(exp)
        
        dq = fi[3]
        if preserve_dq:
            if os.path.exists('./'+exp+'_flt.fits'):
                old = pyfits.open('./'+exp+'_flt.fits')
                dq = old[3]
            
        fi[3] = dq
        fi.writeto('./'+exp+'_flt.fits', clobber=True)
        threedhst.prep_flt_files.apply_best_flat(exp+'_flt.fits', verbose=True)
        
        threedhst.prep_flt_files.apply_persistence_mask(exp+'_flt.fits', limit_sigma=threedhst.options['FLT_PERSISTENCE_THRESH'], filter_size=threedhst.options['FLT_PERSISTENCE_FILTER'], persistence_path=threedhst.options['FLT_PERSISTENCE_PATH'], verbose=True)
        
        #### Make sure most recent IDCTAB is used
        flt = pyfits.open(exp+'_flt.fits','update')
        head = flt[0].header
        if ('INSTRUME' in list(head.keys())) & ('DETECTOR' in list(head.keys())):
            if (head['INSTRUME'].strip() == 'WFC3') & (head['DETECTOR'].strip() == 'IR'):
                #print 'IDCTAB %s -> iref$uab1537ci_idc.fits' %(head['IDCTAB'])
                #head.update('IDCTAB','iref$uab1537ci_idc.fits')
                #flt[0].header.update('IDCTAB','iref$v5r1512fi_idc.fits')
                ### Force use latest IDC tab
                flt[0].header.update('IDCTAB','iref$w3m18525i_idc.fits')
        
                #### Add crosstalk to the pixel uncertainties
                var = flt['ERR'].data**2 
                xtalk = (flt['SCI'].data[::-1,:]/2000.)**2
                time = flt['TIME'].data
                if time is None:
                    time = flt['TIME'].header['PIXVALUE']
                    
                xtalk_mask = (flt['SCI'].data*time)[::-1,:] > 2.e4
                var[xtalk_mask] += xtalk[xtalk_mask]
                flt['ERR'].data = np.sqrt(var)

        #flt[3] = dq
        flt.flush()
        
        #### Apply DQ mask (.mask.reg), if it exists
        threedhst.regions.apply_dq_mask(os.path.basename(fits_file.split('.gz')[0]),
                                        addval=4)

def mask_IR_blobs(flt='ibsyb6q7q_flt.fits'):
    blobs = [[293.2,640],
             [260.5,503.83333],
             [232.16667,418.83333],
             [470.5,375.5],
             [431.83333,287.66667],
             [408,85.8]]
    
    im = pyfits.open(flt, mode='update')
    yc, xc = np.indices((1014,1014))
    
    for blob in blobs:
        r = np.sqrt((xc-blob[0])**2+(yc-blob[1])**2)
        mask = (r < 30) & ((im['DQ'].data & 512) > 0)
        im['DQ'].data[mask] |= 4096
    #
    print('Masked IR blobs: %s.' %(flt))
    im.flush()
    
def swarpOtherBands():
    """
swarpOtherBands()
    
    Align additional bands to the direct image for later use in a fluxcube.
    
    The band information will be read from threedhst.options['OTHER_BANDS']. 
    To define an additional band, add it before running `process_grism`:
    
    options['OTHER_BANDS'].append(['ACS/h_nz*drz_img.fits','F850LP',903.,24.84])
    
    [ search string , filter name,  pivot wavelength (nm),  zeropoint]
    
    The first parameter is a string (including optional wildcards) used 
    to find the "raw" files in the new band that will be swarped to match
    the direct image.
    
    """
    import glob
    
    ROOT_DIRECT = threedhst.options['ROOT_DIRECT']
    ROOT_GRISM = threedhst.options['ROOT_GRISM']
    DIRECT_MOSAIC = threedhst.options['DIRECT_MOSAIC']
    
    other_bands = threedhst.options['OTHER_BANDS']
        
    #### Nothing found
    if not other_bands:
        return None
    
    #### Other bands are defined, use them
    
    #### Initialize SWarp
    sw = threedhst.sex.SWarp()
    sw._aXeDefaults()
    ## get reference image parameters
    sw.swarpMatchImage(DIRECT_MOSAIC,extension=1)  
    ## swarp the reference image to istelf
    sw.swarpImage(DIRECT_MOSAIC+'[1]',mode='wait')    
    ## Refine center position from SWarp's own output
    sw.swarpRecenter()                  
    
    #############################################
    #### Loop through OTHER_BANDS and swarp them 
    #### to match the direct mosaic
    #############################################    
    band_count=0
    for ib, band in enumerate(other_bands):
        #### should be like:
        ####    ['filename', 'Filter', pivot_wavelength, zeropoint]
        ####    ['ACS/h_nz*drz_img.fits','F850LP',903.,24.84]
        if len(band) != 4:
            continue
        else:
            ++band_count
                
        #### SWarp input image list to same pixels as reference (direct) image    
        #otherImages = glob.glob(band[0])
        otherImages = threedhst.shifts.find_align_images_that_overlap(
                  DIRECT_MOSAIC, band[0])
                  
        sw.swarpImage(otherImages,mode='direct')   
        
        #### Put the result from "coadd.fits" into the first extension 
        #### of a copy of the direct drz image.
        direct = pyfits.open(DIRECT_MOSAIC)
        new = pyfits.open('coadd.fits')
        direct[1].data = new[0].data
        
        direct.writeto(ROOT_GRISM+'_'+band[1]+'_drz.fits', clobber=True)
        del(direct)
        try:
            os.remove('coadd.fits')
            os.remove('coadd.weight.fits')
        except:
            pass
        
def prep_name(input_asn):
    """
    make_prep_name(input_asn)
    
    Example: 
        >>> prepFile = make_prep_name("ib3714060_asn.fits")
        >>> print prepFile
        ib3714060_prep.lis
    """
    return input_asn.split('_asn.fits')[0] + '_prep.lis'


def make_aXe_lis(asn_grism_file, asn_direct_file, mode='DIRECT'):
    """
    status = make_aXe_lis(asn_grism_file, asn_direct_file)
    
    Make "inlist" file for aXe routines, with format
    
        grismA_flt.fits directA_flt_1.cat directA_flt.fits 0.0
        grismB_flt.fits directB_flt_1.cat directB_flt.fits 0.0
        ...
    
    If mode != 'DIRECT', then make a lis file formatted using the grism mosaic
    as reference:
        
        grismA_flt.fits grismA_flt_1.cat 0.0
        grismB_flt.fits grismB_flt_1.cat 0.0
        ...
        
    Returns "True" if executes correctly, "False" on an error
    
    """
    from threedhst.utils import ASNFile
    
    asn_grism = ASNFile(asn_grism_file)
    
    if mode.startswith('DIRECT'):
        asn_direct = ASNFile(asn_direct_file)
        
        #### Check that grism and direct ASN tables have same # of entries
        if len(asn_grism.exposures) != len(asn_direct.exposures):
            threedhst.showMessage("""Number of grism exposures (%d) in %s is different from
                       the number of direct images (%d) in %s.
            """ %(len(grism_files), asn_grism, len(direct_files), asn_direct))
            
            return False
        
    #### Make the lis file
    NFILE = len(asn_grism.exposures)
    outfile = prep_name(asn_grism_file)
    fp = open(outfile,'w')
    
    if threedhst.options['ACS_G800L']:
        if mode.startswith('DIRECT'):
            for i in range(NFILE):
                fp.write("%s_flt.fits %s_flt_1.cat,%s_flt_2.cat %s_flt.fits 0.0\n" 
                      %(asn_grism.exposures[i], asn_direct.exposures[i],
                        asn_direct.exposures[i], 
                        asn_direct.exposures[i]))
        else:
            for i in range(NFILE):
                fp.write("%s_flt.fits %s_flt_1.cat,%s_flt_2.cat 0.0\n" 
                      %(asn_grism.exposures[i], asn_grism.exposures[i],
                        asn_grism.exposures[i]))
    else:
        if mode.startswith('DIRECT'):
            for i in range(NFILE):
                fp.write("%s_flt.fits %s_flt_1.cat %s_flt.fits 0.0\n" 
                      %(asn_grism.exposures[i], asn_direct.exposures[i], 
                        asn_direct.exposures[i]))
        else:
            for i in range(NFILE):
                fp.write("%s_flt.fits %s_flt_1.cat 0.0\n" 
                      %(asn_grism.exposures[i], asn_grism.exposures[i]))
        
    fp.close()
    # print "3D-HST / make_aXe_lis: Created %s\n" %outfile
    threedhst.showMessage('Created .list file, %s.' %outfile)
    
    return True

def cleanMultidrizzleOutput():
    """
    cleanMultidrizzleOutput()
    
    Remove *single_[sci/wht].fits, *sci1_blt.fits, *flt*mask1.fits, *coeffs1.dat
    """
    import os,glob
    rmfiles = glob.glob('*single_???.fits')
    rmfiles.extend(glob.glob('*sci[12]_blt.fits'))
    rmfiles.extend(glob.glob('*flt*mask[12].fits'))
    #files.extend(glob.glob('*coeffs1.dat'))
    #rmfiles.extend(glob.glob('ib*.run'))
    rmfiles.extend(glob.glob('*_med.fits'))
    if len(rmfiles) > 0:
        for rmfile in rmfiles:
            os.remove(rmfile)



def flprMulti(n=3):
    """
    flprMulti(n=3)
    
    Run iraf.flpr() `n` times.
    """
    from pyraf import iraf

    if n < 1:
        n=1
    for i in range(n):
        iraf.flpr()

def die():
    """
die()
    
    Print lines to string that will make it easier to restart 
    at a particular place in the reduction script.
    """
    print("""
# step: %s
asn_grism_file = threedhst.currentRun['asn_grism_file']
asn_direct_file = threedhst.currentRun['asn_direct_file']
ROOT_GRISM = threedhst.options['ROOT_GRISM']
ROOT_DIRECT = threedhst.options['ROOT_DIRECT']
sexCat = threedhst.sex.mySexCat(ROOT_GRISM+'_drz.cat')
if 'conf' in threedhst.currentRun.keys():
    conf = threedhst.currentRun['conf']
SPC = threedhst.plotting.SPCFile(ROOT_GRISM+'_2_opt.SPC.fits')
""" %threedhst.currentRun['step'])
    raise IOError

def multidrizzle_defaults(asn_file):
    """
multridrizzle_defaults():

    Set multidrizzle default parameters according to the reference MDZ images
    """
    from pyraf import iraf
    from iraf import stsdas,dither,slitless,axe 
    
    no = iraf.no
    yes = iraf.yes
    INDEF = iraf.INDEF
    
    iraf.unlearn('multidrizzle')
    asn = threedhst.utils.ASNFile(asn_file)
    NEXP = len(asn.exposures)
    header = pyfits.getheader(threedhst.utils.find_fits_gz(asn.exposures[0]+'_flt.fits'))
    INSTRUME = header.get('INSTRUME')
    FILTER = header.get('FILTER')
    
    threedhst.showMessage('Instrument: %s\n    Filter: %s' %(INSTRUME, FILTER))
    
    yes = iraf.yes
    no = iraf.no
    INDEF = iraf.INDEF
    
    iraf.multidrizzle.mdriztab = no
    iraf.multidrizzle.updatewcs = yes
    iraf.multidrizzle.coeffs = 'header'
    iraf.multidrizzle.context = no
    iraf.multidrizzle.clean = no
    iraf.multidrizzle.build = yes
    if NEXP > 1:
        iraf.multidrizzle.static = yes
    else:
        iraf.multidrizzle.static = no
        
    iraf.multidrizzle.static_sig = 4.0
    iraf.multidrizzle.skysub = yes
    iraf.multidrizzle.skywidth = 0.1
    iraf.multidrizzle.skystat = 'mode'
    iraf.multidrizzle.skylower = -100
    iraf.multidrizzle.skyupper = INDEF
    iraf.multidrizzle.skyclip = 5
    iraf.multidrizzle.skylsigma = 4.0
    iraf.multidrizzle.skyusigma = 4.0
    
    iraf.multidrizzle.driz_separate = no
    if (NEXP > 1) & (FILTER != 'G141'):
        iraf.multidrizzle.driz_separate = yes
        
    iraf.multidrizzle.driz_sep_kernel = 'turbo'
    iraf.multidrizzle.driz_sep_wt_scl = 'exptime'
    iraf.multidrizzle.driz_sep_scale = INDEF
    iraf.multidrizzle.driz_sep_pixfrac = 1.0
    iraf.multidrizzle.driz_sep_rot = INDEF
    iraf.multidrizzle.driz_sep_fillval = INDEF
    
    if INSTRUME == 'WFC3':
        iraf.multidrizzle.driz_sep_bits = 576
    if INSTRUME == 'ACS':
        iraf.multidrizzle.driz_sep_bits = 96
    
    iraf.multidrizzle.driz_sep_bits = 0
    
    iraf.multidrizzle.median = no
    iraf.multidrizzle.median_newmasks = no
    if (NEXP > 1) & (FILTER != 'G141'):
        iraf.multidrizzle.median = yes
        iraf.multidrizzle.median_newmasks = yes
    
    iraf.multidrizzle.combine_maskpt = 0.7
    iraf.multidrizzle.combine_type = 'minmed'
    iraf.multidrizzle.combine_nsigma = '4 3'
    iraf.multidrizzle.combine_nlow = 0
    iraf.multidrizzle.combine_nhigh = 1
    iraf.multidrizzle.combine_lthresh = INDEF
    iraf.multidrizzle.combine_hthresh = INDEF
    iraf.multidrizzle.combine_grow = 1
    iraf.multidrizzle.blot = no
    if (NEXP > 1) & (FILTER != 'G141'):
        iraf.multidrizzle.blot = yes
        
    iraf.multidrizzle.blot_interp = 'poly5'
    iraf.multidrizzle.blot_sinscl = 1.0
    iraf.multidrizzle.driz_cr = yes
    iraf.multidrizzle.driz_cr_corr = no
    if INSTRUME == 'WFC3':
        iraf.multidrizzle.driz_cr_snr = '5.0 4.0'
    else:
        iraf.multidrizzle.driz_cr_snr = '3.5 3.0'
        
    iraf.multidrizzle.driz_cr_grow = 1
    iraf.multidrizzle.driz_cr_ctegrow = 0
    iraf.multidrizzle.driz_cr_scale = '1.2 0.7'
    iraf.multidrizzle.driz_combine = yes
    iraf.multidrizzle.final_wht_type = 'EXP'
    iraf.multidrizzle.final_kernel = 'square'
    iraf.multidrizzle.final_wt_scl = 'exptime'
    iraf.multidrizzle.final_scale = INDEF
    iraf.multidrizzle.final_pixfrac = 1.0
    iraf.multidrizzle.final_rot = 0.0
    iraf.multidrizzle.final_fillval = INDEF
    
    #### Set OK bits for output
    if INSTRUME == 'WFC3':
        iraf.multidrizzle.final_bits = 576
    if INSTRUME == 'ACS':
        iraf.multidrizzle.final_bits = 96
    
    #### use only unflagged pixels
    # iraf.multidrizzle.final_bits = 0
    
    iraf.multidrizzle.final_units = 'cps'
    iraf.multidrizzle.gnkeyword = "ATODGNA,ATODGNB,ATODGNC,ATODGND"
    iraf.multidrizzle.rnkeyword = "READNSEA,READNSEB,READNSEC,READNSED"
    iraf.multidrizzle.expkeyword = "EXPTIME"
    # if INSTRUME == 'WFC3':
    #     iraf.multidrizzle.crbit = 0
    # else:
    #     iraf.multidrizzle.crbit = 4096



def fix_g141_sky_image(conf_path='./', verbose=True):
    """
fix_g141_sky_image(conf_path='./', verbose=True)
    
    The G141 average sky image, WFC3.IR.G141.sky.V1.0.fits,  has lots of 
    holes with pixel value = 0,  apparently more than pixels flagged in the BPM.
    
    Fill these in by taking the median of (nonzero) neighboring pixels.
    
    The result is an image, g141_fixed_sky.fits, put in the directory defined
    by `conf_path` (WFC3.IR.G141.sky.V1.0.fits is read from `conf_path`).
    
    """
    
    if verbose:
        threedhst.showMessage("""Filling holes in G141 sky image.
Updated sky image in %s/g141_fixed_sky.fits""" %conf_path)
    
    im = pyfits.open(conf_path+'/'+'WFC3.IR.G141.sky.V1.0.fits')
    sky = im[0].data
    
    x,y = np.where(sky == 0)
    N = len(x)
    pad = 1
    for i in range(N):
        xi = x[i]
        yi = y[i]
        sub = sky[xi-pad:xi+pad+1,yi-pad:yi+pad+1]
        if (np.sum(sub) != 0.0):
            sky[xi,yi] = np.median(sub[sub != 0.0])
    
    im.writeto(conf_path+'/'+'g141_fixed_sky.fits', clobber=True)
  
class Conf(object):
    """
    Conf(infile='WFC3.IR.G141.V1.0.conf')
    
    Read an aXe configuration file for easy manipulation of the parameters.
    
    """            
    def _getPath(self):
        """
        _getPath()
        
        Figure out if we're in the root directory or in DATA
        """
        import os
        if os.getcwd().split('/')[-1] == 'DATA':
            self.path = os.getcwd()+'/../CONF/'
        else:
            self.path = os.getcwd()+'/CONF/'
        
    def _processLines(self):
        """
        _processLines()
        
        Read the lines, extracting PARAMETER VALUE pairs
        """
        self.nlines = len(self.lines)
        self.params = {}
        self._pline = {}
        for i,line in enumerate(self.lines):
            if (line[0] is not '#') & (line.strip() is not  ''):
                spl = line.split(';')[0].split()
                self.params[spl[0]] = ' '.join(spl[1:])
                self._pline[spl[0]] = i
        self.nkeys = list(self.params.keys()).__len__()
        
    def _readLines(self):
        """
        _readlines()
        """
        #self._getPath()
        fp = open(self.path+self.infile,'r')
        self.lines = fp.readlines()
        fp.close()
        
    def __init__(self, infile='WFC3.IR.G141.V1.0.conf', path=None):
        self.infile = infile
        if path is None:
            self._getPath()
        else:
            if not path.startswith('/'):
                path = os.getcwd()+'/'+path
            self.path = path
        self._readLines()
        self._processLines()
        
    def _assignPars(self):
        """
        _assignPars()
        
        Apply parameter values to self.lines
        """
        for key in list(self.params.keys()):
            param = self.params[key]
            if type(param) is not type('xxx'):
                param = np.str(param)
            
            #### New parameter?
            if key in self._pline:
                self.lines[self._pline[key]] = key + ' ' + param +'\n'
            else:
                self.lines.append(key + ' ' + param +'\n')
        
        self._processLines()
    
    def writeto(self, output='tmp.conf'):
        """
        writeto(output='tmp.conf')
        """ 
        #self._getPath()
        self._assignPars()
        fp = open(self.path+output,'w')
        fp.writelines(self.lines)
        fp.close()

def set_ACS_G800L():
    """ Set default parameters for ACS G800L rather than WFC3 """
    import threedhst
    
    #### ACS rather than WFC3 grism
    threedhst.options['ACS_G800L'] = True
    
    threedhst.options['CONFIG_FILE'] = 'ACS.WFC.CHIP2.Cycle13.5.conf,ACS.WFC.CHIP1.Cycle13.5.conf'
    threedhst.options['SKY_BACKGROUND'] = 'ACS.WFC.CHIP2.msky.1.fits,ACS.WFC.CHIP1.msky.1.fits'
    
    threedhst.options['DRZRESOLA'] = '40.0'
    threedhst.options['DRZSCALE'] = '0.05'
    
    threedhst.options['GRISM_NAME'] = 'G800L'
    threedhst.options['MAG_ZEROPOINT'] = 25.94333
    threedhst.options['FILTWAVE'] = 805.6948
    
    
def update_FLX_WCS(path_to_FLT='../RAW/'):
    """
    There seems to be something wrong with the WCS keywords of the FLX images
    produced by fcubeprep.  Force CRVALX and CDX_X keywords to match those 
    in the FLT files.
    """
    import glob
    
    edges = threedhst.options['AXE_EDGES'].split(',')
    padding = float(edges[1])+float(edges[0])
    
    flx_files = glob.glob('ib*FLX.fits')
    for flx_file in flx_files:
        flt_file = threedhst.utils.find_fits_gz('%s/%s.fits' %(path_to_FLT, flx_file.split('_2')[0]))
        flx = pyfits.open(flx_file,'update')
        flt = pyfits.open(flt_file)
        head_flt = flt[1].header
        keys = ['CRVAL1','CRVAL2','CD1_1','CD1_2','CD2_1','CD2_2']
        for ext in flx[1:]:
            head_flx = ext.header
            for key in keys:
                head_flx.update(key,head_flt.get(key))
            
            head_flx.update('CRPIX1',head_flx.get('NAXIS1')/2.-5.*padding/200.)
        #
        print('Updating WCS header keywords: %s\n' %(flx_file))
        flx.flush()
            
def update_catalogs(root='COSMOS-3-G141', HTML_DIR='./HTML/', DATA_DIR='./DATA/', CONT_LAM = 1.4e4, GRISM = 'G141'):
    """
    Add columns to the SExtractor catalogs:
        
        HAS_SPEC = 1/0 if a given object actually has a grism spectrum.
        FCONTAM = fractional contamination at lambda = CONT_LAM 
        FCOVER = fraction of pixels between 1.15 and 1.6 microns that have grism coverage
    (Run in root directory of a given field)
        ACONTAM = average contamination within coverage region
        MCONTAM = maximum contamination within coverage region
    """  
    import os
    import threedhst.catIO as catIO
    
    if GRISM == 'G141':
        COVER = [1.15e4, 1.6e4]
              
    PWD = os.getcwd()
    if not os.path.exists('./DATA/'):
        os.chdir('../')
        if not os.path.exists('./DATA'):
            threedhst.message("Can't find DATA, HTML directories", warn=True)
            return False
    
    #### Read SExtractor catalog
    sexCat = threedhst.sex.mySexCat(DATA_DIR+'/'+root+'_drz.cat')
    
    has_spec = sexCat.id*0
    fcontam = sexCat.id*0.
    acontam = fcontam*0.
    mcontam = fcontam*0.
    fcover = fcontam*0.
    
    #### Loop through objects looking for spectra and contam. fraction
    for idx, id in enumerate(sexCat.id):
        ascii_spec = HTML_DIR+'/ascii/'+root+'_%05d.dat' %(id)
        if os.path.exists(ascii_spec):
            has_spec[idx] = 1
            try:
                spec = catIO.Readfile(ascii_spec)
            except:
                has_spec[idx] = 0
                continue
            
            cint = np.interp(CONT_LAM, spec.lam, spec.contam)
            fint = np.interp(CONT_LAM, spec.lam, spec.flux)
                        
            try:
                fcontam[idx] = cint/fint
            except:
                "Above will have error for 0/0"
                fcontam[idx] = -1
            
            sub = (spec.lam > COVER[0]) & (spec.lam < COVER[1])
            sub_covered = (spec.lam > COVER[0]) & (spec.lam < COVER[1]) & (spec.flux != 0) & np.isfinite(spec.flux)
            
            if len(spec.lam[sub]) > 0:
                fcover[idx] = len(spec.lam[sub_covered])*1./len(spec.lam[sub])
                if fcover[idx] > 0:
                    acontam[idx] = np.mean((spec.contam/spec.flux)[sub_covered])
                    mcontam[idx] = np.max((spec.contam/spec.flux)[sub_covered])
                else:
                    acontam[idx] = -1
                    mcontam[idx] = -1
                
    #### Now add the columns to the catalogs
    if 'HAS_SPEC' not in sexCat.column_names:
        status = sexCat.addColumn(data=has_spec, name='HAS_SPEC', format='%d', comment='Equal to 1 if there is a grism spectrum for a given object')
    
    if 'FCONTAM' not in sexCat.column_names:
        status = sexCat.addColumn(data=fcontam, name='FCONTAM', format='%8.2f', comment='Contamination fraction at %.2e A.' %(CONT_LAM))
    
    if 'ACONTAM' not in sexCat.column_names:
        status = sexCat.addColumn(data=acontam, name='ACONTAM', format='%8.2f', comment='Average contamination at %.2e-%.2e A.' %(COVER[0], COVER[1]))
    
    if 'MCONTAM' not in sexCat.column_names:
        status = sexCat.addColumn(data=mcontam, name='MCONTAM', format='%8.2f', comment='Max contamination at %.2e-%.2e A.' %(COVER[0], COVER[1]))
    
    if 'FCOVER' not in sexCat.column_names:
        status = sexCat.addColumn(data=fcover, name='FCOVER', format='%8.2f', comment='Fraction of pixels with grism coverage at %.2e-%.2e A.' %(COVER[0], COVER[1]))
                
    #### Done
    sexCat.write()
    sexCat.write(HTML_DIR+'/'+os.path.basename(sexCat.filename))
    
    os.chdir(PWD)
