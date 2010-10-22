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
    
CONF can be symlinked from e.g. /research/HST/GRISM/CONF
    
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
    
    if threedhst.options['SKY_BACKGROUND'] is not None:
        if not os.path.exists('CONF/'+threedhst.options['SKY_BACKGROUND']):
            raise IOError("options['SKY_BACKGROUND'] doesn't exist:" +   
                      "CONF/"+threedhst.options['SKY_BACKGROUND'])

def write_to_log(lines, init=False):
    """
write_to_log(lines, init=False)
    
    Write output of IRAF tasks to options['ROOT_DIRECT']+'.iraf.log'
    """
    if init:
        fp = open(threedhst.options['ROOT_DIRECT']+'.iraf.log','w')
    else:
        fp = open(threedhst.options['ROOT_DIRECT']+'.iraf.log','a')
        
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
    import numpy as np
    import aXe2html.sexcat.sextractcat
    import threedhst.dq
    from threedhst.TerminalController import TerminalController
    
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
    ROOT_DIRECT = asn_direct_file.split('_asn.fits')[0]
    
    #### Save parameters for the current run
    threedhst.currentRun['asn_grism_file'] = asn_grism_file
    threedhst.currentRun['asn_direct_file'] = asn_direct_file
    threedhst.options['ROOT_GRISM'] = ROOT_GRISM
    threedhst.options['ROOT_DIRECT'] = ROOT_DIRECT
    
    #### Copy ASN files from RAW, if they don't already exist
    if not os.path.exists(asn_grism_file):
        shutil.copy(threedhst.options['PATH_TO_RAW']+asn_grism_file,'./')
    if not os.path.exists(asn_direct_file):
        shutil.copy(threedhst.options['PATH_TO_RAW']+asn_direct_file,'./')
    
    #print "\n3DHST.process_grism: Fresh ASN and FLT files...\n"
    threedhst.showMessage("Fresh ASN and FLT files...")
    
    #### Read ASN files
    asn_grism  = threedhst.utils.ASNFile(file=asn_grism_file)
    asn_direct = threedhst.utils.ASNFile(file=asn_direct_file)
    ### Force PROD-DTH to be the root of the input ASN file
    asn_grism.product = ROOT_GRISM
    asn_direct.product = ROOT_DIRECT
    asn_grism.writeToFile(asn_grism_file)
    asn_direct.writeToFile(asn_direct_file)
    
    #### Copy fresh FLT files from RAW
    threedhst.process_grism.fresh_flt_files(asn_direct_file, 
                          from_path=threedhst.options['PATH_TO_RAW'])
                                            
    threedhst.process_grism.fresh_flt_files(asn_grism_file, 
                          from_path=threedhst.options['PATH_TO_RAW'])
            
    threedhst.currentRun['step'] = 'COPY_FROM_RAW'
    
    ############################################################################
    ####   Compute shifts
    ############################################################################
    if os.path.exists(ROOT_DIRECT + '_shifts.txt') is False:
        threedhst.options['MAKE_SHIFTFILES'] = True
        
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
        
    ############################################################################
    ####   Make detection image and generate a catalog
    ############################################################################
    
    #### First Multidrizzle run on DIRECT images to create a detection image
    #print "\n3DHST.process_grism: Creating MultiDrizzled detection image\n"
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
        
    status = iraf.multidrizzle(input=asn_direct_file, \
       shiftfile=ROOT_DIRECT + '_shifts.txt', skysub=MDRIZ_SKYSUB, \
       output = '', final_scale = threedhst.options['DRZSCALE'], 
       final_pixfrac = threedhst.options['USED_PIXFRAC'], Stdout=1)
    # status = iraf.multidrizzle(input=asn_direct_file, \
    #    shiftfile=ROOT_DIRECT + '_shifts.txt', skysub=MDRIZ_SKYSUB, \
    #    output = '', final_scale = iraf.INDEF, 
    #    final_pixfrac = 1.0, Stdout=1)
    
    ### Read the '.run' file output by MultiDrizzle
    mdrzRun = MultidrizzleRun(root=ROOT_DIRECT)
    
    cleanMultidrizzleOutput()
    threedhst.currentRun['step'] = 'DETECTION_IMAGE'
    
    #############################################
    #### Align to reference image
    #### e.g. threedhst.options['ALIGN_IMAGE'] = '../ACS/h_sz*drz_img.fits'
    #############################################
    if threedhst.options['ALIGN_IMAGE']:
        
        #print "\n3DHST.process_grism: Aligning WCS\n"
        threedhst.showMessage('Aligning WCS to %s'
                  %threedhst.options['ALIGN_IMAGE'])
        #### Compute the alignment
        xshift, yshift, rot, scale = threedhst.shifts.align_to_reference(
                             ROOT_DIRECT,threedhst.options['ALIGN_IMAGE'],
                             fitgeometry=threedhst.options['ALIGN_GEOMETRY'])
        
        #### shifts measured in DRZ frame.  Translate to FLT frame
        drz = pyfits.open(ROOT_DIRECT+'_drz.fits')
        alpha = (180.-drz[1].header['PA_APER'])/360.*2*np.pi

        xsh = (xshift*np.cos(alpha)-yshift*np.sin(alpha))*np.float(mdrzRun.scl)
        ysh = (xshift*np.sin(alpha)+yshift*np.cos(alpha))*np.float(mdrzRun.scl)
        
        print 'Final shift:', xsh, ysh, drz[1].header['PA_APER']
        fp = open(ROOT_DIRECT+'_align.info','w')
        fp.write('%s %8.3f %8.3f\n' 
                  %(threedhst.options['ALIGN_IMAGE'], xsh, ysh)) 
        fp.close()
        
        #### Read the shiftfile       
        shiftF = threedhst.shifts.ShiftFile(ROOT_DIRECT+'_shifts.txt')
        
        #### Apply the alignment shifts to the shiftfile
        shiftF.xshift = list(np.array(shiftF.xshift)-xsh)
        shiftF.yshift = list(np.array(shiftF.yshift)-ysh)
        shiftF.rotate = list(np.array(shiftF.rotate)+rot)
        shiftF.scale = list(np.array(shiftF.scale)*scale)
        
        shiftF.print_shiftfile(ROOT_DIRECT+'_shifts.txt')
        
        #### Make the grism shiftfile with the same shifts as 
        #### for the direct images
        threedhst.shifts.make_grism_shiftfile(asn_direct_file, asn_grism_file)
        
        #### Update the header WCS to reflect the shifts.
        #### !!! Need to turn updatewcs=no in Multidrizzle to use WCS updates
        ## Direct
        # shiftF = threedhst.shifts.ShiftFile(ROOT_DIRECT+'_shifts.txt')
        # for j in range(shiftF.nrows):
        #     im = pyfits.open(shiftF.images[j], mode='update')
        #     for ext in range(1,6):
        #         im[ext].header['CRPIX1'] += shiftF.xshift[j] 
        #         im[ext].header['CRPIX2'] += shiftF.yshift[j] 
        #     im.flush()
        #     
        # ## Grism
        # shiftF = threedhst.shifts.ShiftFile(ROOT_GRISM+'_shifts.txt')
        # for j in range(shiftF.nrows):
        #     im = pyfits.open(shiftF.images[j], mode='update')
        #     for ext in range(1,6):
        #         im[ext].header['CRPIX1'] += shiftF.xshift[j] 
        #         im[ext].header['CRPIX2'] += shiftF.yshift[j] 
        #     im.flush()
            
        #### Rerun multidrizzle with the new shifts
        # print "\n3DHST.process_grism: MultiDrizzle detection image, refined shifts\n"
        threedhst.showMessage('MultiDrizzle detection image, refined shifts')
        
        status = iraf.multidrizzle(input=asn_direct_file, \
           shiftfile=ROOT_DIRECT + '_shifts.txt', skysub=MDRIZ_SKYSUB, \
           output = '', final_scale = threedhst.options['DRZSCALE'],
           final_pixfrac = threedhst.options['USED_PIXFRAC'],
           Stdout=1)
        # status = iraf.multidrizzle(input=asn_direct_file, \
        #    shiftfile=ROOT_DIRECT + '_shifts.txt', skysub=MDRIZ_SKYSUB, \
        #    output = '', final_scale = iraf.INDEF,
        #    final_pixfrac = 1.0,
        #    Stdout=1)
        
        # ## No shiftfile, shifts contained in FLT WCS
        # status = iraf.multidrizzle(input=asn_direct_file, \
        #    shiftfile='', output = '', final_scale = INDEF, final_pixfrac = 1.0,
        #    updatewcs=no, Stdout=1)
                
        cleanMultidrizzleOutput()
        threedhst.currentRun['step'] = 'ALIGN_TO_REFERENCE'
        
    
    direct_mosaic = ROOT_DIRECT+'_drz.fits'
    
    #############################################
    #### Make a catalog with SExtractor
    #############################################
    
    #print "\n3DHST.process_grism: Making the catalog\n"
    threedhst.showMessage('Making the catalog')
    
    #### Get out SCI and WHT extensions of the Multidrizzle mosaic, 
    #### SExtractor was choking on multi-extension FITS files
    try:
        os.remove(ROOT_DIRECT+'_SCI.fits')
        os.remove(ROOT_DIRECT+'_WHT.fits')
    except:
        pass
    iraf.imcopy(input=direct_mosaic+'[1]',output=ROOT_DIRECT+'_SCI.fits')
    iraf.imcopy(input=direct_mosaic+'[2]',output=ROOT_DIRECT+'_WHT.fits')
    
    #### Run SExtractor on the direct image, with the WHT 
    #### extension as a weight image
    se = threedhst.sex.SExtractor()
    
    ## Set the output parameters required for aXe 
    ## (stored in [threedhst source]/data/aXe.param) 
    se.aXeParams()
    
    ## XXX add test for user-defined .conv file
    se.copyConvFile()
    
    se.overwrite = True
    se.options['CATALOG_NAME']    = ROOT_DIRECT+'_drz.cat'
    se.options['CHECKIMAGE_NAME'] = ROOT_DIRECT+'_seg.fits'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = ROOT_DIRECT+'_WHT.fits'
    se.options['FILTER']    = 'Y'
    
    #### Detect thresholds (default = 1.5)
    se.options['DETECT_THRESH']    = str(threedhst.options['DETECT_THRESH']) 
    se.options['ANALYSIS_THRESH']  = str(threedhst.options['ANALYSIS_THRESH']) 
    se.options['MAG_ZEROPOINT'] = str(threedhst.options['MAG_ZEROPOINT'])
    
    #### Run SExtractor
    status = se.sextractImage(ROOT_DIRECT+'_SCI.fits')
    
    threedhst.currentRun['step'] = 'SEXTRACT_CATALOG'
    
    #### Read catalog to keep around
    sexCat = threedhst.sex.mySexCat(ROOT_DIRECT+'_drz.cat')
    
    #### Trim objects on the edge of the detection image whose
    #### spectra will fall off of the grism image [turned off]
    #threedhst.regions.trim_edge_objects(sexCat)
    
    threedhst.currentRun['sexCat'] = sexCat
    
    #### Replace MAG_AUTO column name to MAG_F1392W
    sexCat.change_MAG_AUTO_for_aXe(filter='F1392W')
    sexCat.writeToFile()
    
    #### Make region file for SExtracted catalog
    threedhst.sex.sexcatRegions(ROOT_DIRECT+'_drz.cat', 
                                ROOT_DIRECT+'_drz.reg', format=2)
    
    #### Make a region file for the pointing itself
    threedhst.regions.asn_region(asn_direct_file)
    
    #### Make zeroth order region file [turned off]
    #threedhst.regions.make_zeroth(sexCat, outfile=ROOT_GRISM+'_zeroth.reg')
        
    threedhst.currentRun['step'] = 'PROCESS_CATALOG'
    
    ############################################################################
    ####   Run aXe scripts
    ############################################################################

    #### Initialize parameters, update the config file in CONF
    conf = Conf(threedhst.options['CONFIG_FILE'])
    threedhst.currentRun['conf'] = conf
    
    #### Need to scale the 0th order sensitivity curve
    if conf.params['SENSITIVITY_B'] == 'wfc3_abscal_IRg141_0th_sens.fits':
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
    conf.params['DRZROOT'] = ROOT_DIRECT
    conf.params['DRZRESOLA'] = threedhst.options['DRZRESOLA']
    conf.params['DRZSCALE'] = threedhst.options['DRZSCALE']
    conf.params['DRZPFRAC'] = threedhst.options['PIXFRAC']
    
    #### Workaround to get 0th order contam. in the right place for the fluxcube
    conf.params['BEAMB'] = '-220 220'    
    conf.writeto(ROOT_DIRECT+'_full.conf')
        
    CONFIG = ROOT_DIRECT+'_full.conf'
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
        
    #### Set the aXe environment variables
    set_aXe_environment(grating=threedhst.options['GRISM_NAME'])
    
    #### Make 'lis' file for input into aXe
    status = make_aXe_lis(asn_grism_file, asn_direct_file)
    shutil.move(prep_name(asn_grism_file),'..')

    #############################################
    #### Run aXe.iolprep to process the catalog as needed [[ in DATA directory]]
    #############################################
    flprMulti()
    status = iraf.iolprep(mdrizzle_ima=ROOT_DIRECT+'_drz.fits',
                 input_cat=ROOT_DIRECT+'_drz.cat', Stdout=1)
    threedhst.currentRun['step'] = 'IOLPREP'
    
    #############################################
    #### In root directory, run aXe.axeprep to subtract the sky background
    #### from the grism images
    #############################################
    os.chdir('../')
    flprMulti()
    status = iraf.axeprep(inlist=prep_name(asn_grism_file), configs=CONFIG,
        backgr=BACKGR, backims=SKY, mfwhm=3.0,norm="NO", Stdout=1)
    threedhst.currentRun['step'] = 'AXEPREP'
    
    #### Multidrizzle the grism images to flag cosmic rays
    os.chdir('./DATA')
    flprMulti()
    status = iraf.multidrizzle(input=asn_grism_file, 
                      shiftfile=ROOT_GRISM +  '_shifts.txt', 
                      output = '', final_scale = threedhst.options['DRZSCALE'], 
                      final_pixfrac = threedhst.options['USED_PIXFRAC'], 
                      skysub = no, Stdout=1)
    # status = iraf.multidrizzle(input=asn_grism_file, 
    #                   shiftfile=ROOT_GRISM +  '_shifts.txt', 
    #                   output = '', final_scale = iraf.INDEF, 
    #                   final_pixfrac = 1.0, 
    #                   skysub = no, Stdout=1)
    
    cleanMultidrizzleOutput()
    threedhst.currentRun['step'] = 'MULTIDRIZZLE_GRISM'
    
    #############################################
    #### Prepare the fluxcube for the contamination model
    #############################################
    swarpOtherBands()
    mag_zeropoint = '%5.2f' %threedhst.options['MAG_ZEROPOINT']
    
    if threedhst.options['OTHER_BANDS']:
        #### If OTHER_BANDS are available, use them for the fluxcube, start
        #### with the F140W image.
        lines = [ROOT_DIRECT+'_drz.fits, 1392.3, %s\n' %mag_zeropoint]
        for band in threedhst.options['OTHER_BANDS']:
            ## 'band' like ['../ACS/h_sz*drz_img.fits','F850LP',903.,24.84]
            if len(band) == 4:
                lines.append(ROOT_DIRECT+'_%s_drz.fits,' %band[1] +
                             ' %6.1f, %5.2f\n' %(band[2],band[3]))
    else:
        #### If no other bands present for the fluxcube, force flat spectrum
        #### in f_nu (constant AB mag)

        ## Make dummy copies of the direct image for the fluxcube
        shutil.copy(ROOT_DIRECT+'_drz.fits','flux1.fits')
        shutil.copy(ROOT_DIRECT+'_drz.fits','flux2.fits')
        lines = ['flux1.fits, 1442.3, %s\n' %mag_zeropoint,
                   ROOT_DIRECT+'_drz.fits, 1392.3, %s\n' %mag_zeropoint,
                   'flux2.fits, 1342.3, %s\n' %mag_zeropoint]
    
    #### Make a 'zeropoints.lis' file needed for fcubeprep   
    fp = open('zeropoints.lis','w')
    fp.writelines(lines)
    fp.close()
    
    #### Run aXe.fcubeprep
    # print "\n3DHST.proces_grism: iraf.FCUBEPREP\n"
    threedhst.showMessage('Prepare FLUXCUBE images with iraf.fcubeprep')
    
    flprMulti()
    status = iraf.fcubeprep(grism_image = ROOT_GRISM+'_drz.fits',
       segm_image = ROOT_DIRECT+'_seg.fits',
       filter_info = 'zeropoints.lis', AB_zero = yes, 
       dimension_info = '0,0,0,0', interpol="poly5", Stdout=1)
    
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
        
    status = iraf.axecore(inlist=prep_name(asn_grism_file), configs=CONFIG,
        back="NO",extrfwhm=5.0, drzfwhm=4.0, backfwhm=0.0,
        slitless_geom=geom_yn, orient=geom_yn, exclude="NO", lambda_mark=1392.0, 
        cont_model="fluxcube", model_scale=4.0, lambda_psf=1392.0,
        inter_type="linear", np=10, interp=-1, smooth_lengt=0, smooth_fwhm=0.0,
        spectr="NO", adj_sens=threedhst.options['AXE_ADJ_SENS'], weights="NO",
        sampling="drizzle", Stdout=1)
    
    threedhst.currentRun['step'] = 'AXECORE'
        
    #############################################
    #### aXe.tdrzprep - prepare for aXe drizzling
    #############################################
    
    # return 
    # 
    # #### Hack to correct for different scale used in MDRZ images
    # BEAMA = np.round( np.cast[float](conf.params['BEAMA'].split())*
    #         0.128254/np.float(threedhst.options['DRZSCALE']))
    # 
    # conf.params['BEAMA'] = '%d %d' %(BEAMA[0], BEAMA[1])
    # conf.writeto(ROOT_DIRECT+'_drz_hack.conf')
    # CONFIG_HACK = ROOT_DIRECT+'_drz_hack.conf'
    
    flprMulti()
    #print "\n3DHST.proces_grism: iraf.DRZPREP\n"
    threedhst.showMessage('Running iraf.drzprep')
    
    status = iraf.drzprep(inlist=prep_name(asn_grism_file), configs=CONFIG,
        opt_extr="YES", back="NO", Stdout=1)
    
    threedhst.currentRun['step'] = 'DRZPREP'
        
    #############################################
    #### aXe.taxedrizzle - drizzle combine the dithers.  This step also 
    #### produces the 1D extraction
    #############################################
    threedhst.options['DRIZZLE_PATH'] = os.environ['AXE_DRIZZLE_PATH']
    #print "\n3DHST.proces_grism: iraf.AXEDRIZZLE (#1)\n"
    threedhst.showMessage('Running iraf.axedrizzle')
    
    flprMulti()
    status = iraf.axedrizzle(inlist=prep_name(asn_grism_file), configs=CONFIG,
                    infwhm=5.0,outfwhm=4.0, back="NO", makespc="YES",
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
        iraf.axedrizzle(inlist=prep_name(asn_grism_file), configs=CONFIG,
                        infwhm=5.0,outfwhm=4.0, back="NO", makespc="YES",
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
        cont = pyfits.open('../OUTPUT_G141/'+expi+'_flt_2.CONT.fits')
        flt[1].data = cont[1].data.copy()
        flt.flush()
    asn_grism.product = None #ROOT_GRISM+'CONT'
    asn_grism.writeToFile(ROOT_GRISM+'CONT_asn.fits')
    flprMulti()
    status = iraf.multidrizzle(input = ROOT_GRISM+'CONT_asn.fits', \
       output = '', shiftfile = ROOT_GRISM + '_shifts.txt', \
       final_scale = threedhst.options['DRZSCALE'], 
       final_pixfrac = threedhst.options['USED_PIXFRAC'], skysub = no,
       static=no, driz_separate=no,median=no, blot=no, driz_cr=no, Stdout=1)
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
        rmfiles = glob.glob('*FLX.fits')
        rmfiles = glob.glob('flux?.fits')
        rmfiles.extend(glob.glob('*_flt.fits'))
        rmfiles.extend(glob.glob('*flt_1.cat'))
        rmfiles.extend(glob.glob(ROOT_DIRECT+'*[SW][CH]?.fits'))
        #rmfiles.extend(glob.glob('*coeffs1.dat'))
        #rmfiles.extend(glob.glob('threedhst_auto.*'))
        rmfiles.extend(glob.glob('zeropoints.lis'))
        rmfiles.extend(glob.glob('default.conv'))
        #for expi in asn_grism.exposures:
        #    rmfiles.extend(glob.glob('../OUTPUT_G141/'+expi+'*'))
        
        if len(rmfiles) > 0:
            for rmfile in rmfiles:
                os.remove(rmfile)
        
        threedhst.currentRun['step'] = 'CLEANUP_DATA'
    
    ############################################################################
    #### Make output webpages with spectra thumbnails    
    ############################################################################
    if threedhst.options['MAKE_WEBPAGE']:
        threedhst.plotting.make_data_products(ROOT_DIRECT, ROOT_GRISM)
                
    #############################################
    #### Write parameter log
    #############################################

    os.chdir('../')

    logfile = ROOT_DIRECT+'.threedhst.info'
    print '\nthreedhst: Parameter log in <%s>.\n' %logfile
    
    threedhst.showOptions(to_file = logfile)
    if threedhst.options['MAKE_WEBPAGE']:
        threedhst.showOptions(to_file = "./HTML/"+logfile)
        
    #############################################
    #### Done!
    #############################################
    
    print '\nthreedhst: cleaned up and Done!\n'
    threedhst.currentRun['step'] = 'DONE'

def fresh_flt_files(asn_filename, from_path="../RAW"):
    """
    """ 
    import threedhst.dq
    
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
        
        print exp
        
        fi.writeto('./'+exp+'_flt.fits', clobber=True)
        
        #### Apply DQ mask (.mask.reg), if it exists
        if threedhst.options['PYSAO_INSTALLED']:
            threedhst.dq.apply_dq_mask(os.path.basename(fits_file), addval=2048)
    
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

    other_bands = threedhst.options['OTHER_BANDS']
        
    #### Nothing found
    if not other_bands:
        return None
    
    #### Other bands are defined, use them
    
    #### Initialize SWarp
    sw = threedhst.sex.SWarp()
    sw._aXeDefaults()
    ## get reference image parameters
    sw.swarpMatchImage(ROOT_DIRECT+'_drz.fits',extension=1)  
    ## swarp the reference image to istelf
    sw.swarpImage(ROOT_DIRECT+'_drz.fits[1]',mode='wait')    
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
                  ROOT_DIRECT,band[0])
                  
        sw.swarpImage(otherImages,mode='direct')   
        
        #### Put the result from "coadd.fits" into the first extension 
        #### of a copy of the direct drz image.
        direct = pyfits.open(ROOT_DIRECT+'_drz.fits')
        new = pyfits.open('coadd.fits')
        direct[1].data = new[0].data
        direct.writeto(ROOT_DIRECT+'_'+band[1]+'_drz.fits', clobber=True)
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
    
    #### Check that grism and direct ASN tables have same # of entries
    if len(asn_grism.exposures) != len(asn_direct.exposures):
#         print """
# 3D-HST / make_aXe_lis: Number of grism exposures (%d) in %s is different from
#                        the number of direct images (%d) in %s.
#               """ %(len(grism_files), asn_grism, len(direct_files), asn_direct)
        
        threedhst.showMessage("""Number of grism exposures (%d) in %s is different from
                       the number of direct images (%d) in %s.
                      """ %(len(grism_files), asn_grism, len(direct_files), asn_direct))
        
        return False
    
    #### Make the lis file
    NFILE = len(asn_grism.exposures)
    outfile = prep_name(asn_grism_file)
    fp = open(outfile,'w')
    for i in range(NFILE):
        fp.write("%s_flt.fits %s_flt_1.cat %s_flt.fits 0.0\n" 
              %(asn_grism.exposures[i], asn_direct.exposures[i], 
                asn_direct.exposures[i]))
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
    rmfiles.extend(glob.glob('*sci1_blt.fits'))
    rmfiles.extend(glob.glob('*flt*mask1.fits'))
    #files.extend(glob.glob('*coeffs1.dat'))
    #rmfiles.extend(glob.glob('ib*.run'))
    rmfiles.extend(glob.glob('ib*_med.fits'))
    if len(rmfiles) > 0:
        for rmfile in rmfiles:
            os.remove(rmfile)
            
class MultidrizzleRun():
    """
MultidrizzleRun(root='IB3728050')
    
    Read a .run output file created by MultiDrizzle.
    """
    def __init__(self, root='IB3728050'):
        import numpy as np
        
        runfile = root+'.run'
        self.root = root
        
        self.flt = []
        self.xsh = []
        self.ysh = []
        self.rot = []
        self.scl = 1.
        
        for line in open(runfile,'r'):
            if line.startswith('drizzle.scale'):
                self.scl = line.split()[2]
                
            if line.startswith('drizzle '):
                spl = line.split()
                self.flt.append(spl[1].split('.fits')[0])
                for tag in spl:
                    if tag.startswith('xsh'):
                        self.xsh.append(np.float(tag.split('=')[1]))
                    if tag.startswith('ysh'):
                        self.ysh.append(np.float(tag.split('=')[1]))
                    if tag.startswith('rot'):
                        self.rot.append(np.float(tag.split('=')[1]))
                        
    def blot_back(self, ii=0, SCI=True, WHT=True, orig_path='./'):
        """
blot_back(ii=0, SCI=True, WHT=True)
        
    `blot` the final multidrizzle DRZ file to the pixels of
    frame #ii as listed in the run file.
        """
        flt_orig = pyfits.open(orig_path+self.flt[ii]+'.fits.gz')
        exptime = flt_orig[0].header.get('EXPTIME')
        
        iraf.delete(self.flt[ii]+'.BLOT.*.fits')
        if SCI:
            iraf.blot(data=self.root+'_drz.fits[1]',
                outdata=self.flt[ii]+'.BLOT.SCI.fits', scale=self.scl,
                coeffs=self.flt[ii]+'_coeffs1.dat', xsh=self.xsh[ii], 
                ysh=self.ysh[ii], 
                rot=self.rot[ii], outnx=1014, outny=1014, align='center', 
                shft_un='input', shft_fr='input', in_un='cps', out_un='cps', 
                interpol='poly5', sinscl='1.0', expout=exptime, 
                expkey='EXPTIME',fillval=0.0)
        
        if WHT:
            iraf.blot(data=self.root+'_drz.fits[2]',
                outdata=self.flt[ii]+'.BLOT.WHT.fits', scale=0.46,
                coeffs=self.flt[ii]+'_coeffs1.dat', xsh=self.xsh[ii], 
                ysh=self.ysh[ii], 
                rot=self.rot[ii], outnx=1014, outny=1014, align='center', 
                shft_un='input', shft_fr='input', in_un='cps', out_un='cps', 
                interpol='poly5', sinscl='1.0', expout=exptime, 
                expkey='EXPTIME',fillval=0.0)

def flprMulti(n=3):
    """
    flprMulti(n=3)
    
    Run iraf.flpr() `n` times.
    """
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
    print """
# step: %s
asn_grism_file = threedhst.currentRun['asn_grism_file']
asn_direct_file = threedhst.currentRun['asn_direct_file']
ROOT_GRISM = threedhst.options['ROOT_GRISM']
ROOT_DIRECT = threedhst.options['ROOT_DIRECT']
sexCat = threedhst.sex.mySexCat(ROOT_DIRECT+'_drz.cat')
if 'conf' in threedhst.currentRun.keys():
    conf = threedhst.currentRun['conf']
SPC = threedhst.plotting.SPCFile(ROOT_DIRECT+'_2_opt.SPC.fits')
""" %threedhst.currentRun['step']
    raise IOError

def multidrizzle_defaults():
    """
multridrizzle_defaults():

    [not used] 
    Set multidrizzle default parameters
    """
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

def fix_g141_sky_image(conf_path='./', verbose=True):
    """
fix_g141_sky_image(conf_path='./', verbose=True)
    
    The G141 average sky image, WFC3.IR.G141.sky.V1.0.fits,  has lots of 
    holes with pixel value = 0,  apparently more than pixels flagged in the BPM.
    
    Fill these in by taking the median of (nonzero) neighboring pixels.
    
    The result is an image, g141_fixed_sky.fits, put in the directory defined
    by `conf_path` (WFC3.IR.G141.sky.V1.0.fits is read from `conf_path`).
    
    """
    import pyfits
    import numpy as np
    
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
            self.path = '../CONF/'
        else:
            self.path = './CONF/'
        
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
        import numpy as np
        for key in self.params.keys():
            param = self.params[key]
            if type(param) is not type('xxx'):
                param = np.str(param)
            
            #### New parameter?
            if self._pline.has_key(key):
                self.lines[self._pline[key]] = key + ' ' + param +'\n'
            else:
                self.lines.append(key + ' ' + param +'\n')
        
        self._processLines()
    
    def writeto(self, output='tmp.conf'):
        """
        writeto(output='tmp.conf')
        """ 
        self._getPath()
        self._assignPars()
        fp = open(self.path+output,'w')
        fp.writelines(self.lines)
        fp.close()
    
    