"""
3DHST

Utilities for processing WFC3 Grism exposures
from the 3D-HST large program.

"""

__version__ = "1.0"


#### Package imports
# from . import utils
# from . import shifts
# from . import sex
# from . import process_grism
# from . import plotting
# from . import plotting
# from . import gmap
# from . import regions
# from . import spec1d   
# from . import TerminalController

options = {}
currentRun = {}
currentRun['step'] = 'INIT'

noNewLine = '\x1b[1A\x1b[1M'

def showMessage(msg, warn=False):
    """
showMessage(msg)
    
    Print a system message formatted like:
    
    ***********************************
    *  THREEDHST.`module`.`function`  *
    ***********************************
    
    `msg`
    
    ***********************************
    """       
    import os
    import sys
    import inspect
    import time
    
    from threedhst.TerminalController import TerminalController
    
    calling_function_name = sys._getframe(1).f_code.co_name    
    module_name =  os.path.basename(inspect.stack()[1][1]).split('.py')[0]
    
    term = TerminalController()
    
    char = '='
    substr = 'THREEDHST.'+module_name+'.'+calling_function_name
    
    substr = char+'  '+substr+'  '+char
    
    NL = len(substr)
    topbar = char*NL
    
    t0 = time.localtime()
    theDate = '%0d/%0d/%0d %0d:%0d' %(t0[0],t0[1],t0[2],t0[3],t0[4])
    N2 = (NL-2-len(theDate))/2
    botbar = char*N2+' '+theDate+' '+char*(NL-N2-2-len(theDate))
    
    if warn:
        text_color = term.WHITE
        bg_color = term.BG_RED
    else:
        text_color = term.BLUE
        bg_color = term.BG_WHITE
        
    print((bg_color+text_color+term.BOLD+'\n'+topbar+
           '\n'+substr+'\n'+topbar+'\n\n'+term.NORMAL+
           msg+'\n\n'+
           bg_color+text_color+term.BOLD+botbar+'\n'+term.NORMAL))

def defaultOptions():
    """
defaultOptions()
    
    Set THREEDHST default options.
    
    To see the defaults, run
    
    >>> threedhst.defaultOptions()
    >>> threedhst.showOptions()
    """    
    #showMessage('Initializing THREEDHST parameters')
    #### Delete all keywords and reset
    for key in list(options.keys()):
        pop = options.popitem()
    
    #### ACS rather than WFC3 grism
    options['ACS_G800L'] = False
    
    #### Optionally supply a ready-made direct image
    #### that you will match the grism image to.  Useful
    #### for when you have an existing mosaic that extends
    #### beyond the coverage of a few simple direct images.
    options['PREFAB_DIRECT_IMAGE'] = None
    options['PREFAB_GRISM_IMAGE'] = None
    
    #### Use a particular catalog rather than creating one with 
    #### SExtractor.
    options['FORCE_CATALOG'] = None
    
    #### Use the direct mosaic as a reference for making the catalog.
    #### Assumes one-to-one grism/direct exposures.  Set to 'GRISM'
    #### if you want to use the grism mosaic as a reference for splitting
    #### the individual exposures, i.e. with PREFAB direct/images
    options['CATALOG_MODE'] = 'DIRECT'
    
    #### Directories
    options['DRIZZLE_PATH'] = './DRIZZLE_G141/'
    options['PATH_TO_RAW'] = '../RAW/'
    
    #### General detection parameters
    options['DETECT_THRESH'] = 5.     ## Default 1.5
    options['ANALYSIS_THRESH']  = 5.  ## Default 1.5
    options['GRISM_NAME'] = 'G141'
    options['MAG_ZEROPOINT'] = 26.46
    options['FILTWAVE'] = 1392.
    
    #### Make shiftfiles (required for multidrizzle)
    options['MAKE_SHIFTFILES'] = True
    #### WCS alignment image
    options['ALIGN_IMAGE'] = None
    options['ALIGN_GEOMETRY'] = "shift"
    
    #### Add other bands to the fluxcube
    options['OTHER_BANDS'] = []
    
    ####### [this doesn't do anything now]
    #### For fluxcube, if these remain 'None', use files 
    #### created from the internal SExtractor run 
    #options['CATALOG_FILE'] = None
    #options['SEGMENTATION_MAP'] = None
    
    #### Config options
    #options['CONFIG_FILE'] = 'WFC3.IR.G141.V1.0.conf'
    options['CONFIG_FILE'] = 'WFC3.IR.G141.V2.0.conf'
    options['SKY_BACKGROUND'] = 'WFC3.IR.G141.sky.V1.0.fits'
    options['DRZRESOLA'] = '46.5'
    options['DRZSCALE'] = '0.128254'
    
    #### Axe extraction geometry to get objects outside the edges
    #### of the grism images.  Currently doesn't work correctly and 
    #### the model is offset w.r.t. the observed spectra.
    options['AXE_EDGES'] = "0,0,0,0"
    
    #### Use higher-resolution drizzle
    options['DRZSCALE'] = '0.06'
    options['PIXFRAC'] = '0.8'
    #### If the number of frames in an ASN table
    #### is < NFRAME_FOR_PIXFRAC, use PIXFRAC=1.0
    options['NFRAME_FOR_PIXFRAC'] = 4
    options['USED_PIXFRAC'] = '1.0'
    
    #### Limiting (direct) magnitude for objects run through the grism 
    #### reduction.  This is useful for cases where you want a very low 
    #### DETECT_THRESH to get good segmentation images but don't want to
    #### extract spectra for all of the faint sources.
    options['LIMITING_MAGNITUDE'] = 99.
    
    #### aXe extraction geometry
    #### currently set slitless_geom=NO, orient=NO in aXecore
    #### to get the 2D spectra to line up with the orientation
    #### of the direct thumbnail.
    options['FULL_EXTRACTION_GEOMETRY'] = False
    #### aXe adjust sensitivity - convolve grism throughput with source profile
    options['AXE_ADJ_SENS'] = "YES"
    #### aXe extract with "optimal weights"
    options['AXE_OPT_EXTR'] = "YES"
    
    #### Use a local installation of taxe rather than the STSDAS aXe
    options['USE_TAXE'] = False
    
    #### Second drizzle to flag and remove cosmic rays
    options['RUN_DRZREJ'] = False
    
    #### Clean intermediate files from MultiDrizzle and aXe
    options['CLEAN_UP'] = True
    
    #### Output webpage with google map for browsing spectra
    options['MAKE_WEBPAGE'] = True
    ## Image format for webpage
    options['WEB_IMAGE_FORMAT'] = 'png'
    
    #### Path to WFC3-IR Persistence files
    options['FLT_PERSISTENCE_PATH'] = '/3DHST/Spectra/Work/PERSISTENCE/All'
    options['FLT_PERSISTENCE_THRESH'] = 0.6
    options['FLT_PERSISTENCE_SUBTRACT'] = False
    options['FLT_PERSISTENCE_FILTER'] = 3

#############################################
#### Set the default options    
#############################################
defaultOptions()

try:
    import pysao
    options['PYSAO_INSTALLED'] = True
except:
    options['PYSAO_INSTALLED'] = False
    print('\nWARNING: No pysao installation found.  `threedhst.dq` won\'t work\n but the reduction scripts should be OK.\n')

def showOptions(to_file=None):
    """
    printOptions()
    
    Show the current THREEDHST option set.
    """
    import time
    
    if to_file is None:
        for key in list(options.keys()):
            print(('%s = %s' %(key,str(options[key]))))
    else:
        fp = open(to_file,"w")
        fp.write('#######################################\n')
        fp.write('###                                 ###\n')
        fp.write('###    threedhst   %s      ###\n' %__version__)
        fp.write('###                                 ###\n')
        fp.write('###    %s     ###\n' %time.asctime())
        fp.write('###                                 ###\n')
        fp.write('#######################################\n')
        for key in list(options.keys()):
            fp.write('%s = %s\n' %(key,str(options[key])))
        fp.close()
    
    
    