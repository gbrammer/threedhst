"""
3DHST

Utilities for processing WFC3 Grism exposures
from the 3D-HST large program.

"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import utils
import shifts
import sex
import process_grism
import plotting
import gmap
import regions
import spec1d   

options = {}
currentRun = {}
currentRun['step'] = 'INIT'

def defaultOptions():
    """
defaultOptions()
    
    Set THREEDHST default options.
    
    To see the defaults, run
    
    >>> threedhst.defaultOptions()
    >>> threedhst.showOptions()
    """
    #### Delete all keywords and reset
    for key in options.keys():
        pop = options.popitem()
        
    options['DETECT_THRESH'] = 5.     ## Default 1.5
    options['ANALYSIS_THRESH']  = 5.  ## Default 1.5
    options['GRISM_NAME'] = 'G141'
    options['MAG_ZEROPOINT'] = 26.46
            
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
    options['CONFIG_FILE'] = 'WFC3.IR.G141.V1.0.conf'
    options['SKY_BACKGROUND'] = 'WFC3.IR.G141.sky.V1.0.fits'
    options['DRZRESOLA'] = '46.5'
    options['DRZSCALE'] = '0.128254'
    
    #### aXe extraction geometry
    #### currently set slitless_geom=NO, orient=NO in aXecore
    #### to get the 2D spectra to line up with the orientation
    #### of the direct thumbnail.
    options['FULL_EXTRACTION_GEOMETRY'] = False
    #### aXe adjust sensitivity - convolve grism throughput with source profile
    options['AXE_ADJ_SENS'] = "YES"
    #### aXe extract with "optimal weights"
    options['AXE_OPT_EXTR'] = "YES"

    #### Second drizzle to flag and remove cosmic rays
    options['RUN_DRZREJ'] = False
    
    #### Clean intermediate files from MultiDrizzle and aXe
    options['CLEAN_UP'] = True
    
    #### Output webpage with google map for browsing spectra
    options['MAKE_WEBPAGE'] = True
    ## Image format for webpage
    options['WEB_IMAGE_FORMAT'] = 'png'
    # options['IMAGE_FORMAT'] = 'svgz'

#############################################
#### Set the default options    
#############################################
defaultOptions()

def showOptions(to_file=None):
    """
    printOptions()
    
    Show the current THREEDHST option set.
    """
    import time
    
    if to_file is None:
        for key in options.keys():
            print '%s = %s' %(key,str(options[key]))
    else:
        fp = open(to_file,"w")
        fp.write('#######################################\n')
        fp.write('###                                 ###\n')
        fp.write('###    threedhst   %s      ###\n' %__version__)
        fp.write('###                                 ###\n')
        fp.write('###    %s     ###\n' %time.asctime())
        fp.write('###                                 ###\n')
        fp.write('#######################################\n')
        for key in options.keys():
            fp.write('%s = %s\n' %(key,str(options[key])))
        fp.close()
    