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

options = {}

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
    options['CONFIG_FILE'] = 'WFC3.IR.G141.V1.0.conf'
    options['SKY_BACKGROUND'] = 'WFC3.IR.G141.sky.V1.0.fits'
    
    options['THIS_RUN'] = 0     
    options['LAST_RUN'] = 0     
    
defaultOptions()

def showOptions():
    """
    printOptions()
    
    Show the current THREEDHST option set.
    """
    for key in options.keys():
        print '%s = %s' %(key,str(options[key]))