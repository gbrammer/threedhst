"""
3DHST

Utilities for processing WFC3 Grism exposures
from the 3D-HST large program.

"""
#   $URL $
#   $Rev $
#   $Author $
#   $Date $

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
    >>> """
    #### Delete all keywords and reset
    for key in options.keys():
        pop = options.popitem()
        
    options['DETECT_THRESH'] = 5.     ## Default 1.5
    options['ANALYSIS_THRESH']  = 5.  ## Default 1.5
    options['GRISM_NAME'] = 'G141'
    
defaultOptions()

def printOptions():
    """
    printOptions()
    
    Show the current THREEDHST option set.
    """
    for key in options.keys():
        print '%s = %s' %(key,str(options[key]))