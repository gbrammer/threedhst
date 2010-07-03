"""
3DHST.plotting

Utilities for plotting grism spectra.

"""
# $URL $
# $Rev $
# $Author $
# $Date $

import pyfits

import matplotlib
import matplotlib.pyplot as pyplot
import pylab

class SPCFile(object):
    """
    SPCFile
    
    Class for reading and plotting spectra from an aXe
    
    """
    