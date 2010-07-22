"""
3DHST.dq

Check DQ (satellite trails, scattered light)
of direct/grism exposures.

"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import pyfits

import numpy as np
import threedhst

def check_data_quality(flt_file):
    """
check_data_quality(flt_file)
    
Display an image and check DQ extension (e.g. for satellite trails)
    
    """
    ##### Start a DS9 instance with pysao
    import pysao
    ds9 = pysao.ds9()
    ##### Open FITS file
    fi = pyfits.open(flt_file)
    ### Display SCI extension [1]
    ds9.frame(1)
    ds9.view(fi[1])
    ds9.scale(-0.1,1)
    ### Display DQ extension [3]
    ds9.frame(2)
    ds9.view(fi[3])
    ds9.scale(0,100)
    ### DS9 Tile and return to Frame 1
    ds9.set('tile yes')
    ds9.set('tile grid')
    ds9.set('zoom to fit')
    ds9.frame(1)
    ds9.set('zoom to fit')
    ##### Ask at prompt if we should define regions and continue
    test = raw_input('3DHST/check_data_quality: Define DQ region y/[n]? ')
    if (test == '') | (test.lower().startswith('n')):
        return False
    else:
        ds9.set('regions shape polygon')
        dummy = raw_input('Define region polygon in DS9 and hit [return] ')
    ##### Define region
    ds9.set('regions system image')
    regions = ds9.get('regions source')
    ##### Initialize DQ image
    dq = np.zeros(fi[1].data.shape,dtype=np.int)
    ##### Loop through user-defined regions
    for region in regions.split('\n'):
        if region.strip().startswith('polygon'):
            #region = 'polygon(375.05333,642.2,465.18667,642.2,751.36,709.8,393.08,326.73333,210.56,461.93333,465.18667,552.06667,375.05333,552.06667,221.82667,509.25333)'
            spl = np.float_(np.array(
                     region[region.find('(')+1:region.find(')')].split(',')
                     ))
            px = spl[0::2]
            py = spl[1::2]
            dq += region_mask(fi[1].data.shape,px,py)
    ##### Set DQ bit
    dq[np.where(dq > 0)] = 2048
    ##### Show region mask in Frame 2        
    ds9.frame(2)
    ds9.view_array(fi[3].data+dq)
    ds9.scale(0,100)
    ##### Save defined regions in output file, [flt_file]+'.mask.reg'
    fp = open(flt_file+'.mask.reg','w')
    fp.write(regions)
    fp.close()
    return True

