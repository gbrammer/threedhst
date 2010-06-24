"""
3DHST.utils

Utilities for reading/writing ASN and DS9/region files.

"""
# $URL$
# $Rev$
# $Author$
# $Date$

import os,pyfits
import numpy as np

def read_asn(asn_file):
    """
asn = read_asn(asn_file)
    
Read the contents of an ASN file.
    
>>> print asn.names
('MEMNAME', 'MEMTYPE', 'MEMPRSNT')
>>> print asn.field('MEMNAME')[0]
'IB3714FUQ'
    
    """
    hdulist = pyfits.open(asn_file)
    asn = hdulist[1].data
    return asn
####### read_asn

def asn_region(asn_file):
    """
asn_region(asn_file)
    
Create a DS9 region file for the exposures defined in an ASN file.
    
    """
    ##### Output file
    output_file = asn_file.split('.fits')[0]+'.pointing.reg'
    fp = open(output_file,'w')
    fp.write('fk5\n') ### WCS coordinates
    ##### Read ASN file
    asn = read_asn(asn_file)
    ##### Find exposures, which have MEMTYPE = EXP-DTH
    expIDX = np.where(asn.field('MEMTYPE').find('EXP-DTH') == 0)[0]
    NEXP = expIDX.shape[0]
    RAcenters = np.arange(NEXP*1.)
    DECcenters = np.arange(NEXP*1.)
    for i in range(NEXP):
        ###### Get FITS header and keywords from extension 1 (SCI)
        flt_file = asn.field('MEMNAME')[expIDX[i]].lower()+'_flt.fits'
        try:
            fptest = open(flt_file,'r')
        except IOError:
            try:
                flt_file += '.gz'
                fptest = open(flt_file,'r')
            except IOError:
                raise IOError('[Errno 2] No such file or directory: \'%s[.gz]\'' %flt_file.split('.gz')[0] )
        fptest.close()
        
        regX, regY = wcs_polygon(flt_file,extension=1)
        line = "polygon(%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f)" \
            %(regX[0],regY[0],regX[1],regY[1],regX[2],regY[2],regX[3],regY[3])
        
        RAcenters[i] = np.mean(regX)
        DECcenters[i] = np.mean(regY)
        fp.write(line+' # color=magenta\n')

    ##### Text label with ASN filename
    fp.write('# text(%10.6f,%10.6f) text={%s} color=magenta' \
        %(np.mean(RAcenters),np.mean(DECcenters),asn_file.split('_asn.fits')[0]))
    fp.close()
    print "3D-HST / ASN_REGION: %s\n" %(output_file)
    
####### asn_region

def wcs_polygon(fits_file, extension=1):
    """    
X, Y = wcs_polygon(fits_file, extension=1)
    
Calculate a DS9/region polygon from WCS header keywords

    """
    from math import cos,pi
    
    ##### Open the FITS file
    hdulist = pyfits.open(fits_file) 
    ##### Get the header
    try:
        sci = hdulist[extension].header
    except IndexError:
        print 'ERROR 3D-HST/wcs_polygon:\n'+\
              'Extension #%d out of range in %s' %(extension, fits_file)
        raise
        
    NAXIS = [sci['NAXIS1'],sci['NAXIS2']]
    CRPIX = [sci['CRPIX1'],sci['CRPIX2']]
    CRVAL = [sci['CRVAL1'],sci['CRVAL2']]
    cosDec = cos(CRVAL[1]/180*pi)
    ##### Make region polygon from WCS keywords
    regX = CRVAL[0] + ( (np.array([0,NAXIS[0],NAXIS[0],0])-CRPIX[0])*sci['CD1_1'] + \
                        (np.array([0,0,NAXIS[1],NAXIS[1]])-CRPIX[1])*sci['CD1_2'] ) / cosDec
    regY = CRVAL[1] + ( (np.array([0,NAXIS[0],NAXIS[0],0])-CRPIX[0])*sci['CD2_1'] + \
                        (np.array([0,0,NAXIS[1],NAXIS[1]])-CRPIX[1])*sci['CD2_2'] )
    return regX, regY
    
    #return line, CRVAL
####### wcs_polygon
    
def region_mask(shape,px,py):
    """
mask = region_mask(image.shape,px,py)

Make a mask image where pixels within the polygon defined by px_i, py_i
are set to 1.

Note: something like this could be used to flag grism 0th order contaminants
"""
    NX=shape[0]
    NY=shape[1]
    y,x = np.mgrid[1:NX+1,1:NY+1]
    ##### Close polygons
    NPOLY = px.shape[0]
    tmp_px = np.append(px,px[0])
    tmp_py = np.append(py,py[0])
    theta = np.zeros((NX,NY),dtype=np.float)
    for i in np.arange(NPOLY):
        ##### Dot, cross products
        X1 = tmp_px[i] - x
        Y1 = tmp_py[i] - y 
        X2 = tmp_px[i+1] - x
        Y2 = tmp_py[i+1] - y
        dp = X1*X2 + Y1*Y2
        cp = X1*Y2 - Y1*X2
        theta += np.arctan(cp/dp) - np.pi*((cp < 0) & (dp < 0)) + np.pi*((cp > 0) & (dp < 0))
    ##### Set up mask
    dq = np.zeros((NX,NY),dtype=np.int)
    flag_idx = np.where(np.abs(theta) > np.pi)
    dq[flag_idx] = 1
    return dq
##### region_mask

def point_in_polygon(x,y,px,py):
    """
test = point_in_polygon(x,y,px,py)

Test if coordinates (x,y) are inside polygon defined by (px, py)

[http://www.dfanning.com/tips/point_in_polygon.html, translated to Python]

"""    
    N = px.shape[0]
    ##### Close polygons
    tmp_px = np.append(px,px[0])
    tmp_py = np.append(py,py[0])
    ##### Counters
    i = np.arange(N)
    ip = np.arange(N)+1
    ##### Dot, cross products
    X1 = tmp_px[i] - x
    Y1 = tmp_py[i] - y 
    X2 = tmp_px[ip] - x
    Y2 = tmp_py[ip] - y
    dp = X1*X2 + Y1*Y2
    cp = X1*Y2 - Y1*X2
    theta = np.arctan(cp/dp) - np.pi*((cp < 0) & (dp < 0)) + np.pi*((cp > 0) & (dp < 0))
    if np.abs(np.sum(theta)) > np.pi:
        return True
    else:
        return False
    
##### inside_polygon

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
            spl = np.float_(np.array(region[region.find('(')+1:region.find(')')].split(',')))
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
##### check_data_quality

def asn_file_info(asn_file, verbose=0):
    """
asn_file_info(asn_file, verbose=0)

Get header information from files defined in an ASN table

"""
    #asn_file = 'ib6o23020_asn.fits'
    fp_asn = pyfits.open(asn_file)
    asn_data = fp_asn[1].data
    roots = asn_data.field('MEMNAME')[np.where(asn_data.field('MEMTYPE') == 'EXP-DTH')]
    NFILES = roots.shape[0]
    lines = ['# %s' %asn_file]
    lines.append('# flt_file  filter  exptime  date_obs  time_obs pos_targ1 pos_targ2')
    ##### Loop through flt files in ASN list
    for i in range(NFILES):
        flt_file = roots[i].lower()+'_flt.fits'
        try:
            fptest = open(flt_file,'r')
        except IOError:
            try:
                flt_file += '.gz'
                fptest = open(flt_file,'r')
            except IOError:
                raise IOError('[Errno 2] No such file or directory: \'%s[.gz]\'' %flt_file.split('.gz')[0] )
        fptest.close()
        #print flt_file
        fp_flt = pyfits.open(flt_file)
        ##### Get general information from extension 0
        fp_header = fp_flt[0].header
        line = '%s %6s %7.1f %s %s %7.2f %7.2f' %(flt_file,fp_header['FILTER'],fp_header['EXPTIME'],
                                              fp_header['DATE-OBS'],fp_header['TIME-OBS'],
                                              fp_header['POSTARG1'],fp_header['POSTARG2'])
        lines.append(line)

    #### Print to stdout
    if verbose > 0:
        for line in lines:
            print line
    
##### asn_file_info