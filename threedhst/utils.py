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

def find_fits_gz(fits_file, hard_break = True):
    """
the_file = find_fits_gz(fits_file, hard_break = True)
    
With ``fits_file`` being some filename with an extension
``.fits`` (ib3713wvq_flt.fits), check to see if the file 
itself or its gzipped counterpart (fits_file+'.gz')
exists.  
    
If neither is found, 
    hard_break = True  : raise an IOError
    hard_break = False : return None
    
    """
    import os
    if os.path.exists(fits_file):
        return fits_file
    if os.path.exists(fits_file+'.gz'):
        return fits_file+'.gz'
    #### File not found.  Either raise an error or return None
    if hard_break:
        raise IOError('File %s[.gz] not found in %s' %(fits_file, os.getcwd()))
    else:
        return None


class ASNFile(object):
    """
ASNFile()
        
    Class for handling ASN fits files.
        
    >>> asn = ASNFile(file='ib3701050_asn.fits')
    >>> asn.exposures
    ['ib3701ryq', 'ib3701sbq', 'ib3701sdq', 'ib3701sqq']
    >>> asn.product
    'IB3701050'
    >>> asn.
     """
    def _read_asn_file(self):
        """
_read_asn_file(self)
        
    Read an ASN FITS file (self.file).
        """
        import numpy as np
        from warnings import warn
        
        self.in_fits = pyfits.open(self.file)
        data = self.in_fits[1].data
        self.header = self.in_fits[0].header
        
        names = data.field('MEMNAME')
        types = data.field('MEMTYPE')
        
        ##### Exposures
        exp_idx  = np.where(types == 'EXP-DTH')
        if exp_idx[0].shape[0] == 0:
            warn ('ASN file %s has no EXP-DTH items')
        else:
            self.exposures = []
            for exp in names[exp_idx]:
                self.exposures.append(exp.lower())
        
        ##### Products
        prod_idx = np.where(types == 'PROD-DTH')
        if prod_idx[0].shape[0] != 1:
            warn ('ASN file %s has N != 1 PROD-DTH items' %self.file )
            self.product = None
        else:
            self.product = names[prod_idx[0]][0].upper()
    
    
    def __init__(self, file=None):
        self.file = file
        self.exposures = []
        self.product = None
        if file:
            self._read_asn_file()
    
    
    def writeToFile(self, out_file=None, clobber=True):
        """
writeToFile(self,out_file=None, clobber=True)
        """
        if not out_file:
            print "USAGE:: writeToFile(self,out_file='output_asn.fits')"
        else:
            nexp  = self.exposures.__len__()
            if self.product:
                nprod = 1
            else:
                nprod = 0
            nrows = nexp + nprod
            #### Primary HDU
            hdu = self.in_fits[0].copy()
            #### BinTable HDU
            tbhdu = pyfits.new_table(self.in_fits[1].columns, nrows=nrows)
            for i in range(nexp):
                tbhdu.data[i] = (self.exposures[i].upper(), 'EXP-DTH', True)
            if nprod > 0:
                tbhdu.data[i+1] = (self.product, 'PROD-DTH', True)
            tbhdu.header = self.in_fits[1].header.copy()
            #### Create HDUList and write it to output file
            self.out_fits = pyfits.HDUList([hdu,tbhdu])
            self.out_fits.writeto(out_file, clobber=clobber)
    
    
    def showContents(self):
        """
showContents()
        
    >>> x = ASNFile(file='ib3702060_asn.fits')
    >>> x.showContents()
    1   ib3703uxq    EXP-DTH      yes
    2   ib3703vaq    EXP-DTH      yes
    3   ib3703vcq    EXP-DTH      yes
    4   ib3703vsq    EXP-DTH      yes
    5   IB3703050    PROD-DTH     yes
        """
        if self.exposures.__len__() > 0:
            for i,exp in enumerate(self.exposures):
                print '%5d   %s    EXP-DTH      yes' %(i+1,exp)
            print '%5d   %s    PROD-DTH     yes' %(i+2,self.product)
    
    


def asn_file_info(asn_file, verbose=1):
    """
asn_file_info(asn_file, verbose=1)
    
Get header information from files defined in an ASN table.
    
    >>> asn_file_info('ib3702060_asn.fits')
    # ib3702060_asn.fits
    # flt_file  filter  exptime  date_obs  time_obs pos_targ1 pos_targ2
    ib3702u4q_flt.fits.gz   G141  1302.9 2010-04-15 19:03:01    0.00    0.00
    ib3702u8q_flt.fits.gz   G141  1302.9 2010-04-15 19:32:22    0.61    0.18
    ib3702ukq_flt.fits.gz   G141  1302.9 2010-04-15 20:37:55    0.27    0.67
    ib3702uoq_flt.fits.gz   G141  1402.9 2010-04-15 21:06:31   -0.34    0.48
    
    """
    #asn_file = 'ib6o23020_asn.fits'
    asn = ASNFile(asn_file)
    lines = ['# %s' %asn_file]
    lines.append('# flt_file  filter  exptime  date_obs  time_obs pos_targ1 pos_targ2')
    ##### Loop through flt files in ASN list
    for exp in asn.exposures:
        flt_file = find_fits_gz(exp.lower()+'_flt.fits')
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
    asn = ASNFile(asn_file)
    NEXP = len(asn.exposures)
    RAcenters  = np.zeros(NEXP)
    DECcenters = np.zeros(NEXP)
    ##### Loop through exposures and get footprints
    for i, exp_root in enumerate(asn.exposures):
        flt_file = find_fits_gz(exp_root.lower()+'_flt.fits', hard_break = True)
        
        regX, regY = wcs_polygon(flt_file,extension=1)
        line = "polygon(%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f)" \
            %(regX[0],regY[0],regX[1],regY[1],regX[2],regY[2],regX[3],regY[3])
                    
        RAcenters[i] = np.mean(regX)
        DECcenters[i] = np.mean(regY)
        fp.write(line+' # color=magenta\n')
        
    ##### Text label with ASN filename
    fp.write('# text(%10.6f,%10.6f) text={%s} color=magenta\n' \
        %(np.mean(RAcenters),np.mean(DECcenters),asn_file.split('_asn.fits')[0]))
    fp.close()
    print "3D-HST / ASN_REGION: %s" %(output_file)
    


def wcs_polygon(fits_file, extension=1):
    """    
X, Y = wcs_polygon(fits_file, extension=1)
    
Calculate a DS9/region polygon from WCS header keywords.  
    
Will try to use pywcs.WCS.calcFootprint if pywcs is installed.  Otherwise
will compute from header directly.
    
    """
    ##### Open the FITS file
    hdulist = pyfits.open(fits_file) 
    ##### Get the header
    try:
        sci = hdulist[extension].header
    except IndexError:
        print 'ERROR 3D-HST/wcs_polygon:\n'+\
              'Extension #%d out of range in %s' %(extension, fits_file)
        raise
    
    #### Try to use pywcs if it is installed
    pywcs_exists = True
    try:
        import pywcs
    except:
        pywcs_exists = False   
    
    if pywcs_exists:
        wcs = pywcs.WCS(sci)
        footprint = wcs.calcFootprint()
        regX = footprint[:,0]    
        regY = footprint[:,1]    
        return regX, regY
    
    #### Do it by hand if no pywcs    
    NAXIS = [sci['NAXIS1'],sci['NAXIS2']]
    CRPIX = [sci['CRPIX1'],sci['CRPIX2']]
    CRVAL = [sci['CRVAL1'],sci['CRVAL2']]
    cosDec = np.cos(CRVAL[1]/180*np.pi)
    ##### Make region polygon from WCS keywords
    regX = CRVAL[0] + ( (np.array([0,NAXIS[0],NAXIS[0],0])-CRPIX[0])*sci['CD1_1'] + \
                        (np.array([0,0,NAXIS[1],NAXIS[1]])-CRPIX[1])*sci['CD1_2'] ) / cosDec
    regY = CRVAL[1] + ( (np.array([0,NAXIS[0],NAXIS[0],0])-CRPIX[0])*sci['CD2_1'] + \
                        (np.array([0,0,NAXIS[1],NAXIS[1]])-CRPIX[1])*sci['CD2_2'] )
    return regX, regY

    
def region_mask(shape,px,py):
    """
mask = region_mask(image.shape,px,py)
    
Make a mask image where pixels within the polygon defined by px_i, py_i
are set to 1.  This is the same algorithm as in :ref:`point_in_polygon`
but with array orders switched around to be much more efficient.
    
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
        #### Make arctan aware of the signs of x,y in arctan(x/y)
        #theta += np.arctan(cp/dp) - np.pi*((cp < 0) & (dp < 0)) + np.pi*((cp > 0) & (dp < 0))
        theta += np.arctan2(cp,dp) 
    ##### Set up mask
    dq = np.zeros((NX,NY),dtype=np.int)
    flag_idx = np.where(np.abs(theta) > np.pi)
    dq[flag_idx] = 1
    return dq


def point_in_polygon(x,y,px,py):
    """
test = point_in_polygon(x,y,px,py)
    
Test if coordinates (x,y) are inside polygon defined by (px, py)
    
<http://www.dfanning.com/tips/point_in_polygon.html>, translated to Python
    
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
    # theta = np.arctan(cp/dp) - np.pi*((cp < 0) & (dp < 0)) + np.pi*((cp > 0) & (dp < 0))
    theta = np.arctan2(cp,dp)
    if np.abs(np.sum(theta)) > np.pi:
        return True
    else:
        return False


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

