"""
3DHST.utils

Utilities for reading/writing ASN and DS9/region files.

"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import os
import string
import time

import pyfits
import numpy as np

def columnFormat(colname):
    """
Set format for common column names.  
    
    id: 'A'
    ra,dec: 'D'
    [everything else]: 'E'
    """
    fmt='E'
    if colname.lower().find('id') >= 0: fmt='A'
    if colname.lower().find('ra') >= 0: fmt='D'
    if colname.lower().find('dec') >= 0: fmt='D'
    return fmt
    
def ASCIItoFITS(infile, comment='#'):
    """
ASCIItoFITS(infile, [comment='#'])

    Read an ASCII file, `infile`, and get column names from first line, 
    which begins with the `comment` character.
    
    Output will be in `infile`+'.FITS'
    """
    ### Get first header line and count commented lines to skip
    file = open(infile,'r')
    line0 = file.readline()
    line=line0
    hskip=0
    while line.startswith(comment):
        line = file.readline()
        hskip +=1
    #file.close()
    
    #### clean up special characters from header
    line0 = line0.replace('.','p')
    line0 = line0.replace('-','_')
    line0 = line0.replace('(','')
    line0 = line0.replace(')','')
    line0 = line0.replace('-','_')
    line0 = line0.replace('[','')
    line0 = line0.replace(']','')
    header=string.split(line0[1:])
    NCOLS = len(header)
    for i in range(NCOLS):
        #### lowercase column names
        header[i] = header[i].lower()
        #### Add an 'n' to column names that start with a digit
        #### so you can address them like column = data.x
        if header[i][0].isdigit():
            header[i] = 'n'+header[i]
            
    if NCOLS == 0:
        print ('No header line found.  I\'m looking for a first line that'+
               'begins with %s followed by column names.' %comment)
        return None
        
    #### Check next NCHECK data lines to look for alphanumeric characters
    formats = ['f4']*NCOLS
    formats_pyfits = ['D']*NCOLS
    NCHECK = 10
    for i in range(NCHECK):
        line = file.readline().split()
        for j, value in enumerate(line):
            nchar = len(value)
            for k in range(nchar):
                if value[k].isalpha():
                    formats[j] = '1S'
                    formats_pyfits[j] = 'A'
    
    file.close()
            
    #### Read data file
    data=np.loadtxt(infile, comments=comment,
                    dtype = {'names': tuple(header), 'formats': tuple(formats)})
                    
    #### Make output FITS table
    # make_struct='str = {'
    go_ColDefs='cols=pyfits.ColDefs(['
    for i in range(NCOLS):
        col_string = 'col%d = pyfits.Column(name=\'%s\',' %(i,header[i]) + \
    ' format=\'%s\', array=data[data.dtype.names[%d]])' %(formats_pyfits[i],i)
        exec(col_string)
        go_ColDefs += 'col%d,' %(i)
        # make_struct += '\'%s\':data[0:,%d],' %(header[i],i)
    
    exec(go_ColDefs[:-1]+'])') # cols=pyfits.ColDefs([col1, col2, ...])
    
    #### Initialize table
    tbhdu = pyfits.new_table(cols)
    
    #### Primary HDU
    hdu = pyfits.PrimaryHDU()
    
    thdulist = pyfits.HDUList([hdu,tbhdu])
    
    #### Add modification time of "infile" to FITS header
    infile_mod_time = time.strftime("%m/%d/%Y %I:%M:%S %p", \
                         time.localtime(os.path.getmtime(infile)))
    thdulist[1].header.update('MODTIME',infile_mod_time)
    
    thdulist.writeto(infile+'.FITS', clobber=True)
    
    #return tbhdu.data, tbhdu.columns
    return tbhdu.data
    
def ReadASCIICat(infile, comment='#', force=False, verbose=False):
    """
data = ReadASCIICat(infile, comment='#', force=False, verbose=False)

    Read FITS table created from an ASCII catalog file.
    
    If ASCIItoFITS output doesn't exist or "force=True", create it.
    """
    
    if os.path.exists(infile) is False:
        print ('File, %s, not found.' %(infile))
        return None
    
    fileExists = False
    if (os.path.exists(infile+'.FITS')):
        fileExists = True
        theFITSFile = infile+'.FITS'
    if (os.path.exists(infile+'.FITS.gz')):
        fileExists = True
        theFITSFile = infile+'.FITS.gz'
    
    if (fileExists is False) or (force):
        if verbose:
            print ('Running ASCIItoFITS: %s' %(infile))
        data = ASCIItoFITS(infile,comment=comment)
        #return data
        #return ASCIItoFITS(infile,comment=comment)
    else:
        if verbose:
            print ('Reading : %s' %(theFITSFile))
        hdulist = pyfits.open(theFITSFile)
        #### Check mod time of 'infile'.  
        #### If changed, then re-read with ASCIItoFITS
        infile_mod_time = time.strftime("%m/%d/%Y %I:%M:%S %p", \
                                time.localtime(os.path.getmtime(infile)))
        if infile_mod_time == hdulist[1].header['MODTIME']:
            data = hdulist[1].data
            #return hdulist[1].data
        else:
            if verbose:
                print('%s has changed.  Re-generating FITS file...' %(infile))
            data = ASCIItoFITS(infile,comment=comment)
            #return data

    data.file = infile
    return data
    
class listArray(list):
    """
listArray(list)
    
    Add + - / * ** operators to Python list objects via Numpy.
    
    >>> x = listArray([1,2,3])
    >>> x * 2
    [2, 4, 6]
    
    """
    import numpy as np

    def __add__(self, addval):
        return list(np.array(self)+addval)
    
    def __sub__(self, addval):
        return list(np.array(self)-addval)

    def __mul__(self, addval):
        return list(np.array(self)*addval)
    
    def __div__(self, addval):
        return list(np.array(self)/addval)
    
    def __pow__(self, addval):
        return list(np.array(self)**addval)
    
    # def __getitem__(self, idx):
    #     return list(np.array(self)[idx])
    
def get_package_data(dataname):
    """
    (taken from astropysics.io)
    Use this function to load data files distributed with astropysics in the 
    astropysics/data directory
    
    `dataname` is the file name of a file in the *threedhst/data* directory, and
    a string with the contents of the file will be returned
    """
    try:
        ### Find the data directory in the root
        ### directory of the threedhst package
        from . import __name__ as rootname
        from . import __file__ as rootfile
        from pkgutil import get_loader
        from os.path import dirname
        path = dirname(rootfile)+'/data/'+dataname
        return get_loader(rootname).get_data(path)
    except:
        ### Hardwired  in case relative import doesn't work
        fp = open('/research/HST/GRISM/3DHST/progs/threedhst/data/'+dataname)
        return fp.read()

def find_fits_gz(fits_file, hard_break = True):
    """
the_file = find_fits_gz(fits_file, hard_break = True)
    
    With ``fits_file`` being some filename with an extension ``.fits``
    (ib3713wvq_flt.fits), check to see if the file itself or its gzipped
    counterpart (fits_file+'.gz') exists.
    
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
        raise IOError('File %s[.gz] not found in %s' 
                              %(fits_file, os.getcwd()))
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
            warn('ASN file %s has no EXP-DTH items')
        else:
            self.exposures = []
            for exp in names[exp_idx]:
                self.exposures.append(exp.lower())
        
        ##### Products
        prod_idx = np.where(types == 'PROD-DTH')
        if prod_idx[0].shape[0] != 1:
            warn('ASN file %s has N != 1 PROD-DTH items' %self.file )
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

    out_file='self' writes to `self.file`.
    
        """
        if not out_file:
            print "USAGE:: writeToFile(out_file='output_asn.fits')"
        else:
            if out_file == 'self':
                out_file = self.file
            
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
    
    def append(self, new):
        """
append(self, new)
        
    `new` must be an instance of ASNFile.
        
    `new.exposures` are added to the `self.exposures` list.
        """
        from warnings import warn
        if not isinstance(new,self.__class__):
            warn("argument is not an instance of ASNFile()")
        else:
            self.exposures.extend(new.exposures)
    


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
    lines.append(
       '# flt_file  filter  exptime  date_obs  time_obs pos_targ1 pos_targ2')
    ##### Loop through flt files in ASN list
    for exp in asn.exposures:
        flt_file = find_fits_gz(exp.lower()+'_flt.fits')
        fp_flt = pyfits.open(flt_file)
        ##### Get general information from extension 0
        fp_header = fp_flt[0].header
        line = '%s %6s %7.1f %s %s %7.2f %7.2f' %(
                      flt_file,fp_header['FILTER'],fp_header['EXPTIME'],    
                      fp_header['DATE-OBS'],fp_header['TIME-OBS'],
                      fp_header['POSTARG1'],fp_header['POSTARG2'])
        lines.append(line)
    
    #### Print to stdout
    if verbose > 0:
        for line in lines:
            print line
    

# def asn_region(asn_file):
#     """
# asn_region(asn_file)
#     
# Create a DS9 region file for the exposures defined in an ASN file.
#     
#     """
#     ##### Output file
#     output_file = asn_file.split('.fits')[0]+'.pointing.reg'
#     fp = open(output_file,'w')
#     fp.write('fk5\n') ### WCS coordinates
#     ##### Read ASN file
#     asn = ASNFile(asn_file)
#     NEXP = len(asn.exposures)
#     RAcenters  = np.zeros(NEXP)
#     DECcenters = np.zeros(NEXP)
#     ##### Loop through exposures and get footprints
#     for i, exp_root in enumerate(asn.exposures):
#         flt_file = find_fits_gz(exp_root.lower()+'_flt.fits', hard_break = True)
#         
#         regX, regY = wcs_polygon(flt_file,extension=1)
#         line = "polygon(%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f)" \
#             %(regX[0],regY[0],regX[1],regY[1],regX[2],regY[2],regX[3],regY[3])
#                     
#         RAcenters[i] = np.mean(regX)
#         DECcenters[i] = np.mean(regY)
#         fp.write(line+' # color=magenta\n')
#         
#     ##### Text label with ASN filename
#     fp.write('# text(%10.6f,%10.6f) text={%s} color=magenta\n' \
#         %(np.mean(RAcenters),np.mean(DECcenters),
#           asn_file.split('_asn.fits')[0]))
#     fp.close()
#     print '3D-HST / ASN_REGION: %s\n' %(output_file)
#     
# 
# 
# def wcs_polygon(fits_file, extension=1):
#     """    
# X, Y = wcs_polygon(fits_file, extension=1)
#     
# Calculate a DS9/region polygon from WCS header keywords.  
#     
# Will try to use pywcs.WCS.calcFootprint if pywcs is installed.  Otherwise
# will compute from header directly.
#     
#     """
#     ##### Open the FITS file
#     hdulist = pyfits.open(fits_file) 
#     ##### Get the header
#     try:
#         sci = hdulist[extension].header
#     except IndexError:
#         print 'ERROR 3D-HST/wcs_polygon:\n'+\
#               'Extension #%d out of range in %s' %(extension, fits_file)
#         raise
#     
#     #### Try to use pywcs if it is installed
#     pywcs_exists = True
#     try:
#         import pywcs
#     except:
#         pywcs_exists = False   
#     
#     if pywcs_exists:
#         wcs = pywcs.WCS(sci)
#         footprint = wcs.calcFootprint()
#         regX = footprint[:,0]    
#         regY = footprint[:,1]    
#         return regX, regY
#     
#     #### Do it by hand if no pywcs    
#     NAXIS = [sci['NAXIS1'],sci['NAXIS2']]
#     CRPIX = [sci['CRPIX1'],sci['CRPIX2']]
#     CRVAL = [sci['CRVAL1'],sci['CRVAL2']]
#     cosDec = np.cos(CRVAL[1]/180*np.pi)
#     ##### Make region polygon from WCS keywords
#     regX = CRVAL[0] + \
#             ((np.array([0,NAXIS[0],NAXIS[0],0])-CRPIX[0])*sci['CD1_1'] +                        
#              (np.array([0,0,NAXIS[1],NAXIS[1]])-CRPIX[1])*sci['CD1_2']) / cosDec
#     
#     regY = CRVAL[1] + \
#             ((np.array([0,NAXIS[0],NAXIS[0],0])-CRPIX[0])*sci['CD2_1'] +         
#              (np.array([0,0,NAXIS[1],NAXIS[1]])-CRPIX[1])*sci['CD2_2'])
#              
#     return regX, regY
#     
# def region_mask(shape,px,py):
#     """
# mask = region_mask(image.shape,px,py)
#     
# Make a mask image where pixels within the polygon defined by px_i, py_i
# are set to 1.  This is the same algorithm as in :ref:`point_in_polygon`
# but with array orders switched around to be much more efficient.
#     
# Note: something like this could be used to flag grism 0th order contaminants
#     """
#     NX=shape[0]
#     NY=shape[1]
#     y,x = np.mgrid[1:NX+1,1:NY+1]
#     ##### Close polygons
#     NPOLY = px.shape[0]
#     tmp_px = np.append(px,px[0])
#     tmp_py = np.append(py,py[0])
#     theta = np.zeros((NX,NY),dtype=np.float)
#     for i in np.arange(NPOLY):
#         ##### Dot, cross products
#         X1 = tmp_px[i] - x
#         Y1 = tmp_py[i] - y 
#         X2 = tmp_px[i+1] - x
#         Y2 = tmp_py[i+1] - y
#         dp = X1*X2 + Y1*Y2
#         cp = X1*Y2 - Y1*X2
#         theta += np.arctan2(cp,dp) 
#     ##### Set up mask
#     dq = np.zeros((NX,NY),dtype=np.int)
#     flag_idx = np.where(np.abs(theta) > np.pi)
#     dq[flag_idx] = 1
#     return dq
# 
# 
# def point_in_polygon(x,y,px,py):
#     """
# test = point_in_polygon(x,y,px,py)
#     
# Test if coordinates (x,y) are inside polygon defined by (px, py)
#     
# <http://www.dfanning.com/tips/point_in_polygon.html>, translated to Python
#     
#     """    
#     N = px.shape[0]
#     ##### Close polygons
#     tmp_px = np.append(px,px[0])
#     tmp_py = np.append(py,py[0])
#     ##### Counters
#     i = np.arange(N)
#     ip = np.arange(N)+1
#     ##### Dot, cross products
#     X1 = tmp_px[i] - x
#     Y1 = tmp_py[i] - y 
#     X2 = tmp_px[ip] - x
#     Y2 = tmp_py[ip] - y
#     dp = X1*X2 + Y1*Y2
#     cp = X1*Y2 - Y1*X2
#     theta = np.arctan2(cp,dp)
#     if np.abs(np.sum(theta)) > np.pi:
#         return True
#     else:
#         return False
# 

