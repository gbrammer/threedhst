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
        output = fp.read()
        fp.close()
        return output

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
        
    >>> asn = ASNFile(file='ib3701050_asn.fits', grow=False)
    >>> asn.exposures
    ['ib3701ryq', 'ib3701sbq', 'ib3701sdq', 'ib3701sqq']
    >>> asn.product
    'IB3701050'
    
    If grow=True, allow file rootnames to be 20 characters rather than 14.
     """
    def _read_asn_file(self, grow=True):
        """
_read_asn_file(self)
        
    Read an ASN FITS file (self.file).
        """
        import numpy as np
        from warnings import warn
        
        self.in_fits = pyfits.open(self.file)
        
        data = self.in_fits[1].data
        
        if grow:
            #### Allow more characters in the MEMNAME column
            memname = pyfits.Column(name='MEMNAME', format='30A', array=self.in_fits[1].columns[0].array.astype('S30'), disp='A30')
            memtype = self.in_fits[1].columns[1]
            memprsnt = self.in_fits[1].columns[2]
            coldefs = pyfits.ColDefs([memname, memtype, memprsnt])
            hdu = pyfits.new_table(coldefs)
            hdu.header = self.in_fits[1].header
            hdu.header['TFORM1'] = '30A'
            hdu.header['TDISP1'] = 'A30'
            hdu.header['NAXIS1'] += 16
            self.in_fits[1] = hdu        
    
        data = self.in_fits[1].data
        #print data
        
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
    
    
    def __init__(self, file=None, grow=True):
        self.file = file
        self.exposures = []
        self.product = None
        if file:
            self._read_asn_file(grow=grow)
    
    
    def write(self, out_file=None, clobber=True):
        """
write(self,out_file=None, clobber=True)

    out_file='self' writes to `self.file`.
    
        """
        if not out_file:
            print "USAGE:: asn.write(out_file='output_asn.fits')"
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
            tbhdu = pyfits.new_table(self.in_fits[1].columns, nrows=nrows, fill=True)
            for i in range(nexp):
                tbhdu.data[i] = (self.exposures[i].upper(), 'EXP-DTH', True)
            if nprod > 0:
                tbhdu.data[i+1] = (self.product, 'PROD-DTH', True)
            
            tbhdu.header = self.in_fits[1].header.copy()
            tbhdu.header.update('ASN_ID',out_file.split('_asn.fits')[0])
            tbhdu.header.update('ASN_TAB',out_file)
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
    

def combine_asn_shifts(asn_list, out_root='combined', path_to_FLT='./', 
                       run_multidrizzle=False):
    """
combine_asn_shifts(asn_list, out_root='combined', path_to_FLT='./', 
                   run_multidrizzle=False)
                   
    Combine a list of ASN tables and their associated shiftfiles into a single 
    output, suitable for making a mosaic across visits with different
    orientations.  The shifts determined by, e.g. `tweakshifts` are all relative
    to the output WCS defined for the exposures in a given visit.  This 
    reference frame will be different for different values of "PA_V3", so to 
    combine the shifts into a single file, the different orientation angles
    need to be taken into account.
    
    This routine combines the ASN tables and the shiftfiles, with the 
    "reference" WCS taken from the first exposure of the first input ASN table.
    Assuming that the shiftfiles all produce images registered to some *fixed*
    reference WCS, this local reference doesn't really matter.
    
    The script looks for FLT images (gzipped or not) in the relative path
    defined by ``path_to_FLT``.
    
    EXAMPLE: 
    
    $ ls *shifts.txt *asn.fits
    ib3725050_asn.fits
    ib3725050_shifts.txt
    ib3726050_asn.fits
    ib3726050_shifts.txt
    ib3727050_asn.fits
    ib3727050_shifts.txt
    ib3728050_asn.fits
    ib3728050_shifts.txt
        
    $ python
    >>> import glob
    >>> asn_files = glob.glob()
    >>> combine_asn_shifts(asn_files, out_root='combined',
                           run_multidrizzle=False)
    """
    import numpy as np
    import pyfits
    import threedhst.utils
    import threedhst.shifts
    
    #### "Reference WCS is the first set of exposures"
    asn_ref = threedhst.utils.ASNFile(asn_list[0])
    shift_ref = threedhst.shifts.ShiftFile(asn_list[0].split('_asn.fits')[0]+'_shifts.txt')
    
    #### Set WCS reference to first image rather than 'tweak.fits'
    shift_ref.headerlines[1] = '# refimage: %s_flt.fits[1]\n' %(asn_ref.exposures[0])
    
    #### Get PA_V3 angle of reference
    fits = threedhst.utils.find_fits_gz('%s/%s_flt.fits' %(path_to_FLT, asn_ref.exposures[0]))
    angle_ref = pyfits.getheader(fits).get('PA_V3')
    
    #### Loop through other ASN files
    for asn_file in asn_list[1:]:
        print asn_file
        
        #### Read the ASN file
        asn = threedhst.utils.ASNFile(asn_file)
        asn_ref.exposures.extend(asn.exposures)
        
        #### Get the PA_V3 angle of the first exposures in this ASN file
        fits = threedhst.utils.find_fits_gz('%s/%s_flt.fits' %(path_to_FLT, asn.exposures[0]))
        angle = pyfits.getheader(fits).get('PA_V3')
        #### Difference angle between current and reference
        alpha = (angle-angle_ref)/360.*2*np.pi
                
        #### Read the shifts and rotate them into reference frame
        shift = threedhst.shifts.ShiftFile(asn_file.split('_asn.fits')[0]+'_shifts.txt')
        for i,im in enumerate(shift.images):
            xshift = shift.xshift[i]
            yshift = shift.yshift[i]
            xsh = (xshift*np.cos(alpha)-yshift*np.sin(alpha))
            ysh = (xshift*np.sin(alpha)+yshift*np.cos(alpha))
            shift_ref.append(im, xshift=xsh, yshift=ysh,
                             rotate=shift.rotate[i], scale=shift.scale[i])
    
    #### Write ASN and shiftfiles
    shift_ref.write('%s_shifts.txt' %(out_root))
    asn_ref.product = out_root
    asn_ref.write('%s_asn.fits' %(out_root))
    
    #### Run Multidrizzle with combined ASN and shifts
    if run_multidrizzle:
        run_asn = out_root+'_asn.fits'
        #threedhst.process_grism.fresh_flt_files(run_asn,
        #    from_path='../../RAW')s
        threedhst.prep_flt_files.startMultidrizzle(root=run_asn, 
            use_shiftfile=True, skysub=False, final_scale=0.06, updatewcs=False, 
            pixfrac=0.8, driz_cr=False, median=False)

def asn_file_info(asn_file, verbose=1, path_to_FLT = './'):
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
        flt_file = find_fits_gz(path_to_FLT+'/'+exp.lower()+'_flt.fits')
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
    
def replace_OrIg():
    """
    Copy the Multidrizzle backup files back to the original 
    flt files.
    """
    import glob
    import shutil
    files=glob.glob('*OrIg_flt.fits')
    for file in files:
        out=file.replace('OrIg_','')
        print file, out
        shutil.move(file, out)

def biweight(xarr, both=False, mean=False):
    """
    Compute the biweight estimator for an input array.
    
    Example:
    
    >>> x = np.random.randn(1000)
    >>> x[0:5] = 5000
    >>> mu, sigma = biweight(x, both=True)
    >>> sig = biweight(x, both=False)            # get just sigma
    >>> print mu, sigma, np.median(x), np.std(x)
    
    """
    bigm = np.median(xarr)
    mad = np.median(np.abs(xarr-bigm))
    
    c = 6.0
    u = (xarr-bigm)/c/mad
    u1 = np.abs(u) < 1
    
    #### biweight mean
    if (u[u1].size > 0):
        cbi = bigm + np.sum((xarr[u1]-bigm)*(1-u[u1]**2)**2)/ \
                  np.sum((1-u[u1]**2)**2)
    else:
        cbi = -99
        
    #### biweight sigma
    c=9.0
    u = (xarr-bigm)/c/mad
    u1 = np.abs(u) < 1

    if (u[u1].size > 0):
        sbi = (np.sqrt(xarr.size)*np.sqrt(np.sum((xarr[u1]-bigm)**2*(1-u[u1]**2)**4))/ \
                                  np.abs(np.sum((1-u[u1]**2)*(1-5*u[u1]**2))))
    else:
        sbi = -99
    
    if mean:
        return cbi
        
    if both:
        return cbi, sbi
    else:
        return sbi

def biweight2(xarr, both=False, mean=False):
    """
    Compute the biweight estimator for an input array.
    
    Example:
    
    >>> x = np.random.randn(1000)
    >>> x[0:5] = 5000
    >>> mu, sigma = biweight(x, both=True)
    >>> sig = biweight(x, both=False)            # get just sigma
    >>> print mu, sigma, np.median(x), np.std(x)
    
    """
    bigm = np.median(xarr, axis=0)
    sh = xarr.shape
    if len(sh) == 2:
        bigm2 = np.dot(np.ones((sh[0],1)), bigm.reshape((1,sh[1])))
    else:
        bigm2 = np.ones((sh[0],))*bigm
            
    mad = np.median(np.abs(xarr-bigm2), axis=0)
    if len(sh) == 2:
        mad2 = np.dot(np.ones((sh[0],1)), mad.reshape((1,sh[1])))
    else:
        mad2 = np.ones((sh[0],))*mad
        
    c = 6.0
    u = (xarr-bigm2)/c/mad2
    u1 = np.abs(u) < 1
    
    #### Force flagged items in sum to be zero
    u1x = (np.abs(u) > 1) & (not np.isfinite(xarr))
    xarr[u1x] = bigm2[u1x]
    u[u1x] = 1
    
    #### biweight mean
    if (u[u1].size > 0):
        cbi = bigm + np.sum((xarr-bigm2)*(1-u**2)**2, axis=0)/ \
                  np.sum((1-u**2)**2, axis=0)
    else:
        cbi = -99
        
    #### biweight sigma
    c=9.0
    u = (xarr-bigm2)/c/mad2
    u1 = np.abs(u) < 1
    
    ####
    u1x = (np.abs(u) > 1) & (not np.isfinite(xarr))
    xarr[u1x] = bigm2[u1x]
    u[u1x] = 1
    
    if (u[u1].size > 0):
        sbi = (np.sqrt(xarr.shape[0])*np.sqrt(np.sum((xarr-bigm2)**2*(1-u**2)**4, axis=0))/ \
                                  np.abs(np.sum((1-u**2)*(1-5*u**2), axis=0)))
    else:
        sbi = -99
    
    if mean:
        return cbi
        
    if both:
        return cbi, sbi
    else:
        return sbi

def nmad(xarr):
    return 1.48*np.median(np.abs(xarr-np.median(xarr)))

def runmed(xi, yi, NBIN=10):
        
    NPER = xi.size/NBIN
    xm = np.arange(NBIN)*1.
    xs = xm*0
    ym = xm*0
    ys = xm*0
    N = np.arange(NBIN)
    
    so = np.argsort(xi)
    idx = np.arange(NPER)
    for i in range(NBIN):
        ym[i], ys[i] = biweight(yi[so][idx+NPER*i], both=True)
        xm[i], xs[i] = biweight(xi[so][idx+NPER*i], both=True)
        N[i] = xi[so][idx+NPER*i].size
    
    return xm, ym, ys, N

def medfilt(xarr, N=3):
    """
    Median filter
    """
    out = xarr*0.
    half = int(N/2)
    for i in range(0,half):
        out[i] = np.median(xarr[i:i+half+1])
        
    for i in range(half, len(xarr)-half,1):
        out[i] = np.median(xarr[i-half:i+half+1])
    
    for i in range(len(xarr)-half,len(xarr)):
        out[i] = np.median(xarr[i-half:i])
    
    return out
    
def xyrot(xin, yin, theta, x0=0., y0=0., radians=False, ccw=False):
    """
    Rotate (xin, yin) coordinates by an angle `theta`
    """
    
    if radians:
        rad = theta           
    else:
        rad=theta*2*np.pi/360.
    
    if ccw:
        rad = 2*np.pi-rad
    
    mat = np.zeros((2,2))
    mat[0,:] = np.array([np.cos(rad),-np.sin(rad)])
    mat[1,:] = np.array([np.sin(rad),np.cos(rad)])
        
    coo = np.zeros((2,len(xin)))
    coo[0,:] = xin-x0
    coo[1,:] = yin-y0

    out = np.dot(mat.transpose(), coo)
    xout = out[0,:]
    yout = out[1,:]
    return xout, yout
