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

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import numpy as np

def skipme(base='_filename', ext='skipme', verbose=True):
    """
    Helper function to make a dummy file with some basename as a placeholder to indicate
    that a script has been started related to the `base` string.  
    
    Returns True if the file `base`.`ext` exists, otherwise "touches" that file and returns False
    """
    skipfile = '%s.%s' %(base, ext)
    if os.path.exists(skipfile):
        print('File %s found' %(skipfile))
        return True
    else:
        os.system('touch %s' %(skipfile))
        return False
        
def test_conserve():
    xfull = np.arange(0,1000001,1)
    yfull = np.sin(xfull/np.pi/2/20)+1
    
    xint = np.arange(0,1000001,10000)
    yint_0 = np.interp(xint, xfull, yfull)
    yint_1 = threedhst.utils.interp_conserve(xint, xfull, yfull)
    yint_2 = threedhst.utils.interp_conserve_c(xint, xfull, yfull)
    np.trapz(yint_0, xint)/np.trapz(yfull,xfull)-1, np.trapz(yint_1, xint)/np.trapz(yfull,xfull)-1, np.trapz(yint_2, xint)/np.trapz(yfull,xfull)-1
    
def interp_conserve_c(x, xp, fp, left=0, right=0):
    """
    Interpolate `xp`,`yp` array to the output x array, conserving flux.  
    `xp` can be irregularly spaced.
    """
    from scipy import weave
    from scipy.weave import converters
    
    templmid = (x[1:]-x[:-1])/2.+x[:-1]
    templmid = np.append(templmid, np.array([x[0], x[-1]]))
    templmid = templmid[np.argsort(templmid)]
    tempfmid = np.interp(templmid, xp, fp, left=left, right=right)
        
    #### Code from eazy
    NTEMPL = len(templmid)
    tf = fp
    tlam = xp
    ntlam = len(xp)
    
    outy = templmid[:-1]
    
    code = """
    long i,k,istart;
    double h, numsum;
    
    ////// Rebin template grid to master wavelength grid, conserving template flux
    i=0;
    for (k=0;k<NTEMPL;++k) {
        numsum=0.;
        
        //// Go to where tlam is greater than the first midpoint
        while ((tlam(i) < templmid(k)) && (i < ntlam)) ++i;
        istart=i;
        
        /////// First point
        if (tlam(i) < templmid(k+1)) {
            h = tlam(i)-templmid(k);
            numsum+=h*(tf(i)+tempfmid(k));
            ++i;
        }
        if (i==0) ++i;
                
        /////// Template points between master grid points
        while ((tlam(i) < templmid(k+1)) && (i < ntlam)) {
            h = tlam(i)-tlam(i-1);
            numsum+=h*(tf(i)+tf(i-1));
            ++i;
        }
        
        //// If no template points between master grid points, then just use interpolated midpoints
        if ( i == istart ) {
            h = templmid(k+1)-templmid(k);
            numsum=h*(tempfmid(k+1)+tempfmid(k));
        } else {  
            ///// Last point              
            --i;
            h = templmid(k+1)-tlam(i);
            numsum+=h*(tempfmid(k+1)+tf(i));
        }
        
        outy(k) = numsum*0.5/(templmid(k+1)-templmid(k));
    }
    return_val = 1;
    
    """
    
    result = weave.inline(code,['templmid','tempfmid','NTEMPL','tlam','tf','ntlam','outy'], type_converters=converters.blitz, compiler = 'gcc', verbose=2)
    return outy
    
def interp_conserve(x, xp, fp, left=0., right=0.):
    """
    Interpolate `xp`,`yp` array to the output x array, conserving flux.  
    `xp` can be irregularly spaced.
    """
    midpoint = (x[1:]-x[:-1])/2.+x[:-1]
    midpoint = np.append(midpoint, np.array([x[0],x[-1]]))
    midpoint = midpoint[np.argsort(midpoint)]
    int_midpoint = np.interp(midpoint, xp, fp, left=left, right=right)
    int_midpoint[midpoint > xp.max()] = 0.
    int_midpoint[midpoint < xp.min()] = 0.
    
    fullx = np.append(xp, midpoint)
    fully = np.append(fp, int_midpoint)
    
    so = np.argsort(fullx)
    fullx, fully = fullx[so], fully[so]
    
    outy = x*0.
    dx = midpoint[1:]-midpoint[:-1]
    for i in range(len(x)):
        bin = (fullx >= midpoint[i]) & (fullx <= midpoint[i+1])
        outy[i] = np.trapz(fully[bin], fullx[bin])/dx[i]
    
    return outy
    
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
        print(('No header line found.  I\'m looking for a first line that'+
               'begins with %s followed by column names.' %comment))
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
        print(('File, %s, not found.' %(infile)))
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
            print(('Running ASCIItoFITS: %s' %(infile)))
        data = ASCIItoFITS(infile,comment=comment)
        #return data
        #return ASCIItoFITS(infile,comment=comment)
    else:
        if verbose:
            print(('Reading : %s' %(theFITSFile)))
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
                print(('%s has changed.  Re-generating FITS file...' %(infile)))
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
    
    ### Prefer the non gzipped version
    if os.path.exists(fits_file):
        return fits_file
    
    ### unzipped not found.  gzipped?
    if os.path.exists(fits_file+'.gz'):
        return fits_file+'.gz'
        
    #### File not found.  Either raise an error or return None
    if hard_break:
        raise IOError('File %s[.gz] not found in %s' 
                              %(fits_file, os.getcwd()))
    else:
        return None

def ASN_footprint(asn_file, color='green', path_to_flt='./', verbose=False):
    """
    Make a region file with footprints for each exposure in an ASN file
    """
    
    try:
        import stwcs
        has_stwcs = True
    except:
        import astropy.wcs as pywcs
        has_stwcs = False
        
    import glob
    
    asn = ASNFile(asn_file)
    fp = open(asn_file.replace('asn.fits', 'asn.reg'),'w')
    fp.write('fk5\n')

    for exp in asn.exposures:
        file=glob.glob('%s/%s_fl?.fi*[tg][sz]' %(path_to_flt, exp))[0]
        flt = pyfits.open(file)
        
        if flt[0].header['DETECTOR'] in ['WFC', 'UVIS']:
            extensions = [1,4]
        else:
            extensions = [1]
        
        for ext in extensions:
            if has_stwcs:
                wcs = stwcs.wcsutil.HSTWCS(flt, ext=ext)
            else:
                wcs = pywcs.WCS(flt[ext].header)
            
            foot = wcs.calc_footprint()
            if verbose:
                print(exp, foot)
                 
            poly_str = ', '.join(['%.6f, %.6f' %(foot[i][0], foot[i][1]) for i in range(foot.shape[0])])
            fp.write('polygon(%s) # color=%s\n' %(poly_str, color))
    
    fp.close()
    
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
            memname = pyfits.Column(name='MEMNAME', format='40A', array=self.in_fits[1].columns[0].array.astype('S40'), disp='A40')
            memtype = self.in_fits[1].columns[1]
            memprsnt = self.in_fits[1].columns[2]
            coldefs = pyfits.ColDefs([memname, memtype, memprsnt])
            try:
                #print 'from_columns'
                hdu = pyfits.BinTableHDU.from_columns(coldefs)
            except:
                #print 'fail pyfits'
                hdu = pyfits.new_table(coldefs)
                
            hdu.header = self.in_fits[1].header
            hdu.header['TFORM1'] = '40A'
            hdu.header['TDISP1'] = 'A40'
            hdu.header['NAXIS1'] += 26
            self.in_fits[1] = hdu        
    
        data = self.in_fits[1].data
        #print data
        
        self.header = self.in_fits[0].header
        
        names = data.field('MEMNAME')
        types = data.field('MEMTYPE')
        
        ##### Exposures
        #exp_idx  = np.where(types == 'EXP-DTH')
        exp_idx = types == 'EXP-DTH'
        #### check if MEMTYPE starts with EXP, have other cases where type is "EXP-RP#"
        for ii, type in enumerate(types):
            if types[ii].startswith('EXP'):
                exp_idx[ii] = True
                
        if exp_idx.sum() == 0:
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
            print("USAGE:: asn.write(out_file='output_asn.fits')")
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
            if 'EXTEND' not in list(hdu.header.keys()):
                hdu.header.update('EXTEND', True, after='NAXIS')
                
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
                print('%5d   %s    EXP-DTH      yes' %(i+1,exp))
            print('%5d   %s    PROD-DTH     yes' %(i+2,self.product))
    
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
        print(asn_file)
        
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
            print(line)
    
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
        print(file, out)
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

def gehrels(Nin,twosig=False,threesig=False):
    """
    Poisson-like confidence intervals for counting data from Gehrels (1986)
    http://adsabs.harvard.edu/cgi-bin/bib_query?1986ApJ...303..336G
    """
    n = Nin
    
    S = 1.0
    bet = 0.0
    gam = 0.
    fixit = 0.
    
    if twosig:
        ## 0.9772
        S = 2.0
        bet = 0.062
        gam = -2.07
        fixit=1.
    
    if threesig:
        ## 0.9987
        S = 3.0
        bet = 0.222
        gam = -1.88
        fixit=1.     
    
    upper = n+S*np.sqrt(n+3./4)+(S**2+3)/4. 
    lower = n*(1-1./9./n-S/3/np.sqrt(n))**3    
    lower = n*(1-1./9./n-S/3/np.sqrt(n)+bet*n**gam*fixit)**3

    if np.isscalar(n):
        if n <= 0:
            upper=0
            lower=0
    else:
        bad = n <= 0
        upper[bad] = 0
        lower[bad] = 0
    
    return (lower, upper)
    
def nmad(xarr):
    """
    result = nmad(arr)

    Get the NMAD statistic of the input array, where
    NMAD = 1.48 * median(ABS(arr) - median(arr)).
    """
    return 1.48*np.median(np.abs(xarr-np.median(xarr)))

def runmed(xi, yi, NBIN=10, use_median=False, use_nmad=False, reverse=False):
    """
    Running median/biweight/nmad
    """
    
    NPER = xi.size/NBIN
    xm = np.arange(NBIN)*1.
    xs = xm*0
    ym = xm*0
    ys = xm*0
    N = np.arange(NBIN)
    
    so = np.argsort(xi)
    if reverse:
        so = so[::-1]
        
    idx = np.arange(NPER)
    for i in range(NBIN):
        ym[i], ys[i] = biweight(yi[so][idx+NPER*i], both=True)
        xm[i], xs[i] = biweight(xi[so][idx+NPER*i], both=True)
        N[i] = xi[so][idx+NPER*i].size
        if use_median:
            xm[i], ym[i] = np.median(xi[so][idx+NPER*i]), np.median(yi[so][idx+NPER*i])
        if use_nmad:
            xs[i], ys[i] = nmad(xi[so][idx+NPER*i]), nmad(yi[so][idx+NPER*i])
            
    return xm, ym, ys, N

def medfilt(xarr, N=3, AVERAGE=False):
    """
    Median filter
    """
    out = xarr*0.
    half = int(N/2)
    
    if AVERAGE:
        for i in range(0,half):
            out[i] = np.mean(xarr[i:i+half+1])     
        for i in range(half, len(xarr)-half,1):
            out[i] = np.mean(xarr[i-half:i+half+1])
        for i in range(len(xarr)-half,len(xarr)):
            out[i] = np.mean(xarr[i-half:i])
    else:
        for i in range(0,half):
            out[i] = np.median(xarr[i:i+half+1])     
        for i in range(half, len(xarr)-half,1):
            out[i] = np.median(xarr[i-half:i+half+1])
        for i in range(len(xarr)-half,len(xarr)):
            out[i] = np.median(xarr[i-half:i])
        
    return out

def diff(xarr):
    """
    Like np diff but make same size as input array filling first element
    with diff[0]
    """    
    d = np.diff(xarr)
    return np.append(d[0], d)
    
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
    xout = out[0,:] + x0
    yout = out[1,:] + y0
    return xout, yout

def color_table(value, table='hsv.rgb', normalized=False, show=False):
    """
    Return (r,g,b) values extracted from a color table.
    """
    import threedhst
    import glob
    
    if show:
        files = glob.glob(os.path.dirname(threedhst.__file__)+'/data/*rgb')
        for file in files:
            print(os.path.basename(file))
            
        return (0,0,0)
        
    try:
        data = np.loadtxt(os.path.dirname(threedhst.__file__)+'/data/'+table)
    except:
        'Color table [%s] not found in `threedhst/data`' %(table)
        return (0,0,0)
        
    idx = np.arange(256)
    if normalized:
        idx /= 256.
        
    ri = np.interp(value, idx, data[:,0])/255.
    gi = np.interp(value, idx, data[:,1])/255.
    bi = np.interp(value, idx, data[:,2])/255.
    
    return (ri, gi, bi)
    
def which_3dhst_pointing(ra, dec, regions_file=None, ancillary=False):
    """
    Compute in which 3D-HST pointing(s) a given object lies
    
    Example:
        
    >>> threedhst.utils.which_3dhst_pointing(34.404790, -5.224868)
    ['UDS-17']
    >>> threedhst.utils.which_3dhst_pointing('02:17:37.21','-05:13:27.96')
    ['UDS-17']
    >>> threedhst.utils.which_3dhst_pointing(34.404790, -5.224868, regions_file = 'match.reg')
    ['UDS-17']
    $ cat match.reg
    fk5
    polygon(34.423760, -5.191379, 34.389569,- 5.191188, 34.389733,  -5.229800, 34.423914, -5.229039) # color=magenta width=2  text={UDS-17}
    
    """
    import threedhst
    try:
        import angles
    except:
        threedhst.showMessage('"angles" package not found.  Install with\n\n    $ pip install angles', warn=True)
        return False
        
    pointings = """polygon(214.890094,52.834909,214.930200,52.858851,214.884783,52.886038,214.845805,52.861430) # color=orange width=2 dash=1 text={AEGIS-1}
polygon(214.966899,52.933864,214.923447,52.912124,214.964807,52.882650,215.007231,52.905113) # color=magenta width=2  text={AEGIS-10}
polygon(214.846564,52.867390,214.886700,52.891331,214.841249,52.918519,214.802242,52.893911) # color=magenta width=2  text={AEGIS-11}
polygon(214.794915,52.704310,214.790330,52.738247,214.726847,52.734719,214.732997,52.700875) # color=magenta width=2  text={AEGIS-12}
polygon(214.798609,52.796166,214.745817,52.807992,214.724236,52.771652,214.777552,52.760725) # color=magenta width=2  text={AEGIS-13}
polygon(214.859481,52.802807,214.854886,52.836744,214.791259,52.833216,214.797423,52.799372) # color=magenta width=2  text={AEGIS-14}
polygon(214.843960,52.739074,214.839372,52.773011,214.775838,52.769483,214.781993,52.735638) # color=magenta width=2  text={AEGIS-15}
polygon(214.928560,52.910372,214.885131,52.888632,214.926468,52.859158,214.968869,52.881622) # color=magenta width=2  text={AEGIS-16}
polygon(214.959157,52.991431,214.915648,52.969691,214.957062,52.940216,214.999542,52.962680) # color=magenta width=2  text={AEGIS-17}
polygon(214.834281,52.750919,214.781543,52.762745,214.759985,52.726405,214.813245,52.715477) # color=magenta width=2  text={AEGIS-18}
polygon(214.927156,52.962106,214.883676,52.940366,214.925062,52.910891,214.967514,52.933355) # color=magenta width=2  text={AEGIS-19}
polygon(214.928703,52.807737,214.968784,52.831678,214.923395,52.858866,214.884442,52.834258) # color=magenta width=2  text={AEGIS-2}
polygon(214.889326,52.937870,214.845870,52.916129,214.887233,52.886655,214.929661,52.909119) # color=magenta width=2  text={AEGIS-20}
polygon(214.962561,52.886956,214.919156,52.865216,214.960471,52.835741,215.002849,52.858205) # color=magenta width=2  text={AEGIS-21}
polygon(214.744714,52.842284,214.733652,52.808896,214.796424,52.801689,214.805936,52.835248) # color=orange width=2  dash=1 text={AEGIS-22}
polygon(214.792188,52.870264,214.782085,52.836765,214.845097,52.830221,214.853644,52.863876) # color=magenta width=2  text={AEGIS-23}
polygon(214.755271,52.774150,214.702505,52.785975,214.680935,52.749635,214.734224,52.738708) # color=magenta width=2  text={AEGIS-24}
polygon(214.783921,52.735885,214.779333,52.769822,214.715804,52.766294,214.721958,52.732450) # color=magenta width=2  text={AEGIS-25}
polygon(215.081955,52.956181,215.038481,52.934441,215.079862,52.904966,215.122308,52.927430) # color=magenta width=2  text={AEGIS-26}
polygon(215.001675,53.009172,214.958148,52.987432,214.999579,52.957958,215.042077,52.980422) # color=magenta width=2  text={AEGIS-27}
polygon(214.887180,52.771637,214.882589,52.805575,214.819007,52.802046,214.825166,52.768202) # color=magenta width=2  text={AEGIS-28}
polygon(215.021249,52.981108,214.997074,52.950329,215.055278,52.934158,215.078012,52.965331) # color=magenta width=2  text={AEGIS-29}
polygon(214.733474,52.710712,214.728889,52.744650,214.665397,52.741121,214.671547,52.707277) # color=magenta width=2  text={AEGIS-3}
polygon(214.889204,52.834725,214.845852,52.812985,214.887117,52.783511,214.929444,52.805974) # color=orange width=2  dash=1 text={AEGIS-30}
polygon(214.687267,52.779628,214.741578,52.770630,214.757844,52.807969,214.703141,52.816044) # color=magenta width=2  text={AEGIS-4}
polygon(214.645468,52.748467,214.699741,52.739469,214.715996,52.776808,214.661332,52.784883) # color=magenta width=2  text={AEGIS-5}
polygon(215.016286,52.905264,214.965872,52.889977,214.995160,52.855646,215.044843,52.871778) # color=magenta width=2  text={AEGIS-6}
polygon(215.038687,52.933396,214.998543,52.909455,215.044003,52.882267,215.083019,52.906875) # color=magenta width=2  text={AEGIS-7}
polygon(214.998182,52.960761,214.954703,52.939021,214.996088,52.909547,215.038538,52.932010) # color=magenta width=2  text={AEGIS-8}
polygon(214.825179,52.770946,214.820587,52.804883,214.757007,52.801355,214.763166,52.767511) # color=magenta width=2  text={AEGIS-9}
polygon(150.099223,2.391097,150.086084,2.422515,150.050573,2.407277,150.064586,2.376241) # color=magenta width=2  text={COSMOS-1}
polygon(150.139837,2.338875,150.126698,2.370293,150.091189,2.355055,150.105202,2.324018) # color=magenta width=2  text={COSMOS-10}
polygon(150.166533,2.457486,150.153392,2.488904,150.117880,2.473666,150.131894,2.442629) # color=magenta width=2  text={COSMOS-11}
polygon(150.132886,2.268319,150.119747,2.299737,150.084240,2.284499,150.098252,2.253463) # color=magenta width=2  text={COSMOS-12}
polygon(150.134031,2.229605,150.124252,2.262223,150.087346,2.250777,150.098034,2.218447) # color=magenta width=2  text={COSMOS-13}
polygon(150.134692,2.195402,150.121554,2.226820,150.086048,2.211583,150.100060,2.180546) # color=magenta width=2  text={COSMOS-14}
polygon(150.117737,2.435305,150.134611,2.405721,150.167999,2.425169,150.150305,2.454268) # color=magenta width=2  text={COSMOS-15}
polygon(150.118987,2.400583,150.135860,2.370999,150.169248,2.390447,150.151553,2.419546) # color=magenta width=2  text={COSMOS-16}
polygon(150.123220,2.369194,150.140093,2.339610,150.173480,2.359058,150.155786,2.388157) # color=magenta width=2  text={COSMOS-17}
polygon(150.198090,2.380605,150.169614,2.399311,150.148703,2.366838,150.177690,2.348937) # color=magenta width=2  text={COSMOS-18}
polygon(150.194210,2.243268,150.160797,2.249953,150.153584,2.212018,150.187169,2.206270) # color=magenta width=2  text={COSMOS-19}
polygon(150.169037,2.210125,150.155899,2.241543,150.120393,2.226305,150.134405,2.195268) # color=magenta width=2  text={COSMOS-2}
polygon(150.121766,2.330861,150.138639,2.301277,150.172025,2.320725,150.154331,2.349824) # color=magenta width=2  text={COSMOS-20}
polygon(150.169844,2.247208,150.156706,2.278626,150.121199,2.263388,150.135211,2.232352) # color=magenta width=2  text={COSMOS-21}
polygon(150.201314,2.280490,150.167900,2.287175,150.160687,2.249240,150.194273,2.243492) # color=magenta width=2  text={COSMOS-22}
polygon(150.168177,2.281930,150.155039,2.313348,150.119531,2.298111,150.133543,2.267074) # color=magenta width=2  text={COSMOS-23}
polygon(150.085642,2.316972,150.102514,2.287388,150.135900,2.306836,150.118206,2.335935) # color=magenta width=2  text={COSMOS-24}
polygon(150.084444,2.388625,150.097583,2.357207,150.133093,2.372445,150.119080,2.403482) # color=magenta width=2  text={COSMOS-25}
polygon(150.049475,2.303083,150.066348,2.273499,150.099733,2.292947,150.082039,2.322046) # color=magenta width=2  text={COSMOS-26}
polygon(150.049059,2.369750,150.065932,2.340166,150.099319,2.359614,150.081625,2.388713) # color=magenta width=2  text={COSMOS-27}
polygon(150.083281,2.419194,150.100154,2.389610,150.133542,2.409058,150.115847,2.438157) # color=magenta width=2  text={COSMOS-28}
polygon(150.104504,2.324430,150.091365,2.355848,150.055856,2.340611,150.069869,2.309574) # color=magenta width=2  text={COSMOS-3}
polygon(150.203510,2.294986,150.190371,2.326404,150.154863,2.311166,150.168875,2.280129) # color=magenta width=2  text={COSMOS-4}
polygon(150.098665,2.252764,150.085527,2.284182,150.050019,2.268944,150.064031,2.237907) # color=magenta width=2  text={COSMOS-5}
polygon(150.100749,2.216930,150.087611,2.248348,150.052105,2.233111,150.066117,2.202074) # color=magenta width=2  text={COSMOS-6}
polygon(150.096442,2.427208,150.083302,2.458626,150.047791,2.443388,150.061805,2.412352) # color=magenta width=2  text={COSMOS-7}
polygon(150.131636,2.442486,150.118496,2.473904,150.082984,2.458666,150.096998,2.427629) # color=magenta width=2  text={COSMOS-8}
polygon(150.198930,2.330264,150.185792,2.361682,150.150282,2.346444,150.164295,2.315407) # color=magenta width=2  text={COSMOS-9}
polygon( 53.047381,-27.696185, 53.088765,-27.685558, 53.099610,-27.718251, 53.058226,-27.728878) # color=magenta text={ERS}
polygon(188.993160, 62.172507,189.030550, 62.201734,188.959032, 62.221146,188.923459, 62.191427) # color=magenta width=2 text={GOODS-N-11}
polygon(189.041309, 62.194782,189.078223, 62.224148,189.006320, 62.243293,188.971234, 62.213441) # color=magenta width=2 text={GOODS-N-12}
polygon(189.089070, 62.217032,189.126012, 62.246398,189.054055, 62.265543,189.018943, 62.235691) # color=magenta width=2 text={GOODS-N-13}
polygon(189.139822, 62.239687,189.173093, 62.269996,189.098828, 62.287176,189.067451, 62.256430) # color=magenta width=2 text={GOODS-N-14}
polygon(189.184466, 62.261505,189.221463, 62.290871,189.149400, 62.310015,189.114236, 62.280163) # color=magenta width=2 text={GOODS-N-15}
polygon(189.224929, 62.334449,189.166375, 62.313976,189.217001, 62.283400,189.274271, 62.304603) # color=magenta width=2 text={GOODS-N-16}
polygon(189.273891, 62.356590,189.214345, 62.336717,189.263561, 62.305634,189.321854, 62.326251) # color=magenta width=2 text={GOODS-N-17}
polygon(189.319660, 62.378937,189.261613, 62.358104,189.313185, 62.327842,189.369928, 62.349398) # color=magenta width=2 text={GOODS-N-18}
polygon(189.048036, 62.147199,189.084893, 62.176565,189.013103, 62.195710,188.978072, 62.165858) # color=magenta width=2 text={GOODS-N-21}
polygon(189.095755, 62.169421,189.132639, 62.198787,189.060796, 62.217932,189.025739, 62.188080) # color=magenta width=2 text={GOODS-N-22}
polygon(189.149002, 62.192488,189.179131, 62.223491,189.103331, 62.239028,189.075140, 62.207630) # color=magenta width=2 text={GOODS-N-23}
polygon(189.191193, 62.213894,189.228131, 62.243260,189.156182, 62.262404,189.121074, 62.232552) # color=magenta width=2 text={GOODS-N-24}
polygon(189.238912, 62.236116,189.275878, 62.265482,189.203875, 62.284626,189.168741, 62.254774) # color=magenta width=2 text={GOODS-N-25}
polygon(189.283487, 62.308738,189.221887, 62.290331,189.267480, 62.258100,189.327909, 62.277279) # color=magenta width=2 text={GOODS-N-26}
polygon(189.332135, 62.330893,189.269821, 62.312978,189.314260, 62.280390,189.375426, 62.299087) # color=magenta width=2 text={GOODS-N-27}
polygon(189.379764, 62.353125,189.317467, 62.335163,189.362052, 62.302608,189.423199, 62.321352) # color=magenta width=2 text={GOODS-N-28}
polygon(189.102483, 62.121810,189.139308, 62.151176,189.067578, 62.170321,189.032576, 62.140469) # color=magenta width=2 text={GOODS-N-31}
polygon(189.150243, 62.144032,189.187095, 62.173398,189.115313, 62.192543,189.080286, 62.162691) # color=magenta width=2 text={GOODS-N-32}
polygon(189.197920, 62.166282,189.234800, 62.195648,189.162965, 62.214793,189.127912, 62.184941) # color=magenta width=2 text={GOODS-N-33}
polygon(189.245681, 62.188505,189.282588, 62.217871,189.210699, 62.237015,189.175621, 62.207163) # color=magenta width=2 text={GOODS-N-34}
polygon(189.286708, 62.261414,189.227863, 62.241212,189.277713, 62.210404,189.335292, 62.231342) # color=magenta width=2 text={GOODS-N-35}
polygon(189.341077, 62.232977,189.378039, 62.262343,189.306044, 62.281487,189.270914, 62.251636) # color=magenta width=2 text={GOODS-N-36}
polygon(189.128256, 62.096387,189.191491, 62.113221,189.149924, 62.146554,189.087682, 62.128882) # color=magenta width=2 text={GOODS-N-41}
polygon(189.204689, 62.118671,189.241511, 62.148037,189.169789, 62.167182,189.134791, 62.137330) # color=magenta width=2 text={GOODS-N-42}
polygon(189.252367, 62.140893,189.289215, 62.170260,189.217440, 62.189404,189.182416, 62.159552) # color=magenta width=2 text={GOODS-N-43}
polygon(189.296421, 62.162711,189.337518, 62.190829,189.268675, 62.212225,189.229320, 62.183564) # color=magenta width=2 text={GOODS-N-44}
polygon(189.336107, 62.184499,189.385510, 62.209552,189.323980, 62.235365,189.276114, 62.209659) # color=magenta width=2 text={GOODS-N-45}
polygon(189.395565, 62.207588,189.432496, 62.236954,189.360561, 62.256099,189.325460, 62.226247) # color=magenta width=2 text={GOODS-N-46}
polygon(53.253844,-27.661061,53.227953,-27.686237,53.260483,-27.711953,53.285567,-27.686147) # color=magenta width=2  text={GOODS-S-1}
polygon(53.146446,-27.865471,53.132847,-27.833614,53.092140,-27.847604,53.106741,-27.879111) # color=magenta width=2  text={GOODS-S-10}
polygon(53.172026,-27.691544,53.133730,-27.688386,53.130113,-27.726866,53.168491,-27.729074) # color=magenta width=2  text={GOODS-S-11}
polygon(53.092665,-27.715852,53.054361,-27.712694,53.050743,-27.751174,53.089130,-27.753383) # color=magenta width=2  text={GOODS-S-12}
polygon(53.160758,-27.896771,53.147155,-27.864914,53.106437,-27.878904,53.121042,-27.910411) # color=magenta width=2  text={GOODS-S-13}
polygon(53.215197,-27.914693,53.201591,-27.882836,53.160866,-27.896827,53.175474,-27.928333) # color=magenta width=2  text={GOODS-S-14}
polygon(53.135409,-27.942252,53.121800,-27.910395,53.081065,-27.924385,53.095676,-27.955891) # color=magenta width=2  text={GOODS-S-15}
polygon(53.175364,-27.928388,53.161756,-27.896531,53.121026,-27.910521,53.135636,-27.942027) # color=magenta width=2  text={GOODS-S-16}
polygon(53.128527,-27.742861,53.090213,-27.739702,53.086594,-27.778182,53.124990,-27.780391) # color=magenta width=2  text={GOODS-S-17}
polygon(53.165707,-27.750025,53.168176,-27.784006,53.211694,-27.781149,53.208153,-27.747245) # color=magenta width=2  text={GOODS-S-18}
polygon(53.043416,-27.824225,53.045885,-27.858206,53.089434,-27.855349,53.085890,-27.821445) # color=magenta width=2  text={GOODS-S-19}
polygon(53.040408,-27.790345,53.042877,-27.824325,53.086412,-27.821469,53.082869,-27.787565) # color=magenta width=2  text={GOODS-S-2}
polygon(53.085097,-27.812309,53.087566,-27.846289,53.131110,-27.843433,53.127566,-27.809529) # color=magenta width=2  text={GOODS-S-20}
polygon(53.200584,-27.883018,53.186983,-27.851161,53.146270,-27.865152,53.160873,-27.896658) # color=magenta width=2  text={GOODS-S-21}
polygon(53.254855,-27.901207,53.241251,-27.869350,53.200531,-27.883340,53.215137,-27.914847) # color=magenta width=2  text={GOODS-S-22}
polygon(53.016930,-27.744426,53.054140,-27.735798,53.042672,-27.698543,53.005752,-27.708089) # color=magenta width=2  text={GOODS-S-23}
polygon(53.028220,-27.783186,53.062602,-27.767899,53.042628,-27.733569,53.008745,-27.749701) # color=magenta width=2  text={GOODS-S-24}
polygon(53.049980,-27.768050,53.015004,-27.782267,53.033618,-27.817193,53.068129,-27.802117) # color=magenta width=2  text={GOODS-S-25}
polygon(53.066424,-27.890689,53.104925,-27.889692,53.103216,-27.851108,53.064764,-27.853057) # color=magenta width=2  text={GOODS-S-26}
polygon(53.224984,-27.840796,53.215228,-27.807857,53.173101,-27.818009,53.183895,-27.850692) # color=magenta width=2  text={GOODS-S-27}
polygon(53.303424,-27.849997,53.271617,-27.830798,53.247348,-27.862901,53.279749,-27.881305) # color=magenta width=2  text={GOODS-S-28}
polygon(53.208351,-27.715639,53.170046,-27.712480,53.166428,-27.750960,53.204815,-27.753169) # color=magenta width=2  text={GOODS-S-29}
polygon(53.127254,-27.669806,53.129720,-27.703786,53.173207,-27.700930,53.169668,-27.667026) # color=magenta width=2  text={GOODS-S-3}
polygon(53.132156,-27.705400,53.093855,-27.702241,53.090237,-27.740721,53.128620,-27.742930) # color=magenta width=2  text={GOODS-S-30}
polygon(53.082352,-27.778275,53.084821,-27.812256,53.128351,-27.809399,53.124809,-27.775495) # color=magenta width=2  text={GOODS-S-31}
polygon(53.123970,-27.763973,53.126438,-27.797953,53.169962,-27.795096,53.166420,-27.761193) # color=magenta width=2  text={GOODS-S-32}
polygon(53.168537,-27.729019,53.130228,-27.725861,53.126609,-27.764341,53.165001,-27.766549) # color=magenta width=2  text={GOODS-S-33}
polygon(53.130655,-27.786407,53.160428,-27.807984,53.187750,-27.777874,53.157308,-27.757043) # color=magenta width=2  text={GOODS-S-34}
polygon(53.085677,-27.683592,53.088144,-27.717572,53.131636,-27.714716,53.128097,-27.680812) # color=magenta width=2  text={GOODS-S-35}
polygon(53.133767,-27.785174,53.163540,-27.806751,53.190862,-27.776640,53.160421,-27.755810) # color=magenta width=2  text={GOODS-S-36}
polygon(53.129712,-27.788174,53.159486,-27.809751,53.186809,-27.779640,53.156366,-27.758810) # color=magenta width=2  text={GOODS-S-37}
polygon(53.132938,-27.787807,53.162711,-27.809384,53.190034,-27.779274,53.159591,-27.758443) # color=magenta width=2  text={GOODS-S-38}
polygon(53.186508,-27.851630,53.172910,-27.819772,53.132209,-27.833763,53.146808,-27.865269) # color=magenta width=2  text={GOODS-S-4}
polygon(53.117080,-27.914693,53.103474,-27.882836,53.062750,-27.896827,53.077357,-27.928333) # color=magenta width=2  text={GOODS-S-5}
polygon(53.240651,-27.869538,53.227051,-27.837681,53.186343,-27.851671,53.200945,-27.883177) # color=magenta width=2  text={GOODS-S-6}
polygon(53.124592,-27.795573,53.127061,-27.829553,53.170598,-27.826696,53.167055,-27.792793) # color=magenta width=2  text={GOODS-S-7}
polygon(53.167125,-27.784323,53.169594,-27.818303,53.213126,-27.815446,53.209584,-27.781543) # color=magenta width=2  text={GOODS-S-8}
polygon(53.089071,-27.753161,53.050753,-27.750002,53.047134,-27.788482,53.085534,-27.790691) # color=magenta width=2  text={GOODS-S-9}
polygon( 53.111766,-27.843283, 53.093896,-27.813137, 53.055472,-27.831405, 53.074280,-27.861099) # color=magenta text={SN-GEORGE}
polygon( 34.455240, -5.243509, 34.421062, -5.243304, 34.421240, -5.281884, 34.455410, -5.281136) # color=magenta text={SN-MARSHALL-0}
polygon( 34.460418, -5.250832, 34.428370, -5.238999, 34.415286, -5.275313, 34.447653, -5.286248) # color=magenta text={SN-MARSHALL-1}
polygon( 53.138492,-27.758488, 53.136643,-27.792488, 53.180229,-27.793926, 53.180988,-27.759895) # color=magenta text={SN-PRIMO}
polygon(34.491140,-5.153326,34.456951,-5.153135,34.457115,-5.191748,34.491294,-5.190987) # color=magenta width=2  text={UDS-1}
polygon(34.376053,-5.248102,34.348154,-5.228414,34.326047,-5.260138,34.354487,-5.279039) # color=magenta width=2  text={UDS-10}
polygon(34.457671,-5.190260,34.423480,-5.190068,34.423644,-5.228681,34.457825,-5.227920) # color=magenta width=2  text={UDS-11}
polygon(34.457233,-5.228174,34.423039,-5.227982,34.423203,-5.266595,34.457387,-5.265834) # color=magenta width=2  text={UDS-12}
polygon(34.331462,-5.189694,34.362367,-5.204258,34.378603,-5.169193,34.347303,-5.155497) # color=magenta width=2  text={UDS-13}
polygon(34.423789,-5.228596,34.389596,-5.228404,34.389760,-5.267017,34.423943,-5.266256) # color=magenta width=2  text={UDS-14}
polygon(34.294046,-5.224082,34.266147,-5.204395,34.244042,-5.236118,34.272480,-5.255020) # color=magenta width=2  text={UDS-15}
polygon(34.319166,-5.244215,34.291267,-5.224528,34.269161,-5.256252,34.297601,-5.275153) # color=magenta width=2  text={UDS-16}
polygon(34.423760,-5.191379,34.389569,-5.191188,34.389733,-5.229800,34.423914,-5.229039) # color=magenta width=2  text={UDS-17}
polygon(34.581579,-5.142743,34.557268,-5.166684,34.584799,-5.193872,34.608427,-5.169264) # color=magenta width=2  text={UDS-18}
polygon(34.295047,-5.188305,34.325952,-5.202869,34.342188,-5.167805,34.310888,-5.154108) # color=magenta width=2  text={UDS-19}
polygon(34.491161,-5.190818,34.456970,-5.190626,34.457134,-5.229239,34.491316,-5.228478) # color=magenta width=2  text={UDS-2}
polygon(34.272577,-5.172115,34.244681,-5.152428,34.222577,-5.184152,34.251013,-5.203053) # color=magenta width=2  text={UDS-20}
polygon(34.264877,-5.208368,34.236980,-5.188681,34.214875,-5.220404,34.243313,-5.239306) # color=magenta width=2  text={UDS-21}
polygon(34.274229,-5.247827,34.246329,-5.228139,34.224223,-5.259863,34.252663,-5.278764) # color=magenta width=2  text={UDS-22}
polygon(34.234691,-5.161227,34.265594,-5.175791,34.281830,-5.140727,34.250532,-5.127031) # color=magenta width=2  text={UDS-23}
polygon(34.423750,-5.152490,34.389561,-5.152299,34.389725,-5.190911,34.423904,-5.190151) # color=magenta width=2  text={UDS-24}
polygon(34.310968,-5.154419,34.341871,-5.168983,34.358107,-5.133918,34.326809,-5.120222) # color=magenta width=2  text={UDS-25}
polygon(34.268621,-5.166497,34.299525,-5.181060,34.315761,-5.145996,34.284462,-5.132300) # color=magenta width=2  text={UDS-26}
polygon(34.251970,-5.201358,34.282876,-5.215922,34.299113,-5.180857,34.267812,-5.167161) # color=magenta width=2  text={UDS-27}
polygon(34.390698,-5.229568,34.356504,-5.229376,34.356668,-5.267989,34.390852,-5.267228) # color=magenta width=2  text={UDS-28}
polygon(34.491141,-5.228313,34.456948,-5.228121,34.457112,-5.266734,34.491295,-5.265973) # color=magenta width=2  text={UDS-3}
polygon(34.457668,-5.152899,34.423478,-5.152707,34.423643,-5.191320,34.457822,-5.190559) # color=magenta width=2  text={UDS-4}
polygon(34.347388,-5.155249,34.378291,-5.169813,34.394527,-5.134749,34.363229,-5.121053) # color=magenta width=2  text={UDS-5}
polygon(34.361166,-5.203999,34.392072,-5.218563,34.408309,-5.183499,34.377008,-5.169803) # color=magenta width=2  text={UDS-6}
polygon(34.396167,-5.214486,34.366652,-5.197295,34.347407,-5.230817,34.377392,-5.247178) # color=magenta width=2  text={UDS-7}
polygon(34.367695,-5.201449,34.339798,-5.181761,34.317693,-5.213485,34.346131,-5.232386) # color=magenta width=2  text={UDS-8}
polygon(34.346372,-5.232274,34.318473,-5.212586,34.296367,-5.244310,34.324806,-5.263211) # color=magenta width=2  text={UDS-9}"""
    
    #### CANDELS SN and COOPER EGS pointings
    if ancillary:
        pointings += """\npolygon( 53.047381,-27.696185, 53.088765,-27.685558, 53.099610,-27.718251, 53.058226,-27.728878) # color=cyan width=2 text={ERS}
polygon( 53.138492,-27.758488, 53.136643,-27.792488, 53.180229,-27.793926, 53.180988,-27.759895) # color=cyan text={SN-PRIMO}
polygon( 53.111766,-27.843283, 53.093896,-27.813137, 53.055472,-27.831405, 53.074280,-27.861099) # color=cyan text={SN-GEORGE}
polygon( 34.455240, -5.243509, 34.421062, -5.243304, 34.421240, -5.281884, 34.455410, -5.281136) # color=cyan text={SN-MARSHALL}
polygon(214.275840,52.475143,214.221331,52.467670,214.235844,52.430083,214.289985,52.438482) # color=cyan width=2 text={COOPER/EGS12004280}
polygon(214.132993,52.411228,214.183039,52.396129,214.210554,52.430905,214.159838,52.445144) # color=cyan width=2 text={COOPER/EGS12004754}
polygon(214.517658,52.546031,214.464967,52.534565,214.486923,52.498337,214.539063,52.510694) # color=cyan width=2 text={COOPER/EGS12007881.12012083}
polygon(214.643818,52.562188,214.591307,52.574014,214.569842,52.537674,214.622873,52.526747) # color=cyan width=2 dash=1 text={COOPER/EGS12011767}
polygon(214.721838,52.607467,214.682415,52.631679,214.637647,52.604260,214.678167,52.580728) # color=cyan width=2 text={COOPER/EGS12015684}
polygon(214.567181,52.670648,214.535231,52.698650,214.483221,52.676369,214.516449,52.648921) # color=cyan width=2 dash=1 text={COOPER/EGS12020405}
polygon(214.549134,52.726809,214.542945,52.692966,214.606352,52.689092,214.610975,52.723025) # color=cyan width=2 dash=1 text={COOPER/EGS12024866}
polygon(214.630347,52.835657,214.607146,52.804629,214.665613,52.789073,214.687368,52.820480) # color=cyan width=2 text={COOPER/EGS13004661}
polygon(214.883627,52.914686,214.834597,52.897827,214.866829,52.864469,214.915055,52.882148) # color=cyan width=2 dash=1 text={COOPER/EGS13011148}
polygon(214.835962,52.882940,214.802256,52.910249,214.751314,52.886881,214.786271,52.860153) # color=cyan width=2 text={COOPER/EGS13011439}
polygon(215.174708,53.026993,215.142495,53.054995,215.090056,53.032713,215.123558,53.005266) # color=cyan width=2 dash=1 text={COOPER/EGS13034445}
polygon(214.994022,53.040119,215.013087,53.008057,215.073327,53.021417,215.052779,53.053144) # color=cyan width=2 dash=1 text={COOPER/EGS13035123}"""
    
    if isinstance(ra, str):
        ra = DMS2decimal(ra, hours=True)
    #
    if isinstance(dec, str):
        dec = DMS2decimal(dec, hours=False)
       
    matched_idx = [] 
    matched_list = []
    
    for i, polystr in enumerate(pointings.split('\n')):
        #print polystr
        poly = polystr[polystr.find('(')+1:polystr.find(')')].split(',')
        pointing = polystr[polystr.find('{')+1:polystr.find('}')]
        px = np.cast[float](poly[0::2])
        py = np.cast[float](poly[1::2])
        if threedhst.regions.point_in_polygon(ra, dec, px, py):
            matched_list.append(pointing)
            matched_idx.append(i)
    
    if regions_file:
        spl = pointings.split('\n')
        fp = open(regions_file, 'w')
        fp.write('fk5\n')
        for idx in matched_idx:
            fp.write(spl[idx])
        
        fp.close()
        
    return matched_list
    
def decimal2HMS(input, hours=False):
    """
    Convert decimal degrees to DD:MM:SS or HH:MM:SS.
    
    >>> ra = 150.06852667
    >>> print threedhst.utils.decimal2HMS(ra, hours=True)
    10:00:16.45
    
    """
    if hours:
        value = input*24/360.
    else:
        value = input
        
    sign = value < 0
    if (hours):
        pm = ''
    else:
        if sign:
            pm = '-'
        else:
            pm = '+'
        
    deg = np.abs(int(value))
    rem = np.abs(value)-deg
    min = int(rem*60)
    sec = (rem-min/60.)*3600
    
    return '%s%02d:%02d:%05.2f' %(pm, deg, min, sec)
    
def DMS2decimal(string, hours=True):
    """
    Convert degree-minute-second string into decimal value.
    
    EXAMPLE:
    
    >>> print threedhst.utils.DMS2decimal('100127.5', hours=True)
    >>> print threedhst.utils.DMS2decimal('10:01:27.5', hours=True)
    >>> print threedhst.utils.DMS2decimal('+02:00:00')
    
    """
    
    if string.startswith('+'):
        string = string[1:]
        
    split_DMS = string.split(':')
    #### If no colons in the string, assume DDMMSS.S
    if len(split_DMS) != 3:
        repl = string.replace('-','')
        split_DMS = [repl[0:2], repl[2:4], repl[4:]]
        #print split_DMS
        
    dms = np.abs(np.cast[float](split_DMS))
    if string.startswith('-'):
        neg = True
    else:
        neg = False
            
    decimal = dms[0]+dms[1]/60.+dms[2]/3600.
    if hours:
        decimal *= 360./24
    
    if neg:
        decimal *= -1
    
    return decimal
    
def image_size_and_wcs(input="test.fits", extension=1):
    """
    Get image size and center world cordinate for a North-up image (e.g. drz.fits)
    """
    header = pyfits.getheader(input,extension)
    NX, NY = header['NAXIS1'], header['NAXIS2']
    center_ra = header['CRVAL1']+(NX/2.-header['CRPIX1'])*header['CD1_1']
    center_dec = header['CRVAL2']+(NY/2.-header['CRPIX2'])*header['CD2_2']
    return (NX, NY, center_ra, center_dec, np.abs(header['CD2_2']*3600.))
    
def subimage(input="test.fits", output="sibimage.fits", ra=0, dec=0, size=10, ext=0, verbose=True):
    """
    Extract a thumbnail from a larger image, centered on physical 
    coordinate (ra,dec) with `size` arcsec on a side
    """    
    # ra, dec = 34.404794, -5.2248747
    # size=40
    # input, ext = 'UDS-17-F140W_drz.fits', 1
    #
    im = pyfits.open(input)
    wcs = im[ext].header
    #### for now, assume aligned pixels N-up / E-left, square pixels
    xpix = (ra-wcs['CRVAL1'])/wcs['CD1_1']+wcs['CRPIX1']
    ypix = (dec-wcs['CRVAL2'])/wcs['CD2_2']+wcs['CRPIX2']
    #
    xr, yr = np.round(xpix), np.round(ypix)
    size_pix = np.round(size / np.abs(wcs['CD1_1']) / 3600. / 2)
    #
    subim = im[ext].data[yr-size_pix:yr+size_pix, xr-size_pix:xr+size_pix]
    #
    crval1 = (xr-size_pix-wcs['CRPIX1'])*wcs['CD1_1']+wcs['CRVAL1']
    crval2 = (yr-size_pix-wcs['CRPIX2'])*wcs['CD2_2']+wcs['CRVAL2']
    wcs['CRVAL1'] = crval1
    wcs['CRVAL2'] = crval2
    wcs['CRPIX1'] = 0
    wcs['CRPIX2'] = 0
    wcs['NAXIS1'] = size_pix*2
    wcs['NAXIS2'] = size_pix*2
    #
    out = pyfits.PrimaryHDU(data=subim, header=wcs)
    out.writeto(output, clobber=True)
    #
    if verbose:
        print('Wrote: %s' %(output))
        
def gen_tempname(root='tmp'):
    """
    Generate a temporary filename from the current time.time().
    
    Default is like 'tmp12345'
    """
    import time
    t0 = int(time.time()*100 % 1.e5)
    return '%s%05d' %(root, t0)
    
def contiguous_extent(array, x0, y0):
    """
    Find extent of contigous region of a segmentation image.
    
    compare to scipy.ndimage.label --> should be better and faster
    
    im = np.zeros((2048,2048))
    yi, xi = np.indices(im.shape)
    r = np.sqrt((xi-1024)**2+(yi-1024)**2)
    object = (r < 10)*1.
    
    xmi, xma, ymi, yma = threedhst.utils.contiguous_extent(object, 1024, 1024)
    
    import scipy.ndimage as nd
    s = [[0,1,0],
         [1,1,1],
         [0,1,0]]
    
    labeled_array, num_features = nd.label(object, structure=s)
    
    """
    import scipy.ndimage as nd
    s = [[0,1,0],
         [1,1,1],
         [0,1,0]]
    
    labeled_array, num_features = nd.label(array, structure=s)
    matched = labeled_array == labeled_array[y0, x0]
    yi, xi = np.indices(array.shape)
    return xi[matched].min(), xi[matched].max(), yi[matched].min(), yi[matched].max()
    
    # mask = array == array[y0, x0]    
    # 
    # ### Grow the mask by shifting 1 pix in each direction to join a region 
    # ### like
    # #             o
    # #              o    oo
    # #              oooooo
    # #               ooooo            
    # mask = mask | np.roll(mask, 1, axis=1) | np.roll(mask, -1, axis=1)
    # mask = mask | np.roll(mask, 1, axis=0) | np.roll(mask, -1, axis=0)
    # 
    # ##### First pass
    # xmi, xma = x0, x0
    # ymi, yma = y0, y0
    # 
    # ## compare current row with row above and below to find if any overlap
    # while (mask[:, xma] & mask[:, xma+1]).sum() > 0:
    #     xma += 1
    #     
    # while (mask[:, xmi] & mask[:, xmi-1]).sum() > 0:
    #     xmi -= 1
    # 
    # ## compare columns
    # while (mask[yma, xmi:xma] & mask[yma+1, xmi:xma]).sum() > 0:
    #     yma += 1
    # 
    # while (mask[ymi, xmi:xma] & mask[ymi-1, xmi:xma]).sum() > 0:
    #     ymi -= 1
    # 
    # ##### Second pass on x after finding y extent
    # xmi, xma = x0, x0
    # while (mask[ymi:yma, xma] & mask[ymi:yma, xma+1]).sum() > 0:
    #     xma += 1
    # 
    # while (mask[ymi:yma, xmi] & mask[ymi:yma, xmi-1]).sum() > 0:
    #     xmi -= 1
    # 
    # return xmi+1, xma, ymi+1, yma
    #return xmi+1, xma-1, ymi+1, yma-1
    
def calc_mag(wavelength, flux, xfilt, yfilt, fnu_units=False, CCD=True):
    """
    Integrate a spectrum through a filter passband.  
    
    If CCD=True, assume filter is for a photon-counting detector
    (see Maiz-Apellaniz 2006)
    
    #### Photon-counting
    flux = int(flam * R * lam, lam) / int(R * f_ref_lam * lam, lam)
    
      ## AB mags
      f_ref_lam = 3631 Jy * 3.e18 / lam**2 = 3631.e-23 * 3.e18 / lam**2 (flam)
      flux = int(flam * R * lam, lam) / int(R * 3631 Jy * 3.e18 / lam, lam)
    
    #### Energy-counting (photometer???)
    flux = int(flam * R, lam) / int(R * f_ref_lam, lam)
    
    """
    
    yf_int = np.interp(wavelength, xfilt, yfilt, left=0, right=0)
    
    if not fnu_units:
        fnu = flux*wavelength**2/3.e18
    else:
        fnu = flux
    
    if CCD:
        detector_scale = wavelength
    else:
        detector_scale = np.ones(len(wavelength))
    #
    filter_flux = np.trapz(fnu*3.e18/wavelength**2*yf_int*detector_scale, wavelength)/np.trapz(yf_int*detector_scale*3631.e-23*3.e18/wavelength**2, wavelength)
    filter_mag = -2.5*np.log10(filter_flux) #-48.6
    return filter_mag
    
def survey_area(ra_in, dec_in):
    """
    Compute survey area sampled by a list of ra/dec points
    
    Requires Shapely: http://toblerity.github.com/shapely/manual.html
    
    Returns: (x, y, area, Shapely)
        x, y = Boundary arrays
        area = area in square arcsec
        Shapely = shapely polygon, etc. object for the convex hull around the 
                  points.  Note that the point coordinates are scaled to the 
                  median positions of the input ra/dec coordinates
                      
    Example:
    
    xfoot, yfoot, area, geom = survey_area(ra, dec)
    plt.scatter(ra, dec)
    plt.plot(xfoot, yfoot, color='red')
    
    """
    from shapely.geometry import Point, MultiPoint

    d0 = np.median(dec_in)
    ra = (ra_in-np.median(ra_in))*np.cos(dec_in/360.*2*np.pi)*3600.
    de = (dec_in-np.median(dec_in))*3600.
    N = ra.shape[0]
    
    points = []
    for i in range(N):
        points.append((ra[i], de[i]))
    
    points_list = MultiPoint(points)
    points_geom = points_list.convex_hull
    xb, yb = points_geom.boundary.xy
    
    xout = np.asarray(xb)/3600./np.cos(d0/360.*2*np.pi) + np.median(ra_in)
    yout = np.asarray(yb)/3600. + d0
    
    return xout, yout, points_geom.area, points_geom
    ##
    # from shapely.ops import cascaded_union
    # polygons = [Point(ra[i], de[i]).buffer(pointsize) for i in range(N)]
    # geom = cascaded_union(polygons)
    # x, y = geom.boundary.xy
    # plt.scatter(ra, de)
    # plt.plot(xg, yg, color='red')
    # 
    # The function is particularly useful in dissolving MultiPolygons.
    # 
    # m = MultiPolygon(polygons)
    # geom = cascaded_union(m)
    
def gethead(image, ext=0, keys=['EXPTIME'], parse_dtype=True, make_dict=False):
    """
    Shell wrapper around wcstools "gethead" to extract header keywords much
    faster than `pyfits.getheader`.
    
    `ext`: image extension
    `keys`: list of header keywords to extract
    `make_dict`: return a dictionary with the keywords as keys
    `parse_dtype`: parse keyword value strings to float/int if possible
    """

    from subprocess import Popen,PIPE
    stdout, stderr = Popen('gethead -a -x %d %s %s' %(ext, image, ' '.join(keys)), shell=True, stdout=PIPE).communicate()
    result = []
    if not parse_dtype:
        return stdout.split()[1:]
        
    for key in stdout.split()[1:]:
        try:
            val = float(key)
            if val.is_integer():
                result.append(int(val))
            else:
                result.append(val)
        except:
            result.append(key)
    
    if make_dict:
        d = {}
        for i in range(len(keys)):
            d[keys[i]] = result[i]
        
        return d
    else:
        return result
#
def pointing_region(ra=177.400395109, dec=22.4054119338, orient=291, name='', optimal_grism=False):
    """
    Make region file for arbitrary coords / ORIENT
    """
    
    #### GOODS-S 03 as reference
    # import stwcs.wcsutil
    # wfc3 = pyfits.open('ibhj03xmq_flt.fits')
    # wcs = stwcs.wcsutil.HSTWCS(wfc3, ext=('SCI',1))
    # footprint = wcs.calc_footprint()
    # ### Optimal grism offset to the left from the center pixel
    # rd0 = wcs.wcs_pix2world([507-91], [507],1)
    
    footprint = np.array([[ 53.12761393, -27.67010485],
           [ 53.13000295, -27.70407857],
           [ 53.17346751, -27.70126514],
           [ 53.16999254, -27.66736672]])
    
    wx, wy = footprint[:,0], footprint[:,1]
    x0, y0 = 53.15022231714418, -27.68537223320745  ## RA_APER / DEC_APER
    if optimal_grism:
        x0, y0 = 53.14636648, -27.68562166
        
    ## ACS
    # import stwcs.wcsutil
    # acs = pyfits.open('jbhj03xnq_flc.fits')
    # wcs = stwcs.wcsutil.HSTWCS(acs, ext=('SCI',2))
    #footprint = wcs.calc_footprint()
    
    footprint = np.array([[ 53.2032301 , -27.78275401],
           [ 53.1781686 , -27.7653409 ],
           [ 53.13469503, -27.80643922],
           [ 53.15993549, -27.82475767]])
    
    ax1, ay1 = footprint[:,0], footprint[:,1]
    
    footprint = np.array([[ 53.17757919, -27.76486217],
           [ 53.15311246, -27.7482801 ],
           [ 53.10962477, -27.78829072],
           [ 53.13406647, -27.80594152]])
    
    ax2, ay2 = footprint[:,0], footprint[:,1]
    
    wx = (wx-x0)*np.cos(y0/180*np.pi); wy = (wy-y0)
    ax1 = (ax1-x0)*np.cos(y0/180*np.pi); ay1 = (ay1-y0)
    ax2 = (ax2-x0)*np.cos(y0/180*np.pi); ay2 = (ay2-y0)
    x0, y0 = 0, 0
    
    # ### Reference is taken from COSMOS-28 actual images
    # stri = '150.083305,  2.419194,150.100164,  2.389610,150.133522,  2.409058,150.115843,  2.438157'
    # wx = np.cast[float](stri.split(',')[0::2])                
    # wy = np.cast[float](stri.split(',')[1::2])
    # 
    # 
    # stri = '150.193030,  2.347347,150.165445,  2.353287,150.148879,  2.299487,150.176995,  2.292809'
    # ax1 = np.cast[float](stri.split(',')[0::2])                
    # ay1 = np.cast[float](stri.split(',')[1::2])
    # 
    # stri = '150.16557,2.3535789,150.13799,2.3595173,150.12142,2.3057169,150.14954,2.2990394'
    # ax2 = np.cast[float](stri.split(',')[0::2])                
    # ay2 = np.cast[float](stri.split(',')[1::2])
    
    # plt.plot(wx, wy)
    # plt.plot(ax1,ay1)
    # plt.plot(ax2,ay2)

    #x0 = np.mean(wx); y0 = np.mean(wy)
    #x0 = 150.110890833; y0 = 2.41361111
    ax1 -= x0; ay1 -= y0
    ax2 -= x0; ay2 -= y0
    wx -= x0; wy -= y0
    
    wx = np.append(wx, wx[0]); wy = np.append(wy, wy[0])
    ax1 = np.append(ax1, ax1[0]); ay1 = np.append(ay1, ay1[0])
    ax2 = np.append(ax2, ax2[0]); ay2 = np.append(ay2, ay2[0])
    
    angle = orient-311
    
    wx, wy = xyrot(wx, wy, angle)
    ax1, ay1 = xyrot(ax1, ay1, angle)
    ax2, ay2 = xyrot(ax2, ay2, angle)
    
    wx = wx/np.cos(dec/360*2*np.pi)+ra; wy += dec
    ax1 = ax1/np.cos(dec/360*2*np.pi)+ra; ay1 += dec
    ax2 = ax2/np.cos(dec/360*2*np.pi)+ra; ay2 += dec
    
    wfc = 'polygon('
    a1 = 'polygon('
    a2 = 'polygon('
    for i in range(4):
        if i < 3:
            com = ','
        else:
            com = ')'
        
        wfc += '%f,%f' %(wx[i],wy[i]) + com
        a1 += '%f,%f' %(ax1[i],ay1[i]) + com
        a2 += '%f,%f' %(ax2[i],ay2[i]) + com
    
    #if (pi.status == 'Archived') | (pi.status == 'Executed'):
    dash = ''
    font='helvetica 9 bold roman'
    acs_color="cyan"
    wfc_color="magenta"
                
    out = """fk5
%s # color=%s width=2 %s
%s # color=%s width=2 %s
%s # color=%s width=2 %s
""" %(wfc, wfc_color, dash, a1, acs_color, dash, a2, acs_color, dash)
    
    if name != '':
        out += "# text(%f,%f) %s font=\"%s\" color=%s" %( ra, dec, name_label, font, wfc_color)
        
    return out
            
def roll2(array, dx, dy):
    """
    `numpy.roll` on both x & y axes
    """
    return np.roll(np.roll(array, dx, axis=1), dy, axis=0)
    