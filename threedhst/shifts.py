"""
3DHST.shifts

Utilities for computing shifts and processing shiftfiles.
"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import os
import pyfits
import pyraf
from pyraf import iraf
from iraf import stsdas,dither

import threedhst

no = iraf.no
yes = iraf.yes
INDEF = iraf.INDEF

def compute_shifts(asn_direct):
    """compute_shifts(asn_direct)"""
        
    root = asn_direct.split('_asn.fits')[0].lower()
    ### Get shifts with tweakshifts
    iraf.flpr()
    iraf.flpr()
    iraf.flpr()
    iraf.tweakshifts(input=asn_direct, shiftfile='',
                     reference=root+'_tweak.fits',
                     output = root+'_shifts.txt', findmode = 'catalog',
                     gencatalog = 'daofind', sextractpars = '', 
                     undistort = yes, computesig = yes, idckey = 'idctab', \
       clean = yes, verbose = no, catfile = '', xcol = 1, ycol = 2, \
       fluxcol = 3, fluxmax = INDEF, fluxmin = INDEF, fluxunits = 'counts', \
       nbright = INDEF, refcat = '', refxcol = 1, refycol = 2, rfluxcol = 3, \
       rfluxmax = INDEF, rfluxmin = INDEF, rfluxunits = 'counts', \
       refnbright = INDEF, minobj = 15, nmatch = 30, matching = 'tolerance', \
       xyxin = INDEF, xyyin = INDEF, tolerance = 1.0, fwhmpsf = 2.5, \
       sigma = 0.0, datamin = INDEF, datamax = INDEF, threshold = 4.0, \
       nsigma = 1.5, fitgeometry = 'rscale', function = 'polynomial', \
       maxiter = 3, reject = 3.0, crossref = '', margin = 50, tapersz = 50, \
       pad = no, fwhm = 7.0, ellip = 0.05, pa = 45.0, fitbox = 7, \
    )
    ### !!! Need to add steps to match WCS to astrometric references (e.g. ACS)
    
def align_to_reference():
    """
xshift, yshift = align_to_reference()
    """        
    import glob
    import os
    
    #### Clean slate
    rmfiles = ['SCI.fits','align.cat',
               'align.map','align.match','align.reg','align.xy',
               'direct.cat','direct.reg','direct.xy']
    for file in rmfiles:
        try:
            os.remove(file)
        except:
            pass
            
    root_direct = threedhst.currentRun['root_direct']
    matchImagePixels(input=glob.glob(threedhst.options['ALIGN_IMAGE']),
                     matchImage=root_direct+'_drz.fits',
                     output='align.fits', match_extension = 1)
                     
    #### Run SExtractor on the direct image, with the WHT 
    #### extension as a weight image
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CHECKIMAGE_TYPE'] = 'NONE'
    se.options['WEIGHT_TYPE']     = 'NONE'
    se.options['FILTER']    = 'Y'
    ## Detect thresholds (default = 1.5)
    se.options['DETECT_THRESH']    = '10' 
    se.options['ANALYSIS_THRESH']  = '10' 
    se.options['MAG_ZEROPOINT'] = str(threedhst.options['MAG_ZEROPOINT'])

    ## Run SExtractor on direct and alignment images
    se.options['CATALOG_NAME']    = 'direct.cat'
    iraf.imcopy(root_direct+'_drz.fits[1]',"SCI.fits")
    status = se.sextractImage('SCI.fits')
    
    se.options['CATALOG_NAME']    = 'align.cat'
    status = se.sextractImage('align.fits')
    
    directCat = threedhst.sex.mySexCat('direct.cat')
    alignCat = threedhst.sex.mySexCat('align.cat')
    
    #### Get x,y coordinates of detected objects
    directCat.x = directCat.columns[directCat.searchcol('X_IMAGE')].entry
    directCat.y = directCat.columns[directCat.searchcol('Y_IMAGE')].entry
    fp = open('direct.xy','w')
    for i in range(len(directCat.x)):
        fp.write('%s  %s\n' %(directCat.x[i],directCat.y[i]))
    fp.close()
    
    alignCat.x = alignCat.columns[alignCat.searchcol('X_IMAGE')].entry
    alignCat.y = alignCat.columns[alignCat.searchcol('Y_IMAGE')].entry
    fp = open('align.xy','w')
    for i in range(len(alignCat.x)):
        fp.write('%s  %s\n' %(alignCat.x[i],alignCat.y[i]))
    fp.close()
     
    #### use IRAF to compute shifts
    iraf.xyxymatch(input="direct.xy", reference="align.xy",
                   output="align.match",
                   tolerance=8, separation=0, verbose=yes)
    iraf.geomap(input="align.match", database="align.map",
                fitgeometry="shift", interactive=no)
    
    #### Parse geomap.output 
    fp = open("align.map","r")
    for line in fp.readlines():
        spl = line.split()
        if spl[0] == 'xshift':
            xshift = float(spl[1])    
        if spl[0] == 'yshift':
            yshift = float(spl[1])    
    fp.close()
    
    #### Cleanup
    rmfiles = ['SCI.fits','align.cat',
               'align.map','align.match','align.reg','align.xy',
               'direct.cat','direct.reg','direct.xy']
    for file in rmfiles:
        try:
            os.remove(file)
        except:
            pass
        
    return xshift, yshift
    
def matchImagePixels(input=None,matchImage=None,output=None,
                     match_extension=0):
    """
    matchImagePixels(input=None,matchImage=None,output=None,
                     match_extension=0)
    
    SWarp input image to same size/scale as matchIMage.
    Output is input_root+'.match.fits'
    
    """
    from threedhst.sex import SWarp
    
    #input = '/research/HST/GRISM/WEINER/ACS/h_nz_sect33_v2.0_drz_img.fits'
    #matchImage = 'IB3721050_SCI.fits'
    
    sw = threedhst.sex.SWarp()
    sw._aXeDefaults()
    sw.overwrite = True
    ### Get first guess coordinates
    sw.swarpMatchImage(matchImage,extension=match_extension)
    status = sw.swarpImage(matchImage+'[%d]' %match_extension,mode='wait')
    os.remove('coadd.fits')
    os.remove('coadd.weight.fits')
    ### Recenter
    sw.swarpRecenter()
    
    if not output:
        output = os.path.basename(input).split('.fits')[0]+'.match.fits'
    
    sw.options['IMAGEOUT_NAME'] = output
    #sw.options['WEIGHTOUT_NAME'] = base+'.match.weight.fits'
    status = sw.swarpImage(input,mode='direct')
    os.remove('coadd.weight.fits')


def make_grism_shiftfile(asn_direct, asn_grism):
    """make_grism_shiftfile(asn_direct, grism_direct)"""
    from threedhst.utils import ASNFile
    root_direct = asn_direct.split('_asn.fits')[0].lower()
    root_grism  =  asn_grism.split('_asn.fits')[0].lower()
    sf = ShiftFile(root_direct+'_shifts.txt')
    asn = ASNFile(asn_grism)
    for i,exp in enumerate(asn.exposures):
        sf.images[i] = exp+'_flt.fits'
    #for i in range(sf.nrows):
    #    sf.images[i] = asn.field('MEMNAME')[i].lower()+'_flt.fits'
        
    sf.print_shiftfile(root_grism+'_shifts.txt')
    print "3d-HST / make_grism_shiftfile: %s_shifts.txt" %root_grism


def checkShiftfile(asn_direct):
    """
    checkShiftfile(asn_direct)
    
    Make sure that there is a line in the shiftfile for each exposure 
    in the ASN table.
    """
    from threedhst.utils import ASNFile
    asn = ASNFile(asn_direct)
    
    sf_file = asn_direct.split('_asn.fits')[0]+'_shifts.txt'
    sf = ShiftFile(sf_file)
    for exp in asn.exposures:
        if exp+'_flt.fits' not in sf.images:
            raise NameError('Exposure, %s, not in %s' %(exp,sf_file))
    print '3DHST.shifts.checkShiftfile: %s looks OK.' %sf_file
    
class ShiftFile():
    """
    ShiftFile(infile)
    
    Based on aXe2html.sextractcat
    """
    def __init__(self, filename):
        linelist = self.opencat(filename)
        self.headerlines = self.extractheader(linelist)
        rowlines    = self.extractrows(linelist)
        self.processrows(rowlines)
        self.nrows = len(rowlines)
    
    def opencat(self, filename):
        """
        Input:
            filename - the name of the sextractor ascii catalog
            
        Return:
                linelist - a list wit all lines of the sextractor
                       catalog
                       
        Description:
            Opens the file and stores all the lines in a list.

        """
        listfile = open(filename,'r')
        linelist = listfile.readlines()
        return linelist
    
    
    def extractheader(self, linelist):
        """
    Input:
        linelist - all lines of a sextractor catalog

    Return:
        headerlines - the lines which contain header
                      information

    Description:
        Extracts the header lines from the list of lines
        Header lines have a '#' as the first digits and are
        longer than two characters. This allows for small
        cosmetical changes to the original Sextractor tables
        """
        headerlines = []
        for index in range(len(linelist)):
            oneline = linelist[index]
            if (oneline[0] == "#") and (len(oneline)) > 2:
                headerlines.append(oneline)
        return headerlines
    
    
    def extractrows(self, linelist):
        """
    Input:
        linelist - all lines of a sextractor catalog
        
    Return:
        rowlines - the content lines of a sextractor
                   catalog

    Description:
        Extracts the content lines from the list of lines.
        Content lines have a all lines except the ones which start with a 
        '#' and are longer than one character. This allows for
        e.g. whitespaces as the last 'rows', which often remain 
        after editing an ascii-file
        
        """
        rowlines = []
        for index in range(len(linelist)):
            oneline = linelist[index]
            if oneline[0] != "#" and len(oneline) > 1:
                rowlines.append(oneline)        
        return rowlines
    
    
    def processrows(self, rowlines):
        """
        processrows(self, rowlines)
        
        Read: image xshift yshift rotate scale from shiftfile.
            
        """
        nlines = len(rowlines)
        self.images = []
        self.xshift = []
        self.yshift = []
        self.rotate = []
        self.scale = []
        for i in range(nlines):
            line = rowlines[i].split()
            self.images.append(line[0])
            self.xshift.append(float(line[1]))
            self.yshift.append(float(line[2]))
            if len(line) > 3:
                self.rotate.append(float(line[3]))
            else:
                self.rotate.append(0.)
                self.scale.append(1.0)
            if len(line) > 4:
                self.scale.append(float(line[4]))
    
    
    def print_shiftfile(self, outfile):
        """print_shiftfile(outfile)"""
        fp = open(outfile,'w')
        fp.writelines(self.headerlines)
        for i in range(self.nrows):
            line = '%-20s %8.4f %8.4f %8.3f %8.3f\n' %(self.images[i],
                self.xshift[i], self.yshift[i], self.rotate[i], self.scale[i])
            fp.write(line)
        fp.close()
    
    def pop(self,idx):
        """
pop(self,idx)
        """
        out = self.images.pop(idx) 
        out = self.xshift.pop(idx) 
        out = self.yshift.pop(idx) 
        out = self.rotate.pop(idx) 
        out = self.scale.pop(idx) 
        self.nrows -= 1
    
    def append(self,image, xshift=0., yxhift=0.,
                    rotate=0.0, scale=1.0):
        """
append(self,image, xshift=0., yxhift=0.,
                   rotate=0.0, scale=1.0)
        """
        self.images.extend([image])
        self.xshift.extend([xshift])
        self.yshift.extend([yshift])
        self.rotate.extend([rotate])
        self.scale.extend([scale])
        self.nrows += 1