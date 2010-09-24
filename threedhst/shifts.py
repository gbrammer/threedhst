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

def run_tweakshifts(asn_direct):
    """
run_tweakshifts(asn_direct)
    
    asn_direct - filename of ASN table of direct images [...]_asn.fits
    
    This routine only uses dither.tweakshifts to compute the relative shifts of 
    the direct images
    """
    
    root = asn_direct.split('_asn.fits')[0].lower()

    try:
        os.remove(root+'_tweak.fits')
    except:
        pass        

    iraf.flpr()
    iraf.flpr()
    iraf.flpr()
    status = iraf.tweakshifts(input=asn_direct, shiftfile='',
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
       nsigma = 1.5, fitgeometry = 'rxyscale', function = 'polynomial', \
       maxiter = 3, reject = 3.0, crossref = '', margin = 50, tapersz = 50, \
       pad = no, fwhm = 7.0, ellip = 0.05, pa = 45.0, fitbox = 7, \
    Stdout=1)
    
def find_align_images_that_overlap(ROOT_DIRECT, ALIGN_IMAGE):
    """
align_img_list = find_align_images_that_overlap()
    
    Look through the images defined in threedhst.options['ALIGN_IMAGE'] and 
    return a list of only those that overlap.  
    
    This was written for the GOODS fields, where the multiple GOODS ACS tiles
    can be used to align the F140W images.  ALIGN_IMAGE will be something like 
    h_nz_sect*, but you don't want to waste time swarping large images that 
    don't overlap with the target image.
    """
    import glob
    
    #ROOT_DIRECT = threedhst.currentRun['ROOT_DIRECT']
    align_images = glob.glob(ALIGN_IMAGE)
    
    #### Get polygon of the direct mosaic edges
    px, py = threedhst.regions.wcs_polygon(ROOT_DIRECT+'_drz.fits', extension=1)
    
    #### Loop through align_images and check if they overlap with the 
    #### direct mosaic
    align_img_list = []
    for align_image in align_images:
        qx, qy = threedhst.regions.wcs_polygon(align_image, extension=0)
        if threedhst.regions.polygons_intersect(px, py, qx, qy):
            align_img_list.append(align_image)
    
    return align_img_list
    
def align_to_reference(ROOT_DIRECT, ALIGN_IMAGE, fitgeometry="shift"):
    """
xshift, yshift = align_to_reference()
    """        
    import os
    import glob
    import shutil
    
    #### Clean slate    
    rmfiles = ['SCI.fits','align.cat',
               'align.map','align.match','align.reg','align.xy', 
               ROOT_DIRECT+'_align.fits',
               'direct.cat','direct.reg','direct.xy']
    
    for file in rmfiles:
        try:
            os.remove(file)
        except:
            pass
                
    #### Get only images that overlap from the ALIGN_IMAGE list    
    align_img_list = find_align_images_that_overlap(ROOT_DIRECT, ALIGN_IMAGE)
    if not align_img_list:
        print 'threedhst.shifts.align_to_reference: no alignment images overlap.'
        return 0,0
    
    #### Use swarp to combine the alignment images to the same image 
    #### dimensions as the direct mosaic
    matchImagePixels(input=align_img_list,
                     matchImage=ROOT_DIRECT+'_drz.fits',
                     output=ROOT_DIRECT+'_align.fits', match_extension = 1)
                     
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

    #### Run SExtractor on direct and alignment images
    ## direct image
    se.options['CATALOG_NAME']    = 'direct.cat'
    iraf.imcopy(ROOT_DIRECT+'_drz.fits[1]',"SCI.fits")
    status = se.sextractImage('SCI.fits')
    
    ## alignment image
    se.options['CATALOG_NAME']    = 'align.cat'
    status = se.sextractImage(ROOT_DIRECT+'_align.fits')
    
    ## Read the catalogs
    directCat = threedhst.sex.mySexCat('direct.cat')
    alignCat = threedhst.sex.mySexCat('align.cat')
    
    #### Get x,y coordinates of detected objects
    ## direct image
    fp = open('direct.xy','w')
    for i in range(len(directCat.X_IMAGE)):
        fp.write('%s  %s\n' %(directCat.X_IMAGE[i],directCat.Y_IMAGE[i]))
    fp.close()
    
    ## alignment image
    fp = open('align.xy','w')
    for i in range(len(alignCat.X_IMAGE)):
        fp.write('%s  %s\n' %(alignCat.X_IMAGE[i],alignCat.Y_IMAGE[i]))
    fp.close()
     
    #### iraf.xyxymatch to find matches between the two catalogs
    status1 = iraf.xyxymatch(input="direct.xy", reference="align.xy",
                   output="align.match",
                   tolerance=8, separation=0, verbose=yes, Stdout=1)
    
    #### Compute shifts with iraf.geomap
    status2 = iraf.geomap(input="align.match", database="align.map",
                fitgeometry=fitgeometry, interactive=no, Stdout=1)
    
    #fp = open(root+'.iraf.log','a')
    #fp.writelines(status1)
    #fp.writelines(status2)
    #fp.close()
    
    xshift = 0
    yshift = 0
    rot = 0
    scale = 0
    
    #### Parse geomap.output 
    fp = open("align.map","r")
    for line in fp.readlines():
        spl = line.split()
        if spl[0].startswith('xshift'):
            xshift = float(spl[1])    
        if spl[0].startswith('yshift'):
            yshift = float(spl[1])    
        if spl[0].startswith('xrotation'):
            rot = float(spl[1])    
        if spl[0].startswith('xmag'):
            scale = float(spl[1])    
        
    fp.close()
    
    shutil.copy('align.map',ROOT_DIRECT+'_align.map')
    
    #### Cleanup
    rmfiles = ['SCI.fits','align.cat',
               'align.map','align.match','align.reg','align.xy',
               'direct.cat','direct.reg','direct.xy']
    
    for file in rmfiles:
        try:
            os.remove(file)
        except:
            pass
        
    return xshift, yshift, rot, scale
    
def matchImagePixels(input=None,matchImage=None,output=None,
                     match_extension=0):
    """
matchImagePixels(input=None,matchImage=None,output=None,
                 match_extension=0)
    
    SWarp input image to same size/scale as matchIMage.
    Output default is input_root+'.match.fits'
    
    """
    from threedhst.sex import SWarp
    
    #input = '/research/HST/GRISM/WEINER/ACS/h_nz_sect33_v2.0_drz_img.fits'
    #matchImage = 'IB3721050_SCI.fits'
    
    if not input:
        return False
    
    if not output:
        output = os.path.basename(input).split('.fits')[0]+'.match.fits'
    
    #### initialize SWarp
    sw = threedhst.sex.SWarp()
    sw._aXeDefaults()
    sw.overwrite = True
    
    #### Get first guess coordinates
    sw.swarpMatchImage(matchImage,extension=match_extension)
    status = sw.swarpImage(matchImage+'[%d]' %match_extension,mode='wait')
    os.remove('coadd.fits')
    os.remove('coadd.weight.fits')
    
    #### Recenter
    sw.swarpRecenter()
    
    #### Make the final output image
    sw.options['IMAGEOUT_NAME'] = output
    #sw.options['WEIGHTOUT_NAME'] = base+'.match.weight.fits'
    status = sw.swarpImage(input,mode='direct')
    os.remove('coadd.weight.fits')


def make_grism_shiftfile(asn_direct, asn_grism):
    """
make_grism_shiftfile(asn_direct, grism_direct)
    
    Make a shiftfile for grism exposures to match
    corresponding direct images
    """
    from threedhst.utils import ASNFile
    ROOT_DIRECT = asn_direct.split('_asn.fits')[0].lower()
    ROOT_GRISM  =  asn_grism.split('_asn.fits')[0].lower()
    #### Read shiftfile and ASN table
    sf = ShiftFile(ROOT_DIRECT+'_shifts.txt')
    asn = ASNFile(asn_grism)
    
    #### Change the image names in the shiftfile to the grism exposures
    for i,exp in enumerate(asn.exposures):
        sf.images[i] = exp+'_flt.fits'
    
    #### Write the new shiftfile
    sf.print_shiftfile(ROOT_GRISM+'_shifts.txt')
    
    print "\n3DHST.shifts.make_grism_shiftfile: %s_shifts.txt\n" %ROOT_GRISM


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
    print "\n3DHST.shifts.checkShiftfile: %s looks OK.\n" %sf_file
    
class ShiftFile():
    """
ShiftFile(infile)
    
    Read and manipulate a shiftfile produced by, e.g., tweakshifts
    
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