"""
3DHST.compute_shifts

Utilities for computing shifts and processing shiftfiles.
"""
#   $URL$
#   $Rev$
#   $Author$
#   $Date$

import os
import pyfits
from pyraf import iraf
no = iraf.no
yes = iraf.yes

def compute_shifts(asn_direct):
    """compute_shifts(asn_direct)"""
    root = asn_direct.split('_asn.fits')[0].lower()
    ### Get shifts with tweakshifts
    iraf.tweakshifts ( input = asn_direct, shiftfile = '', reference = root+'_tweak.fits', \
       output = root+'_shifts.txt', findmode = 'catalog', gencatalog = 'daofind', \
       sextractpars = '', undistort = yes, computesig = yes, idckey = 'idctab', \
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
    ### !!! Need to add steps to match WCS to astrometric references (ACS)
    
def make_grism_shiftfile(asn_direct, asn_grism):
    """make_grism_shiftfile(asn_direct, grism_direct)"""
    direct_root = asn_direct.split('_asn.fits')[0].lower()
    grism_root  =  asn_grism.split('_asn.fits')[0].lower()
    sf = ShiftFile(direct_root+'_shifts.txt')
    asn = read_asn(asn_grism)
    for i in range(sf.nrows):
        sf.images[i] = asn.field('MEMNAME')[i].lower()+'_flt.fits'
        
    sf.print_shiftfile(grism_root+'_shifts.txt')
    print "3d-HST / make_grism_shiftfile: %s_shifts.txt" %grism_root

    
class ShiftFile():
    """ShiftFile(infile)
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
        """processrows(self, rowlines)
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
            line = "%s %f %f %f %f\n" %(self.images[i], self.xshift[i], self.yshift[i], self.rotate[i], self.scale[i])
            fp.write(line)
        fp.close()
