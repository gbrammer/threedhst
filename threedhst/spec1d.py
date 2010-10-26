"""
3DHST.spec1d

Process a 1-D spectrum and search for emission/absorption lines.
"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import matplotlib.pyplot as pyplot
import numpy as np
import numpy.random as nprand
import threedhst

### Random number seed so results are repeatable!!!
SEED = 67

def estimateRedshift(lines):
    """
estimateRedshift(lines)
    
    `lines` is a list of `spLineNew` objects, found with `findLines`.
    
    IDEA: for starters, just make a plot showing where all of the other lines
    would be if an observed line were associated with a particular line from a
    predefined linelist.
    """
    lines = readLinelist()
    idHa = np.where(np.array(lines.species) == 'H a')[0]
    ido3 = np.where(np.array(lines.species) == 'O III')[0]
    ido2 = np.where(np.array(lines.species) == 'O II')[0]
    
    
def readLinelist():
    """
lines = readLinelist()

    Read SDSS-like linelist from *threedhst/data/linelist.txt*
    """
    linelist = threedhst.utils.get_package_data("linelist.txt")
    lines = linelist.split('\n')
    species = lineList()
    for line in lines:
        if not line.startswith("#"):
            spl = line.split()
            species.append(lineSpecies(spl[0], spl[1], spl[2], 
            ' '.join(spl[3:])))
    
    return species
    
class lineList():
    """
lineList()
    """
    def __init__(self):
        """
__init__()
        """
        self.wave = threedhst.utils.listArray([])
        self.gal_weight = threedhst.utils.listArray([])
        self.qso_weight = threedhst.utils.listArray([])
        self.species = threedhst.utils.listArray([])
        self.lines = []
    
    def append(self, line):
        """
append(lineSpecies)
        """
        if not isinstance(line, lineSpecies):
            print 'Input line is not a `lineSpecies` object'
            return None
        
        self.wave.append(line.wave)
        self.gal_weight.append(line.gal_weight)
        self.qso_weight.append(line.qso_weight)
        self.species.append(line.species)
        self.lines.append(line)
    
    def pop(self, idx):
        """
pop(idx)
        """
        val = self.wave.pop(idx)
        val = self.gal_weight.pop(idx)
        val = self.qso_weight.pop(idx)
        val = self.species.pop(idx)
        val = self.lines.pop(idx)
        

class lineSpecies():
    """
lineSpecies(wave=0., gal_weight=0, qso_weight=0, species="")

    e.g. halpha = lineSpecies(6563.,8,8,"H a")
    """
    def __init__(self, wave=0., gal_weight=0, qso_weight=0, species=""):
        self.wave = np.float(wave)
        self.gal_weight = np.float(gal_weight)
        self.qso_weight = np.float(qso_weight)
        self.species = species
    
    def show(self):
        print '%7.1f %d %d %s' %(self.wave, self.gal_weight, self.qso_weight, 
                                self.species)
                                
def findLines(SPCFile, idx=195, show=False, verbose=False):
    """
lines = findLines(SPCFile, idx=195, show=False, verbose=False)
    
    Find emission lines and filter 0th order contamination that looks like
    emission lines.
    """
    
    lines  = spWFindLines(SPCFile, idx=idx, show=show, check_contam=False)
    
    NLINE = len(lines)
    
    if NLINE == 0:
        return None
    
    if verbose:
        print '\nBefore: \n'
        for line in lines:
            print line.wave, line.type, line.flag
    
    waves = np.zeros(NLINE)
    for i in range(NLINE): waves[i] = lines[i].wave
    
    #### find duplicates and remove
    kill = np.zeros(NLINE)
    # minimum separation between lines
    dup_toler = 100
    for i in range(NLINE-1):
        for j in range(i+1,NLINE):
            if np.abs(waves[j]-waves[i]) < dup_toler:
                kill[j] = 1
    for i in reversed(range(NLINE)):
        if kill[i] > 0:
            out = lines.pop(i)
            NLINE-=1
    
    #### Flag absorption lines near emission lines that are 
    #### probably artifacts of the line-finder
    abs_toler = 750
    for i in range(NLINE):
        if lines[i].type == 'abs':
            for j in range(NLINE):
                if (i != j) & (np.abs(lines[i].wave-lines[j].wave) < abs_toler):
                    lines[i].flag = 'artifact'
    
    ## Check if emission lines are in contamination spectrum
    contam = spWFindLines(SPCFile, idx=idx, check_contam=True, show=False)
    NCONTAM = len(contam)
    if NCONTAM > 0:
        if verbose:
            print '\nContam: \n'
            for line in contam:
                print line.wave, line.type, line.flag
        contam_toler = 70
        for i in range(NLINE):
            for j in range(NCONTAM):
                if ( (contam[j].type == 'em') & 
                     (np.abs(contam[j].wave-lines[i].wave) < contam_toler) ):
                     lines[i].flag = 'contam'
    
    if verbose:
        print '\nAfter: \n'
        for line in lines:
            print line.wave, line.type, line.flag
    
    return lines

  
def spWFindLines(SPCFile, idx=195, show=True, check_contam=False):
    """
lines = spWFindLines(SPCFile, idx=195, show=True, check_contam=False)
    
    >>> print lines[0].wave, lines[1].scale, lines[0].type
    
    Find emission lines with a wavelet transform. Translated from the SDSS DR7
    algorithm.
    
    """
    
    debug = 0
    iplot = 0
    
    ### input parameters
    wavemin = 1.1e4
    wavemax = 1.68e4
    nfilt = 3
    gthresh = 5
    wthresh = 3
    sigmaMin = 1.0
    ewMin = 2.0
    
    gthresh = 1.5
    wthresh = 1.5
    
    spec = SPCFile.getSpec(idx)
    lam  = spec.field('LAMBDA')*1.
    ok = np.where( (lam > wavemin) & (lam < wavemax) )
    
    flux = (spec.field('FLUX')*1.)[ok]
    ferr = (spec.field('FERROR')*1.)[ok]
    contam = (spec.field('CONTAM')*1.)[ok]
    lam = (lam*1.)[ok]
    npix = len(ok)
    
    corr = flux-contam
            
    if show: 
        pl = pyplot.plot(lam,corr)
    
    mask = np.array([1.0, 4.0, 6.0, 4.0, 1.0])
    
    ## Allocate space ##
    npix = len(lam)
    x = np.zeros(npix)
    t = np.zeros(npix)
    value = np.zeros(npix)
    origspec = np.zeros(npix)
    err = np.zeros(npix)
    ssmooth = np.zeros(npix)
    esmooth = np.zeros(npix)
    wave = np.zeros(npix)
    gauss = np.zeros(npix)
    gsmooth = np.zeros(npix)
    
    #### Search for "lines" in the contamination image
    if check_contam:
        nprand.seed(SEED)
        corr = contam+nprand.normal(0.,1.,(npix))*ferr/20.
        ferr /= 20.
    
    xmin = lam[0]
    xmax = lam[-1]

    x = lam*1.
    value = corr*1.
    origspec = corr*1.
    err = ferr*1.
    nprand.seed(SEED)
    gauss = nprand.normal(0.,1.,(npix))
    smin = np.min(value)
    smax = np.max(value)
    ok = np.where((err > 0) & (lam > wavemin) & (lam < wavemax))
    sn = np.median(value/err)
    
    if (debug > 0):
        print ("spWFindLines: mean s/n = %f\n" %sn)
    
    emLines = []
    
    ## Loop over smoothing scales ##
    scale = 1
    ## Smooth at this scale ##
    for ifilt in range(nfilt):
        gmean = 0.0
        grms = 0.0
        wmean = 0.0
        wrms = 0.0
        wmin = 0
        wmax = 0
        for i in range(npix):
            ssum = 0
            esum = 0
            gsum = 0
            for j in range(5):
                k = i + scale*(j-2)
                if (k < 0): k = 0
                if (k > npix-1): k = npix-1
                ssum += value[k]*mask[j]
                esum += err[k]*mask[j]
                gsum += gauss[k]*mask[j]
            
            ssmooth[i] = ssum/16.0
            esmooth[i] = esum/16.0
            gsmooth[i] = gsum/16.0
            
            ## Wavelet at this scale is simply difference between previous 
            ## and current smoothed spectra
            wave[i] = value[i] - ssmooth[i]
            gwave = gauss[i] - gsmooth[i]
            gmean +=  gwave
            grms +=  gwave*gwave
            wmean +=  wave[i]
            wrms +=  wave[i]*wave[i]
            if (wave[i] < wmin): wmin = wave[i]
            if (wave[i] > wmax): wmax = wave[i]
        
        gmean /= npix
        grms = np.sqrt(grms/npix - gmean*gmean)
        wmean /= npix
        wrms = np.sqrt(wrms/npix - wmean*wmean)
        if (debug > 0):
            print ("mean, rms gaussian noise = %f %f\n" %(gmean, grms))
            print ("mean, rms wavelet signal = %f %f\n" %(wmean, wrms))
        
        ## Copy smoothed spec etc and set threshold ##
        value = ssmooth*1.
        err = esmooth*1.
        gauss = gsmooth*1.
        t = gthresh*grms*esmooth
        idx = np.where(wthresh*wrms > t)
        if (idx[0].shape[0] > 0):
            t[idx] = wthresh*wrms
        
        # #for (i = 0 i < npix i++) 
        # for i in range(npix):
        #     value[i] = ssmooth[i]
        #     err[i] = esmooth[i]
        #     gauss[i] = gsmooth[i]
        #     t[i] = gthresh*grms*esmooth[i]
        #     if (wthresh*wrms > t[i]): t[i] = wthresh*wrms
        
        ## Plot current wavelet ##
        if (iplot > 0):
            cpgsvp(xvmin, xvmax, yvmin + (ifilt+1)*dy, yvmin + (ifilt+2)*dy)
            cpgswin(xmin, xmax, wmin, wmax)
            cpgbox("BCST", 0.0, 0, "BCNSTV", 0.0, 0)
            cpgline(npix, x, wave)
            cpgsci(3)
            cpgline(npix, x, t)
            cpgsci(1)
            #sprintf(label, "\\gs = %d", (int) pow(2,ifilt))
            cpgmtxt("T", -1.5, 0.85, 0.0, label)
        
        if show: 
            opl = pyplot.plot(x,wave)
            pyplot.plot(x,t,color=opl[0]._color)
            
        # Do the line finding. Locate pixels above local threshold (t[i]). We
        # require that adjacent pixels are also above threshold and that pixel
        # higher than surrounding four. We could also require that the central
        # pixel has mask value SP_MASK_OK, but this condition excludes too many
        # genuine lines.
        #for (i = 2 i < npix-2 ++i) {
        for i in range(2,npix-2):
            ### Has to be within range
            testRange = (x[i] > wavemin) & (x[i] < wavemax)
            ### Emission
            testPos = ( (wave[i] > t[i]) & (wave[i-1] > t[i-1]) & 
                        (wave[i+1] > t[i+1]) & 
                        (wave[i] >= wave[i-1]) & (wave[i] >= wave[i+1]) & 
                        (wave[i] >= wave[i-2]) & (wave[i] >= wave[i+2]) )
            
            ### Absorption
            nwave = -1*wave            
            testNeg = ( (nwave[i] > t[i]) & (nwave[i-1] > t[i-1]) & 
                        (nwave[i+1] > t[i+1]) & 
                        (nwave[i] >= nwave[i-1]) & (nwave[i] >= nwave[i+1]) & 
                        (nwave[i] >= nwave[i-2]) & (nwave[i] >= nwave[i+2]) )
            
            # if ( (x[i] > wavemin) & (x[i] < wavemax) & 
            #     (wave[i] > t[i]) & (wave[i-1] > t[i-1]) & (wave[i+1] > t[i+1]) & 
            #     (wave[i] >= wave[i-1]) & (wave[i] >= wave[i+1]) & 
            #     (wave[i] >= wave[i-2]) & (wave[i] >= wave[i+2]) ):
                # & spLineFound(spec->foundLines, x[i]) == 0
            if (testRange & (testPos | testNeg)):
                
                ## Find zero-crossings to estimate line width ##
                ilo = i
                while( (ilo > 0) & (wave[ilo] > 0.0) ):
                    ilo -= 1
                
                ihi = i
                while( (ihi < npix-2) & (wave[ihi] > 0.0) ):
                    ihi += 1
                    
                ## Create a new line ##
                line = spLineNew()
                line.wave = x[i]
                line.waveMin = x[ilo]
                line.waveMax = x[ihi]
                line.height = corr[i] - ssmooth[i]
                line.sigma = 0.75*(x[ihi] - x[ilo])
                line.sigmaMin = sigmaMin
                line.ew = 0
                line.ewMin = ewMin
                line.continuum = 0
                if (t[i] > 0):
                  line.nsigma = wave[i]/t[i]
                else:
                  line.nsigma = 0.0
                
                line.restWave = 0
                if (debug > 0):
                    print ("line centre %f, limits %f %f, height %f\n" 
                              %(x[i], x[ilo], x[ihi], line.height))
                
                if show: 
                    pyplot.plot(line.wave*np.array([1,1]),np.array([-1,1]),
                        color=opl[0]._color)
                
                line.scale = scale
                if testPos:
                    line.type = 'em'
                else:
                    line.type = 'abs'
                    
                emLines.append(line)
        
        scale *= 2
    
    if show:
        pyplot.xlim(wavemin,wavemax)
        pyplot.ylim(-1.e-19,5.e-19)
    
    # for line in emLines:
    #     print line.wave, line.scale, line.sigma, line.type
        
    return emLines
    
class spLineNew:
    wave = 0.
    waveMin = 0.
    waveMax = 0.
    height = 0.
    sigma = 0.
    sigmaMin = 0.
    ew = 0.
    ewMin = 0.
    continuum = 0.
    nsigma = 0.
    sn = 0.
    restwave = 0.
    scale = 0.
    type = ''
    flag = 'ok'