"""
3DHST.spec1d

Process a 1-D spectrum and search for emission/absorption lines.
"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import numpy as np

import matplotlib.pyplot as pyplot
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
            print('Input line is not a `lineSpecies` object')
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
        print('%7.1f %d %d %s' %(self.wave, self.gal_weight, self.qso_weight, 
                                self.species))
                                
def findLines(SPCFile, idx=195, show=False, verbose=False, trim_abs=False):
    """
lines = findLines(SPCFile, idx=195, show=False, verbose=False, trim_abs=False)
    
    Find emission lines and filter 0th order contamination that looks like
    emission lines.
    
    If `trim_abs` is set, only return emission lines.
    """
    
    lines  = spWFindLines(SPCFile, idx=idx, show=show, check_contam=False)
    
    NLINE = len(lines)
    
    if NLINE == 0:
        return None
    
    if verbose:
        print('\nBefore: \n')
        for line in lines:
            print(line.wave, line.type, line.flag)
    
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
    for i in reversed(list(range(NLINE))):
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
            print('\nContam: \n')
            for line in contam:
                print(line.wave, line.type, line.flag)
        contam_toler = 150
        for i in range(NLINE):
            for j in range(NCONTAM):
                if ( (contam[j].type == 'em') & 
                     (np.abs(contam[j].wave-lines[i].wave) < contam_toler) ):
                     lines[i].flag = 'contam'
    
    if trim_abs:
        lsave = []
        for line in lines:
            if line.type == 'em':
                lsave.append(line)
        lines = lsave
        
    if verbose:
        print('\nAfter: \n')
        for line in lines:
            print(line.wave, line.type, line.flag)
    
    return lines

  
def spWFindLines(SPCFile, idx=195, show=True, check_contam=False):
    """
lines = spWFindLines(SPCFile, idx=195, show=True, check_contam=False)
    
    >>> print lines[0].wave, lines[1].scale, lines[0].type
    
    Find emission lines with a wavelet transform. Translated from the SDSS DR7
    algorithm.
    
    """
    import threedhst.plotting
    
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
    
    if isinstance(SPCFile, threedhst.plotting.SPCFile):
        spec = SPCFile.getSpec(idx)
        lam  = spec.field('LAMBDA')*1.
        flux = spec.field('FLUX')*1.
        ferr = spec.field('FERROR')*1.
        contam = spec.field('CONTAM')*1.
    else:
        #### Alternative is an ascii spectra file
        lam = SPCFile.lam*1.
        flux = SPCFile.flux*1.
        ferr = SPCFile.error*1.
        contam = SPCFile.contam*1.
        
    ok = (lam > wavemin) & (lam < wavemax) 
    if ok.sum() < 3:
        return []
        
    flux = flux[ok]
    ferr = ferr[ok]
    contam = contam[ok]
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
        print(("spWFindLines: mean s/n = %f\n" %sn))
    
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
            print(("mean, rms gaussian noise = %f %f\n" %(gmean, grms)))
            print(("mean, rms wavelet signal = %f %f\n" %(wmean, wrms)))
        
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
                line.height = corr[i] #- ssmooth[i]
                #line.sigma = 0.75*(x[ihi] - x[ilo])
                line.sigmaMin = sigmaMin
                line.ew = 0
                line.ewMin = ewMin
                line.continuum = np.median(corr)
                line.sn = (line.height-line.continuum) / np.median(ferr)
                line.fcontam = contam[i]/flux[i]
                
                if (t[i] > 0):
                  line.nsigma = wave[i]/t[i]
                else:
                  line.nsigma = 0.0
                
                line.restWave = 0
                if (debug > 0):
                    print(("line centre %f, limits %f %f, height %f\n" 
                              %(x[i], x[ilo], x[ihi], line.height)))
                
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
    
class readLinesDat():
    def __init__(self, infile):
        wave = []
        sigma = []
        eqw = []
        sn = []
        id = []
        ndet = []
        
        fp = open(infile,'r')
        lines = fp.readlines()
        fp.close()
        
        for line in lines:
            if not line.startswith('#'):
                spl = line.split()
                N = len(spl)
                if N > 1:
                    for i in range(N/4):
                        id.append(spl[0])
                        wave.append(spl[i*4+1])
                        sigma.append(spl[i*4+2])
                        eqw.append(spl[i*4+3])
                        sn.append(spl[i*4+4])
                        ndet.append(N/4)
                        
        self.id = np.cast[int](id)
        self.ndet = np.cast[int](ndet)
        self.wave = np.cast[float](wave)
        self.sigma = np.cast[float](sigma)
        self.eqw = np.cast[float](eqw)
        self.sn = np.cast[float](sn)
        
def test1D():
    import threedhst
    
    root='AEGIS_F140W'
    path='/research/HST/GRISM/3DHST/AEGIS/HTML/'
    ID=418
    
    threedhst.spec1d.extract1D(ID, root=root, path=path)
    
def extract1D(ID, root='orient1', path='../HTML', show=False, out2d=False):
    """ 
    Own extraction of 1D spectra.
    
    Errors are estimated directly from the 2D extractions.
    
    Include background estimate.
    
    """
    import matplotlib.pyplot as plt

    import threedhst
    import threedhst.catIO as catIO
    #from threedhst.dq import myDS9
    
    from scipy import polyfit, polyval
    
    twod = pyfits.open('%s/images/%s_%05d_2D.fits.gz' %(path, root, ID))
    
    head = twod[1].header
    data = twod[1].data
    cont = twod[4].data
    model = twod[5].data
    weight = twod[7].data
    weight = weight*0+1
    twod.close()
    
    flux_limit = model.max()*5.e-3    
    ### Mask where the 2D flux is less than 1% of the maximum
    mask = (model < flux_limit) & (cont < flux_limit) & (data != 0)
    mask_model = (model < flux_limit) & (data != 0)
    #mask = mask & (data < 0.5*model.max())
    
    # ma = model
    # ma[cont > model] = cont[cont > model]
    # flux_limit = ma.max()*1.e-2    
    # mask = (ma < flux_limit) 
    # mask = model < flux_limit
    
    wht = data*1.
    wht[mask == False] = 0
    
    ### S/N
    rms = threedhst.utils.biweight(wht[mask])
    #ds9.view(data/rms)
    
    ### 1D profile
    contam_mask = (cont > 1.e-3) & (data != 0)
    tt = data*1.
    tt[contam_mask] = 0
    N = data*0+1
    N[contam_mask] = 0
    
    N *= weight
    
    profile = np.sum(tt*N, axis=1) / np.sum(N, axis=1)
    profile[np.sum(N, axis=1) == 0] = 0.
    profile /= np.sum(profile)
    
    model_profile = np.sum(model, axis=1)
    model_profile /= np.sum(model_profile)
    
    #### Extract where within 0.2 of the peak of the profile
    xpix = np.arange(profile.size)
    
    use = xpix[profile > 0.2*profile.max()]
    use = xpix[model_profile > 0.2*model_profile.max()]
    if use.size == 0:
        use = xpix
        
    # xmax = list(xpix[profile == profile.max()])[0]
    # use = [xpix[xmax]]
    # i=1
    # while (profile[xmax+i] >= 0.2*profile.max()):
    #     use.append(xmax+i)
    #     i+=1
    # 
    # i= -1
    # while (profile[xmax+i] >= 0.2*profile.max()):
    #     use.append(xmax+i)
    #     i-=1
    # 
    # use.sort()
    # use = np.array(use)

    yi, xi = np.indices(data.shape)
    profile_mask = (yi >= use.min()) & (yi <= use.max())
    masked = data*1
    masked[profile_mask == False] = 0
    #masked[mask_model == True] = 0
    
    masked_cont = cont*1
    masked_cont[profile_mask == False] = 0
    #masked_cont[mask_model == True] = 0
    
    background = data*1
    background[mask == False] = 0
    Nbg = data*0
    Nbg[mask == True] = 1
    
    #### Do the sum
    N = data*0
    N[profile_mask] = 1
    #N[mask_model == True] = 0
    
    GAIN = 1 # WFC3 FLT images already in e-
    sha = data.shape
    twod_model_profile = np.dot(model_profile.reshape(sha[0],-1), np.ones((1,sha[1]))) * (data != 0)
    #N *= twod_model_profile #> 0# weight > 0
    
    oned_dn = np.sum(masked, axis=0) / np.sum(N, axis=0)
    oned_dn_cont = np.sum(masked_cont*N, axis=0) / np.sum(N, axis=0)
    poisson_var = oned_dn*head['EXPTIME']*GAIN ## GAIN = 2.5 e- / DN
    poisson_var[poisson_var < 0] = 0

    oned_dn_var = (rms*GAIN)**2 / np.sum(N, axis=0)*head['EXPTIME']**2 + poisson_var
    #oned_dn_var = rms**2 / np.sum(N, axis=0)*head['EXPTIME']**2 
    #oned_dn_var = rms**2 * head['EXPTIME']**2 + poisson_var
    oned_dn_err = np.sqrt(oned_dn_var)/head['EXPTIME']

    oned_dn_bg = np.sum(background, axis=0) / np.sum(Nbg, axis=0)
    
    #oned_dn -= oned_dn_bg
    
    lam = (np.arange(oned_dn.size)-head['CRPIX1']+1)*head['CDELT1']+head['CRVAL1']
    
    #### Fit a polygon to the background
    xfit = np.arange(oned_dn_bg.size*1.)/oned_dn_bg.size
    xfit_use = (lam > 1.08e4) & (lam < 1.68e4) & (np.isfinite(oned_dn_bg)) & (oned_dn_bg != 0)
    oned_dn_bg_fit = oned_dn*0.
    if xfit[xfit_use].size > 0:
        #print xfit_use.shape
        #print lam[xfit_use].min()
        if (lam[xfit_use].min() < 1.1e4) & (lam[xfit_use].max() > 1.6e4):
            polycoeffs = polyfit(xfit[xfit_use], oned_dn_bg[xfit_use], 5)
            oned_dn_bg_fit = polyval(polycoeffs, xfit)

    #print 'Fit ORDER: %d' %(ORDER)
    
    ### sensitivity
    if threedhst.options['GRISM_NAME'] == 'G800L':
        sens = pyfits.open('../CONF/ACS.WFC.1st.sens.7.fits')[1].data

    if threedhst.options['GRISM_NAME'] == 'G141':
        sens = pyfits.open('../CONF/WFC3.IR.G141.1st.sens.2.fits')[1].data
    
    yint = np.interp(lam, sens.WAVELENGTH, sens.SENSITIVITY)
    yint_err = np.interp(lam, sens.WAVELENGTH, sens.ERROR)
    
    ### Final fluxed spectra
    oned_flux = oned_dn * GAIN / yint
    oned_flux_cont = oned_dn_cont * GAIN / yint
    oned_flux_var = (oned_dn_err / yint)**2 + (oned_flux * yint_err / yint)**2
    oned_flux_err = np.sqrt(oned_flux_var)
    oned_flux_bg = oned_dn_bg * GAIN / yint
    oned_flux_bg_fit = oned_dn_bg_fit * GAIN / yint
    
    #### Compare to aXe
    lint = 1.3e4
    
    #### Show results
    if show:
        spec = catIO.ReadASCIICat('%s/ascii/%s_%05d.dat' %(path, root, ID))    
        spec.lam = spec.field('lam')
        
        fig = plt.figure(figsize=[6,4.1],dpi=100)
        fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.16,
                            bottom=0.15,right=0.98,top=0.98)
        ax = fig.add_subplot(111)

        plt.plot(lam, oned_flux / np.interp(lint, lam, oned_flux), color='red', linestyle='-')
        plt.plot(lam, (oned_flux - oned_flux_bg_fit) / np.interp(lint, lam, oned_flux), color='orange', linestyle='-')
        plt.plot(lam, (oned_flux_bg) / np.interp(lint, lam, oned_flux), color='green', linestyle='-')
        plt.plot(lam, (oned_flux_bg_fit) / np.interp(lint, lam, oned_flux), color='green', linestyle='-', alpha=0.5)
        plt.plot(lam, oned_flux_cont / np.interp(lint, lam, oned_flux), color='red', alpha=0.2, linestyle='-')
        plt.plot(lam, oned_flux_err / np.interp(lint, lam, oned_flux), color='red')

        plt.plot(spec.lam, spec.flux / np.interp(lint, spec.lam, spec.flux), color='blue')
        plt.plot(spec.lam, spec.contam / np.interp(lint, spec.lam, spec.flux), color='blue', alpha=0.2)
        plt.plot(spec.lam, spec.error / np.interp(lint, spec.lam, spec.flux), color='blue')

        setlabel()

        # print 'Norm: %6.3f' %(np.interp(lint, lam, oned_flux)/np.interp(lint, spec.lam, spec.flux))
        #print lam.size
        
    if out2d:
        ## return the masked 2D thumbnails for the background and object extractions
        return background, masked, model, cont, data
    else:
        ## return a structure that looks like the 1D spec FITS files
        out = {}
        out['lam'] = lam
        out['flux'] = oned_flux
        out['error'] = oned_flux_err
        out['contam'] = oned_flux_cont
        out['background'] = oned_flux_bg_fit
        return out

def show_extract1D():
    import os
    
    import threedhst
    import matplotlib.pyplot as plt
    
    os.chdir('/research/HST/GRISM/3DHST/COSMOS/DATA')
    
    ID=110
    root='COSMOS-26-G141'
    background, masked, model, cont, data = threedhst.spec1d.extract1D(ID, root=root, path='../HTML/', out2d=True)
    
    plt.rcParams['font.family'] = 'serif'
    plt.rcParams['font.serif'] = ['Times']
    plt.rcParams['font.size'] = 9
    
    fig = plt.figure(figsize=[4.7,7.2],dpi=100)
    fig.subplots_adjust(wspace=0.02,hspace=0.02,left=0.02,
                        bottom=0.01,right=0.98,top=0.97)
    
    vmi, vma = -0.02, 0.07
    inter = 'nearest'
    
    ax = fig.add_subplot(5,1,1)
    ax.imshow(data, vmin=vmi, vmax=vma, interpolation=inter)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    plt.title('2D spectrum')

    ax = fig.add_subplot(5,1,2)
    ax.imshow(cont, vmin=vmi, vmax=vma, interpolation=inter)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    plt.title('Contamination model')

    ax = fig.add_subplot(5,1,3)
    ax.imshow(model, vmin=vmi, vmax=vma, interpolation=inter)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    plt.title('Object model')

    ax = fig.add_subplot(5,1,4)
    ax.imshow(masked-cont, vmin=vmi, vmax=vma, interpolation=inter)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    plt.title('Extraction region')
    
    ax = fig.add_subplot(5,1,5)
    ax.imshow(background, vmin=vmi, vmax=vma, interpolation=inter)
    ax.set_yticklabels([])
    ax.set_xticklabels([])
    plt.title('Background')
    
    plt.savefig('/tmp/example_extract1D.pdf')
    
    SPC = threedhst.plotting.SPCFile(root+'_2_opt.SPC.fits',
                    axe_drizzle_dir=os.environ['AXE_DRIZZLE_PATH'])
    
    threedhst.plotting.plot1Dspec(SPC, ID, outfile='/tmp/own_extraction.pdf',
                  close_window=True, show_test_lines=False, own_extraction=True)
    #
    threedhst.plotting.plot1Dspec(SPC, ID, outfile='/tmp/axe_extraction.pdf',
                  close_window=True, show_test_lines=False, own_extraction=False)
    
def setlabel():
    import matplotlib.pyplot as plt

    plt.xlim(1.e4,1.7e4)
    ymax = 2.5
    plt.ylim(-0.3*ymax, ymax)
    plt.xlabel(r'$\lambda$')
    plt.ylabel(r'$f_\lambda$')
