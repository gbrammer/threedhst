"""
eazyPy: routines for reading and plotting Eazy output

    EazyParam
    readEazyBinary
    getEazySED
    getEazyPz
    plotExampleSED
    nMAD
    zPhot_zSpec
    
"""

import os
import string

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc
import pylab

import threedhst.catIO as catIO

class FilterDefinition:
    def __init__(self, bp=None):
        """
        Placeholder for the filter definition information.
        """
        self.name = None
        self.wavelength = None
        self.transmission = None
        
        if bp is not None:
            self.wavelength = np.cast[np.double](bp.wave)
            self.transmission =  np.cast[np.double](bp.throughput)
            self.name = bp.name
            
    def extinction_correction(self, EBV, Rv=3.1, mag=True, source_lam=None, source_flux = None):
        """
        Get the MW extinction correction within the filter.  
        
        Optionally supply a source spectrum.
        """
        
        if self.wavelength is None:
            print 'Filter not defined.'
            return False
        
        if source_flux is None:
            source_flux = self.transmission*0.+1
        else:
            source_flux = np.interp(self.wavelength, source_lam, source_flux, left=0, right=0)
            
        Av = EBV*Rv
        Alambda = milkyway_extinction(lamb = self.wavelength, Rv=Rv)
        delta = np.trapz(self.transmission*source_flux*10**(-0.4*Alambda*Av), self.wavelength) / np.trapz(self.transmission*source_flux, self.wavelength)
        
        if mag:
            return 2.5*np.log10(delta)
        else:
            return 1./delta
    
    def ABVega(self):
        """
        Compute AB-Vega conversion
        """
        try:
            import pysynphot as S
        except:
            print 'Failed to import "pysynphot"'
            return False
        
        vega=S.FileSpectrum(S.locations.VegaFile)
        abmag=S.FlatSpectrum(0,fluxunits='abmag')
        #xy, yy = np.loadtxt('hawki_y_ETC.dat', unpack=True)
        bp = S.ArrayBandpass(wave=self.wavelength, throughput=self.transmission, name='')
        ovega = S.Observation(vega, bp)
        oab = S.Observation(abmag, bp)
        return -2.5*np.log10(ovega.integrate()/oab.integrate())
    
    def pivot(self):
        """
        PySynphot pivot wavelength
        """
        try:
            import pysynphot as S
        except:
            print 'Failed to import "pysynphot"'
            return False
            
        self.bp = S.ArrayBandpass(wave=self.wavelength, throughput=self.transmission, name='')
        return self.bp.pivot()
        
    def rectwidth(self):
        """
        Synphot filter rectangular width
        """
        try:
            import pysynphot as S
        except:
            print 'Failed to import "pysynphot"'
            return False
            
        self.bp = S.ArrayBandpass(wave=self.wavelength, throughput=self.transmission, name='')
        return self.bp.rectwidth()

    #
    def ctw95(self):
        """
        95% cumulative throughput width
        http://www.stsci.edu/hst/acs/analysis/bandwidths/#keywords
        
        """
        
        dl = np.diff(self.wavelength)
        filt = np.cumsum((self.wavelength*self.transmission)[1:]*dl)
        ctw95 = np.interp([0.025, 0.975], filt/filt.max(), self.wavelength[1:])
        return np.diff(ctw95)
            
        
class FilterFile:
    def __init__(self, file='FILTER.RES.v8.R300'):
        """
        Read a EAZY (HYPERZ) filter file.
        """
        fp = open(file)
        lines = fp.readlines()
        fp.close()
        
        filters = []
        wave = []
        for line in lines:
            if 'lambda_c' in line:
                if len(wave) > 0:
                    new_filter = FilterDefinition()
                    new_filter.name = header
                    new_filter.wavelength = np.cast[float](wave)
                    new_filter.transmission = np.cast[float](trans)
                    filters.append(new_filter)
                    
                header = ' '.join(line.split()[1:])
                wave = []
                trans = []
            else:
                lspl = np.cast[float](line.split())
                wave.append(lspl[1])
                trans.append(lspl[2])
        # last one
        new_filter = FilterDefinition()
        new_filter.name = header
        new_filter.wavelength = np.cast[float](wave)
        new_filter.transmission = np.cast[float](trans)
        filters.append(new_filter)
           
        self.filters = filters
        self.NFILT = len(filters)
    
    def names(self, verbose=True):
        """
        Print the filter names.
        """
        if verbose:
            for i in range(len(self.filters)):
                print '%5d %s' %(i+1, self.filters[i].name)            
        else:
            string_list = ['%5d %s\n' %(i+1, self.filters[i].name) for i in range(len(self.filters))]
            return string_list
            
    def write(self, file='xxx.res', verbose=True):
        """
        Dump the filter information to a filter file.
        """
        fp = open(file,'w')
        for filter in self.filters:
            fp.write('%6d %s\n' %(len(filter.wavelength), filter.name))
            for i in range(len(filter.wavelength)):
                fp.write('%-6d %.5e %.5e\n' %(i+1, filter.wavelength[i], filter.transmission[i]))
        
        fp.close()
        
        string_list = self.names(verbose=False)
        fp = open(file+'.info', 'w')
        fp.writelines(string_list)
        fp.close()
        
        if verbose:
            print 'Wrote <%s[.info]>' %(file)
            
    def search(self, search_string, case=True, verbose=True):
        """ 
        Search filter names for `search_string`.  If `case` is True, then
        match case.
        """
        import re
        
        if not case:
            search_string = search_string.upper()
        
        matched = []
        
        for i in range(len(self.filters)):
            filt_name = self.filters[i].name
            if not case:
                filt_name = filt_name.upper()
                
            if re.search(search_string, filt_name) is not None:
                if verbose:
                    print '%5d %s' %(i+1, self.filters[i].name)
                matched.append(i)
        
        return np.array(matched)
        
class ParamFilter(FilterDefinition):
    def __init__(self, line='#  Filter #20, RES#78: COSMOS/SUBARU_filter_B.txt - lambda_c=4458.276253'):
        
        self.lambda_c = float(line.split('lambda_c=')[1])
        self.name = line.split()[4]
        self.fnumber = int(line.split('RES#')[1].split(':')[0])
        self.cnumber = int(line.split('Filter #')[1].split(',')[0])
        
class EazyParam():
    """
    Read an Eazy zphot.param file.
    
    Example: 
    
    >>> params = EazyParam(PARAM_FILE='zphot.param')
    >>> params['Z_STEP']
    '0.010'

    """    
    def __init__(self, PARAM_FILE='zphot.param', READ_FILTERS=False):
        self.filename = PARAM_FILE
        self.param_path = os.path.dirname(PARAM_FILE)
        
        f = open(PARAM_FILE,'r')
        self.lines = f.readlines()
        f.close()
        
        self._process_params()
        
        filters = []
        templates = []
        for line in self.lines:
            if line.startswith('#  Filter'):
                filters.append(ParamFilter(line))
            if line.startswith('#  Template'):
                templates.append(line.split()[3])
                
        self.NFILT = len(filters)
        self.filters = filters
        self.templates = templates
        
        if READ_FILTERS:
            RES = FilterFile(self.params['FILTERS_RES'])
            for i in range(self.NFILT):
                filters[i].wavelength = RES.filters[filters[i].fnumber-1].wavelength
                filters[i].transmission = RES.filters[filters[i].fnumber-1].transmission
                
    def _process_params(self):
        params = {}
        formats = {}
        self.param_names = []
        for line in self.lines:
            if line.startswith('#') is False:
                lsplit = line.split()
                if lsplit.__len__() >= 2:
                    params[lsplit[0]] = lsplit[1]
                    self.param_names.append(lsplit[0])
                    try:
                        flt = float(lsplit[1])
                        formats[lsplit[0]] = 'f'
                        params[lsplit[0]] = flt
                    except:
                        formats[lsplit[0]] = 's'
                    
        self.params = params
        #self.param_names = params.keys()
        self.formats = formats
    
    def show_filters(self):
        for filter in self.filters:
            print ' F%d, %s, lc=%f' %(filter.fnumber, filter.name, filter.lambda_c)
    
    def write(self, file=None):
        if file == None:
            print 'No output file specified...'
        else:
            fp = open(file,'w')
            for param in self.param_names:
                if isinstance(self.params[param], np.str):
                    fp.write('%-25s %s\n' %(param, self.params[param]))
                else:
                    fp.write('%-25s %f\n' %(param, self.params[param]))
                    #str = '%-25s %'+self.formats[param]+'\n'
            #
            fp.close()
            
    #
    def __getitem__(self, param_name):
        """
    __getitem__(param_name)

        >>> cat = mySexCat('drz.cat')
        >>> print cat['NUMBER']

        """
        if param_name not in self.param_names:
            print ('Column %s not found.  Check `column_names` attribute.'
                    %column_name)
            return None
        else:
            #str = 'out = self.%s*1' %column_name
            #exec(str)
            return self.params[param_name]
    
    def __setitem__(self, param_name, value):
        self.params[param_name] = value

class TranslateFile():
    def __init__(self, file='zphot.translate'):
        self.file=file
        self.ordered_keys = []
        lines = open(file).readlines()
        self.trans = {}
        self.error = {}
        for line in lines:
            spl = line.split()
            key = spl[0]
            self.ordered_keys.append(key)
            self.trans[key] = spl[1]
            if len(spl) == 3:
                self.error[key] = float(spl[2])
            else:
                self.error[key] = 1.
            #
            
    def change_error(self, filter=88, value=1.e8):
        
        if isinstance(filter, str):
            if 'f_' in filter:
                err_filt = filter.replace('f_','e_')
            else:
                err_filt = 'e'+filter

            if err_filt in self.ordered_keys:
                self.error[err_filt] = value
                return True
        
        if isinstance(filter, int):
            for key in self.trans.keys():
                if self.trans[key] == 'E%0d' %(filter):
                    self.error[key] = value
                    return True
        
        print 'Filter %s not found in list.' %(str(filter))
    
    def write(self, file=None, show_ones=False):

        lines = []
        for key in self.ordered_keys:
            line = '%s  %s' %(key, self.trans[key])
            if self.trans[key].startswith('E') & ((self.error[key] != 1.0) | show_ones):
                line += '  %.1f' %(self.error[key])

            lines.append(line+'\n')

        if file is None:
            file = self.file
        
        if file:
            fp = open(file,'w')
            fp.writelines(lines)
            fp.close()
        else:
            for line in lines:
                print line[:-1]

def readRFBinary(file='OUTPUT/test.153-155.coeff'):
    """
    Read Rest-frame coefficients file
    """
    f = open(file, 'rb')
    NOBJ, NFILT, NTEMP = np.fromfile(file=f,dtype=np.int32, count=3)
    rftempfilt = np.fromfile(file=f,dtype=np.double,count=NFILT*NTEMP).reshape((NTEMP,NFILT)).transpose()
    rfcoeff = np.fromfile(file=f,dtype=np.double,count=NOBJ*NTEMP).reshape((NOBJ, NTEMP)).transpose()
    if NFILT == 1: 
        rftempfilt = rftempfilt.flatten()
    
    f.close()
    
    d = {'NOBJ':NOBJ, 'NFILT':NFILT, 'NTEMP':NTEMP, 'tempfilt':rftempfilt, 'coeffs':rfcoeff}
    return d
    
def readEazyBinary(MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same'):
    """
tempfilt, coeffs, temp_sed, pz = readEazyBinary(MAIN_OUTPUT_FILE='photz', \
                                                OUTPUT_DIRECTORY='./OUTPUT', \
                                                CACHE_FILE = 'Same')

    Read Eazy BINARY_OUTPUTS files into structure data.
    
    If the BINARY_OUTPUTS files are not in './OUTPUT', provide either a relative or absolute path
    in the OUTPUT_DIRECTORY keyword.
    
    By default assumes that CACHE_FILE is MAIN_OUTPUT_FILE+'.tempfilt'.
    Specify the full filename if otherwise. 
    """
    
    #root='COSMOS/OUTPUT/cat3.4_default_lines_zp33sspNoU'
    
    root = OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE
    
    ###### .tempfilt
    if CACHE_FILE == 'Same':
        CACHE_FILE = root+'.tempfilt'
    
    if os.path.exists(CACHE_FILE) is False:
        print ('File, %s, not found.' %(CACHE_FILE))
        return -1,-1,-1,-1
    
    f = open(CACHE_FILE,'rb')
    
    s = np.fromfile(file=f,dtype=np.int32, count=4)
    NFILT=s[0]
    NTEMP=s[1]
    NZ=s[2]
    NOBJ=s[3]
    tempfilt = np.fromfile(file=f,dtype=np.double,count=NFILT*NTEMP*NZ).reshape((NZ,NTEMP,NFILT)).transpose()
    lc = np.fromfile(file=f,dtype=np.double,count=NFILT)
    zgrid = np.fromfile(file=f,dtype=np.double,count=NZ)
    fnu = np.fromfile(file=f,dtype=np.double,count=NFILT*NOBJ).reshape((NOBJ,NFILT)).transpose()
    efnu = np.fromfile(file=f,dtype=np.double,count=NFILT*NOBJ).reshape((NOBJ,NFILT)).transpose()
    
    f.close()
    
    tempfilt  = {'NFILT':NFILT,'NTEMP':NTEMP,'NZ':NZ,'NOBJ':NOBJ,\
                 'tempfilt':tempfilt,'lc':lc,'zgrid':zgrid,'fnu':fnu,'efnu':efnu}
    
    ###### .coeff
    f = open(root+'.coeff','rb')
    
    s = np.fromfile(file=f,dtype=np.int32, count=4)
    NFILT=s[0]
    NTEMP=s[1]
    NZ=s[2]
    NOBJ=s[3]
    coeffs = np.fromfile(file=f,dtype=np.double,count=NTEMP*NOBJ).reshape((NOBJ,NTEMP)).transpose()
    izbest = np.fromfile(file=f,dtype=np.int32,count=NOBJ)
    tnorm = np.fromfile(file=f,dtype=np.double,count=NTEMP)
    
    f.close()
    
    coeffs = {'NFILT':NFILT,'NTEMP':NTEMP,'NZ':NZ,'NOBJ':NOBJ,\
              'coeffs':coeffs,'izbest':izbest,'tnorm':tnorm}
              
    ###### .temp_sed
    f = open(root+'.temp_sed','rb')
    s = np.fromfile(file=f,dtype=np.int32, count=3)
    NTEMP=s[0]
    NTEMPL=s[1]
    NZ=s[2]
    templam = np.fromfile(file=f,dtype=np.double,count=NTEMPL)
    temp_seds = np.fromfile(file=f,dtype=np.double,count=NTEMPL*NTEMP).reshape((NTEMP,NTEMPL)).transpose()
    da = np.fromfile(file=f,dtype=np.double,count=NZ)
    db = np.fromfile(file=f,dtype=np.double,count=NZ)
    
    f.close()
    
    temp_sed = {'NTEMP':NTEMP,'NTEMPL':NTEMPL,'NZ':NZ,\
              'templam':templam,'temp_seds':temp_seds,'da':da,'db':db}
              
    ###### .pz
    if os.path.exists(root+'.pz'):
        f = open(root+'.pz','rb')
        s = np.fromfile(file=f,dtype=np.int32, count=2)
        NZ=s[0]
        NOBJ=s[1]
        chi2fit = np.fromfile(file=f,dtype=np.double,count=NZ*NOBJ).reshape((NOBJ,NZ)).transpose()

        ### This will break if APPLY_PRIOR No
        s = np.fromfile(file=f,dtype=np.int32, count=1)
        
        if len(s) > 0:
            NK = s[0]
            kbins = np.fromfile(file=f,dtype=np.double,count=NK)
            priorzk = np.fromfile(file=f, dtype=np.double, count=NZ*NK).reshape((NK,NZ)).transpose()
            kidx = np.fromfile(file=f,dtype=np.int32,count=NOBJ)
            pz = {'NZ':NZ,'NOBJ':NOBJ,'NK':NK, 'chi2fit':chi2fit, 'kbins':kbins, 'priorzk':priorzk,'kidx':kidx}
        else:
            pz = None
        
        f.close()
        
    else:
        pz = None
    
    if False:
        f = open(root+'.zbin','rb')
        s = np.fromfile(file=f,dtype=np.int32, count=1)
        NOBJ=s[0]
        z_a = np.fromfile(file=f,dtype=np.double,count=NOBJ)
        z_p = np.fromfile(file=f,dtype=np.double,count=NOBJ)
        z_m1 = np.fromfile(file=f,dtype=np.double,count=NOBJ)
        z_m2 = np.fromfile(file=f,dtype=np.double,count=NOBJ)
        z_peak = np.fromfile(file=f,dtype=np.double,count=NOBJ)
        f.close()
        
    ###### Done.    
    return tempfilt, coeffs, temp_sed, pz

        
def getEazySED(idx, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same', scale_flambda=True, verbose=False, individual_templates=False):
    """
lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
     getEazySED(idx, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same')
    
    Get best-fit Eazy template for object number 'idx' from the specified Eazy output files. 

    Output variables are as follows:
        
        lambdaz: full best-fit template (observed) wavelength, interpolated at WAVELENGTH_GRID
        temp_sed:          "        "              flux (F_lambda)
        lci: filter pivot wavelengths
        fobs: observed fluxes, including zeropoint offsets if used, F_lambda
        efobs: observed flux errors,    "            "        "        "
    """
    tempfilt, coeffs, temp_seds, pz = readEazyBinary(MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, CACHE_FILE = CACHE_FILE)
    
    ##### Apply zeropoint factors
    param = EazyParam(PARAM_FILE=OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.param')
    fnumbers = np.zeros(len(param.filters), dtype=np.int)
    for i in range(len(fnumbers)):
        fnumbers[i] = int(param.filters[i].fnumber)
    
    zpfile = OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zeropoint'
    if os.path.exists(zpfile):
        zpfilts, zpf_file = np.loadtxt(zpfile, unpack=True, dtype=np.str)                                    
        zpf = np.ones(tempfilt['NFILT'])
        for i in range(len(zpfilts)):
            match = fnumbers == int(zpfilts[i][1:])
            zpf[match] = np.float(zpf_file[i])
    else:
        zpf = np.ones(tempfilt['NFILT'])

    zpfactors = np.dot(zpf.reshape(tempfilt['NFILT'],1),\
                       np.ones(tempfilt['NOBJ']).reshape(1,tempfilt['NOBJ']))

    if verbose:
        print zpf
        
    tempfilt['fnu'] *= zpfactors
    tempfilt['efnu'] *= zpfactors
    
    lci = tempfilt['lc'].copy()
    
    params = EazyParam(PARAM_FILE=OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.param')
    abzp = np.float(params['PRIOR_ABZP'])
        
    # fobs = tempfilt['fnu'][:,idx]/(lci/5500.)**2*flam_factor
    # efobs = tempfilt['efnu'][:,idx]/(lci/5500.)**2*flam_factor
    ### Physical f_lambda fluxes, 10**-17 ergs / s / cm2 / A
    if scale_flambda:
        flam_factor = 10**(-0.4*(params['PRIOR_ABZP']+48.6))*3.e18/1.e-17
    else:
        flam_factor = 5500.**2
    
    missing = (tempfilt['fnu'][:,idx] < -99) | (tempfilt['efnu'][:,idx] < 0)
    fobs = tempfilt['fnu'][:,idx]/lci**2*flam_factor
    efobs = tempfilt['efnu'][:,idx]/lci**2*flam_factor
    fobs[missing] = -99
    efobs[missing] = -99
    #print lci, tempfilt['fnu'][:,idx], tempfilt['efnu'][:,idx]
    
    ##### Broad-band SED
    obs_sed = np.dot(tempfilt['tempfilt'][:,:,coeffs['izbest'][idx]],\
                     coeffs['coeffs'][:,idx])/(lci)**2*flam_factor
    
    zi = tempfilt['zgrid'][coeffs['izbest'][idx]]
    
    ###### Full template SED, observed frame
    lambdaz = temp_seds['templam']*(1+zi)
    temp_sed = np.dot(temp_seds['temp_seds'],coeffs['coeffs'][:,idx])
    if individual_templates:
        temp_sed = temp_seds['temp_seds']*coeffs['coeffs'][:,idx]
    
    temp_sed /= (1+zi)**2
    
    temp_sed *= (1/5500.)**2*flam_factor
    
    ###### IGM absorption
    lim1 = np.where(temp_seds['templam'] < 912)
    lim2 = np.where((temp_seds['templam'] >= 912) & (temp_seds['templam'] < 1026))
    lim3 = np.where((temp_seds['templam'] >= 1026) & (temp_seds['templam'] < 1216))
    
    if lim1[0].size > 0: temp_sed[lim1] *= 0.
    if lim2[0].size > 0: temp_sed[lim2] *= 1.-temp_seds['db'][coeffs['izbest'][idx]]
    if lim3[0].size > 0: temp_sed[lim3] *= 1.-temp_seds['da'][coeffs['izbest'][idx]]
        
    ###### Done
    return lambdaz, temp_sed, lci, obs_sed, fobs, efobs

def getAllPz(MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same'):
    """
    Return a matrix with *all* normalized p(z) for a given catalog
    """
    import threedhst.eazyPy as eazy
    
    tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, CACHE_FILE = CACHE_FILE)
    
    full_pz = np.zeros((pz['NZ'], pz['NOBJ']))
    for i in xrange(pz['NOBJ']):
        #print i
        zz, pzi = eazy.getEazyPz(i, MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, CACHE_FILE=CACHE_FILE, binaries=(tempfilt, pz))
        full_pz[:,i] = pzi
    
    return zz, full_pz
            
def getEazyPz(idx, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same', binaries=None, get_prior=False):
    """
zgrid, pz = getEazyPz(idx, \
                      MAIN_OUTPUT_FILE='photz', \
                      OUTPUT_DIRECTORY='./OUTPUT', \
                      CACHE_FILE='Same', binaries=None)
                      
    Get Eazy p(z) for object #idx.
    
    To avoid re-reading the binary files, supply binaries = (tempfilt, pz)
    
    """
    if binaries is None:
        tempfilt, coeffs, temp_seds, pz = readEazyBinary(MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                                                    OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                                                    CACHE_FILE = CACHE_FILE)
    else:
        tempfilt, pz = binaries
        
    if pz is None:
        return None, None
    
    ###### Get p(z|m) from prior grid
    kidx = pz['kidx'][idx]
    #print kidx, pz['priorzk'].shape
    if (kidx > 0) & (kidx < pz['priorzk'].shape[1]):
        prior = pz['priorzk'][:,kidx]
    else:
        prior = np.ones(pz['NZ'])
        
    ###### Convert Chi2 to p(z)
    pzi = np.exp(-0.5*(pz['chi2fit'][:,idx]-min(pz['chi2fit'][:,idx])))*prior#*(1+tempfilt['zgrid'])
    
    if np.sum(pzi) > 0:
        pzi/=np.trapz(pzi, tempfilt['zgrid'])
    
    ###### Done
    if get_prior:
        return tempfilt['zgrid'], pzi, p
    else:
        return tempfilt['zgrid'], pzi
        
class TemplateError():
    """
    Make an easy (spline) interpolator for the template error function
    """
    def __init__(self, file='templates/TEMPLATE_ERROR.eazy_v1.0'):
        from scipy import interpolate
        self.te_x, self.te_y = np.loadtxt(file, unpack=True)
        self._spline = interpolate.InterpolatedUnivariateSpline(self.te_x, self.te_y)
        
    def interpolate(self, filter_wavelength, z):
        """
        observed_wavelength is observed wavelength of photometric filters.  But 
        these sample the *rest* wavelength of the template error function at lam/(1+z)
        """
        return self._spline(filter_wavelength/(1+z))
        
class Template():
    def __init__(self, sp=None, file=None):
        self.wavelength = None
        self.flux = None
        self.flux_fnu = None
        
        if sp is not None:
            self.wavelength = np.cast[np.double](sp.wave)
            self.flux = np.cast[np.double](sp.flux)
            self.flux_fnu = self.flux
            
        if file is not None:
            self.wavelength, self.flux = np.loadtxt(file, unpack=True)
            self.set_fnu()
    
    def set_fnu(self):
        self.flux_fnu = self.flux * (self.wavelength/5500.)**2
        
    def integrate_filter(self, filter, z=0):
        """
        Integrate the template through a `FilterDefinition` filter object.
        
        "unicorn" is part of a yet non-public module.  np.interp can 
        stand in for the meantime.
        """
        try:
            import unicorn
            temp_filter = unicorn.utils_c.interp_conserve_c(filter.wavelength, 
                                     self.wavelength*(1+z), self.flux_fnu)
        
        except ImportError:
            temp_filter = np.interp(filter.wavelength, self.wavelength*(1+z),
                                    self.flux_fnu)
            
        temp_int = np.trapz(filter.transmission*temp_filter/filter.wavelength, filter.wavelength) / np.trapz(filter.transmission/filter.wavelength, filter.wavelength)
        #temp_int = np.trapz(filter.transmission*temp_filter, filter.wavelength) / np.trapz(filter.transmission, 1./filter.wavelength)
        return temp_int
        
class TemplateInterpolator():
    """
    Class to use scipy spline interpolator to interpolate pre-computed eazy template 
    photometry at arbitrary redshift(s).
    """
    def __init__(self, bands=None, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same', zout=None, f_lambda=True):
        from scipy import interpolate
        import threedhst.eazyPy as eazy
        
        #### Read the files from the specified output
        tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, CACHE_FILE = CACHE_FILE)
        
        if bands is None:
            self.bands = np.arange(tempfilt['NFILT'])
        else:
            self.bands = np.array(bands)
        
        self.band_names = ['' for b in self.bands]
        
        if zout is not None:
            param = eazy.EazyParam(PARAM_FILE=zout.filename.replace('.zout','.param'))
            self.band_names = [f.name for f in param.filters]
            self.bands = np.array([f.fnumber-1 for f in param.filters])
                        
        self.NFILT = len(self.bands)
        self.NTEMP = tempfilt['NTEMP']
        self.lc = tempfilt['lc'][self.bands]
        self.sed = temp_seds
        self.templam = self.sed['templam']
        self.temp_seds = self.sed['temp_seds']
        
        # if True:
        #     import threedhst
        #     import unicorn
        #     threedhst.showMessage('Conroy model', warn=True)
        #     cvd12 = np.loadtxt(unicorn.GRISM_HOME+'/templates/cvd12_t11_solar_Chabrier.dat')
        #     self.temp_seds[:,0] = np.interp(self.templam, cvd12[:,0], cvd12[:,1])
        
        self.in_zgrid = tempfilt['zgrid']
        self.tempfilt = tempfilt['tempfilt'][self.bands, :, :]
        if f_lambda:
            for i in range(self.NFILT):
                self.tempfilt[i,:,:] /= (self.lc[i]/5500.)**2
                
        ###### IGM absorption
        self.igm_wave = []
        self.igm_wave.append(self.templam < 912)
        self.igm_wave.append((self.templam >= 912) & (self.templam < 1026))
        self.igm_wave.append((self.templam >= 1026) & (self.templam < 1216))
        
        self._spline_da = interpolate.InterpolatedUnivariateSpline(self.in_zgrid, temp_seds['da'])
        self._spline_db = interpolate.InterpolatedUnivariateSpline(self.in_zgrid, temp_seds['db'])
        
        #### Make a 2D list of the spline interpolators
        self._interpolators = [range(self.NTEMP) for i in range(self.NFILT)]                
        for i in range(self.NFILT):
            for j in range(self.NTEMP):
                self._interpolators[i][j] = interpolate.InterpolatedUnivariateSpline(self.in_zgrid, self.tempfilt[i, j, :])
        #
        self.output = None
        self.zout = None
    
    def interpolate_photometry(self, zout):
        """
        Interpolate the EAZY template photometry at `zout`, which can be a number or an 
        array.
        
        The result is returned from the function and also stored in `self.output`.
        """               
        output = [range(self.NTEMP) for i in range(self.NFILT)]                
        for i in range(self.NFILT):
            for j in range(self.NTEMP):
                output[i][j] = self._interpolators[i][j](zout)
        
        self.zgrid = np.array(zout)
        self.output = np.array(output)
        return self.output
        
    def check_extrapolate(self):
        """
        Check if any interpolated values are extrapolated from the original redshift grid
        
        Result is both returned and stored in `self.extrapolated`
        """
        if self.zout is None:
            return False
        
        self.extrapolated = np.zeros(self.output.shape, dtype=np.bool) ## False

        bad = (self.zgrid < self.in_zgrid.min()) | (self.zgrid > self.in_zgrid.max())
        self.extrapolated[:, :, bad] = True
        
        return self.extrapolated
        
    def get_IGM(self, z, matrix=False, silent=False):
        """
        Retrieve the full SEDs with IGM absorption
        """
        ###### IGM absorption
        # lim1 = self.templam < 912
        # lim2 = (self.templam >= 912) & (self.templam < 1026)
        # lim3 = (self.templam >= 1026) & (self.templam < 1216)
        
        igm_factor = np.ones(self.templam.shape[0])
        igm_factor[self.igm_wave[0]] = 0.
        igm_factor[self.igm_wave[1]] = 1. - self._spline_db(z)
        igm_factor[self.igm_wave[1]] = 1. - self._spline_da(z)
        
        if matrix:
            self.igm_factor = np.dot(igm_factor.reshape(-1,1), np.ones((1, self.NTEMP)))
        else:
            self.igm_factor = igm_factor
            
        self.igm_z = z
        self.igm_lambda = self.templam*(1+z)
        
        if not silent:
            return self.igm_lambda, self.igm_factor
#
def interpolate_tempfilt_loop(tempfilt, zgrid, zi, output):
    """    
    Linear interpolate an Eazy "tempfilt" grid at z=zi.  
    
    `tempfilt` is [NFILT, NTEMP, NZ] integrated flux matrix
    `zgrid` is [NZ] redshift grid
    `output` is empty [NFILT, NTEMP] grid to speed up execution
    """
    sh = tempfilt.shape
    NF, NT, NZ = sh[0], sh[1], sh[2]
    #output = np.zeros((NF, NT))
    for iz in range(NZ-1):
        dz = zgrid[iz+1]-zgrid[iz]
        fint = 1 - (zi-zgrid[iz])/dz
        if (fint > 0) & (fint <= 1):
            fint2 = 1 - (zgrid[iz+1]-zi)/dz
            # print iz, zgrid[iz], fint, fint2
            for ifilt in range(NF):
                for itemp in range(NT):
                    #print ifilt, itemp
                    output[ifilt, itemp] = tempfilt[ifilt, itemp, iz]*fint + tempfilt[ifilt, itemp, iz+1]*fint2
            break              
    
    return output

try:
    from numba import double, jit
    faster_interpolate_tempfilt = jit(double[:,:](double[:,:,:], double[:], double, double[:,:]))(interpolate_tempfilt_loop)
except:
    faster_interpolate_tempfilt = interpolate_tempfilt_loop
    #pass
    
def convert_chi_to_pdf(tempfilt, pz):
    """
    Convert the `chi2fit` array in the `pz` structure to probability densities.
    
    The `tempfilt` structure is needed to get the redshift grid.
    """
    pdf = pz['chi2fit']*0.

    for idx in range(pz['NOBJ']):
        ###### Get p(z|m) from prior grid
        kidx = pz['kidx'][idx]
        #print kidx, pz['priorzk'].shape
        if (kidx > 0) & (kidx < pz['priorzk'].shape[1]):
            prior = pz['priorzk'][:,kidx]
        else:
            prior = np.ones(pz['NZ'])

        ###### Convert Chi2 to p(z)
        pzi = np.exp(-0.5*(pz['chi2fit'][:,idx]-min(pz['chi2fit'][:,idx])))*prior
        if np.sum(pzi) > 0:
            pzi/=np.trapz(pzi, tempfilt['zgrid'])
            pdf[:,idx] = pzi
            
    return tempfilt['zgrid']*1., pdf
    
def plotExampleSED(idx=20, writePNG=True, MAIN_OUTPUT_FILE = 'photz', OUTPUT_DIRECTORY = 'OUTPUT', CACHE_FILE = 'Same', lrange=[3000,8.e4], axes=None, individual_templates=False, fnu=False):
    """
PlotSEDExample(idx=20)

    Plot an example Eazy best-fit SED.
    """

    #zout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
    zout = catIO.Readfile(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
    #qz = np.where(zout.z_spec > 0)[0]
    print zout.filename
    qz = np.arange(len(zout.id))
    
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
        getEazySED(qz[idx], MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                          OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                          CACHE_FILE = CACHE_FILE, individual_templates=individual_templates, scale_flambda=True)
    
    zgrid, pz = getEazyPz(qz[idx], MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                                   OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                                   CACHE_FILE = CACHE_FILE)
    ##### plot defaults
    #rc('font',**{'family':'serif','serif':['Times']})
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.serif'] = ['Helvetica']
    #plt.rcParams['ps.useafm'] = True
    plt.rcParams['patch.linewidth'] = 0.
    plt.rcParams['patch.edgecolor'] = 'black'
    #plt.rcParams['text.usetex'] = True
    plt.rcParams['text.usetex'] = False
    #plt.rcParams['text.latex.preamble'] = ''

    ##### start plot
    if axes is None:
        fig = plt.figure(figsize=[8,4],dpi=100)
        fig.subplots_adjust(wspace=0.18, hspace=0.0,left=0.09,bottom=0.15,right=0.98,top=0.98)
    
    #### Plot parameters
    plotsize=35
    alph=0.9
    
    if fnu:
        temp_sed *= (lambdaz / 5500.)**2
        fobs *= (lci/5500.)**2
        efobs *= (lci/5500.)**2
        obs_sed *= (lci/5500.)**2

    #### Full best-fit template
    if axes is None:
        ax = fig.add_subplot(121)
        axp = fig.add_subplot(122)
    else:
        ax = axes[0]
        axp = axes[1]
        
    if individual_templates:
        ax.plot(lambdaz, temp_sed, linewidth=1.0, color='blue',alpha=0.4)
        ax.plot(lambdaz, temp_sed.sum(axis=1), linewidth=1.0, color='blue',alpha=alph)
    else:
        ax.plot(lambdaz, temp_sed, linewidth=1.5, color='blue',alpha=alph*0.8, zorder=-3)
    
    #### template fluxes integrated through the filters
    ax.scatter(lci, obs_sed,
               c='red',marker='o',s=plotsize,alpha=alph, zorder=-1)

    #### Observed fluxes w/ errors
    #ax.errorbar(lci,fobs,yerr=efobs,ecolor=None,
    #           color='black',fmt='o',alpha=alph)
    #
    # ax.errorbar(lci, fobs, yerr=efobs, ecolor='black',
    #            color='black',fmt='o',alpha=alph, markeredgecolor='black', markerfacecolor='None', markeredgewidth=1.5, ms=8, zorder=1)

    highsn = fobs/efobs > 2
    ax.errorbar(lci[highsn], fobs[highsn], yerr=efobs[highsn], ecolor='black',
               color='black',fmt='o',alpha=alph, markeredgecolor='black', markerfacecolor='None', markeredgewidth=1.5, ms=8, zorder=2)
    #
    ax.errorbar(lci[~highsn], fobs[~highsn], yerr=efobs[~highsn], ecolor='0.7',
               color='black',fmt='o',alpha=alph, markeredgecolor='0.7', markerfacecolor='None', markeredgewidth=1.5, ms=8, zorder=1)
    
    for i in range(len(lci)):
        print '%f %e %e %e' %(lci[i], obs_sed[i], fobs[i], efobs[i])
        
    #### Set axis range and titles
    ax.semilogx()
    ax.set_xlim(lrange[0],lrange[1])
    ax.set_ylim(-0.05*max(obs_sed),1.1*max(fobs))
    ax.set_xlabel(r'$\lambda$ [$\AA$]')
    ax.set_ylabel(r'$f_\lambda$')
    
    ##### P(z)
    if pz is not None:            
        axp.plot(zgrid, pz, linewidth=1.0, color='orange',alpha=alph)
        axp.fill_between(zgrid,pz,np.zeros(zgrid.size),color='yellow')

        if zout.z_spec[qz[idx]] > 0:
            axp.plot(zout.z_spec[qz[idx]]*np.ones(2), np.array([0,1e6]),color='red',alpha=0.4)

        #### Set axis range and titles
        axp.set_xlim(0,np.ceil(np.max(zgrid)))
        axp.set_ylim(0,1.1*max(pz))
        axp.set_xlabel(r'$z$')
        axp.set_ylabel(r'$p(z)$')
        
    if (writePNG is not False) & (axes is None):
        if isinstance(writePNG, str):
            out=writePNG
        else:
            out='/tmp/test.pdf'
            
        fig.savefig(out,dpi=100)

    if axes is None:
        return [ax, axp]
    

def nMAD(arr):
    """
result = nMAD(arr)

    Get the NMAD statistic of the input array, where
    NMAD = 1.48 * median(ABS(arr) - median(arr)).
    """
    med = np.median(arr)
    return 1.48*np.median(np.abs(arr-med))
    
def _add_ticks():
    """
    NAME:
       _add_ticks
    PURPOSE:
       add minor axis ticks to a plot
    INPUT:
       (none; works on the current axes)
    OUTPUT:
       (none; works on the current axes)
    HISTORY:
       2009-12-23 - Written - Bovy (NYU)
    """
    ax=plt.gca()
    xstep= ax.xaxis.get_majorticklocs()
    xstep= xstep[1]-xstep[0]
    ax.xaxis.set_minor_locator(ticker.MultipleLocator(xstep/2.))
    ystep= ax.yaxis.get_majorticklocs()
    ystep= ystep[1]-ystep[0]
    ax.yaxis.set_minor_locator(ticker.MultipleLocator(ystep/2.))

def zPhot_zSpec(zoutFile='../COSMOS/OUTPUT/cat3.4_default_lines_zp33sspNoU.zout', \
                zmax=4, marker='o', color='black', alpha=0.2, ticks=None):
    """
zPhot_zSpec(zoutfile="./OUTPUT/photz.zout', zmax=4)

    Make a nice zphot-zspec plot for an Eazy output file.
    """
    #zout = '../COSMOS/OUTPUT/cat3.4_default_lines_zp33sspNoU.zout'
    if os.path.exists(zoutFile) is False:
        print ('File, %s, not found.' %(zoutFile))
        return
        
    zout = catIO.Readfile(zoutFile)
    
    ##### plot defaults
    #rc('font',**{'family':'serif','serif':['Times']})
    #rc('font',**{'family':'sans-serif','sans-serif':['Helvetica']})    
    #plt.rcParams['ps.useafm'] = True
    plt.rcParams['font.family'] = 'sans-serif'
    plt.rcParams['font.serif'] = ['Times']    
    plt.rcParams['patch.linewidth'] = 0.25
    plt.rcParams['patch.edgecolor'] = 'black'
    plt.rcParams['text.usetex'] = False
    #plt.rcParams['text.latex.preamble'] = ''
    plt.rcParams['font.size'] = 12
    plt.rcParams['xtick.major.size'] = 6
    plt.rcParams['xtick.minor.size'] = 3
    plt.rcParams['ytick.major.size'] = plt.rcParams['xtick.major.size']
    plt.rcParams['ytick.minor.size'] = plt.rcParams['xtick.minor.size']

    ##### start plot
    fig = plt.figure(figsize=[5,5],dpi=100)
    fig.subplots_adjust(wspace=0.13,hspace=0.0,left=0.12,bottom=0.10,right=0.97,top=0.98)
    
    #### Plot parameters
    plotsize=20
    alph=alpha
    
    #### zphot - zspec
    qz = np.where(zout.z_spec > 0)
    ax = fig.add_subplot(111)
    ax.scatter(zout.z_spec[qz],zout.z_peak[qz],
               c=color, marker=marker,s=plotsize,alpha=alph)
    ax.plot(np.array([0,zmax]),np.array([0,zmax]),alpha=0.2,color='white',linestyle='-',linewidth=2)
    ax.plot(np.array([0,zmax]),np.array([0,zmax]),alpha=0.7,color='red',linestyle='-',linewidth=1)

    #### Add labels
    dz = (zout.z_peak[qz]-zout.z_spec[qz])/(1+zout.z_spec[qz])

    ax.text(0.1*zmax, 0.95*zmax, zoutFile.replace('_','\_'), \
            ha="left", va="center", size=6)

    ax.text(0.1*zmax, 0.9*zmax, r"$\sigma_\mathrm{NMAD} = %8.3f$" %(nMAD(dz)), \
            ha="left", va="center", size=16)
    
    ax.text(0.1*zmax, 0.82*zmax, r"$N = %0d$" %(len(dz)), \
            ha="left", va="center", size=16)
    
    #### Set axis range and titles
    ax.set_xlim(0,zmax)
    ax.set_ylim(0,zmax)
    if ticks is None:
        ticks = np.arange(0,zmax,1)
    
    ax.set_xticks(ticks)
    ax.set_yticks(ticks)
    
    _add_ticks()
    
    plt.xlabel(r"$z_\mathrm{spec}$")
    plt.ylabel(r"$z_\mathrm{peak}$")
    
    #fig.savefig('/tmp/test.pdf',dpi=100)
    
def show_fit_residuals(root='photz_v1.7.fullz', PATH='./OUTPUT/', savefig=None, adjust_zeropoints='zphot.zeropoint', fix_filter=None, ref_filter=None, get_resid=False, wclip=[1200, 3.e4]):
    """
    Plot the EAZY fit residuals to evaluate zeropoint updates
    """
    import threedhst
    import threedhst.eazyPy as eazy
    
    if not PATH.endswith('/'):
        PATH += '/'
    
    ##### Read the param file
    param = eazy.EazyParam('%s%s.param' %(PATH, root))
    
    ##### Read template fluxes and coefficients
    tempfilt, coeffs, temp_seds, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE=root, OUTPUT_DIRECTORY=PATH,CACHE_FILE = 'Same')
    
    if coeffs['izbest'].max() == 0:
        STAR_FIT = True
    else:
        STAR_FIT = False
        
    param = eazy.EazyParam(PARAM_FILE=PATH+'/'+root+'.param')
    fnumbers = np.zeros(len(param.filters), dtype=np.int)
    for i in range(len(fnumbers)):
        fnumbers[i] = int(param.filters[i].fnumber)
        
    zpfile = PATH+'/'+root+'.zeropoint'
    if os.path.exists(zpfile):
        zpfilts, zpf_file = np.loadtxt(zpfile, unpack=True, dtype=np.str)                                    
        zpf = np.ones(tempfilt['NFILT'])
        for i in range(len(zpfilts)):
            match = fnumbers == int(zpfilts[i][1:])
            zpf[match] = np.float(zpf_file[i])
    else:
        zpf = np.ones(tempfilt['NFILT'])
        
    zpfactors = np.dot(zpf.reshape(tempfilt['NFILT'],1), np.ones(tempfilt['NOBJ']).reshape(1,tempfilt['NOBJ']))
    
    ok = (tempfilt['fnu'] > -90) & (tempfilt['efnu'] > 0)
    tempfilt['fnu'][ok] *= zpfactors[ok]
    tempfilt['efnu'][ok] *= zpfactors[ok]
    
    obs_sed = np.zeros((tempfilt['NFILT'], tempfilt['NOBJ']), dtype=np.float)
    for i in xrange(tempfilt['NOBJ']):
        obs_sed[:,i] = np.dot(tempfilt['tempfilt'][:,:,coeffs['izbest'][i]], coeffs['coeffs'][:,i])
    
    zi = tempfilt['zgrid'][coeffs['izbest']]
    lc = tempfilt['lc']
    offsets = lc*0.
    
    lc_rest = np.dot(lc.reshape(-1,1), 1./(1+zi.reshape(1,-1)))
    
    resid = (obs_sed-tempfilt['fnu']) / obs_sed + 1
    signoise = tempfilt['fnu']/np.sqrt(tempfilt['efnu']**2+(0.01*tempfilt['fnu'])**2)
    
    if get_resid:
        return lc_rest, obs_sed, tempfilt['fnu'], signoise
        
    #### Plot colors
    colors = range(tempfilt['NFILT'])
    for ci, i in enumerate(np.argsort(lc)):
        colors[i] = threedhst.utils.color_table((ci+1.)/tempfilt['NFILT']*250, table='rainbow.rgb')
    
    fig = plt.figure(figsize=(12,4.8))
    fig.subplots_adjust(wspace=0.0, hspace=0.0, left=0.09, bottom=0.10, right=0.98, top=0.98)
    
    ax = fig.add_axes((0.06, 0.12, 0.6, 0.86))
    
    #### Plot the residuals
    xx, yy, ss = [], [], []
    stats = range(len(lc))
    
    keep = np.isfinite(signoise[0,:])
    for i in np.argsort(lc):
        keep &= (tempfilt['efnu'][i,:] > 0) #& (tempfilt['fnu'][i,:] > 0)
        # if signoise[i,:].max() > 3:
        #     keep & (signoise[i,:] > 3)
        #keep &= (signoise[i,:] > 3) #| (signoise[i,:] < 0.0001)#& (np.abs(resid[i,:]-1)/tempfilt['efnu'][i,:] < 5)
    
    nfilt = ((tempfilt['efnu'] > 0) & ((tempfilt['fnu']/zpfactors) > -90)).sum(axis=0)
    keep = nfilt > (nfilt.max()-5)
    
    for i in np.argsort(lc):
        #keep = signoise[i,:] > 3
        ok = keep & (resid[i,:] > 0) & (tempfilt['fnu'][i,:] > 0) #& (signoise[i,:] > 3)
        if np.std(zi[keep]) == 0:
            rnd = np.random.normal(size=keep.sum())*0.01*lc[i]
        else:
            rnd = 0.
        #
        sc = ax.plot(lc[i]/(1+zi[ok])+rnd, resid[i,ok], marker='.', alpha=0.05, linestyle='None', color=colors[i])
        #xm, ym, ys, nn = threedhst.utils.runmed(lc[i]/(1+zi[keep]), resid[i,keep], NBIN=int(len(keep)/1000.))
        xm, ym, ys, nn = threedhst.utils.runmed(lc[i]/(1+zi[ok]), resid[i,ok], NBIN=np.maximum(int(ok.sum()/1000.), 8))
        xx.append(xm)
        yy.append(ym)
        ss.append(ys)
        val = np.sum(resid[i,ok]*signoise[i,ok]**2)/np.sum(signoise[i,ok]**2)
        #print lc[i], i, ok.sum()
        if ok.sum() == 0:
            stats[i] = {'mean':0,
                      'median':0,
                      'std':1,
                      'stdmean':1,
                      'p':[0,0,0,0,0],
                      'pstd':0,
                      'val':0}            
        else:
            stats[i] = {'mean':np.mean(resid[i,ok]),
                      'median':np.median(resid[i,ok]),
                      'std':np.std(resid[i,ok]),
                      'stdmean':np.std(resid[i,ok])/np.sqrt(ok.sum()),
                      'p':np.percentile(resid[i,ok], [2.5,16,50,84,97.5]),
                      'pstd':(np.percentile(resid[i,ok], 84)-np.percentile(resid[i,ok], 16))/2/np.sqrt(ok.sum()),
                      'val':val}
        #
        # offsets[i] = stats[i]['median']
        # #
        # print '%s %.3f %.3f %.3f %.3f %d %.3f' %(param.filters[i].name, stats[i]['median'], stats[i]['pstd'], stats[i]['p'][1], stats[i]['p'][-2], ok.sum(), val)
        
    for ci, i in enumerate(np.argsort(lc)):
        ax.plot(xx[ci], yy[ci], alpha=0.2, color='black', linewidth=4)
        ax.plot(xx[ci], yy[ci], alpha=0.8, color=colors[i], linewidth=2)
    
    #### Adjustments to zeropoints
    lcfull = []
    residfull = []
    for i in np.argsort(lc):
        ok = keep & (resid[i,:] > 0) & (tempfilt['fnu'][i,:] > 0) & (signoise[i,:] > 3)
        lcfull.extend(lc[i]/(1+zi[keep]))
        residfull.extend(resid[i,keep]/stats[i]['median'])
        
    xmfull, ymfull, ysfull, nnfull = threedhst.utils.runmed(np.array(lcfull), np.array(residfull), NBIN=np.maximum(int(len(residfull)/2000.), 10))
    
    ymfull[xmfull > wclip[1]] = 1.
    ymfull[xmfull < wclip[0]] = 1.
    
    #### XXX gaussians
    try:
        from astroML.sum_of_norms import sum_of_norms, norm    
        n_gaussians = 30
        w_best, rms, locs, widths = sum_of_norms(np.log10(xmfull), ymfull-1, n_gaussians, spacing='linear', full_output=True)
        norms = w_best * norm(np.log10(xmfull)[:, None], locs, widths)
        ymfull = norms.sum(1)+1
    except:
        pass
        
    #plt.plot(xmfull, norms.sum(1)+1, color='orange', linewidth=2, zorder=1000)
        
    ax.plot(xmfull, ymfull, color='black', alpha=0.75, linewidth=2)
    
    if os.path.exists('tweak_%s/tweak.dat' %(root)):
        tx, ty = np.loadtxt('tweak_%s/tweak.dat' %(root), unpack=True)
        ty_int = np.interp(xmfull, tx, ty, left=1, right=1)
    else:
        ty_int = xmfull*0.+1
    
    ax.plot(xmfull, ymfull*ty_int, color='black', alpha=0.3, linewidth=2)
    np.savetxt('tweak_%s/tweak.dat' %(root), np.array([xmfull, ymfull*ty_int]).T, fmt='%.5e')
    
    NT = len(param.templates)
    fp = open('tweak_%s/spectra.param' %(root),'w')
    for i in range(NT):
        wt, ft = np.loadtxt(param.templates[i], unpack=True)
        ft /= np.interp(wt, xmfull, ymfull, left=1, right=1)
        np.savetxt('tweak_%s/%s' %(root, os.path.basename(param.templates[i])), np.array([wt, ft]).T, fmt='%.5e')
        fp.write('%d tweak_%s/%s 1.0 0 1.0\n' %(i+1, root, os.path.basename(param.templates[i])))
    #
    fp.close()
    
    ### account for overall wiggles
    for i in np.argsort(lc):
        lcz = lc[i]/(1+zi)
        yint = np.interp(lcz, xmfull, ymfull, left=-99, right=-99)
        ok = keep & (yint > 0) & (resid[i,:] > 0) & (tempfilt['fnu'][i,:] > 0) #& (signoise[i,:] > 3)
        #### ignore Lyman series absorption and IR
        ok = ok & (lcz > 1500) & (lcz < 2.5e4)
        #
        rfix = resid[i,:]/yint
        val = np.sum(rfix[ok]*signoise[i,ok]**2)/np.sum(signoise[i,ok]**2)
        if ok.sum() == 0:
            stats[i] = {'mean':0,
                      'median':0,
                      'std':1,
                      'stdmean':1,
                      'p':[0,0,0,0,0],
                      'pstd':0,
                      'val':0}            
        else:
            stats[i] = {'mean':np.mean(rfix[ok]),
                      'median':np.median(rfix[ok]),
                      'std':np.std(rfix[ok]),
                      'stdmean':np.std(rfix[ok])/np.sqrt(ok.sum()),
                      'p':np.percentile(rfix[ok], [2.5,16,50,84,97.5]),
                      'pstd':(np.percentile(rfix[ok], 84)-np.percentile(rfix[ok], 16))/2/np.sqrt(ok.sum()),
                      'val':val}
        #
        offsets[i] = stats[i]['median']
        #
        print '%s %d %.3f %.3f %.3f %.3f %d %.3f' %(param.filters[i].name, ok.sum(), stats[i]['median'], stats[i]['pstd'], stats[i]['p'][1], stats[i]['p'][-2], keep.sum(), val)
        
    if not os.path.exists(adjust_zeropoints):
        fp = open(adjust_zeropoints,'w')
        for filt in param.filters:
            fp.write('F%0d  1.0\n' %(filt.fnumber))
        
        fp.close()
    
    zpfilt, zpval = np.loadtxt(adjust_zeropoints, dtype=np.str, unpack=True)
    zpval = np.cast[float](zpval)
    
    # for ci, i in enumerate(np.argsort(lc)):
    #     if not STAR_FIT:
    #         keep_i = keep & (signoise[i,:] > 3) & (resid[i,:] > 0.2) & (np.abs(resid[i,:]-1)/tempfilt['efnu'][i,:] < 3)
    #         lcz = lc[i]/(1+zi[keep_i])
    #         yint = np.interp(lcz, xmfull, ymfull)
    #         err = np.sqrt(tempfilt['efnu'][i,keep_i]**2+(0.02*tempfilt['fnu'][i,keep_i])**2)
    #         offsets[i] = 1./(np.sum(resid[i,keep_i]*yint/err**2)/np.sum(resid[i,keep_i]**2/err**2))
    #         print i, keep_i.sum(), offsets[i]
            
    ## Normalize to first filter
    #offsets /= offsets[0]
    if ref_filter is not None:
        for ci, i in enumerate(np.argsort(lc)):
            if param.filters[i].fnumber == ref_filter:
                offsets /= offsets[i]
                print 'Norm to %s.' %(param.filters[i].name)
                break
            
    ref_offset = 1.
    if isinstance(fix_filter, dict):
        for ci, i in enumerate(np.argsort(lc)):
            if param.filters[i].fnumber in fix_filter.keys():
                offsets[i] = fix_filter[param.filters[i].fnumber]
                fstr = 'F%0d' %(param.filters[i].fnumber)
                if fstr in zpfilt:
                    zpval[zpfilt == fstr] = 1.
                # print 'Filt %d: %f' %(param.filters[i].fnumber, offsets[i])
    
    ## Write a new zphot.translate file
    for ci, i in enumerate(np.argsort(lc)):
        mat = zpfilt == 'F%0d' %(param.filters[i].fnumber)
        if (len(mat[mat]) > 0): # & (param.filters[i].lambda_c < 3.e4):
            zpval[mat] *= offsets[i]/ref_offset
    
    fp = open(adjust_zeropoints,'w')
    for i in range(len(zpval)):
        fp.write('%s  %.4f\n' %(zpfilt[i], zpval[i]))
    
    fp.close()
            
    #### Plot things
    ax.set_xlim(800,1.e5)
    ax.semilogx()
    ax.set_ylim(0.5,1.5)
    ax.set_xlabel(r'$\lambda_\mathrm{rest}\ [\mu\mathrm{m}]$')
    ax.set_ylabel(r'(temp - phot) / temp')
    ax.set_xticklabels([0.1,0.5,1,5])
    ytick = ax.set_xticks([1000,5000,1.e4,5e4])
    
    #### Add labels for filters
    NF = len(lc)
    font = np.minimum(9*28./NF, 14)
    
    for ci, i in enumerate(np.argsort(lc)):
        ax.text(0.99, 0.99-0.04*ci*28./NF, '%s %3d %.3f' %(os.path.basename(param.filters[i].name).replace('_','_'), param.filters[i].fnumber, offsets[i]/ref_offset), transform = ax.transAxes, color='0.5', horizontalalignment='right', va='top', fontsize=font)
        ax.text(0.99, 0.99-0.04*ci*28./NF, '%s %3d %.3f' %(os.path.basename(param.filters[i].name).replace('_','_'), param.filters[i].fnumber, offsets[i]/ref_offset), transform = ax.transAxes, color=colors[i], alpha=0.8, horizontalalignment='right', va='top', fontsize=font)

    #### Add zphot zspec plot
    ax = fig.add_axes((0.67, 0.12, 0.32, 0.86))
    zout = catIO.Readfile('%s/%s.zout' %(PATH, root))
    if 'z_peak' not in zout.columns:
        zout.z_peak = zout.z_spec
    
    dz = (zout.z_peak-zout.z_spec)/(1+zout.z_spec)
    keep = (zout.z_spec > 0)
    sigma = threedhst.utils.nmad(dz[keep])
    outlier = np.abs(dz[keep]) > 0.15
    
    ax.plot(np.log10(1+zout.z_spec), np.log10(1+zout.z_peak), marker='.', alpha=0.1, linestyle='None')
    ax.plot([0,10],[0,10], color='white', alpha=0.4, linewidth=3)
    ax.plot([0,10],[0,10], color='black', alpha=0.4, linewidth=1)
    ax.set_xticklabels([0,1,2,3,4])
    xtick = ax.set_xticks(np.log10(1+np.array([0,1,2,3,4])))
    ax.set_yticklabels([])
    ytick = ax.set_yticks(np.log10(1+np.array([0,1,2,3,4])))
    
    ### histograms
    if 'z_mc' in zout.columns:
        zz = zout.z_mc
    else:
        zz = zout.z_a
        
    yh, xh = np.histogram(np.log10(1+zz), range=[np.log10(1), np.log10(5)], bins=100)
    ax.fill_between(xh[1:], xh[1:]*0., (np.log10((yh+1.)/yh.max())/2.+1)*np.log10(2), color='black', alpha=0.1, zorder=-100)
    yh, xh = np.histogram(np.log10(1+zout.z_spec), range=[np.log10(1), np.log10(5)], bins=100)
    ax.fill_between(xh[1:], xh[1:]*0., (np.log10((yh+1.)/yh.max())/2.+1)*np.log10(2), color='blue', alpha=0.1, zorder=-100)
    
    ax.text(0.5, 0.9, r'$\sigma_\mathrm{nmad}=%.3f$, $f_\mathrm{>0.15}=%.3f$' %(sigma, outlier.sum()*1./keep.sum()), transform = ax.transAxes, color='black', horizontalalignment='center', fontsize=10)
    
    ax.set_xlim(0,np.log10(1+4))
    ax.set_ylim(0,np.log10(1+4))
    
    #### Save an output image
    if savefig is not None:
        fig.savefig(savefig)
        plt.close()
    
    return fnumbers, lc, offsets

def x_get_template_error(root='cosmos', PATH='OUTPUT/', apply=False):

    lc_rest, obs_sed, fnu, signoise = eazy.show_fit_residuals(root=root, PATH=PATH, savefig=None, adjust_zeropoints='zphot.zeropoint.XX', fix_filter=None, ref_filter=None, get_resid=True, wclip=[1200, 3.e4])
    
    resid = (obs_sed-fnu) / obs_sed + 1
    err = np.sqrt((fnu/signoise)**2-(0.01*fnu)**2)    
    full_err = np.sqrt(err**2) # + (0.05*fnu)**2) #
    
    deviate = (obs_sed - fnu) / full_err
    ok = (fnu > 0) & (signoise > 20) & (obs_sed > 0) & (np.abs(deviate) < 10) & (err > 0)
    xm, ym, ys, nn = threedhst.utils.runmed(lc_rest[ok], deviate[ok], NBIN=np.maximum(int(ok.sum()/500.), 8), use_nmad=True, use_median=True)

    if False:
        #plt.plot(xm, ym)
        plt.plot(xm, ys, color='black')
        for i in range(lc_rest.shape[0]):
            ok_i = (fnu[i,:] > 0) & (signoise[i,:] > 20) & (obs_sed[i,:] > 0) & (np.abs(deviate[i,:]) < 10) & (err[i,:] > 0)
            xmi, ymi, ysi, nni = threedhst.utils.runmed(lc_rest[i,ok_i], deviate[i,ok_i], NBIN=np.maximum(int(ok_i.sum()/500.), 8), use_nmad=True, use_median=True)
            plt.plot(xmi, ysi, alpha=0.5)
    
    #### Try to do it analytically
    g2 = 1./obs_sed**2*((fnu-obs_sed)**2/0.455-err**2)
    ok = ok & (g2 > 0)
    xg, yg, ygs, nn = threedhst.utils.runmed(lc_rest[ok], np.sqrt(g2[ok]), NBIN=np.maximum(int(ok.sum()/500.), 8), use_nmad=True, use_median=False)
    g_int = np.interp(lc_rest, xg, yg)
    full_err = np.sqrt(err**2+ (g_int*obs_sed)**2) #
    deviate = (obs_sed - fnu) / full_err
    
    xg = np.append(xg, 6.e4)
    yg = np.append(yg, 0.3)

    xg = np.append(xg, 800)
    yg = np.append(yg, yg[0])
    
    so = np.argsort(xg)
    xg, yg = xg[so], yg[so]
    
    from matplotlib import pyplot as plt
    from astroML.datasets import fetch_vega_spectrum
    from astroML.sum_of_norms import sum_of_norms, norm
    
    n_gaussians = 15
    w_best, rms, locs, widths = sum_of_norms(np.log10(xg), yg, n_gaussians,
                                             spacing='linear',
                                             full_output=True)

    tx, ty = np.loadtxt('templates/TEMPLATE_ERROR.eazy_v1.0', unpack=True)
    xx = np.log10(tx)
    #xx = np.log10(np.logspace(np.log10(800), np.log10(8.e4), 200))
    norms = w_best * norm(xx[:, None], locs, widths)
    plt.plot(xg, yg, color='black', linewidth=2, alpha=0.7)
    plt.plot(10**xx, norms, ls='-', c='#FFAAAA', alpha=0.2)
    #plt.plot(10**xx, norms.sum(1), '-r', label='sum of gaussians', alpha=0.8, linewidth=2)
    plt.plot(tx, ty*0.5, color='green', linewidth=1, alpha=0.5, label='TE v1.0')
    
    plt.xlim(400, 9.e4); plt.semilogx()
    
    new_te = norms.sum(1)
    new_te[tx >= 6.e4] = 0.3
    new_te[tx <= 1000] = yg[0]
    plt.plot(10**xx, new_te, '-r', label='sum of gaussians', alpha=0.8, linewidth=2)
    
    ttx, tty = np.loadtxt('templates/TEMPLATE_ERROR.v2.0.zfourge', unpack=True)
    plt.plot(ttx, tty, color='orange')
    
    np.savetxt('TEMPLATE_ERROR.v2.x', np.array([tx, new_te]).T, fmt='%.3e')
    
    #plt.scatter(lc_rest[ok], deviate[ok], alpha=0.02)
    # xm, ym, ys, nn = threedhst.utils.runmed(lc_rest[ok], deviate[ok], NBIN=np.maximum(int(ok.sum()/500.), 8), use_nmad=True, use_median=True)
    # plt.plot(xm, ys)
    # 
    # tx, ty = np.loadtxt('templates/TEMPLATE_ERROR.eazy_v1.0', unpack=True)
    # plt.plot(tx, ty/2.*5)
    
    
def log_hist(z, r0=(0,4), dz=0.002, plot=False, *args, **kwargs):
    """
    Make a histogram of redshifts with uniform spacing in dz/(1+z)
    """
    ran = [np.log10(1+r0[0]), np.log10(1+r0[1])]
    yh, xh = np.histogram(np.log10(1+z), range=ran, bins=int((ran[1]-ran[0])*1./dz*np.log(10)))
    lx = 10**xh-1
    
    if plot:
        plt.plot(lx[1:], np.maximum(yh, 0.0001), linestyle='steps-post', *args, **kwargs)
    
    return lx, yh
    
    
    
def loop_zeropoints(root='cosmos', tfile='zphot.translate.cosmos',  zfile='zphot.zeropoint.cosmos', fix_filter={}, ref_filter=None, init_filter={}, ignore_initial=[], ignore_all=[], toler=0.005, PATH='./OUTPUT/', fix_zspec=False, check_uvj=False, use_tweaked_templates=True, MAXITER=15, MIN_ITER=2, wclip=[1200, 3.e4]):
    """
    Wrapper around `show_fit_residuals` to allow iterative fitting
    of the zeropoint offsets.
    
    There is a bit of a black art in deciding which filters to ignore for the
    first iteration, which should typically be those with the most discrepant
    zeropoint offsets. By "ignoring" bands with `ignore_initial`, the function
    sets a very large error for those bands by adding a very large scaling to
    the error column in the translate file. That way a residual will still be
    measured but the band won't influence the initial photo-z fit.
    
     `ref_filter` is A reference filter used to normalize the offets so that
    the iteration doesn't go off and change the overall normalization too
    much. This is typically the detection band used for the catalog, say Ks or
    WFC3/F160W.
    
     The `init_filter` parameter contains offsets applied at the first step,
    which otherwise defaults to unity (1) for all bands.
    
     `fix_filter` is a dictionary containing fixed offsets used at all
    iterations.
    
     If `use_tweaked_templates` is set, after the second iteration the
    "tweaked" template set will be used, which contains the
    wavelength-dependent multiplicative term to all of the templates.
    
     If `fix_zspec` is used, the EAZY FIX_ZSPEC parameter will be set. Note
    that this requires only objects with z_spec > 0 in the catalog, a "bug"
    that may fix in the future.
    
     `toler` is the tolerance used to evaluate convergence of the iterations.
    The default value of 0.005 means that the iterations will continue while
    the offsets in any band change by more than 0.5% from one iteration to
    another. Note that the offsets at wavelengths < 4500 A and > 2.5 microns
    are not considered in the tolerance calculation, as those "wagging tail"
    bands at the wavelength extremes are typically less well constrained.
    
    Example: 
    ========
    
    #### Initial guess to get the very discrepant bands close
    init_filter = {88:1.18, 81: 0.85, 82:1.2, 79:0.8, 190: 0.8, 192:0.7}
    
    #### Run once to generate the necessary files in OUTPUT/
    os.system('eazy -p zphot.param.cosmos -t zphot.translate.cosmos')

    #### Run the iteration
    threedhst.eazyPy.loop_zeropoints(root='cosmos', tfile='zphot.translate.cosmos',  zfile='zphot.zeropoint.cosmos', ref_filter=205, fix_filter={205:1}, ignore_initial=['f_Zp', 'f_Ip', 'f_Rp'], toler=0.005, init_filter=init_filter) 
    
    """
    import threedhst
    import threedhst.eazyPy as eazy
    import os
    
    try:
        import unicorn.zp
    except ImportError:
        if check_uvj:
            print "`unicorn` module not found.  Can't do check_uvj."
            check_uvj = False
        
    if not PATH.endswith('/'):
        PATH += '/'
    
    ##### Read the param file
    param = eazy.EazyParam('%s.param' %(os.path.join(PATH, root)))
    tf = eazy.TranslateFile(tfile)
    
    ##### Make an empty zeropoint file
    fp = open(zfile,'w')
    for filter in param.filters:
        if filter.fnumber in fix_filter.keys():
            fp.write('F%d  %s\n' %(filter.fnumber, fix_filter[filter.fnumber]))
        else:
            if filter.fnumber in init_filter.keys():
                fp.write('F%d  %s\n' %(filter.fnumber, init_filter[filter.fnumber]))
            else:
                fp.write('F%d  1.0\n' %(filter.fnumber))
    #
    fp.close()
    
    try:
        os.mkdir('tweak_%s' %(root))
    except:
        pass
        
    clean_files = [os.path.join(PATH, root)+'.zeropoint', 'tweak_%s/tweak.dat' %(root)]
    for file in clean_files:
        if os.path.exists(file):
            os.remove(file)
    #
    NT = len(param.templates)
    fp = open('tweak_%s/spectra.param' %(root),'w')
    for i in range(NT):
        wt, ft = np.loadtxt(param.templates[i], unpack=True)
        np.savetxt('tweak_%s/%s' %(root, os.path.basename(param.templates[i])), np.array([wt, ft]).T, fmt='%.5e')
        fp.write('%d tweak_%s/%s 1.0 0 1.0\n' %(i+1, root, os.path.basename(param.templates[i])))
    
    fp.close()
    
    param.params['GET_ZP_OFFSETS'] = True
    param.params['FIX_ZSPEC'] = fix_zspec
    
    param.write('zphot.param.iter.%s' %(root))
    
    #### Scale huge errors for first step
    for filter in ignore_initial:
        tf.change_error(filter, 1.e6)
    
    for filter in ignore_all:
        tf.change_error(filter, 1.e6)
    
    tf.write('zphot.translate.iter.%s' %(root))
    
    os.system('eazy -p zphot.param.iter.%s -t zphot.translate.iter.%s -z %s' %(root, root, zfile))
    
    fnumbers, lc_i, delta_i = eazy.show_fit_residuals(root=root, fix_filter=fix_filter, ref_filter=ref_filter, adjust_zeropoints=zfile, savefig='%s_iter_%03d.png' %(root, 0), wclip=wclip)
    
    #### Extract UVJ
    if check_uvj:
        for color in ['153,155', '155,161']:
            param.params['REST_FILTERS'] = color
            param.params['READ_ZBIN'] = 'y'
            param.write('zphot.param.iter.%s' %(root))
            os.system('eazy -p zphot.param.iter.%s -t zphot.translate.iter.%s -z %s' %(root, root, zfile))
        #
        param.params['READ_ZBIN'] = 'n'
        param.write('zphot.param.iter.%s' %(root))
        unicorn.zp.diagnostic_UVJ(root=root, ext='UVJ_%04d' %(0))
        
    fp = open('%s_iter.log' %(root),'w')
    fp.write('\n\nIter #%d\n======\n' %(0))
    eazy.log_offsets(fp, fnumbers, lc_i, delta_i, toler)
    
    ### Now loop
    for filter in ignore_initial:
        tf.change_error(filter, 1.)
    #
    for filter in ignore_all:
        tf.change_error(filter, 1.e6)
    
    tf.write('zphot.translate.iter.%s' %(root))
    
    param.params['GET_ZP_OFFSETS'] = True
    param.write('zphot.param.iter.%s' %(root))
    
    for i in range(MAXITER):
        #### Use tweaked templates after two iterations
        if i == 0:
            try:
                os.remove('tweak_%s/tweak.dat' %(root))
            except:
                pass
            #
        
        if (i == 1) & (use_tweaked_templates):
            param.params['TEMPLATES_FILE'] = 'tweak_%s/spectra.param' %(root)
            param.write('zphot.param.iter.%s' %(root))
        
        os.system('eazy -p zphot.param.iter.%s -t zphot.translate.iter.%s -z %s' %(root, root, zfile))
        fnumbers, lc_i, delta_i = eazy.show_fit_residuals(root=root, fix_filter=fix_filter, ref_filter=ref_filter, adjust_zeropoints=zfile, savefig='%s_iter_%03d.png' %(root, i+1), wclip=wclip)
        fp.write('\n\nIter #%d\n======\n' %(i))
        eazy.log_offsets(fp, fnumbers, lc_i, delta_i, toler)
        
        #
        # #### Extract UVJ
        # if check_uvj:
        #     for color in ['153,155', '155,161']:
        #         param.params['REST_FILTERS'] = color
        #         param.params['READ_ZBIN'] = 'y'
        #         param.write('zphot.param.iter')
        #         os.system('eazy -p zphot.param.iter -t zphot.translate.iter -z %s' %(zfile))
        #     #
        #     param.params['READ_ZBIN'] = 'n'
        #     param.write('zphot.param.iter')
        #     unicorn.zp.diagnostic_UVJ(root=root, ext='UVJ_%04d' %(i+1))
        
        fixed_bands = lc_i < 0
        for fi in fix_filter.keys():
            fixed_bands = fixed_bands | (fnumbers == fi)
        
        use_bands = (lc_i > 4500) & (lc_i < 2.5e4) & ~fixed_bands
        
        if (np.abs(delta_i[use_bands]-1).max() < toler) & (i >= MIN_ITER):
            break
    
    fp.close()
    
    pass

#
def log_offsets(fp, fnumbers, lc_i, delta_i, toler):
    so = np.argsort(lc_i)
    for j in so:
        if np.abs(delta_i[j]-1) > toler:
            log = '* F%d  %.4f' %(fnumbers[j], delta_i[j])
        else:
            log = '  F%d  %.4f' %(fnumbers[j], delta_i[j])
        #
        print log
        fp.write(log+'\n')

def compute_taylor_mass(root='photz', PATH='OUTPUT', gi=[157,159], ABZP=25):
    """
    Compute stellar masses from Ned Taylor's g-i / Mi relation
    
    log M/Msun = 1.15 + 0.70 (g-i) - 0.4 Mi
    
    (Taylor 2011, Taylor 2014)
    
    """
    
    rfg = catIO.Readfile('%s/%s.%0d.rf' %(PATH, root, gi[0]))
    mg = ABZP - 2.5*np.log10(rfg['l%0d' %(gi[0])])
    Mg = mg-rfg.dm

    rfi = catIO.Readfile('%s/%s.%0d.rf' %(PATH, root, gi[1]))
    mi = ABZP - 2.5*np.log10(rfi['l%0d' %(gi[1])])
    Mi = mi-rfi.dm
    
    logM = 1.15 + 0.70*(mg-mi) - 0.4*Mi
    return (mg-mi), Mi, logM
    
def compute_eazy_mass(root='photz', PATH='OUTPUT', rf_file='153-155', V_filter='155', ABZP=25, MLv_templates = [4.301, 0.059, 0.292, 0.918, 2.787, 0.940, 4.302], line_correction_V=[ 0.997,  0.976,  0.983,  0.999,  0.999,  0.951]):
    """
    Estimate stellar masses simply from the estimate M/Lv of the templates.  
    Requires a rest-frame color file that samples the v-band (e.g., #155).
    
    The default MLv_templates is determined for the eazy v1.1_lines template set fit with
    Conroy & Gunn (2009) templates.
    
    Includes a correction for emission line fluxes computed for the 
    EAZY v1.1 line set in `line_correction_V`, which is the ratio of rest
    V flux with and without the template emission lines.
    
    returns:
        Lv = V-band luminosity in solar units
        M/Lv = Mass-to-light in V-band (solar units)
        Mass = Stellar mass
        
    """
    import threedhst
    from threedhst import eazyPy as eazy
    from threedhst import catIO
    
    tempfilt, coeffs, temp_sed, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE=root, OUTPUT_DIRECTORY=PATH, CACHE_FILE='Same')
    
    rf = catIO.Readfile('%s/%s.%s.rf' %(PATH, root, rf_file))
    mv = ABZP - 2.5*np.log10(rf['l%s' %(V_filter)])
    Mv = mv-rf.dm
    fnu = 10**(-0.4*(Mv+48.6))            # erg / s / cm2 / Hz, D=10pc
    Lnu = fnu*4*np.pi*(10*3.09e18)**2     # erg / s / Hz
    nuLnu = (3.e8/5500.e-10)*Lnu/3.839e33 # Lsun
    
    MLv_template = np.array(MLv_templates) #[4.301, 0.059, 0.292, 0.918, 2.787, 0.940, 4.302])
    
    ### Check that have enough parameters in the line correction array.  
    ### Fill with ones otherwise
    if len(MLv_templates) != len(line_correction_V):
        threedhst.showMessage(' '*57+'\nlen(MLv_templates) [%d] != len(line_correction_V) [%d]     ' %(len(MLv_templates), len(line_correction_V)), warn=True)
        line_correction_V = np.append(np.array(line_correction_V), np.ones(len(MLv_templates) - len(line_correction_V)))
    
    MLv_template *= np.array(line_correction_V)
    
    scl = np.mean(coeffs['tnorm'][0:5])
    MLv = np.dot(MLv_template.reshape(1,-1), coeffs['coeffs']/scl)/np.sum(coeffs['coeffs'], axis=0)
    MLv = MLv[0,:]
    eazyMass = (MLv * nuLnu)
    
    return nuLnu, MLv, np.log10(eazyMass)
    
def compute_eazy_lineflux(root='photz', PATH='OUTPUT', rf_file='153-155', V_filter='155', ABZP=25):
    """
    Estimate emission line fluxes from hard-coded line / rest V flux 
    ratios
    """
    from threedhst import eazyPy as eazy
    
    tempfilt, coeffs, temp_sed, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE=root, OUTPUT_DIRECTORY=PATH, CACHE_FILE='Same')
    
    rf = catIO.Readfile('%s/%s.%s.rf' %(PATH, root, rf_file))
    mv = ABZP - 2.5*np.log10(rf['l%s' %(V_filter)])
    fnu_V = 10**(-0.4*(mv+48.6)) # erg / s / cm2 / Hz
    flam_V = fnu_V * 3.e18 / 5479.35**2 # erg / s / cm2 / A
    
    ### EAZY v1.1 templates.  
    ### Integrated line flux divided by rest-frame V (#155) flux
    corr = {'Ha':np.array([17.3, 122.1,108.8,2.444,0.8685,346.2]),
            'O3':np.array([2.474, 31.74,  20.9, 0.6329, 0.2167, 62.05]), 
            'O2':np.array([2.386,  60.86,  26.98, 1.419, 0.4854, 70.04])}
    
    line_flux = {}
    
    for key in corr.keys():
        line_ratio = corr[key]
        ### fill with zeros for missing templates
        if len(line_ratio) != tempfilt['NTEMP']:
            line_ratio = np.append(line_ratio, np.zeros(tempfilt['NTEMP']-len(line_ratio)))
        #
        scl = np.mean(coeffs['tnorm'][0:5])
        line_flux[key] = np.dot(line_ratio.reshape(1,-1), coeffs['coeffs']/scl)/np.sum(coeffs['coeffs'], axis=0)*flam_V
    
    return line_flux
        
        
def compute_template_line_fluxes():
    """
    Compute two parameters for emission line fluxes in the EAZY templates:
    
        1) Contribution of lines to the rest-frame V flux for fitting masses
        
        2) Coefficients for predicting [OII], [OIII], Halpha line fluxes
           from the templates:  line_i / f_V  (scale to flam with APZP)
    """
    
    import threedhst.eazyPy as eazy
    import matplotlib.pyplot as plt
    import research.v4 as cat2
    
    #os.chdir("/Users/brammer/3DHST/Spectra/Release/v4.0/Eazy/EazyRerun")
    
    res = eazy.FilterFile('FILTER.RES.latest')
    
    files = glob.glob('templates/EAZY_v1.0_lines/e*nolines.dat')
    no_lines = []
    for file in files:
        no_lines.append(eazy.Template(file))
    #
    files = glob.glob('templates/EAZY_v1.0_lines/e*sed?.dat')
    files = glob.glob('templates/EAZY_v1.1_lines/e*sed?.dat')

    with_lines = []
    for file in files:
        with_lines.append(eazy.Template(file))
    
    # files = glob.glob('templates/EAZY_v1.1_lines/e*sed?.dat')
    # v11_lines = []
    # for file in files:
    #     v11_lines.append(eazy.Template(file))
    # for i in range(6):
    #     with_V = with_lines[i].integrate_filter(res.filters[161-1])
    #     plt.plot(with_lines[i].wavelength, with_lines[i].flux/with_V)
    #     plt.plot(v11_lines[i].wavelength, v11_lines[i].flux/with_V)
    
    x = eazy.readRFBinary(file='OUTPUT/aegis.dusty3.155.coeff')
    tempfilt, coeffs, temp_sed, pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE='aegis.dusty3', OUTPUT_DIRECTORY='OUTPUT', CACHE_FILE='Same')
    
    xx = []
    print '# file UmV  lineV  HaV O3V O2V'
    for i in range(6):
        no_V = no_lines[i].integrate_filter(res.filters[155-1])
        with_V = with_lines[i].integrate_filter(res.filters[155-1])
        with_U = with_lines[i].integrate_filter(res.filters[153-1])
        print '%s %.3f' %(files[i], with_V)
        xx = np.append(xx, with_V)
        #
        #print '%s %.3f  %.3f' %(files[i], -2.5*np.log10(with_U/with_V), no_V / with_V)
        ### Integrate H-alpha flux
        line_only = with_lines[i].flux - np.interp(with_lines[i].wavelength, no_lines[i].wavelength, no_lines[i].flux)
        limits = {'Ha':[6530,6590], 'O3':[4900, 5040], 'O2':[3690, 3760]}
        line_flux_ratio = {}
        for key in limits.keys():
            wrange = (with_lines[i].wavelength >= limits[key][0]) & (with_lines[i].wavelength < limits[key][1])
            line_flux = np.trapz(line_only[wrange], with_lines[i].wavelength[wrange])
            line_flux_ratio[key] =  line_flux / with_V
        #
        print '%s %.4f  %.3f  %.3e %.3e %.3e' %(files[i], -2.5*np.log10(with_U/with_V), no_V / with_V, line_flux_ratio['Ha'], line_flux_ratio['O3'], line_flux_ratio['O2'])

def init_nmf(obj, iz, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same', verbose=True):
    import threedhst.eazyPy as eazy
    
    MAIN_OUTPUT_FILE = 'cosmos-1.deblend.v5.1'
    OUTPUT_DIRECTORY = './cosmos-1.deblend.redshifts'
    CACHE_FILE='Same'
    obj = 2407
    iz = coeffs['izbest'][obj]
    
    eazy.param = eazy.EazyParam(PARAM_FILE=OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.param')
    
    try:
        eazy.NTEMP = eazy.tempfilt['NTEMP']
        eazy.NFILT = eazy.tempfilt['NFILT']
        eazy.NZ = eazy.tempfilt['NZ']
    except:
        if verbose:
            print 'Read EAZY binary files (%s/%s)....' %(OUTPUT_DIRECTORY, MAIN_OUTPUT_FILE)
        eazy.tempfilt, eazy.coeffs, eazy.temp_sed, eazy.pz = eazy.readEazyBinary(MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, CACHE_FILE=CACHE_FILE)
        eazy.werr, eazy.terr = np.loadtxt('templates/TEMPLATE_ERROR.eazy_v1.0', unpack=True)
    
    #
    
    fnu = tempfilt['fnu'][:,obj]
    efnu = tempfilt['efnu'][:,obj]
    
    lc_rest = tempfilt['lc']/(1+eazy.tempfilt['zgrid'][iz])
    template_error = np.interp(lc_rest, eazy.werr, eazy.terr)

    var = (fnu*eazy.param['SYS_ERR'])**2++(fnu*template_error)**2+efnu**2
    sig = np.sqrt(var)
    
    mask = (fnu > -99) & (efnu > 0)
    maskNFILT = len(mask[mask])
    
    #### A matrix for NMF
    aa = tempfilt['tempfilt'][:,:,iz]/np.dot(sig.reshape(eazy.NFILT, 1), np.ones((1, eazy.NTEMP)))
    eazy.amatrix = np.dot(aa[mask,:].transpose(), aa[mask,:])
    
    eazy.bvector = -np.dot((fnu[mask]/var[mask]).reshape(1,maskNFILT), tempfilt['tempfilt'][mask,:,iz]).reshape(eazy.NTEMP)
    
    nit, coeffs_fit = eazy.nonneg_fact()
    
    #### Test
    obs_fit = np.dot(tempfilt['tempfilt'][:,:,iz], coeffs_fit)
    plt.plot(lc[mask], fnu[mask]/(lc[mask]/5500.)**2)
    plt.plot(lc[mask], obs_fit[mask]/(lc[mask]/5500.)**2)
    
def nonneg_fact(toler = 1.e-4, verbose=False):
    import threedhst.eazyPy as eazy
    from scipy import weave
    from scipy.weave import converters
    
    out_coeffs = np.ones(eazy.NTEMP)
    
    NTEMP = eazy.NTEMP
    
    MAXITER=100000 # //// 1.e5
    tol = 100
    itcount=0
    
    while (tol > toler) & (itcount < MAXITER):
        tolnum=0.
        toldenom = 0.
        tol=0
        for i in range(NTEMP):
            vold = out_coeffs[i]
            av = np.sum(eazy.amatrix[i,:]*out_coeffs)
            #////// Update coeffs in place      
            out_coeffs[i] *= -1.*eazy.bvector[i]/av
            tolnum += np.abs(out_coeffs[i]-vold)
            toldenom += vold
        #
        tol = tolnum/toldenom
        itcount += 1
    
    return itcount, out_coeffs
    
    c_code = """
    long i,j;
    long itcount,MAXITER,NTEMP;
    double tolnum,toldenom,tol;
    double vold,av;
    //
    NTEMP = 7;
    MAXITER=100000; //// 1.e5
    tol = 100;
    itcount=0;
    while (tol>toler && itcount<MAXITER) {
        tolnum=0.;
        toldenom = 0.;
        tol=0;
        for (i=0;i<NTEMP;++i) {
            vold = coeffs[i];
            av = 0.;
            for (j=0;j<NTEMP;++j) av += amatrix[i][j]*coeffs[j];
            ////// Update coeffs in place      
            coeffs[i]*=-1.*bvector[i]/av;
            tolnum+=fabs(coeffs[i]-vold);
            toldenom+=vold;
        }
        tol = tolnum/toldenom;
        ++itcount;
    }
    """ 
    
    #### Doesn't work with the matrix....
    #result = weave.inline(c_code,['amatrix','bvector','coeffs','toler'], compiler = 'gcc', verbose=2)
    
def milkyway_extinction(lamb=None, Rv=3.1):
    """
    Return the average milky way extinction curve A(lambda)/A(V~5500 A) from
    Cardelli, Clayton & Mathis (1989).  
    
    Functional form of the MW extinction curve taken from the `astropysics` library: 
    http://packages.python.org/Astropysics/coremods/obstools.html
    
    Input `lamb` is a scalar or array of wavelength values, in Angstroms.
    """
    scalar=np.isscalar(lamb)
    x=1e4/np.array(lamb,ndmin=1) #CCM x is 1/microns
    a,b=np.ndarray(x.shape,x.dtype),np.ndarray(x.shape,x.dtype)

    #### Allow extrapolation
    #if any((x<0.3)|(10<x)):
    #    raise ValueError('some wavelengths outside CCM 89 extinction curve range')
    #
    if any((10<x)):
        print "\nWARNING: MW extinction curve extrapolated at lam < 1000 A\n"
    
    if any((0.3>x)):
        print "\nWARNING: MW extinction curve extrapolated at lam > 3.3 micron\n"
    
    
    #irs=(0.3 <= x) & (x <= 1.1)
    irs=(x <= 1.1)
    opts = (1.1 <= x) & (x <= 3.3)
    nuv1s = (3.3 <= x) & (x <= 5.9)
    nuv2s = (5.9 <= x) & (x <= 8)
    fuvs = (8 <= x) #& (x <= 10)

    #TODO:pre-compute polys

    #CCM Infrared
    a[irs]=.574*x[irs]**1.61
    b[irs]=-0.527*x[irs]**1.61

    #CCM NIR/optical
    a[opts]=np.polyval((.32999,-.7753,.01979,.72085,-.02427,-.50447,.17699,1),x[opts]-1.82)
    b[opts]=np.polyval((-2.09002,5.3026,-.62251,-5.38434,1.07233,2.28305,1.41338,0),x[opts]-1.82)

    #CCM NUV
    y=x[nuv1s]-5.9
    Fa=-.04473*y**2-.009779*y**3
    Fb=-.2130*y**2-.1207*y**3
    a[nuv1s]=1.752-.316*x[nuv1s]-0.104/((x[nuv1s]-4.67)**2+.341)+Fa
    b[nuv1s]=-3.09+1.825*x[nuv1s]+1.206/((x[nuv1s]-4.62)**2+.263)+Fb

    a[nuv2s]=1.752-.316*x[nuv2s]-0.104/((x[nuv2s]-4.67)**2+.341)
    b[nuv2s]=-3.09+1.825*x[nuv2s]+1.206/((x[nuv2s]-4.62)**2+.263)

    #CCM FUV
    a[fuvs]=np.polyval((-.070,.137,-.628,-1.073),x[fuvs]-8)
    b[fuvs]=np.polyval((.374,-.42,4.257,13.67),x[fuvs]-8)

    AloAv = a+b/Rv

    if scalar:
        return AloAv[0]
    else:
        return AloAv
     
#
from . import __file__ as rootfile
igm_z, igm_da, igm_db = np.loadtxt(os.path.join( os.path.dirname(rootfile), 'data', 'igm_factors.txt'), unpack=True)

def igm_factor(wavelength, z):
    
    wz = wavelength/(1+z)
    lim1 = (wz < 912)
    lim2 = (wz >= 912) & (wz < 1026)
    lim3 = (wz >= 1026) & (wz < 1216)
    
    da = np.interp(z, igm_z, igm_da)
    db = np.interp(z, igm_z, igm_db)
    
    factor = wavelength*0.+1
    
    if lim1.sum() > 0: factor[lim1] *= 0.
    if lim2.sum() > 0: factor[lim2] *= 1.-db
    if lim3.sum() > 0: factor[lim3] *= 1.-da
    
    return factor
#
def add_filters():
    """
    Script to add new filters, including calculating AB-Vega offsets and central wavelengths
    """
    import pysynphot as S

    from threedhst import eazyPy as eazy
    from threedhst import catIO
    import unicorn.utils_c

    os.chdir('/usr/local/share/eazy-filters')
    
    files, scale, apply_atm = np.loadtxt('add.list', unpack=True, skiprows=2, dtype=str)
    res = eazy.FilterFile('FILTER.RES.latest.mod')
    atm = catIO.Readfile('mktrans_zm_10_10.dat')
    sp_atm = S.ArraySpectrum(wave=atm.wave*1.e4, flux=atm.transmission)
    
    N = len(files)
    for i in range(N):
        #### Already in filter file?
        s = res.search(files[i], verbose=False)
        if len(s) > 0:
            continue
        #
        filt = eazy.FilterDefinition()
        filter_name = files[i]
        #
        wf, tf = np.loadtxt(files[i], unpack=True)
        wf *= float(scale[i])
        bp = S.ArrayBandpass(wave=wf, throughput=tf)
        R = 500.
        dl_R = bp.pivot()/R
        x_resamp = np.arange(wf.min()-2.*dl_R, wf.max()+3.1*dl_R, dl_R)
        bp_resamp = unicorn.utils_c.interp_conserve_c(x_resamp, bp.wave, bp.throughput)
        if int(apply_atm[i]):
            filter_name += ' +atm '
            sp_atm_resamp = unicorn.utils_c.interp_conserve_c(x_resamp, sp_atm.wave, sp_atm.flux)
            bp_resamp *= sp_atm_resamp
        #
        bp_resamp[x_resamp <= wf.min()] = 0.
        bp_resamp[x_resamp >= wf.max()] = 0.
        filt.wavelength = x_resamp*1
        filt.transmission = bp_resamp*1
        filt.name = '%s lambda_c= %.4e AB-Vega=%.3f' %(filter_name, filt.pivot(), filt.ABVega())
        print filt.name
        res.filters.append(filt)
        res.NFILT += 1
    #
    res.write('FILTER.RES.latest.mod')

def quadri_pairs(zoutfile='OUTPUT/cdfs.zout', catfile=''):
    pass
    c = catIO.Readfile('../Cat2.2/uds.v0.0.15.a.cat.bright', force_lowercase=False)    
    z = catIO.Readfile('OUTPUT/uds.zout')

    c = catIO.Readfile('Catalog/cosmos.v0.10.7.a.cat.bright', force_lowercase=False)    
    z = catIO.Readfile('OUTPUT/cosmos.zout')

    c = catIO.Readfile('Catalog/cdfs.v0.4.8.a.cat.bright', force_lowercase=False)    
    z = catIO.Readfile('OUTPUT/cdfs.zout')
    kmag = 25-2.5*np.log10(c.Kstot)
    
    c = catIO.Readfile('../Catalogs/goodss_3dhst.v4.0.nzpcat.fixconv.bright', force_lowercase=False)    
    z = catIO.Readfile('OUTPUT/goodss.zout')

    c = catIO.Readfile('../Catalogs/uds_3dhst.v4.0.nzpcat.HAWKI.bright', force_lowercase=False)    
    z = catIO.Readfile('OUTPUT/uds.zout')

    kmag = 25-2.5*np.log10(c.f_F160W)
    
    m = catIO.CoordinateMatcher(c)
    ## find pairs
    first, next, drs = [], [], []
    for i in range(c.N):
        print unicorn.noNewLine+'%d/%d' %(i,c.N)
        dr, idx = m.find_nearest(c.ra[i], c.dec[i], N=60, distance_upper_bound=30./3600.)
        ok = (dr > 0) & (dr < 30.)
        first.extend([i]*ok.sum())
        next.extend(list(idx[ok]))
        drs.extend(list(dr[ok]))
        
    first = np.array(first)
    next = np.array(next)
    drs = np.array(drs)
    
    mlim = (18,24)
    zlim = (0.2,10)
    ok = (drs < 25) & (kmag[first] > mlim[0]) & (kmag[next] > mlim[0]) & (kmag[first] < mlim[1]) & (kmag[next] < mlim[1]) & (z.z_peak[first] > zlim[0]) & (z.z_peak[next] > zlim[0]) & (z.z_peak[first] < zlim[1]) & (z.z_peak[next] < zlim[1])
    
    dz = (z.z_peak[first]-z.z_peak[next])/(1+z.z_peak[first])

    #plt.scatter(drs, dz, alpha=0.03)

    yh, xh = np.histogram(dz[ok], bins=200, range=(-0.3,0.3))
    
    NEXTRA = 20
    z_rnd1 = z.z_peak[first][ok][np.cast[int](np.random.rand(ok.sum()*NEXTRA)*ok.sum())]
    z_rnd2 = z.z_peak[first][ok][np.cast[int](np.random.rand(ok.sum()*NEXTRA)*ok.sum())]
    #next_rnd = np.cast[int](np.random.rand(len(next))*z.N)
    dz_rnd = (z_rnd1-z_rnd2)/(1+z_rnd1)
    yhr, xhr = np.histogram(dz_rnd, bins=200, range=(-0.3,0.3))
    
    #plt.plot(xh[:-1], yh, linestyle='steps-mid', alpha=0.5)
    #plt.plot(xhr[:-1], yhr*1./NEXTRA, linestyle='steps-mid', alpha=0.5)
    plt.plot(xhr[1:], yh-yhr*1./NEXTRA, alpha=0.5, color='blue') # , linestyle='steps')
    err = np.sqrt(yhr)
    plt.fill_between(xhr[1:], yh-yhr*1./NEXTRA+err, yh-yhr*1./NEXTRA-err, color='blue', alpha=0.1)
    
def spatial_offset(root='cosmos', PATH='OUTPUT/', apply=False, candels=False):
    """
    Make a figure showing the *spatial* zeropoint residuals for each filter in a catalog
    """
    import threedhst.eazyPy as eazy
    from threedhst import catIO
    
    lc_rest, obs_sed, fnu, signoise = eazy.show_fit_residuals(root=root, PATH=PATH, savefig=None, adjust_zeropoints='zphot.zeropoint.XX', fix_filter=None, ref_filter=None, get_resid=True, wclip=[1200, 3.e4])
    
    param = eazy.EazyParam(PARAM_FILE=PATH+'/'+root+'.param')
    fnumbers = np.zeros(len(param.filters), dtype=np.int)
    for i in range(len(fnumbers)):
        fnumbers[i] = int(param.filters[i].fnumber)
    
    c = catIO.Table(param['CATALOG_FILE'])
    NF = len(fnumbers)
    NX = int(np.round(np.sqrt(NF)))
    NY = int(np.ceil(NF*1./NX))

    plt.ioff()
    #plt.gray()
    plt.set_cmap('RdYlBu')
    fig = plt.figure(figsize=[12,12])
    
    ###
    trans = [line.split()[0] for line in open(PATH+'/'+root+'.translate').readlines()]
    
    cnames = []
    for col in c.columns:
        if candels:
            if col.endswith('FLUX') & (col in trans):
                cnames.append(col)
        else:
            if col.startswith('f_') & (col in trans):
                cnames.append(col)
    
    if candels:
        c.rename_column('RA','ra')
        c.rename_column('DEC','dec')
        c.add_column(catIO.Column(name='x', data=(c['ra']-np.median(c['ra'])*3600.)))
        c.add_column(catIO.Column(name='y', data=(c['dec']-np.median(c['dec'])*3600.)))
            
    print cnames
        
    fitters = {}
    translate = {}
    
    for i in range(len(fnumbers)):
        #print i, len(fnumbers), len(cnames)
        
        ax = fig.add_subplot(NY, NX, i+1)
        print param.filters[i].name, cnames[i]
        dmag = -2.5*np.log10((fnu/obs_sed)[i,:])
        ok = np.isfinite(dmag) & (signoise[i,:] > 3) & (np.abs(dmag) < 0.1)
        ratio = (fnu/obs_sed)[i,:]
        
        if ok.sum() == 0:
            ax.set_xticklabels([]); ax.set_yticklabels([])
            continue
            
        ### Interpolation
        if param.filters[i].lambda_c < 3.e4:
            from astropy.modeling import models, fitting
            p_init = models.Polynomial2D(degree=8)
            fit_p = fitting.LevMarLSQFitter()
            #fit_p = fitting.LinearLSQFitter()
            p = fit_p(p_init, c['x'][ok], c['y'][ok], ratio[ok]) #, weights=signoise[i,ok]**2)
            dmag = -2.5*np.log10((fnu/p(c['x'], c['y'])/obs_sed)[i,:])
        
            ok = np.isfinite(dmag) & (signoise[i,:] > 3) & (np.abs(dmag) < 0.1)
            if apply:
                ok = c[cnames[i]] > -90
                c[cnames[i]][ok] /= p(c['x'], c['y'])[ok]
                c[cnames[i].replace('f_', 'e_')][ok] /= p(c['x'], c['y'])[ok]
            
            fitters[cnames[i]] = p
            translate[cnames[i]] = param.filters[i].fnumber
            
        ax.scatter(c['ra'], c['dec'], c='black', vmin=-0.08, vmax=0.08, alpha=0.1, s=1, marker='.', edgecolor='None')
        ax.scatter(c['ra'][ok], c['dec'][ok], c=dmag[ok], vmin=-0.08, vmax=0.08, alpha=0.2, s=10, marker='s', edgecolor='None')
        sc = ax.scatter(c['ra'][ok][0], c['dec'][ok][0], c=dmag[ok][0], vmin=-0.08, vmax=0.08, alpha=1, s=10, marker='s', edgecolor='None')
        #ax.set_xlabel('RA'); ax.set_ylabel('Dec')
        ax.text(0.5, 0.95, '%s\n(%d)' %(param.filters[i].name, param.filters[i].fnumber), ha='center', va='top', fontsize=8, transform=ax.transAxes, color='red')
        
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.set_xlim(ax.get_xlim()[::-1])
    
    if apply:
        c.write('%s.spatial' %(c.filename), format='ascii.commented_header')
                
    fig.tight_layout(pad=0.1)
    cb = plt.colorbar(sc)
    cb.set_label(r'$\Delta$mag')
    fig.savefig('%s.spatial.png' %(c.filename))
    plt.close()
    plt.ion()
    
    import pickle
    fp = open('%s.spatial.pkl' %(c.filename),'wb')
    pickle.dump(fitters, fp)
    pickle.dump(translate, fp)
    fp.close()
   
def apply_spatial_offsets(cat='cdfs.v1.3.1.nzp.cat', pfile='cdfs.v1.3.1.nzp.cat.AB25.spatial.pkl'):
    """
    Apply 2D "zeropoint" residuals to the photometric catalogs
    """
    import pickle
    import astropy.table
    
    #### Read the catalog
    c = astropy.table.Table.read(cat, format='ascii.commented_header')
    
    #### Read the 2D scalings
    fp = open(pfile,'rb')
    fitter = pickle.load(fp)
    translate = pickle.load(fp)
    fp.close()
    
    #### Fitter is a dictionary with keys equal to the column headings
    #### and whose values are astropy.modeling.Polynomial2D models with 
    #### coordinates 'x' and 'y' from the catalog
    for col in fitter.keys():
        print col, col in c.columns
        ok = c[col] > -90
        c[col][ok] *= 1./fitter[col](c['x'], c['y'])[ok]
        c[col.replace('f_', 'e_')][ok] *= 1./fitter[col](c['x'], c['y'])[ok]
        
        c[col].format='%.5e'
        c[col.replace('f_', 'e_')].format='%.5e'
        
    c.write('%s.spatial' %(cat), format='ascii.commented_header') 

def template_contributions(root='photz', wave=15418.99, PATH='OUTPUT', CACHE_FILE='Same'):
    """
    Get template contributions to the full fit, normalized at the filter closest to the
    specified rest-frame wavelength (default is F160W)
    """
    from collections import OrderedDict
    tempfilt, coeffs, temp_sed, pz = readEazyBinary(MAIN_OUTPUT_FILE='%s' %(root), OUTPUT_DIRECTORY=PATH, CACHE_FILE=CACHE_FILE)
    
    dl = np.abs(tempfilt['lc']-wave)
    ix = np.argmin(dl)
    
    renorm = tempfilt['tempfilt'][ix,:,0]    
    coeffs_h = coeffs['coeffs']*np.dot(renorm.reshape(-1,1), np.ones((1,coeffs['NOBJ'])))
    coeffs_total = np.sum(coeffs_h, axis=0)
    coeffs_h /= coeffs_total
        
    out = OrderedDict()
    out['quiescent'] = coeffs_h[[0,3,4,6], :].sum(axis=0)
    
    if tempfilt['NTEMP'] > 7:
        out['dusty_old'] = coeffs_h[7, :]
    else:
        out['dusty_old'] = out['quiescent']*0.
        
    out['dusty_young'] = coeffs_h[5, :]
    out['dusty_total'] = out['dusty_young'] + out['dusty_old']
    out['young'] = coeffs_h[1, :]
    out['mid'] = coeffs_h[2, :]
    
    ### Additional Erb 2010 template
    if tempfilt['NTEMP'] == 9:
        out['young'] += coeffs_h[8, :]
    
    return out
    
    
def anneal_pz(zgrid, pz, level=0.68):
    """
    "anneal" method for computing confidence intervals.
    
    Start at peak probability and step down in probability density until you get the 
    desired integrated probability
    
    """
    
    pznorm = pz/np.trapz(pz, zgrid)
    zmax = np.argmax(pznorm)
    dz = np.diff(zgrid)
    pzmax = pznorm.max()
    probs = np.sort(np.unique(pznorm))[::-1]
    out = np.zeros((len(pz)-1, len(probs)))
    for i in range(len(probs)):
        ok = pznorm[1:] >= probs[i]
        out[:,i] = pznorm[1:]*dz*ok
        
    ix_level = int(np.round(np.interp(level, out.sum(axis=0), np.arange(len(probs)))))
    z_level = zgrid[1:]*(out[:,ix_level] > 0) 
    
    return out, z_level, np.min(z_level[z_level > 0]), np.max(z_level[z_level > 0])