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

import numpy as np

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib import rc
import pylab

import threedhst.catIO as catIO

class FilterDefinition:
    def __init__(self):
        """
        Placeholder for the filter definition information.
        """
        self.name = None
        self.wavelength = None
        self.transmission = None
        
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
        
    ###### Done.    
    return tempfilt, coeffs, temp_sed, pz

        
def getEazySED(idx, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same', scale_flambda=True, verbose=False):
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
        
    fobs = tempfilt['fnu'][:,idx]/lci**2*flam_factor
    efobs = tempfilt['efnu'][:,idx]/lci**2*flam_factor
    
    ##### Broad-band SED
    obs_sed = np.dot(tempfilt['tempfilt'][:,:,coeffs['izbest'][idx]],\
                     coeffs['coeffs'][:,idx])/(lci)**2*flam_factor
    
    zi = tempfilt['zgrid'][coeffs['izbest'][idx]]
    
    ###### Full template SED, observed frame
    lambdaz = temp_seds['templam']*(1+zi)
    temp_sed = np.dot(temp_seds['temp_seds'],coeffs['coeffs'][:,idx])
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
    
def getEazyPz(idx, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same'):
    """
zgrid, pz = getEazyPz(idx, \
                      MAIN_OUTPUT_FILE='photz', \
                      OUTPUT_DIRECTORY='./OUTPUT', \
                      CACHE_FILE='Same')
                      
    Get Eazy p(z) for object #idx.
    """
    tempfilt, coeffs, temp_seds, pz = readEazyBinary(MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                                                    OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                                                    CACHE_FILE = CACHE_FILE)
    
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
    pzi = np.exp(-0.5*(pz['chi2fit'][:,idx]-min(pz['chi2fit'][:,idx])))*prior
    if np.sum(pzi) > 0:
        pzi/=np.trapz(pzi, tempfilt['zgrid'])
    
    ###### Done
    return tempfilt['zgrid'],pzi
    
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
    
def plotExampleSED(idx=20, writePNG=True, MAIN_OUTPUT_FILE = 'photz', OUTPUT_DIRECTORY = 'OUTPUT', CACHE_FILE = 'Same', lrange=[3000,8.e4], axes=None):
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
                          CACHE_FILE = CACHE_FILE)
    
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
    plotsize=20
    alph=0.7
    
    #### Full best-fit template
    if axes is None:
        ax = fig.add_subplot(121)
    else:
        ax = axes[0]
        
    ax.plot(lambdaz, temp_sed, linewidth=1.0, color='blue',alpha=alph)
    
    #### template fluxes integrated through the filters
    ax.scatter(lci,obs_sed,
               c='red',marker='o',s=plotsize,alpha=alph)

    #### Observed fluxes w/ errors
    ax.errorbar(lci,fobs,yerr=efobs,ecolor=None,
               color='black',fmt='o',alpha=alph)
    
    #### Set axis range and titles
    ax.semilogx()
    ax.set_xlim(lrange[0],lrange[1])
    ax.set_ylim(-0.05*max(obs_sed),1.1*max(obs_sed))
    ax.set_xlabel(r'$\lambda$ [$\AA$]')
    ax.set_ylabel(r'$f_\lambda$')
    
    ##### P(z)
    if pz is not None:
        if axes is None:
            axp = fig.add_subplot(122)
        else:
            axp = axes[1]
            
        axp.plot(zgrid, pz, linewidth=1.0, color='orange',alpha=alph)
        axp.fill_between(zgrid,pz,np.zeros(zgrid.size),color='yellow')

        if zout.z_spec[qz[idx]] > 0:
            axp.plot(zout.z_spec[qz[idx]]*np.ones(2), np.array([0,1e6]),color='red',alpha=0.4)

        #### Set axis range and titles
        axp.set_xlim(0,np.ceil(np.max(zgrid)))
        axp.set_ylim(0,1.1*max(pz))
        axp.set_xlabel(r'$z$')
        axp.set_ylabel(r'$p(z)$')
        
    if writePNG & (axes is None):
        fig.savefig('/tmp/test.pdf',dpi=100)

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
    
def show_fit_residuals(root='photz_v1.7.fullz', PATH='./OUTPUT/', savefig=None, adjust_zeropoints='zphot.zeropoint', fix_filter=None, ref_filter=None):
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
    
    tempfilt['fnu'] *= zpfactors
    tempfilt['efnu'] *= zpfactors
    
    obs_sed = np.zeros((tempfilt['NFILT'], tempfilt['NOBJ']), dtype=np.float)
    for i in xrange(tempfilt['NOBJ']):
        obs_sed[:,i] = np.dot(tempfilt['tempfilt'][:,:,coeffs['izbest'][i]], coeffs['coeffs'][:,i])
    
    zi = tempfilt['zgrid'][coeffs['izbest']]
    lc = tempfilt['lc']
    offsets = lc*0.
    
    resid = (obs_sed-tempfilt['fnu']) / obs_sed + 1
    signoise = tempfilt['fnu']/np.sqrt(tempfilt['efnu']**2+(0.01*tempfilt['fnu'])**2)
    
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
        keep &= tempfilt['efnu'][i,:] > 0
        if signoise[i,:].max() > 3:
            keep & (signoise[i,:] > 3)
        #keep &= (signoise[i,:] > 3) #| (signoise[i,:] < 0.0001)#& (np.abs(resid[i,:]-1)/tempfilt['efnu'][i,:] < 5)
    
    for i in np.argsort(lc):
        #keep = signoise[i,:] > 3
        if np.std(zi[keep]) == 0:
            rnd = np.random.normal(size=keep.sum())*0.01*lc[i]
        else:
            rnd = 0.
        #
        sc = ax.plot(lc[i]/(1+zi[keep])+rnd, resid[i,keep], marker='.', alpha=0.05, linestyle='None', color=colors[i])
        #xm, ym, ys, nn = threedhst.utils.runmed(lc[i]/(1+zi[keep]), resid[i,keep], NBIN=int(len(keep)/1000.))
        xm, ym, ys, nn = threedhst.utils.runmed(lc[i]/(1+zi[keep]), resid[i,keep], NBIN=np.maximum(int(keep.sum()/1000.), 8))
        xx.append(xm)
        yy.append(ym)
        ss.append(ys)
        val = np.sum(resid[i,keep]*signoise[i,keep]**2)/np.sum(signoise[i,keep]**2)
        stats[i] = {'mean':np.mean(resid[i,keep]),
                      'median':np.median(resid[i,keep]),
                      'std':np.std(resid[i,keep]),
                      'stdmean':np.std(resid[i,keep])/np.sqrt(keep.sum()),
                      'p':np.percentile(resid[i,keep], [2.5,16,50,84,97.5]),
                      'pstd':(np.percentile(resid[i,keep], 84)-np.percentile(resid[i,keep], 16))/2/np.sqrt(keep.sum()),
                      'val':val}
        #
        # offsets[i] = stats[i]['median']
        # #
        # print '%s %.3f %.3f %.3f %.3f %d %.3f' %(param.filters[i].name, stats[i]['median'], stats[i]['pstd'], stats[i]['p'][1], stats[i]['p'][-2], keep.sum(), val)
        
    for ci, i in enumerate(np.argsort(lc)):
        ax.plot(xx[ci], yy[ci], alpha=0.2, color='black', linewidth=4)
        ax.plot(xx[ci], yy[ci], alpha=0.8, color=colors[i], linewidth=2)
    
    #### Adjustments to zeropoints
    lcfull = []
    residfull = []
    for i in np.argsort(lc):
        #keep = signoise[i,:] > 3
        #keep = signoise[i,:] > 0
        lcfull.extend(lc[i]/(1+zi[keep]))
        residfull.extend(resid[i,keep]/stats[i]['median'])
        
    xmfull, ymfull, ysfull, nnfull = threedhst.utils.runmed(np.array(lcfull), np.array(residfull), NBIN=np.maximum(int(len(residfull)/2000.), 10))
    #ymfull[xmfull > 3.e4] = 1.
    #ymfull[xmfull < 1200] = 1.
    
    ax.plot(xmfull, ymfull, color='black', alpha=0.75, linewidth=2)
    
    if os.path.exists('tweak/tweak.dat'):
        tx, ty = np.loadtxt('tweak/tweak.dat', unpack=True)
        ty_int = np.interp(xmfull, tx, ty, left=1, right=1)
    else:
        ty_int = xmfull*0.+1
    
    ax.plot(xmfull, ymfull*ty_int, color='black', alpha=0.3, linewidth=2)
    np.savetxt('tweak/tweak.dat', np.array([xmfull, ymfull*ty_int]).T, fmt='%.5e')
    
    NT = len(param.templates)
    fp = open('tweak/spectra.param','w')
    for i in range(NT):
        wt, ft = np.loadtxt(param.templates[i], unpack=True)
        ft /= np.interp(wt, xmfull, ymfull, left=1, right=1)
        np.savetxt('tweak/%s' %(os.path.basename(param.templates[i])), np.array([wt, ft]).T, fmt='%.5e')
        fp.write('%d tweak/%s 1.0 0 1.0\n' %(i+1, os.path.basename(param.templates[i])))
    #
    fp.close()
    
    ### account for overall wiggles
    for i in np.argsort(lc):
        lcz = lc[i]/(1+zi)
        yint = np.interp(lcz, xmfull, ymfull, left=-99, right=-99)
        ok = keep & (yint > 0)
        rfix = resid[i,:]/yint
        val = np.sum(rfix[ok]*signoise[i,ok]**2)/np.sum(signoise[i,ok]**2)
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
        print '%s %.3f %.3f %.3f %.3f %d %.3f' %(param.filters[i].name, stats[i]['median'], stats[i]['pstd'], stats[i]['p'][1], stats[i]['p'][-2], keep.sum(), val)
        
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
        ax.text(0.99, 0.99-0.04*ci*28./NF, '%s %.3f' %(os.path.basename(param.filters[i].name).replace('_','_'), offsets[i]/ref_offset), transform = ax.transAxes, color='0.5', horizontalalignment='right', va='top', fontsize=font)
        ax.text(0.99, 0.99-0.04*ci*28./NF, '%s %.3f' %(os.path.basename(param.filters[i].name).replace('_','_'), offsets[i]/ref_offset), transform = ax.transAxes, color=colors[i], alpha=0.8, horizontalalignment='right', va='top', fontsize=font)

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
    yh, xh = np.histogram(np.log10(1+zout.z_peak), range=[np.log10(1), np.log10(5)], bins=100)
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
    
    