"""
eazyPy: routines for reading and plotting Eazy output

    EazyParam
    readEazyBinary
    getEazySED
    trapz
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
        self.transmision = None
        
class FilterFile:
    def __init__(self, file='FILTER.RES.v8.R300'):
        """
        Read a filter file.
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
    
    def names(self):
        """
        Print the filter names.
        """
        for i in range(len(self.filters)):
            print '%5d %s' %(i+1, self.filters[i].name)
    
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
        if verbose:
            print 'Wrote <%s>.' %(file)
            
class ParamFilter:
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
    def __init__(self, PARAM_FILE='zphot.param'):
        self.filename = PARAM_FILE
        
        f = open(PARAM_FILE,'r')
        self.lines = f.readlines()
        f.close()
        
        self._process_params()
        
        filters = []
        for line in self.lines:
            if line.startswith('#  Filter'):
                filters.append(ParamFilter(line))
        
        self.NFILT = len(filters)
        self.filters = filters
        
    def _process_params(self):
        params = {}
        formats = {}
        for line in self.lines:
            if line.startswith('#') is False:
                lsplit = line.split()
                if lsplit.__len__() >= 2:
                    params[lsplit[0]] = lsplit[1]
                    try:
                        flt = float(lsplit[1])
                        formats[lsplit[0]] = 'f'
                        params[lsplit[0]] = flt
                    except:
                        formats[lsplit[0]] = 's'
                    
        self.params = params
        self.param_names = params.keys()
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
                str = '%-25s %'+self.formats[param]+'\n'
                fp.write(str %(param, self.params[param]))
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

        
def getEazySED(idx, MAIN_OUTPUT_FILE='photz', OUTPUT_DIRECTORY='./OUTPUT', CACHE_FILE='Same'):
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
    tempfilt, coeffs, temp_seds, pz = readEazyBinary(MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                                                    OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                                                    CACHE_FILE = CACHE_FILE)
    
    ##### Apply zeropoint factors
    zpfile = OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zeropoint'
    if os.path.exists(zpfile):
        zpfilts,zpf = np.genfromtxt(zpfile, unpack=True)                                    
    else:
        zpf = np.ones(tempfilt['NFILT'])

    zpfactors = np.dot(zpf.reshape(tempfilt['NFILT'],1),\
                       np.ones(tempfilt['NOBJ']).reshape(1,tempfilt['NOBJ']))
    
    tempfilt['fnu'] *= zpfactors
    tempfilt['efnu'] *= zpfactors
    
    lci = tempfilt['lc'].copy()
    
    params = EazyParam(PARAM_FILE=OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.param')
    abzp = np.float(params['PRIOR_ABZP'])
    
    ##### Broad-band SED
    obs_sed = np.dot(tempfilt['tempfilt'][:,:,coeffs['izbest'][idx]],\
                     coeffs['coeffs'][:,idx])/(lci/5500.)**2
                     
    fobs = tempfilt['fnu'][:,idx]/(lci/5500.)**2
    efobs = tempfilt['efnu'][:,idx]/(lci/5500.)**2
    
    zi = tempfilt['zgrid'][coeffs['izbest'][idx]]
    
    ###### Full template SED, observed frame
    lambdaz = temp_seds['templam']*(1+zi)
    temp_sed = np.dot(temp_seds['temp_seds'],coeffs['coeffs'][:,idx])
    temp_sed /= (1+zi)**2
    
    ###### IGM absorption
    lim1 = np.where(temp_seds['templam'] < 912)
    lim2 = np.where((temp_seds['templam'] >= 912) & (temp_seds['templam'] < 1026))
    lim3 = np.where((temp_seds['templam'] >= 1026) & (temp_seds['templam'] < 1216))
    
    if lim1[0].size > 0: temp_sed[lim1] *= 0.
    if lim2[0].size > 0: temp_sed[lim2] *= 1.-temp_seds['db'][coeffs['izbest'][idx]]
    if lim3[0].size > 0: temp_sed[lim3] *= 1.-temp_seds['da'][coeffs['izbest'][idx]]
        
    ###### Done
    return lambdaz, temp_sed, lci, obs_sed, fobs, efobs

def trapz(x, y):
    """
result = trapz(x, y)

    Integrate y(x) with the trapezoid rule.
    """
    if x.size != y.size:
       print ('N(x) != N(y)')
       return None   
    h = x[1:]-x[:-1]
    trap = y[1:]+y[:-1]
    return np.sum(h*trap/2.)

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
        pzi/=trapz(tempfilt['zgrid'],pzi)
    
    ###### Done
    return tempfilt['zgrid'],pzi
    
def plotExampleSED(idx=20, writePNG=True,
    MAIN_OUTPUT_FILE = 'photz',
    OUTPUT_DIRECTORY = 'OUTPUT',
    CACHE_FILE = 'Same'):
    """
PlotSEDExample(idx=20)

    Plot an example Eazy best-fit SED.
    """

    zout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
    #qz = np.where(zout.z_spec > 0)[0]
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
    plt.rcParams['text.usetex'] = True
    #plt.rcParams['text.latex.preamble'] = ''

    ##### start plot
    fig = plt.figure(figsize=[8,4],dpi=100)
    fig.subplots_adjust(wspace=0.15,hspace=0.0,left=0.08,bottom=0.15,right=0.98,top=0.99)
    
    #### Plot parameters
    plotsize=20
    alph=0.7
    
    #### Full best-fit template
    ax = fig.add_subplot(121)
    ax.plot(lambdaz, temp_sed, linewidth=1.0, color='blue',alpha=alph)
    
    #### template fluxes integrated through the filters
    ax.scatter(lci,obs_sed,
               c='red',marker='o',s=plotsize,alpha=alph)

    #### Observed fluxes w/ errors
    ax.errorbar(lci,fobs,yerr=efobs,ecolor=None,
               color='black',fmt='o',alpha=alph)
    
    #### Set axis range and titles
    ax.semilogx()
    ax.set_xlim(3000,9.e4)
    ax.set_ylim(-0.05*max(obs_sed),1.1*max(obs_sed))
    plt.xlabel(r'$\lambda$ [$\AA$]')
    plt.ylabel(r'$f_\lambda$')
    
    ##### P(z)
    if pz is not None:
        ax = fig.add_subplot(122)
        ax.plot(zgrid, pz, linewidth=1.0, color='orange',alpha=alph)
        ax.fill_between(zgrid,pz,np.zeros(zgrid.size),color='yellow')

        if zout.z_spec[qz[idx]] > 0:
            ax.plot(zout.z_spec[qz[idx]]*np.ones(2),np.array([0,1e6]),color='red',alpha=0.4)

        #### Set axis range and titles
        ax.set_xlim(0,np.ceil(np.max(zgrid)))
        ax.set_ylim(0,1.1*max(pz))
        plt.xlabel(r'$z$')
        plt.ylabel(r'$p(z)$')
        
    if writePNG:
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
    
    
    