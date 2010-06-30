"""
SExtractor class from astropysics.

Modified by GBB to read the output from "sex -dd" and "sex -dp" from files, rather
from a piped output because the "pparm" line dind't run for some reason.
"""

import numpy as np

def _get_package_data(dataname):
    """
    Use this function to load data files distributed with astropysics in the 
    astropysics/data directory
    
    dataname is the file name of a file in the data directory, and a string
    with the contents of the file will be returned
    """
    from . import __name__ as rootname
    from . import __file__ as rootfile
    from pkgutil import get_loader
    from os.path import dirname
    path = dirname(rootfile)+'/data/'+dataname
    return get_loader(rootname).get_data(path)

class SExtractor(object):
    """
    This class is an adaptor to the Sextractor program 
    (Bertin & Arnouts 1996, http://astromatic.iap.fr/software/sextractor/).
    Sextractor must be installed and on the path for this class to function.
    
    options are set by changing values in the options dictionary
    
    output parameters are chosen by setting True/False values in the 
    params dictionary

    """
    
    @staticmethod
    def _getSexDefaults():
        from subprocess import Popen,PIPE
        
        optinfo = {}
        opts = {}
        optorder = []
        parinfo = {}
        parorder = []
        
        try:
            pconf = Popen('sex -dd'.split(),executable='sex',stdout=PIPE,stderr=PIPE)
            pparm = Popen('sex -dp'.split(),executable='sex',stdout=PIPE,stderr=PIPE)
            pconf.wait()
            pparm.wait()
            confstr = pconf.communicate()[0]
            parmstr = pparm.communicate()[0]
        except OSError:
            raise OSError('Sextractor not found on system path')
        
        comm = ''
        newk = k = None
        for l in confstr.split('\n'):
            commenti = l.find('#')
            if commenti>0:
                ls = l[:commenti].split()
                if len(ls)>1:
                    newk = ls[0]
                    newoptval = ' '.join(ls[1:])
                elif len(ls)>0:
                    newk = ls[0]
                    newoptval = ''
                newcomm = l[commenti+1:].strip()
            else:
                newcomm = ''
                
            if newk:
                if k:
                    opts[k] = optval
                    optinfo[k] = comm
                    optorder.append(k)
                k = newk
                optval = newoptval
                
                newk = None
                comm = ''
            comm+=newcomm
              
        for l in parmstr.split('\n'):
            ls = l.split()
            if len(ls) > 1:
                k = ls[0].strip().replace('#','')
                unit = ls[-1].replace('[','').replace(']','') if '[' in ls[-1] else None
                info = ' '.join(ls[1:(-1 if unit else None)])
                parinfo[k] = (info,unit if unit else '')
                parorder.append(k)
        
        SExtractor._optinfo = optinfo
        SExtractor._defaultopts = opts
        SExtractor._optorder = optorder #TODO:OrderedDict for 2.7
        SExtractor._parinfo = parinfo
        SExtractor._parorder = parorder #TODO:OrderedDict for 2.7
    
    @staticmethod   
    def getOptInfo(aslist=False):
        """
        returns the  dictionary of input options and the associated information
        
        if aslist is True, returns an list of 
        """
        if aslist:
            return [(k,SExtractor._optinfo[k]) for k in SExtractor._optorder]
        else:
            return dict(SExtractor._optinfo)
    
    @staticmethod    
    def getParamInfo(aslist=False):
        """
        returns the dictionary of parameters and the associated information
        """
        if aslist:
            return [(k,SExtractor._parinfo[k]) for k in SExtractor._parorder]
        else:
            return dict(SExtractor._parinfo)
        
        
    def __init__(self,sexfile=None,parfile=None):
        
        opts = dict(SExtractor._defaultopts)
        pars = dict([(k,False) for k in  SExtractor._parinfo])
                    
        if sexfile:
            #with open(sexfile) as f:
            f =  open(sexfile)
            for l in f:
                commenti = l.find('#')
                if l > -1:
                    l = l[:commenti]
                ls = l.split()
                if len(ls) > 1:
                    k = ls[0].strip()
                    if k not in opts:
                        raise ValueError('sexfile has invalid option %s'%k)
                    opts[k] = ls[1].strip()
            self.name = sexfile.replace('.sex','')
        else:
            self.name = 'astropysics_auto'
                        
                        
        pfile = opts['PARAMETERS_NAME'] if parfile is None else parfile
        if pfile != 'default.param':
            f = open(pfile)
            for l in f:
                if not (l.strip().startswith('#') or l.strip()==''):
                    k = l.split()[0].strip()
                    if k not in pars:
                        raise ValueError('param file has invalid parameter %s'%k)
                    pars[k] = True
        if not np.any(pars.values()):
            #if no outputs provided, use defaults (from v 2.8.6)
            defs286 =['NUMBER','FLUX_ISO', 'FLUXERR_ISO', 'FLUX_AUTO', 'FLUXERR_AUTO', 'X_IMAGE', 'Y_IMAGE', 'FLAGS']
            for p in defs286:
                if p in pars:
                    pars[p] = True
        
        self.options = opts
        self.params = pars
        
        self.overwrite = False
        
            
    def getParamList(self):
        """
        returns a list of all selected parameters
        """    
        return [k for k in self._parorder if self.params[k]]
    
    def _saveFiles(self,fnbase):
        import os
        
        fnbase = fnbase.replace('.sex','')
        
        dir = os.path.split(fnbase)[0]
        dir = '' if dir=='' else dir+os.sep
        
        self.options['PARAMETERS_NAME'] = dir+fnbase+'.param'
        
        ostr = self._makeOptionStr()
        pstr = self._makeParamStr()
        
        # with open(fnbase+'.sex','w') as f:
        #     f.write(ostr)
        f = open(fnbase+'.sex','w')
        f.write(ostr)
            
        # with open(self.options['PARAMETERS_NAME'],'w') as f:
        #     f.write(pstr)
        f = open(self.options['PARAMETERS_NAME'],'w')
        f.write(pstr)
                
    
    def _makeOptionStr(self):
        ostr = [o+'\t'+str(self.options[o])+'\t # '+SExtractor._optinfo[o] for o in self._optorder]
        ostr = '\n'.join(ostr)
        return ostr
        
    def _makeParamStr(self):
        pstr = [p+'\t # '+str(SExtractor._parinfo[p]) for p in self._parorder if self.params[p]]
        pstr = '\n'.join(pstr)
        return pstr
    
    
    def getOptionList(self,incval=False):
        """
        returns a list of all options.  If incval is True, returns a list
        of tuples of the form (optionname,optionvalue)
        """
        if incval:
            return [(k,self.options[k]) for k in self._optorder]
        else:
            return [k for k in self._optorder]
        
    def sextractImage(self,detimfn,analysisimfn=None,mode='waiterror'):
        """
        writes configuration files and runs sextractor on the input image
        
        mode can be:
        
        * 'waiterror': waits for sextractor to finish, and raises an 
          SExtractorError if it does not complete sucessfully. stdout 
          and sterr are saved to self.lastout and self.lasterr (returns 0)
        * 'wait': waits for sextractor to finish and returns the return code
          stdout and sterr are saved to self.lastout and self.lasterr
        * 'proc': stars the processes but does not wait - returns the Popen 
          instance of the processes
        """
        from subprocess import Popen,PIPE
        from os.path import exists
        
        fnbase = self.name
        if not self.overwrite:
            fnbase = fnbase.replace('.sex','')
            if exists(fnbase+'.sex'):
                fns = fnbase.split('-')
                try:
                    i = int(fns[-1])
                    i+=1
                except ValueError:
                    i = 2
                if len(fns)<2:
                    fns.append(str(i))
                else:
                    fns[-1] = str(i)
                fnbase = '-'.join(fns)
            self.name = fnbase
                
        self._saveFiles(fnbase)
        if analysisimfn:
            #clstr = 'sex {0} {1} -c {2}'.format(detimfn,analysisimfn,self.name+'.sex')
            clstr = 'sex %s %s -c %s' %(detimfn,analysisimfn,self.name+'.sex')
        else:
            # clstr = 'sex {0} -c {1}'.format(detimfn,self.name+'.sex')
            clstr = 'sex %s -c %s' %(detimfn,self.name+'.sex')
        proc = Popen(clstr.split(),executable='sex',stdout=PIPE,stderr=PIPE)
        
        if mode == 'waiterror' or mode =='wait':
            res = proc.wait()
            sout,serr = proc.communicate()
            
            self.lastout = sout
            self.lasterr = serr
            
            if res!=0 and mode == 'waiterror' :
                raise SExtractorError(serr,sout)
            return res
        elif mode == 'proc':
            return proc
        else:
            raise ValueError('unrecognized mode argument '+str(mode))
# try:
#     SExtractor._getSexDefaults()
# except OSError:
#     from warnings import warn
#     warn('SExtractor not found on system - phot.SExtractor class will not function')
    
class SExtractorError(Exception):
    def __init__(self,*args):
        super(SExtractorError,self).__init__(*args)
