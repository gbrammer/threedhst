"""
SExtractor class from astropysics, written by Erik Tollerud.
http://packages.python.org/Astropysics/

Modified by GBB to read the output from "sex -dd" and "sex -dp" from files, 
rather than from a piped output because the "pparm" line didn't run for some 
reason.

Also added similar SWarp class wrapper around SWarp.
"""

import numpy as np

#### This should probably go in utils.py
def _get_package_data(dataname):
    """
    (taken from astropysics.io)
    Use this function to load data files distributed with astropysics in the 
    astropysics/data directory
    
    dataname is the file name of a file in the data directory, and a string
    with the contents of the file will be returned
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
        return fp.read()


class SExtractor(object):
    """
    This class is an adaptor to the Sextractor program 
    (Bertin & Arnouts 1996, http://astromatic.iap.fr/software/sextractor/).
    Sextractor must be installed and on the path for this class to function.
    
    options are set by changing values in the options dictionary
    
    output parameters are chosen by setting True/False values in the 
    params dictionary
    
    Workflow is something like
    
        s = sex.SExtractor([sexfile=''])
        s.aXeParams()
        s.sextractImage('image.fits')
        
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
            pconf.wait()
            confstr = pconf.communicate()[0]
            # pparm = Popen('sex -dp'.split(),executable='sex',stdout=PIPE,stderr=PIPE)
            # pparm.wait()
            # parmstr = pparm.communicate()[0]
            parmstr = _get_package_data('sexdp') #gbb
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
    
    def copyConvFile(self):
        """
        copyConvFile()
        
        Copy default.conv from threedhst/data to ./
        """
        self.conv = _get_package_data('default.conv')
        fp = open('default.conv','w')
        fp.write(self.conv)
        fp.close()
        self.conv = self.conv.split('\n')
        
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
        from warnings import warn
        
        self._getSexDefaults()
        opts = dict(SExtractor._defaultopts)
        pars = dict([(k,False) for k in  SExtractor._parinfo])
            
        if sexfile:
            #with open(sexfile) as f:
            f =  open(sexfile) #gbb
            for l in f:
                commenti = l.find('#')
                if commenti > -1:
                    l = l[:commenti]
                ls = l.split()
                if len(ls) > 1:
                    k = ls[0].strip()
                    if k not in opts:
                        # raise ValueError('sexfile has invalid option %s'%k) 
                        warn('sexfile \'%s\' has invalid option %s' %(sexfile,k)) #gbb
                    # opts[k] = ls[1].strip()
                    if len(ls) > 2:
                        #print ls
                        opts[k] = ' '.join(ls[1:]) #.strip() #gbb
                    else:
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
    
    
    def aXeParams(self):
        """
        Set the columns needed for input to aXe.
        (gbb)
        """
        for k in self._parorder:
            self.params[k] = False
        useParams = _get_package_data('aXe.param').split('\n')[:-1]
        for par in useParams:
            self.params[par] = True
    
    
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
        f = open(fnbase+'.sex','w') #gbb
        f.write(ostr)
            
        # with open(self.options['PARAMETERS_NAME'],'w') as f:
        #     f.write(pstr)
        f = open(self.options['PARAMETERS_NAME'],'w') #gbb
        f.write(pstr)
                
    
    
    def _makeOptionStr(self, maxlen=14):
        #### get longest number of characters of parameter names
        parlen = 0
        for o in self._optorder:
            parlen = np.max([parlen,o.__len__()])
        ostr = []
        for o in self._optorder:
            optlen = np.max([maxlen,self.options[o].__len__()])
            fmt = '%-'+str(parlen)+'s  %-'+str(optlen)+'s  # %s'
            ostr.append(fmt %(o, self.options[o], self._optinfo[o]))
        ostr = '\n'.join(ostr)
        return ostr    
    
    
    def _makeParamStr(self):
        parlen = 0
        for p in self._parorder:
            if self.params[p]:
                parlen = np.max([parlen,p.__len__()])
        fmt = '%-'+str(parlen)+'s   #  %s'
        
        pstr = []
        for p in self._parorder:
            if self.params[p]:
                pstr.append(fmt %(p,self._parinfo[p]))        
        # pstr = [p+'\t # '+str(SExtractor._parinfo[p]) for p in self._parorder if self.params[p]]
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
    
    
    def sextractImage(self,detectionImage,analysisImage=None,mode='waiterror'):
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
        if analysisImage:
            #clstr = 'sex {0} {1} -c {2}'.format(detectionImage,analysisImage,self.name+'.sex')
            clstr = 'sex %s %s -c %s' %(detectionImage,analysisImage,self.name+'.sex')
        else:
            # clstr = 'sex {0} -c {1}'.format(detectionImage,self.name+'.sex')
            clstr = 'sex %s -c %s' %(detectionImage,self.name+'.sex')
        proc = Popen(clstr.split(),executable='sex',stdout=PIPE,stderr=PIPE)
        
        print 'THREEDHST/sex: %s' %clstr
        
        if mode == 'waiterror' or mode =='wait':
            res = proc.wait()
            sout,serr = proc.communicate()
            
            self.lastout = sout
            self.lasterr = serr
            
            if res!=0 and mode == 'waiterror' :
                raise SError(serr,sout)
            return res
        elif mode == 'proc':
            return proc
        else:
            raise ValueError('unrecognized mode argument '+str(mode))


class SError(Exception):
    def __init__(self,*args):
        super(SError,self).__init__(*args)
    

def sexcatRegions(sexcat, regfile, format=1):
	"""
sexcat_regions(sexcat, regfile, format=1)
	
Make DS9 region file from SExtractor catalog.  The coordinate system 
is determined by the format argument, with
	
format = 1
	image coordinates x,y (X_IMAGE, Y_IMAGE)
format = 2
	world coordinates ra,dec (X_WORLD, Y_WORLD)
		
If A, B, THETA columns are present, will make elliptical regions.
	"""
	import os,sys
	import aXe2html.sexcat.sextractcat
    
    #sexcat = 'F140W_SCI.cat'
	#regfile = 'test.reg'
	#format = 1
	if os.access(sexcat, os.R_OK) is False:
		print "SExtractor catalog, %s, not found." %(sexcat)
		return False
	cat = aXe2html.sexcat.sextractcat.SexCat(sexcat)
	## Force format=1 if out of range
	if format < 1 or format > 2:
		format = 1
	if format == 1:
		header = 'image'
		ext = '_IMAGE'
		asec = 1.
		pp = ''
		theta_sign = 1
	else:
		header = 'fk5'
		ext = '_WORLD'
		asec = 3600.
		pp = '"'
		theta_sign = -1
	useEllipse = (cat.searchcol('A'+ext) > -1) and \
				 (cat.searchcol('B'+ext) > -1) and \
				 (cat.searchcol('THETA'+ext) > -1)
	## X,Y columns
	x_col = cat.columns[cat.searchcol('X'+ext)].entry
	y_col = cat.columns[cat.searchcol('Y'+ext)].entry
	## NUMBER column
	num_col = cat.searchcol('NUMBER')
	if num_col > -1:
		num = cat.columns[num_col].entry
	else:
		num = srange(1,cat.nrows+1)
	## Write output file
	fp = open(regfile,'w')
	if useEllipse:
		#print "Ellipse"
		fp.write(header+'\n')
		a_col = cat.columns[cat.searchcol('A'+ext)].entry
		b_col = cat.columns[cat.searchcol('B'+ext)].entry
		theta_col = cat.columns[cat.searchcol('THETA'+ext)].entry
		for i in range(cat.nrows):
			line = "ellipse(%s, %s, %6.2f%s, %6.2f%s, %6.2f) # text={%s}\n" \
				%(x_col[i], y_col[i], float(a_col[i])*asec, pp, \
				  float(b_col[i])*asec, pp, float(theta_col[i])*theta_sign, \
				  str(num[i]))
			fp.write(line)
	else:
		#print "Circle"
		fp.write(header+'\n')
		for i in range(cat.nrows):
			line = "circle(%s, %s, 1\") # text={%s}\n" \
					   %(x_col[i],y_col[i],str(num[i]))
			fp.write(line)
	fp.close()
	
	print "3D-HST / make_region_file: %s." %regfile



class SWarp(object):
    """
    SWarp(object)
    
    This is a class to wrap around SWARP, modeled after the 
    SExtractor class
    
    A workflow might be something like:
    
    >>> sw = threedhst.sex.SWarp()
        >>> sw._aXeDefaults()
    
        >>> sw.swarpMatchImage('IB3714050_drz.fits')  # get reference image parameters from IB3714050_drz.fits
        SWarp.swarpMatchImage: PIXEL_SCALE= 0.128250047918
                                IMAGE_SIZE= 1426,1380
                                    CENTER=  12:36:44.98,  62:08:36.93
        
        >>> sw.swarpImage('IB3714050_drz.fits[1]')    # swarp the reference image to istelf
        THREEDHST/swarp: swarp IB3714050_drz.fits[1] -c auto_default.swarp
        
        >>> sw.swarpRecenter()                        # Refine center position from SWarp's own output
        >>> sw.swarpImage('xxx.fits')                 # SWarp xxx.fits to same pixels as ref image 
    """
    #@staticmethod
    def _getSwarpDefaults(self):
        from subprocess import Popen,PIPE
        
        optinfo = {}
        opts = {}
        optorder = []
        
        try:
            pconf = Popen('swarp -dd'.split(),executable='swarp',stdout=PIPE,stderr=PIPE)
            pconf.wait()
            confstr = pconf.communicate()[0]
        except OSError:
            raise OSError('SWarp not found on system path')
        
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
                      
        self._optinfo = optinfo
        self._defaultopts = opts
        self._optorder = optorder #TODO:OrderedDict for 2.7
    
    
    def _aXeDefaults(self):
        """
        _aXeDefaults
        
            Set SUBTRACT_BACK, WRITE_XML to 'N'
        """
        self.options['SUBTRACT_BACK'] = 'N'
        self.options['WRITE_XML'] = 'N'
    
    
    #@staticmethod   
    def getOptInfo(self,aslist=False):
        """
        returns the  dictionary of input options and the associated information
        
        if aslist is True, returns an list of 
        """
        if aslist:
            return [(k,self._optinfo[k]) for k in self._optorder]
        else:
            return dict(self._optinfo)        
    
    
    def __init__(self,swarpfile=None):
        from warnings import warn
        
        self.lasterr = None
        self.lastout = None
        
        self._getSwarpDefaults()
        opts = dict(self._defaultopts)
            
        if swarpfile:
            #with open(swarpfile) as f:
            f =  open(swarpfile) #gbb
            for l in f:
                commenti = l.find('#')
                if commenti > -1:
                    l = l[:commenti]
                ls = l.split()
                if len(ls) > 1:
                    k = ls[0].strip()
                    if k not in opts:
                        # raise ValueError('swarpfile has invalid option %s'%k) 
                        warn('swarpfile \'%s\' has invalid option %s' %(swarpfile,k)) #gbb
                    # opts[k] = ls[1].strip()
                    if len(ls) > 2:
                        #print ls
                        opts[k] = ' '.join(ls[1:]) #.strip() #gbb
                    else:
                        opts[k] = ls[1].strip()
            self.name = swarpfile.replace('.swarp','')
        else:
            self.name = 'auto_default'
                                
        self.options = opts
        
        self.overwrite = False
    
    
    def _saveFiles(self,fnbase):
        import os
        
        fnbase = fnbase.replace('.swarp','')
        
        dir = os.path.split(fnbase)[0]
        dir = '' if dir=='' else dir+os.sep
                
        ostr = self._makeOptionStr()
        f = open(fnbase+'.swarp','w') #gbb
        f.write(ostr)                
    
    
    def _makeOptionStr(self, maxlen=14):
        #### get longest number of characters of parameter names
        parlen = 0
        for o in self._optorder:
            parlen = np.max([parlen,o.__len__()])
        ostr = []
        for o in self._optorder:
            optlen = np.max([maxlen,self.options[o].__len__()])
            fmt = '%-'+str(parlen)+'s  %-'+str(optlen)+'s  # %s'
            ostr.append(fmt %(o, self.options[o], self._optinfo[o]))
        ostr = '\n'.join(ostr)
        return ostr    
    
    
    def getOptionList(self,incval=False):
        """
        returns a list of all options.  If incval is True, returns a list
        of tuples of the form (optionname,optionvalue)
        """
        if incval:
            return [(k,self.options[k]) for k in self._optorder]
        else:
            return [k for k in self._optorder]
    
    
    @staticmethod
    def decimalToDMS(degrees, hours=False):
        """
        decimalToDMS(self, degrees, hours=False)
        
        Convert decimal degrees to DD:MM:SS.SS.
        
        If hours==True, then return HH:MM:SS.SS
        """
        convert = 1.
        if hours:
            convert = 24/360.
        d = np.abs(degrees*convert)
        
        deg = np.floor(d)
        min = (d-deg)*60.
        sec = (min-np.floor(min))*60
        decimal_sec = (sec-np.floor(sec))*100
        
        if type(deg).__name__ == 'float64':
            if deg < 0:
                out_str = '-'
            else:
                out_str = ' '
            out_str += '%s:%s:%s.%d' %(str(np.int(deg)).rjust(2,'0'),
                                       str(np.int(min)).rjust(2,'0'),
                                       str(np.int(sec)).rjust(2,'0'),
                                       decimal_sec)
        else:
            out_str = []
            print deg.shape, deg
            for i in range(deg.shape[0]):
                if deg[i] < 0:
                    out_i = '-'
                else:
                    out_i = ' '
                out_i += '%s:%s:%s.%d' %(str(np.int(deg[i])).rjust(2,'0'),
                                           str(np.int(min[i])).rjust(2,'0'),
                                           str(np.int(sec[i])).rjust(2,'0'),
                                           decimal_sec[i])
                out_str.append(out_i)
        return out_str
    
    
    def swarpMatchImage(self, matchImage, extension=1, verbose=True):
        """
        swarpMatchImage(self, matchImage, extension=1, verbose=True)
        
        Get WCS image from matchImage[extension] and set swarp parameters so that 
        the output image will have the same size/position.
        """
        import pyfits, pywcs
        from pyraf import iraf
        
        im = pyfits.open(matchImage)
        head = im[extension].header
        wcs = pywcs.WCS(head)
        coord = wcs.all_pix2sky([[head['NAXIS1']/2.,head['NAXIS1']/2.]],0)
        print coord
        ra0 = self.decimalToDMS(coord[0][0],hours=True)
        de0 = self.decimalToDMS(coord[0][1],hours=False)
        #iraf.xy2rd(infile=matchImage+'['+str(extension)+']', x=head['NAXIS1']/2., y = head['NAXIS2']/2.)
        #ra0 = iraf.xy2rd.ra
        #de0 = iraf.xy2rd.dec
        #print ra0,de0
        self.options['CENTER'] = ra0+', '+de0
        if 'IDCSCALE' in head.keys():
            self.options['PIXEL_SCALE'] = str(head['IDCSCALE'])
        else:
            self.options['PIXEL_SCALE'] = str(np.abs(head['CD1_1']*3600.))
        self.options['IMAGE_SIZE']  = '%s,%s' %(head['NAXIS1'],head['NAXIS2'])
        self.options['CENTER_TYPE'] = 'MANUAL'
        self.options['PIXELSCALE_TYPE'] = 'MANUAL'
        
        if verbose:
            print """
SWarp.swarpMatchImage: PIXEL_SCALE= %s
                        IMAGE_SIZE= %s
                            CENTER= %s""" %(self.options['PIXEL_SCALE'],self.options['IMAGE_SIZE'],self.options['CENTER'])
    
    
    def swarpImage(self,inputImage,mode='waiterror'):
        """
        swarpImage(self,inputImage,mode='waiterror')
        
        Writes configuration files and runs swarp on the input image
        
        mode can be:
        
        * 'waiterror': waits for swarp to finish, and raises an 
          SExtractorError if it does not complete sucessfully. stdout 
          and sterr are saved to self.lastout and self.lasterr (returns 0)
        * 'wait': waits for swarp to finish and returns the return code
          stdout and sterr are saved to self.lastout and self.lasterr
        * 'proc': stars the processes but does not wait - returns the Popen 
          instance of the processes
        """
        from subprocess import Popen,PIPE
        from os.path import exists
        
        self.swarpInputImage = inputImage
        
        fnbase = self.name
        if not self.overwrite:
            fnbase = fnbase.replace('.swarp','')
            if exists(fnbase+'.swarp'):
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
        clstr = 'swarp %s -c %s' %(inputImage,self.name+'.swarp')
        proc = Popen(clstr.split(),executable='swarp',stdout=PIPE,stderr=PIPE)
        
        print 'THREEDHST/SWarp.swarpImage: %s' %clstr
        
        if mode == 'waiterror' or mode =='wait':
            res = proc.wait()
            sout,serr = proc.communicate()
            
            self.lastout = sout
            self.lasterr = serr
            
            if res!=0 and mode == 'waiterror' :
                raise SError(serr,sout)
            return res
        elif mode == 'proc':
            return proc
        else:
            raise ValueError('unrecognized mode argument '+str(mode))
            
    
    
    def swarpRecenter(self):
        """
        swarpRecenter(self)
        
        Rerun swarp getting the exact center coordinates from the previous SWarp run.
        This is required to get the pixels in the input and output images to
        coincide exactly.
        
        For best results, edit the ``degtosexal`` and ``degtosexde`` functions in
        ``swarp/src/fitswcs.c`` to print out 4 decimal places in the coordinates:
    	    
    	    sprintf(str,"%c%02d:%02d:%07.4f", sign, dd, dm, ds);
        
        The first swarp run tries to compute the center coordinates directly using the WCS
        information of the center pixel (NAXIS1/2, NAXIS2/2), but this doesn't match
        swarp's internal WCS computation close enough.
        
        """
        if self.lasterr:
            #### Find the second time 'Center' appears 
            idx = self.lasterr.find('Center')
            idx = self.lasterr.find('Center',idx+5)
            idx2 = self.lasterr.find('\n',idx)
            coords = self.lasterr[idx:idx2].split()
            ra0 = coords[1] 
            de0 = coords[2]
            self.options['CENTER'] = ra0+', '+de0
            print 'THREEDHST/Swarp.recenter: %s' %(self.options['CENTER'])
            self.overwrite = True
            #self.swarpImage(self.swarpInputImage)
        else:
            print 'THREEDHST/Swarp.recenter: No SWarp output found'
    


# End