"""
SExtractor class from astropysics, written by Erik Tollerud.
http://packages.python.org/Astropysics/

Modified by GBB to read the output from "sex -dd" and "sex -dp" from files, 
rather than from a piped output because the "pparm" line didn't run for some 
reason.

Also added similar SWarp class wrapper around SWarp.
"""

__version__ = "$Rev: 449 $"
# $URL: https://threedhst.googlecode.com/svn/threedhst/sex.py $
# $Author: gbrammer@gmail.com $
# $Date: 2015-01-19 14:13:55 -0500 (Mon, 19 Jan 2015) $

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import numpy as np

import aXe2html.sexcat.sextractcat
import threedhst

RUN_MODE = 'waiterror'
USE_CONVFILE = 'default.conv'
USE_CONVFILE = 'gauss_4.0_7x7.conv'

class SExtractor(object):
    """
SExtractor()
    
    This class is an adaptor to the Sextractor program 
    (Bertin & Arnouts 1996, http://astromatic.iap.fr/software/sextractor/).
    Sextractor must be installed and in your $PATH for this class to function.
    
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
            pconf = Popen('sex -dd'.split(),
                          executable='sex',stdout=PIPE,stderr=PIPE)
            pconf.wait()
            confstr = pconf.communicate()[0]
            # pparm = Popen('sex -dp'.split(),
            #               executable='sex',stdout=PIPE,stderr=PIPE)
            # pparm.wait()
            # parmstr = pparm.communicate()[0]
            parmstr = threedhst.utils.get_package_data('sexdp') #gbb
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
                if '[' in ls[-1]:
                    unit = ls[-1].replace('[','').replace(']','')                 
                else:
                    unit = None
                info = ' '.join(ls[1:(-1 if unit else None)])
                parinfo[k] = (info,unit if unit else '')
                parorder.append(k)
        
        SExtractor._optinfo = optinfo
        SExtractor._defaultopts = opts
        SExtractor._optorder = optorder #TODO:OrderedDict for 2.7
        SExtractor._parinfo = parinfo
        SExtractor._parorder = parorder #TODO:OrderedDict for 2.7
    
    def copyConvFile(self, grism=False):
        """
        copyConvFile()
        
        Copy a convolution kernel from threedhst/data to ./.
        
        if `grism` is False:
            the kernel is specified by threedhst.sex.USE_CONVFILE ['default.conv']
        
        if `grism` is True:
            Use a special kernel elongated along the spectral axis,
            "threedhst/data/grism.conv".
        """
        if grism:
            self.conv = threedhst.utils.get_package_data('grism.conv')
            fp = open('grism.conv','w')
            fp.write(self.conv)
            fp.close()
            self.conv = self.conv.split('\n')
            self.options['FILTER_NAME'] = 'grism.conv'
        else:
            self.conv = threedhst.utils.get_package_data(USE_CONVFILE)
            fp = open(USE_CONVFILE,'w')
            fp.write(self.conv)
            fp.close()
            self.conv = self.conv.split('\n')
            self.options['FILTER_NAME'] = USE_CONVFILE
        
        #### NNW file for CLASS_STAR
        nnw = threedhst.utils.get_package_data('default.nnw')
        fp = open('default.nnw','w')
        fp.write(nnw)
        fp.close()
        
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
            
        self.executable = 'sex'
        
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
                        warn('sexfile \'%s\' has invalid option %s'
                             %(sexfile,k)) #gbb
                    # opts[k] = ls[1].strip()
                    if len(ls) > 2:
                        #print ls
                        opts[k] = ' '.join(ls[1:]) #.strip() #gbb
                    else:
                        opts[k] = ls[1].strip()
            self.name = sexfile.replace('.sex','')
        else:
            self.name = 'threedhst_auto'
                        
                        
        pfile = opts['PARAMETERS_NAME'] if parfile is None else parfile
        if pfile != 'default.param':
            f = open(pfile)
            for l in f:
                if not (l.strip().startswith('#') or l.strip()==''):
                    k = l.split()[0].strip()
                    if k not in pars:
                        raise ValueError('param file has invalid parameter %s' 
                                         %k)
                    pars[k] = True
        if not np.any(list(pars.values())):
            #if no outputs provided, use defaults (from v 2.8.6)
            defs286 =['NUMBER','FLUX_ISO', 'FLUXERR_ISO', 'FLUX_AUTO',
                      'FLUXERR_AUTO', 'X_IMAGE', 'Y_IMAGE', 'FLAGS']
            for p in defs286:
                if p in pars:
                    pars[p] = True
        
        self.options = opts
        self.params = pars
        
        self.overwrite = False
        
        self.MULTIOPT=False
            
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
        useParams = threedhst.utils.get_package_data('aXe.param').split('\n')[:-1]
        for par in useParams:
            self.params[par] = True
        
        useInputParams = threedhst.utils.get_package_data('sexdd').split('\n')[:-1]
        for line in useInputParams:
            spl = line.strip().split('#')
            if spl[0] is not '':
                spl2 = spl[0].split()
                if len(spl2) >= 2:
                    par = spl2[0]
                    value = ''.join(spl2[1:])
                    #print '%s %s' %(par, value)
                    self.options[par] = value
                    
    def _saveFiles(self,fnbase):
        """
        
        make threedhst_auto.[sex/param] files
        
        """
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
        
        multi_dict = {'FLUX_RADIUS':'PHOT_FLUXFRAC', 'FLUX_APER':'PHOT_APERTURES', 'FLUXERR_APER':'PHOT_APERTURES', 'MAG_APER':'PHOT_APERTURES', 'MAGERR_APER':'PHOT_APERTURES'}
        self.MULTIOPT = False
        
        pstr = []
        for p in self._parorder:
            if self.params[p]:
                ### Test if one of the parameters with optionally multiple inputs, 
                ### like PHOT_FLUXFRAC or PHOT_APERTURES, has multiple values.  If so,
                ### append the necessary (N) to the parameter, like FLUX_RADIUS
                if p in list(multi_dict.keys()):
                    optvalue = self.options[multi_dict[p]]
                    NOPT = len(optvalue.split(','))
                    if NOPT > 1:
                        self.MULTIOPT=True
                        for i in range(1,NOPT+1):
                            pp = '%s(%0d)' %(p, i)
                            pstr.append(fmt %(pp,self._parinfo[p]))        
                    else:
                        pstr.append(fmt %(p,self._parinfo[p]))        
                else:
                    pstr.append(fmt %(p,self._parinfo[p]))
                            
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
    
    
    def sextractImage(self,detectionImage,analysisImage=None,mode=None):
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
        
        if mode is None:
            mode = RUN_MODE
            
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
        
        ### write .param and .sex file
        self._saveFiles(fnbase)
        
        if analysisImage:
            clstr = '%s %s %s -c %s' %(self.executable, detectionImage,
                     analysisImage,self.name+'.sex')
        else:
            # clstr = 'sex {0} -c {1}'.format(detectionImage,self.name+'.sex')
            clstr = '%s %s -c %s' %(self.executable, detectionImage,self.name+'.sex')
        
        print('THREEDHST/sex: %s' %clstr)
        
        if mode == 'waiterror' or mode =='wait':
            
            fp = open('sex_stderr','w')
            proc = Popen(clstr.split(),executable=self.executable, stdout=PIPE,stderr=fp)
            res = proc.wait()
            fp.close()
            
            sout, serr = proc.communicate()
            
            ## Read stderr output
            fp = open('sex_stderr','r')
            serr = ' '.join(fp.readlines())
            fp.close()
                        
            self.lastout = sout
            self.lasterr = serr
            
            if res!=0 and mode == 'waiterror' :
                print(serr, sout)
                raise SError(serr,sout)
            
            self._fix_ascii_head()
            return res
        elif mode == 'proc':
            proc = Popen(clstr.split(),executable=self.executable,stdout=PIPE,stderr=PIPE)

            self._fix_ascii_head()
            return proc
        elif mode == 'direct':
            proc = Popen(clstr.split())
            res = proc.wait()
            self._fix_ascii_head()
        else:
            raise ValueError('unrecognized mode argument '+str(mode))
        
    def _fix_ascii_head(self):
        """
    Fix header for multiple columns like FLUX_APER(2).  The default SExtractor 
    output skips putting these columns in the header, which makes it awkward for 
    reading the catalog later.  
    
    This function is run automatically at the end of `_sextractImage()` if CATALOG_TYPE == 'ASCII_HEAD' and inserts a header line for columns like FLUX_APER(N), calling the value FLUX_APERN for N>1.
        """
        
        if (self.options['CATALOG_TYPE'] == 'ASCII_HEAD') & (self.MULTIOPT):
            #print 'Fix catalog for multi-input parameters...'
            fp = open(self.options['CATALOG_NAME'])
            cat_lines = fp.readlines()
            fp.close()
            NHEADER=0
            while cat_lines[NHEADER].startswith('#'):
                NHEADER+=1
            
            fp = open(self.options['PARAMETERS_NAME'])
            output_params = fp.readlines()
            fp.close()
            
            multi_dict = {'FLUX_RADIUS':'PHOT_FLUXFRAC', 'FLUX_APER':'PHOT_APERTURES', 'FLUXERR_APER':'PHOT_APERTURES', 'MAG_APER':'PHOT_APERTURES', 'MAGERR_APER':'PHOT_APERTURES'}
            
            #### Loop through parameters.  Look for multiple output columns, like
            #### FLUX_RADIUS(1) and also put the option values in the header line
            for ii,output_param in enumerate(output_params):
                par, comment = output_param.strip().split('#')
                
                ## First parameter (parN=1)
                
                if '(' in par:
                    parN = np.int(par.split('(')[1].split(')')[0])
                    parsp = '%s%0d' %(par.split('(')[0], parN)
                    if parN > 1:
                        newline = '# %3d %-22s %-58s [%s]\n' %(ii+1, parsp, comment.split('\'')[1], comment.split('\'')[3])
                        cat_lines.insert(ii, newline)
                
                #### Put a bit more info in the header, like the radius for
                #### aperture photometry
                pari = par.split('(')[0]
                if pari in list(multi_dict.keys()):
                    optval = self.options[multi_dict[pari]]
                    line = cat_lines[ii]
                    sp1 = line.split('[')
                    if len(sp1) == 2:
                        unit = '['+sp1[1]
                    else:
                        unit = '\n'
                    
                    sp2 = sp1[0].strip().split()
                    ix = np.int(sp2[1])
                    parsp = sp2[2]
                    comment = ' '.join(sp2[3:])+' '+optval
                    
                    newline = '# %3d %-22s %-58s %s' %(ix, parsp, comment, unit)
                    cat_lines[ii] = newline
                    
            fp = open(self.options['CATALOG_NAME'],'w')
            fp.writelines(cat_lines)
            fp.close()
            
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
    
    if os.access(sexcat, os.R_OK) is False:
        print("SExtractor catalog, %s, not found." %(sexcat))
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
    
    print('3D-HST / make_region_file: %s.\n' %regfile)

class mySexCat(aXe2html.sexcat.sextractcat.SexCat):
    """
    Extend aXe2html.sexcat.sextractcat.SexCat Class to include option to 
    pop lines from the catalog.
    """
    def __init__(self, filename):
        """
        __init__(self, filename)
        
        Same as base __init__ but save ``filename`` object.
        """
        self.filename = filename
        self.linelist    = self.opencat(filename)
        self.add_missing_header_lines()
        self.headerlines = self.extractheader(self.linelist)
        self.rowlines    = self.extractrows(self.linelist)
        allheads    = self.makeheads(self.headerlines)
        self.ncols  = len(allheads)
        self.nrows  = self.makecols(allheads, self.rowlines)
        success     = self.makeorder()
        
        #### populate columns
        self._easy_columns()
    
    def add_missing_header_lines(self, verbose=False):
        """
        There is only one line in the SExtractor catalog header for some cases where there 
        may be multiple columns, like if multiple apertures are specified for FLUX_APER. 
        
        Find these lines and add them as needed to the header of the SExtractor catalog.
        """
        indexes = []
        names = []
        comments = []
        for line in self.linelist:
            if not line.startswith('#'):
                break
            else:
                sp = line[:-1].split()
                indexes.append(sp[1])
                names.append(sp[2])
                comments.append(' '.join(sp[3:]))
        
        names = np.array(names)
        comments = np.array(comments)
        indexes = np.cast[int](indexes)
        delta = indexes[1:]-indexes[:-1]
        
        skip = delta > 1
        count = 0
        for nskip, idx, name, comment in zip(delta[skip], indexes[skip], names[skip], comments[skip]):
            if verbose:
                print('# %d %s %s' %(idx, name, comment))
            for j in range(1,nskip):
                if verbose:
                    print('# %d %s%d  %s' %(idx+j, name, j+1, comment))
                newname = '%s%d' %(name, j+1)
                self.linelist.insert(idx+j-1, '#  %d %-21s  %s \n' %(idx+j, newname, comment) )
            
    def popItem(self, number_out, verbose=False):
        """
        popItem(self, NUMBER[s], verbose=False)
        
        Pop item id#NUMBER from a SExtractor catalog.
        """
        
        #### search for line with object NUMBER itself and pop it. 
        #### Return message if number not found)
        numbers = list(np.cast[int](self.columns[self.searchcol('NUMBER')].entry))
        
        for number in number_out:
            if number not in numbers:
                continue
        
            idx = np.where(np.array(numbers) == number)[0]
            
            lineOut = self.rowlines.pop(idx)
            num = numbers.pop(idx)
            ll = self.linelist.pop(idx+self.ncols)
            if verbose:
                if verbose > 1:
                    print(lineOut)
                else:
                    print(number)
                    
        allheads    = self.makeheads(self.headerlines)
        self.ncols  = len(allheads)
        self.nrows  = self.makecols(allheads, self.rowlines)
        success     = self.makeorder()
        #### repopulate columns
        self._easy_columns()
        
    def write(self, outfile=None, reformat_header=False):
        """
        write(self, outfile=None)
        
        Write catalog lines to file.  Default overwrites the initial file 
        (``self.filename``).
        """
        if not outfile:
            outfile = self.filename
        fp = open(outfile,'w')
        if reformat_header:
            """ 
            Make a header like
            
            # id ra dec ....
            
            rather than the SExtractor format.
            """
            head = '# '
            for col in self.column_names:
                head+=' '+col
            fp.write(head+'\n')
        else:
            fp.writelines(self.headerlines)
        
        fp.writelines(self.rowlines)
        fp.close()
    
    def change_MAG_AUTO_for_aXe(self, filter='F1392W'):
        """
        change_MAG_AUTO(self, filter='F1392W')
        
        Change the MAG_AUTO column in the catalog to be MAG_{filter} for aXe.
        """
        from warnings import warn
        found_match = False
        for i, line in enumerate(self.headerlines):
            if ' MAG_AUTO ' in line:
                found_match = True
                break
        
        if found_match:
            spl = line.split(' MAG_AUTO ')
            self.headerlines[i] = spl[0]+' MAG_'+filter.upper()+' '+spl[1]
            warn('change_MAG_AUTO_for_aXe: MAG_AUTO -> MAG_'+\
                 filter.upper()+'\n')
            allheads    = self.makeheads(self.headerlines)
            self.ncols  = len(allheads)
            self.nrows  = self.makecols(allheads, self.rowlines)
            success     = self.makeorder()
            #### repopulate columns
            self._easy_columns()
        else:
            warn('change_MAG_AUTO_for_aXe: No MAG_AUTO column found\n')
    
    def _easy_columns(self):
        """
easy_columns()
        
        Populate self.column_names and add make column data easier to get out, 
        like:
        
        >>> id = sexCat.ID   (memory alias, changes to id also in sexCat.ID)
        >>> id = sexCat.ID+0 (make a copy of the array)
        """
        self.column_names = []
        for col in self.columns:
            self.column_names.append(col.getname())
            str = 'self.%s = col.entry*1' %col.getname() # *1 makes a copy
            exec(str)
        
        self.makeRaDec()
        
    def __getitem__(self, column_name):
        """
__getitem__(column_name)
    
    >>> cat = mySexCat('drz.cat')
    >>> print cat['NUMBER']
    
        """
        
        if column_name.upper() not in self.column_names:
            print(('Column %s not found.  Check `column_names` attribute.'
                    %column_name))
            return None
        else:
            str_exec = 'first_item = self.%s[0]' %(column_name.upper())
            exec(str_exec)
            
            if first_item.isdigit():
                type='int'
            else:
                type='float'
            
            try:
                str = 'out = np.cast[%s](self.%s)' %(type, column_name.upper())
                exec(str)
            except:
                str = 'out = self.%s*1' %column_name.upper()
                exec(str)
            
            return out
    
    def makeRaDec(self):
        """
makeRaDec()
    
    id = int(NUMBER)
    ra = float(X_WORLD)
    dec = float(Y_WORLD)
        """
        if 'NUMBER' in self.column_names:
            self.id = np.cast[int](np.array(self.NUMBER))
        if 'X_WORLD' in self.column_names:
            self.ra = np.cast[float](np.array(self.X_WORLD))
        if 'Y_WORLD' in self.column_names:
            self.dec = np.cast[float](np.array(self.Y_WORLD))
    
    def addColumn(self, data=np.arange(2), format='%f', name='NEWDATA', comment='', verbose=False):
        """
        Add a column to a SExtractor catalog
        """
        if not isinstance(data, np.array(1).__class__):
            print("ERROR: `data` is not a numpy array.")
            return False
                
        if data.shape != (self.nrows, ):
            print("ERROR: `data` must have shape (%0d,); has" %(self.nrows), data.shape)
            return False
        
        #### Data array checks out.
        newheader = '#%4d %-22s %s\n' %(self.ncols+1, name, comment)
        self.headerlines.append(newheader)
        self.linelist.insert(self.ncols, newheader)
        self.ncols += 1
        
        self.rowlines    = self.extractrows(self.linelist)
        for i,line in enumerate(self.rowlines):
            newline = line.split('\n')[0]+' '+format %(data[i])+'\n'
            self.linelist[i+self.ncols] = newline
            self.rowlines[i] = newline
            
        ## Reprocess lines
        allheads    = self.makeheads(self.headerlines)
        self.nrows  = self.makecols(allheads, self.rowlines)
        success     = self.makeorder()
        self._easy_columns()
        
        if verbose:
            print('Added column, %s, with format, %s' %(name, format))
            
        return success == 0
    
    def renameColumn(self, original='X_IMAGE', new='X_NEW', verbose=True):
        if original not in self.column_names:
            print('Column %s not in the current catalog.' %(original))
            return False
        
        for i, line in enumerate(self.headerlines):
            if ' %s ' %(original) in line:
                #print line
                new_header = self.headerlines[i].replace(' %s ' %(original), ' %s ' %(new))
                self.headerlines[i] = new_header
                self.linelist[i] = new_header
                
        ### Reprocess the header with the new column
        allheads    = self.makeheads(self.headerlines)
        self.nrows  = self.makecols(allheads, self.rowlines)
        success     = self.makeorder()
        self._easy_columns()
        
        if verbose:
            print('Renamed column %s -> %s' %(original, new))
            
        return success == 0
        
class SWarp(object):
    """
    SWarp()
    
    This is a class to wrap around SWARP, modeled after the 
    SExtractor class.
    
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
    def _getSwarpDefaults(self):
        from subprocess import Popen,PIPE
        
        optinfo = {}
        opts = {}
        optorder = []
        
        try:
            pconf = Popen('swarp -dd'.split(),
                          executable='swarp',stdout=PIPE,stderr=PIPE)
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
        self.options['VMEM_MAX'] = '8192'
        self.options['MEM_MAX'] = '8192'
        self.options['COMBINE_BUFSIZE'] = '8192'
    
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
                        warn('swarpfile \'%s\' has invalid option %s'
                             %(swarpfile,k)) #gbb
                    # opts[k] = ls[1].strip()
                    if len(ls) > 2:
                        #print ls
                        opts[k] = ' '.join(ls[1:]) #.strip() #gbb
                    else:
                        opts[k] = ls[1].strip()
            self.name = swarpfile.replace('.swarp','')
            f.close()
        else:
            self.name = 'threedhst_auto'
                                
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
        f.close()
        
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
            if degrees < 0:
                out_str = '-'
            else:
                out_str = ' '
            out_str += '%s:%s:%s.%d' %(str(np.int(deg)).rjust(2,'0'),
                                       str(np.int(min)).rjust(2,'0'),
                                       str(np.int(sec)).rjust(2,'0'),
                                       decimal_sec)
        else:
            out_str = []
            # print deg.shape, deg
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
        
        Get WCS image from matchImage[extension] and set swarp parameters
        so that the output image will have the same size/position.
        """
        try:
            import astropy.wcs as pywcs
        except:
            import pywcs
        
        from pyraf import iraf
        
        im = pyfits.open(matchImage)
        head = im[extension].header
        
        #### pywcs incompatible with astrodrizzle lookup table distortion?
        if 'CPDIS1' in list(head.keys()):
            head['CPDIS1'] = 'SIP'
            head['CPDIS2'] = 'SIP'
            
        wcs = pywcs.WCS(head)
        # coord = wcs.all_pix2sky([[head['NAXIS1']/2.,head['NAXIS1']/2.]],0)
        try:
            coord = wcs.all_pix2sky([[head['CRPIX1'], head['CRPIX2']]],1)
            coord = wcs.all_pix2sky([[head['NAXIS1']/2., head['NAXIS2']/2.]],1)
        except:
            coord = wcs.all_pix2world([[head['CRPIX1'], head['CRPIX2']]],1)
            coord = wcs.all_pix2world([[head['NAXIS1']/2., head['NAXIS2']/2.]],1)
            
        # print coord
        ra0 = self.decimalToDMS(coord[0][0],hours=True)
        de0 = self.decimalToDMS(coord[0][1],hours=False)
        #iraf.xy2rd(infile=matchImage+'['+str(extension)+']', x=head['NAXIS1']/2., y = head['NAXIS2']/2.)
        #ra0 = iraf.xy2rd.ra
        #de0 = iraf.xy2rd.dec
        #print ra0,de0
        self.options['CENTER'] = ra0+', '+de0
        #if 'IDCSCALE' in head.keys():
        #    self.options['PIXEL_SCALE'] = str(head['IDCSCALE'])
        #else:
        #    #self.options['PIXEL_SCALE'] = str(np.abs(head['CD1_1']*3600.))
        if 'CD1_2' in list(head.keys()):
            CD1_2 = head['CD1_2']
        else:
            CD1_2 = 0.
        self.options['PIXEL_SCALE'] = str(np.sqrt(head['CD1_1']**2+CD1_2**2)*3600.)
        
        self.options['IMAGE_SIZE']  = '%s,%s' %(head['NAXIS1'],head['NAXIS2'])
        self.options['CENTER_TYPE'] = 'MANUAL'
        self.options['PIXELSCALE_TYPE'] = 'MANUAL'
        
        self.NX = head['NAXIS1']
        self.NY = head['NAXIS2']

        if verbose:
            print("""
SWarp.swarpMatchImage: PIXEL_SCALE=  %s
                        IMAGE_SIZE=  %s
                            CENTER= %s\n""" %(self.options['PIXEL_SCALE'],
                            self.options['IMAGE_SIZE'],self.options['CENTER']))
    
    def swarpImage(self,inputImage, mode='direct', verbose=True):
        """
        swarpImage(self,inputImage,mode='waiterror', verbose=True)
        
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
        from subprocess import Popen, PIPE
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
        if isinstance(inputImage,list):
            imgList = ' '.join(inputImage)
        else:
            imgList = inputImage
        
        clstr = 'swarp %s -c %s' %(imgList,self.name+'.swarp')
        
        #print "\n3DHST.sex.swarp.swarpImage:\n\n %s\n" %clstr
        if verbose:
            threedhst.showMessage('Running swarp: %s' %clstr)
        
        if mode == 'waiterror' or mode =='wait':
            # proc = Popen(clstr.split(),
            #              executable='swarp',stdout=PIPE,stderr=PIPE)
            #### Send STDERR output to temporary file because 
            #### swarp seems to spawn additional processed for large files
            #### and proc.wait() never finishes
            fp = open('sex_stderr','w')
            proc = Popen(clstr.split(),executable='swarp', stdout=PIPE,
                         stderr=fp)
            res = proc.wait()
            fp.close()
            
            if verbose: 
                print('Done.\n')
            
            sout, serr = proc.communicate()
            
            ## Read stderr output
            fp = open('sex_stderr','r')
            serr = ' '.join(fp.readlines())
            fp.close()
            
            self.lastout = sout
            self.lasterr = serr
            
            if res!=0 and mode == 'waiterror' :
                raise SError(serr,sout)
            return res
        elif mode == 'proc':
            return proc
        elif mode == 'direct':
            proc = Popen(clstr.split()) #,executable='swarp' #,stdout=PIPE,stderr=PIPE)
            res = proc.wait()
        else:
            raise ValueError('unrecognized mode argument '+str(mode))
    
    def swarpRecenter(self):
        """
        swarpRecenter(self)
        
        Rerun swarp getting the exact center coordinates from the previous 
        SWarp run.  This is required to get the pixels in the input and output
        images to coincide exactly.
        
        For best results, edit the ``degtosexal`` and ``degtosexde`` 
        functions in ``swarp/src/fitswcs.c`` to print out 4 decimal places 
        in the coordinates:
            
            sprintf(str,"%c%02d:%02d:%07.4f", sign, dd, dm, ds);
        
        The first swarp run tries to compute the center coordinates directly
        using the WCS information of the center pixel (NAXIS1/2, NAXIS2/2), 
        but this doesn't match swarp's internal WCS computation close enough.
        
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
            print('THREEDHST/Swarp.recenter: %s\n' %(self.options['CENTER']))
            self.overwrite = True
            #self.swarpImage(self.swarpInputImage)
        else:
            print('THREEDHST/Swarp.recenter: No SWarp output found\n')

def grow_segmentation_map(seg_file='MACS1149-IR_drz_seg.fits', cat_file='MACS1149-IR.cat', size=30, mag=20):
    import scipy.ndimage as nd
    import pyregion
    import os
    #import stsci.convolve
    
    #seg_file='MACS1149-IR_drz_seg.fits'; cat_file='MACS1149-IR.cat'; size=30; mag=20.5
    
    seg_im = pyfits.open(seg_file)
    seg = seg_im[0].data*1
    
    cat = threedhst.catIO.Table(cat_file)
    
    ok = cat['MAG_APER'] < mag
    so = np.argsort(cat['MAG_APER'][ok])
    idx = np.arange(len(cat))[ok][so]
    
    print('Grow masks:')
    for i in idx:
        id = cat['NUMBER'][i]
        print('%5d (%.2f)' %(id, cat['MAG_APER'][i]))
        mask = (seg == id)*1
        grow_mask = nd.maximum_filter(mask, size=size*2)
        #grow_mask = nd.convolve(mask, np.ones((size, size)))
        #grow_mask = stsci.convolve.convolve2d(mask, np.ones((size, size)), output=None, mode='nearest', cval=0.0, fft=1)
        seg[seg == 0] += grow_mask[seg == 0]*id
     
    reg_file = seg_file.replace('.fits', '_mask.reg')
    if os.path.exists(reg_file):
        print('Region mask: %s' %(reg_file))
        reg = pyregion.open(reg_file).as_imagecoord(seg_im[0].header)
        N = len(reg)
        for i in range(N):
            if 'text=' not in reg[i].comment:
                continue
        
            id = int(reg[i].comment.strip('text={')[:-1])
            print('%4d' %(id))
            mask = reg[i:i+1].get_mask(seg_im[0])
            seg[seg == 0] += mask[seg == 0]*id
    
    pyfits.writeto(seg_file.replace('seg.fits','seg_grow.fits'), data=seg, header=seg_im[0].header, clobber=True)
          
# End