import string
import os
import time

import numpy as np 
import pyfits

import threedhst

def columnFormat(colname):
    """
Set format for common column names.  
    
    id: 'A'
    ra,dec: 'D'
    [everything else]: 'E'
    """
    fmt='E'
    if colname.lower().find('id') >= 0: fmt='D'
    if colname.lower().find('ra') >= 0: fmt='D'
    if colname.lower().find('dec') >= 0: fmt='D'
    return fmt
    
def ASCIItoFITS(infile, comment='#'):
    """
ASCIItoFITS(infile, [comment='#'])

    Read an ASCII file, infile, and get column names from first line, 
    which begins with the 'comment' character.
    
    Output will be in infile+'.FITS'
    """
    ### Get first header line and count commented lines to skip
    file = open(infile,'r')
    line0 = file.readline()
    line=line0
    hskip=0
    while line.startswith(comment):
        line = file.readline()
        hskip +=1
    file.close()
    
    #### Read data file
    data=np.loadtxt(infile, comments=comment)
    
    #### clean up special characters from header
    line0.replace('.','p')
    line0.replace('-','_')
    line0.replace('(','')
    line0.replace(')','')
    line0.replace('-','_')
    
    #### Make output FITS table
    header=string.split(line0[1:-1])
    # make_struct='str = {'
    go_ColDefs='cols=pyfits.ColDefs(['
    for i in range(header.__len__()):
        col_string = 'col%d = pyfits.Column(name=\'%s\',' %(i,header[i]) + \
            ' format=\'%s\', array=data[0:,%d])' %(columnFormat(header[i]),i)
        exec(col_string)
        go_ColDefs += 'col%d,' %(i)
        # make_struct += '\'%s\':data[0:,%d],' %(header[i],i)
    
    exec(go_ColDefs[:-1]+'])') # cols=pyfits.ColDefs([col1, col2, ...])
    
    #### Initialize table
    tbhdu = pyfits.new_table(cols)
    
    #### Primary HDU
    hdu = pyfits.PrimaryHDU()
    
    thdulist = pyfits.HDUList([hdu,tbhdu])
    
    #### Add modification time of "infile" to FITS header
    infile_mod_time = time.strftime("%m/%d/%Y %I:%M:%S %p", \
                         time.localtime(os.path.getmtime(infile)))
    thdulist[1].header.update('MODTIME',infile_mod_time)
    
    thdulist.writeto(infile+'.FITS', clobber=True)
    
    return tbhdu.data, tbhdu.columns
    
def ReadASCIICat(infile, comment='#', force=False, verbose=False, getColumns=False):
    """
data, columns = ReadASCIICat(infile,comment='#', force=False, verbose=False)

    Read FITS table created from an ASCII catalog file.
    
    If ASCIItoFITS output doesn't exist or "force=True", create it.
    """
    
    if os.path.exists(infile) is False:
        print ('File, %s, not found.' %(infile))
        return -1
    
    fileExists = False
    if (os.path.exists(infile+'.FITS')):
        fileExists = True
        theFITSFile = infile+'.FITS'
    if (os.path.exists(infile+'.FITS.gz')):
        fileExists = True
        theFITSFile = infile+'.FITS.gz'
    
    if (fileExists is False) or (force):
        if verbose:
            print ('Running ASCIItoFITS: %s' %(infile))
        data, columns = ASCIItoFITS(infile,comment=comment)
        if getColumns:
            return data, columns
        else:
            return data
        #return ASCIItoFITS(infile,comment=comment)
    else:
        if verbose:
            print ('Reading : %s' %(theFITSFile))
        hdulist = pyfits.open(theFITSFile)
        data, header, columns = hdulist[1].data, hdulist[1].header, hdulist[1].columns
        hdulist.close()
        
        #### Check mod time of 'infile'.  
        #### If changed, then re-read with ASCIItoFITS
        infile_mod_time = time.strftime("%m/%d/%Y %I:%M:%S %p", \
                                time.localtime(os.path.getmtime(infile)))
        if infile_mod_time == header['MODTIME']:
            if getColumns:
                return data, columns
            else:
                return data
        else:
            if verbose:
                print('%s has changed.  Re-generating FITS file...' %(infile))
            data, columns = ASCIItoFITS(infile,comment=comment)
            if getColumns:
                return data, columns
            else:
                return data

def easyCat(cat):
    for field in cat.names:
        str = 'cat.%s = cat.field(field)' %(field.lower())
        try:
            exec(str)
        except:
            pass
        
    return cat
    
class EmptyCat(dict):
    """
    Make a simple class to make a dict behave rather like the Readfile class
    with attributes rather than dict elements for the catalog columns.
    
    Note: 
        If you set a column first as an attribute, it isn't added to the 
        columns until you try to retrieve it as an item of the dict:
        
        >>> x = EmptyCat()
        >>> x.ra = np.arange(5)
        >>> x.columns
        []
        >>> x['ra']
        array([0, 1, 2, 3, 4])
        >>> x.columns
        ['ra']
        
    """
    def __init__(self):
        self.d = {}
        #self = self.d
        self.N = None
        self.filename = None
        
    def __setitem__(self, key, value):
        self.d[key] = value
        self.__setattr__(key, self.d[key])
        if self.N is None:
            self.N = len(self.d[self.d.keys()[0]])
            
    def __getitem__(self, key):
        if key not in self.d.keys():
            self.__setitem__(key, self.__getattribute__(key))
        
        return self.d[key]
    
    def __repr__(self):
        return self.d.__repr__()
    
    def __str__(self):
        return self.d.__str__()
        
    @property
    def columns( self ):
        return self.d.keys()
                    
#infile='AEGIS/OUTPUT/cat1.0_default_lines.rf'
#data = ReadASCIIFile('AEGIS/OUTPUT/cat1.0_default_lines.rf', verbose=True)
#data = ReadASCIIFile('AEGIS/aegis-n2.v3.3.cat', verbose=True)
class Readfile():
    """
    Read column data from an ASCII file where the data contain a mix of 
    strings, floats, and ints.  
    
    This function now seems to be both reasonably fast and robust.  It
    supports storing FITS versions of the catalogs and is to be preferred
    over "ReadASCIICat" because it is able to handle files with string
    columns.
    """
    def __init__(self, infile='files.info', force_lowercase = True,
                 comment_char='#', verbose=False, save_fits = True):
        
        self.filename = infile
        self.verbose = verbose
        
        #### Load the FITS version of the catalog, if it exists
        status = self.load_fits()
        if status:
            return None
            
        self.comment_char = comment_char
        
        #### read the lines of the file
        fp = open(infile,'r')
        lines = fp.readlines()
        fp.close()
        
        if len(lines) < 2:
            threedhst.showMessage('Only %d lines in %s.' %(len(lines), infile), warn=True)
            self.status = None
            return None
        
        if not lines[0].startswith(comment_char):
            threedhst.showMessage('First line of %s doesn\'t start with \'%s\':\n%s' %(infile,
                                   comment_char, lines[0]), warn=True)
            self.status = None
            return None
            
        #### get the column names from the first line
        header = lines[0]
        columns = header.replace(comment_char,'').split()
        NCOLUMNS = len(columns)
        
        #### parse column names, fixing characters.
        dict = {}
        for i in range(NCOLUMNS):
            if verbose > 1:
                print columns[i]
            col = columns[i].replace('-','_').replace('.','p')
            if force_lowercase:
                col = col.lower()
            for str in '()[]':
                col = col.replace(str,'')
            #
            if col[0].isdigit():
                col = '_'+col
            #    
            columns[i] = col
            dict[col] = []
        
        #### skip header lines
        ix=0
        line = lines[ix]
        while line.startswith(comment_char) & (ix < len(lines)-1):
            ix+=1
            line = lines[ix]
        
        if ix == len(lines):
            self.status = None
            return None
            
        #### Parse the lines into the data columns
        N=0
        for line in lines[ix:]:
            spl = line.split()
            if len(spl) == NCOLUMNS:
                N+=1
                for i in range(NCOLUMNS):
                    dict[columns[i]].append(spl[i])
        
        if verbose > 3:
            print dict
            
        #### Convert data to numpy arrays and change data types from
        #### strings to int and float as necessary
        for col in columns:
            dict[col] = np.array(dict[col])
            item = dict[col][0]
            try:
                fl = float(item)
                isNumber = True
            except:
                isNumber = False
            
            if isNumber:
                try:
                    dict[col] = np.cast[float](dict[col])
                except:
                    pass
            #int
            if item.isdigit():
                try:
                    dict[col] = np.cast[int](dict[col])
                except:
                    pass
                    
            str = 'self.%s = dict[col]' %(col)
            #print str
            exec(str)
            
        self.NCOLUMNS = NCOLUMNS
        self.columns = columns
        self.N = N
        self.status = True
        
        if save_fits:
            self.write_fits()
            
    def __getitem__(self, key):
        """
    Allow you to address the column names as strings, e.g.
    
    >>> data = Readfile('junk.cat')
    >>> column = 'id'
    >>> print data[column]    
    >>> print data.id         
    
    """
        str = 'x = self.%s' %(key)
        exec(str)
        return x
    
    #
    def __setitem__(self, key, value):
        """
    Allow you to set the data in the columns, e.g.
    
    >>> data = Readfile('junk.cat')
    >>> column = 'id'
    >>> data[column] = 1 
    >>> print data.id         
    
    """
        str = 'self.%s = value' %(key)
        exec(str)
        
    def keys(self):
        return self.columns
        
    def field(self, key):
        """
        field(self, key)
        
        Replicate the PYFITS `field` way of getting the data
        """
        return self.__getitem__(key)
        
    def addColumn(self, name, value):
        """ 
        addColumn(self, name, value)
        """
        str = 'self.%s = value' %name
        exec(str)
    
    def write_text(self, filename='junk.cat', columns=None, select=None, comment='#'):
        """
        Write an ascii file with a boolean selection `select` on the input
        columns
        """
        if select is None:
            select = np.ones(self.N) > 0
        
        format_codes = {'int64':'%d','float64':'%.5e','>i8':'%d', '>f8':'%.5e'}
        
        data = []
        formats = []
        header = comment
        if columns is None:
            columns = self.columns
        
        for column in columns:
            header = '%s %s' %(header, column)
            data.append(self.__getitem__(column)[select])
            dtype = str(data[-1].dtype)
            if (column == 'ra') | (column == 'dec'):
                formats.append('%.6f')
            else:
                formats.append(format_codes[dtype])
        
        fp = open(filename,'w')
        fp.write(header+'\n')
        np.savetxt(fp, np.array(data).T, fmt=tuple(formats))
        fp.close()
        
    def write_fits(self):
        """
        Save the ascii catalog data into a FITS bintable.
        
        The modification date of the ascii catalog is saved in the 'MODTIME'
        keyword of the FITS file
        """
        import time
        
        formats = {}
        formats['bool'] = 'L'
        formats['int16'] = 'I'
        formats['int32'] = 'J'
        formats['int64'] = 'K'
        formats['float32'] = 'E'
        formats['float64'] = 'D'
        
        formats['>i8'] = 'K'
        formats['>f8'] = 'D'
        
        #### Make the table columns, translating numpy data types to "TFORM"
        coldefs = []
        for column in self.columns:
            dtype = str(self.__getitem__(column).dtype)
            #print column, dtype
            if dtype in formats.keys():
                TFORM=formats[dtype]
            else:
                if 'S' not in dtype:
                    threedhst.showMessage('Unrecognized data type in %s: %s' %(self.filename, dtype), warn=True)
                    return False
                #
                TFORM = 'A'+dtype.split('S')[1]
            #
            data = self.__getitem__(column)
            if '>' in dtype:
                cast_types = {'>i8':np.int64, '>f8':np.float64}
                data = np.cast[cast_types[dtype]](data)
            #
            coldefs.append(pyfits.Column(name=column, array=data, format=TFORM))
        
        #### Done, now make the binary table
        tbhdu = pyfits.new_table(coldefs)

        #### Primary HDU
        hdu = pyfits.PrimaryHDU()
        thdulist = pyfits.HDUList([hdu,tbhdu])

        #### Add modification time of "infile" to FITS header
        infile_mod_time = time.strftime("%m/%d/%Y %I:%M:%S %p",
                            time.localtime(os.path.getmtime(self.filename)))
        
        thdulist[1].header.update('MODTIME',infile_mod_time)

        thdulist.writeto(self.filename+'.FITS', clobber=True)
        return True
        
    def load_fits(self):
        """
        Read the FITS bintable version of the catalog.   
        
        If the modification date of the ascii file is different than that
        founc in the FITS file, return a status of False and re-read the
        file generating the FITS file again, if "save_fits" is set during
        __init__
        """
        import time
        if not os.path.exists(self.filename+'.FITS'):
            return False
        
        infile_mod_time = time.strftime("%m/%d/%Y %I:%M:%S %p",
                            time.localtime(os.path.getmtime(self.filename)))
        
        im = pyfits.open(self.filename+'.FITS')[1]
        
        if infile_mod_time != im.header['MODTIME']:
            print('%s has changed.  Re-generating FITS file...' %(self.filename))
            return False
        
        if self.verbose:
            print 'Read from: %s.FITS' %(self.filename)
                   
        self.NCOLUMNS = len(im.data.names)
        self.columns = im.data.names
        self.N = len(im.data[im.data.names[0]])
        for column in self.columns:
            run_str = 'self.%s = im.data[\'%s\']' %(column, column)
            exec(run_str)
            
        return True
        
class markerXML():
    def __init__(self, ra, dec, mag):
        self.ra = np.float(ra)
        self.dec = np.float(dec)
        self.mag = np.float(mag)

class CoordinateMatcher():
    """
    Class for providing automatic methods for retrieving 
    nearest neighbors at a given position from a Readfile object
    
    Example:
        >>> cat = catIO.Readfile('fireworks.cat')
        >>> cat_coords = catIO.CoordinateMatcher(cat)
        >>> ra, dec = cat.ra[100], cat.dec[100]
        >>> cat_coords.find_nearest(ra, dec, N=5)
        (array([  0., 3.21631753, 9.56851659, 9.57823153, 10.48355936]),
         array([100,  90,  95, 157, 174], dtype=int32))
    
    """
    def __init__(self, cat, ra_column = 'ra', dec_column = 'dec', USE_WORLD=False):
        
        if USE_WORLD:
            ra_column, dec_column = 'x_world', 'y_world'
                
        for test in [ra_column, dec_column]:
            if test not in cat.columns:
                print 'Column %s not found in the input catalog' %(test)
                self.status = False
                return None
                
        else:
            self.ra_column = ra_column
            self.dec_column = dec_column
            self.cat = cat
            self.init_tree()
    
        self.status = True
    
    def init_tree(self):
        """
        Set up the tree for retrieving nearest zSpec objects from
        test coordinate positions.
        
        The match is done using the scipy.spatial.cKDTree function, which is
        apparently some 80 times faster than doing it "by hand" with numpy.
        """
        import scipy.spatial
        
        cosd = self.cat[self.ra_column] * np.cos(self.cat[self.dec_column]/360*2*np.pi)
        self.xy = np.array([cosd, self.cat[self.dec_column]]).T
        self.tree = scipy.spatial.cKDTree(self.xy, 10)
    
    def find_nearest(self, ra, dec, N=1):
        """
        Find N nearest neighbors to (ra, dec) in the zSpec catalogs.  
        
        Example: 
            
            >>> dist, ids = zsp.find_nearest(zsp.ra[100], zsp.dec[100], N=5)
            >>> print dist, ids
            (array([  0.        ,  12.96253365,  17.63697491,  29.72497372,  31.16232403]), array([100,  86, 119, 116,  80], dtype=int32))
            
        """
        if self.tree is None:
            self.init_tree()
            
        xy_test = [ra*np.cos(dec/360.*2*np.pi), dec]
        dist, ids = self.tree.query(xy_test, k=N)
        return dist*3600, ids
    
    def match_list(self, ra=[], dec=[], N=1, MATCH_SELF=False, verbose=True):
        """
        Make a full matched list, input 'ra' and 'dec' are
        arrays
        
        If MATCH_SELF, find nearest matches *within* the self catalog
        """
        noNewLine = '\x1b[1A\x1b[1M'
        
        if MATCH_SELF:
            ra = self.cat[self.ra_column]
            dec = self.cat[self.dec_column]
        
        Nlist = len(ra)
        dr_match = ra*0.
        id_match = np.cast[int](ra*0)
        
        for i in range(Nlist):
            if verbose:
                print noNewLine+'%d of %d' %(i+1, Nlist)
            
            dist, ids = self.find_nearest(ra[i], dec[i], N=1+N)
            dr_match[i] = dist[N-1+MATCH_SELF]
            id_match[i] = ids[N-1+MATCH_SELF]
        
        return dr_match, id_match
        
def readMarkerXML(xml_file):
    """
    Read a threedhst XML file, which provides ra, dec, mag
    for each object in a catalog.  Output is a dictionary
    where the tags are the float ID numbers.
    """
    fp = open(xml_file)
    lines = fp.readlines()[0].split('<marker id')
    fp.close()
    
    dict = {}
    for line in lines:
        sp = line.split('"')
        if len(sp) > 3:
            dict[np.int(sp[1])] = markerXML(sp[3], sp[5], sp[7])
    
    return dict