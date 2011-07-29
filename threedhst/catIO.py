import string
import os
import time

import numpy as np 
import pyfits

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
        hdulist.close()
        #### Check mod time of 'infile'.  
        #### If changed, then re-read with ASCIItoFITS
        infile_mod_time = time.strftime("%m/%d/%Y %I:%M:%S %p", \
                                time.localtime(os.path.getmtime(infile)))
        if infile_mod_time == hdulist[1].header['MODTIME']:
            if getColumns:
                return hdulist[1].data, hdulist[1].columns
            else:
                return hdulist[1].data
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
    
#infile='AEGIS/OUTPUT/cat1.0_default_lines.rf'
#data = ReadASCIIFile('AEGIS/OUTPUT/cat1.0_default_lines.rf', verbose=True)
#data = ReadASCIIFile('AEGIS/aegis-n2.v3.3.cat', verbose=True)
class Readfile():
    """
    Read column data from an ASCII file where the data contain a mix of 
    strings, floats, and ints.  Should rewrite using numpy.loadtxt with
    format keywords.
    """
    def __init__(self, infile='files.info', force_lowercase = True,
                 verbose=False):
        
        #### read the lines of the file
        fp = open(infile,'r')
        lines = fp.readlines()
        fp.close()
        
        #### get the column names from the first line
        header = lines[0]
        columns = header.replace('#','').split()
        NCOLUMNS = len(columns)
        
        #### parse column names, fixing characters.
        dict = {}
        for i in range(NCOLUMNS):
            col = columns[i].replace('-','_').replace('.','p')
            if force_lowercase:
                col = col.lower()
            for str in '()[]':
                col = col.replace(str,'')
            
            if col[0].isdigit():
                col = '_'+col
                
            columns[i] = col
            dict[col] = []
        
        #### skip header lines
        ix=0
        line = lines[ix]
        while line.startswith('#'):
            ix+=1
            line = lines[ix]
        
        #### Parse the lines into the data columns
        N=0
        for line in lines[ix:]:
            spl = line.split()
            if len(spl) == NCOLUMNS:
                N+=1
                for i in range(NCOLUMNS):
                    dict[columns[i]].append(spl[i])
        
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

class markerXML():
    def __init__(self, ra, dec, mag):
        self.ra = np.float(ra)
        self.dec = np.float(dec)
        self.mag = np.float(mag)

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