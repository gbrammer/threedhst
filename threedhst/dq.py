"""
3DHST.dq

Check DQ (satellite trails, scattered light)
of direct/grism exposures.

"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import pyfits

import numpy as np
import threedhst
import pysao

class myDS9(pysao.ds9):
    """
myDS9(pysao.ds9)
    
    Add shortcut methods to extend `pysao.ds9`
    """
    def scale(self, min, max):
        """
scale(self, min, max)

    ds9.set('scale limits min max')
        """
        self.set('scale limits %f %f' %(min,max))
    #
    def match_all(self, alignType='wcs'):
        """
match_all(self, alignType='wcs')

    ds9.set('match frames `alignType`')
    ds9.set('match colorbars')
    ds9.set('match scales')
        """
        self.set('match frames %s' %(alignType))
        self.set('match colorbars')
        self.set('match scales')
    #
    def delete_all_frames(self):
        """
delete_all_frames(self)

    Delete all frames in open window.
        """
        frames = ds9.get('frame all').split()
        for frame in frames:
            self.set('frame %s' %frame)
            self.set('frame delete')

def flag_dq_polygon(ds9, flt_file):
    """
check_data_quality(flt_file)
    
Display an image and check DQ extension (e.g. for satellite trails)
    
    """
    ##### Start a DS9 instance
    #ds9 = myDS9()
    ##### Open FITS file
    fi = pyfits.open(flt_file)
    ### Display SCI extension [1]
    ds9.frame(1)
    ds9.view(fi[1])
    ds9.scale(-0.1,1)
    ### Display DQ extension [3]
    ds9.frame(2)
    ds9.view(fi[3])
    ds9.scale(0,100)
    ### DS9 Tile and return to Frame 1
    ds9.set('tile yes')
    ds9.set('tile grid')
    ds9.set('zoom to fit')
    ds9.frame(1)
    ds9.set('zoom to fit')
    ##### Ask at prompt if we should define regions and continue
    test = raw_input('3DHST/check_data_quality: Define DQ region y/[n]? ')
    if (test == '') | (test.lower().startswith('n')):
        return None
    else:
        ds9.set('regions shape polygon')
        dummy = raw_input('Define region polygon in DS9 and hit [return] ')
    ##### Define region
    ds9.set('regions system image')
    regions = ds9.get('regions source')
    ##### Initialize DQ image
    dqflag = np.zeros(fi[1].data.shape,dtype=np.int)
    ##### Loop through user-defined regions
    for region in regions.split('\n'):
        if region.strip().startswith('polygon'):
            #region = 'polygon(375.05333,642.2,465.18667,642.2,751.36,
            # 709.8,393.08,326.73333,210.56,461.93333,465.18667,
            # 552.06667,375.05333,552.06667,221.82667,509.25333)'
            spl = np.float_(np.array(
                     region[region.find('(')+1:region.find(')')].split(',')
                     ))
            px = spl[0::2]
            py = spl[1::2]
            dqflag += threedhst.utils.region_mask(fi[1].data.shape,px,py)
            
    ##### Set DQ bit
    dqflag[np.where(dqflag > 0)] = 2048
    ##### Show region mask in Frame 2        
    ds9.frame(2)
    ds9.view_array(fi[3].data+dqflag)
    ds9.scale(0,100)
    ##### Save defined regions in output file, [flt_file]+'.mask.reg'
    fp = open(flt_file.split('.gz')[0]+'.mask.reg','w')
    fp.write(regions)
    fp.close()
    return True

import Tkinter as tk

class QueryWindow:
    """
    Widget application.
    """
    def __init__(self, master, asn_grism, asn_direct):
        self.direct = ''
        self.grism = ''
        frame = tk.Frame(master)
        frame.pack()
        
        ds9 = myDS9()
        self.idx = 0
        
        #### need to add method for displaying grism/direct SCI/DQ images
        
        #
        self.button = tk.Button(frame, 
            text="QUIT", fg="red", command=frame.quit)
        self.button.pack(side=tk.LEFT)
        #
        self.button_next = tk.Button(frame,
            text="NEXT", bg="green", command=self.dummy)
        self.button_next.pack(side=tk.LEFT)
        #
        self.button_kill = tk.Button(frame,
            text="Kill exposure", bg="red", command=self.dummy)
        self.button_kill.pack(side=tk.LEFT)
        #
        self.button_dqdirect = tk.Button(frame,
            text="Flag DIRECT", bg="white", fg="blue", command=self.dummy)
        self.button_dqdirect.pack(side=tk.LEFT)
        #
        self.button_dqgrism = tk.Button(frame,
            text="Flag GRISM", bg="blue", fg="white", command=self.dummy)
        self.button_dqgrism.pack(side=tk.LEFT)
    #
    def dummy(self):
        print 'Placeholder'
        
def testWidget():  
    root = tk.Tk()
    app = QueryWindow(root,'x')
    root.mainloop()

def apply_dq_mask(flt_file, addval=2048):
    """
apply_dq_mask(flt_file, addval=2048)

    Read mask polygons from `flt_file`+'.mask.reg', if available,
    and apply to the DQ extension of `flt_file`.
    
    DQnew = DQold + `addval` within the polygon.
    """
    try:
        fp = open(flt_file.split('.gz')[0]+'.mask.reg','r')
    except:
        return None
    #
    print 'Applying mask from %s.mask_reg' %(flt_file.split('.gz')[0])
    regions = ''.join(fp.readlines())
    fi = pyfits.open(flt_file,mode='update')
    dqflag = np.zeros(fi[3].data.shape,dtype=np.int)
    ##### Loop through user-defined regions
    for region in regions.split('\n'):
        if region.strip().startswith('polygon'):
            #region = 'polygon(375.05333,642.2,465.18667,642.2,751.36,
            # 709.8,393.08,326.73333,210.56,461.93333,465.18667,
            # 552.06667,375.05333,552.06667,221.82667,509.25333)'
            spl = np.float_(np.array(
                     region[region.find('(')+1:region.find(')')].split(',')
                     ))
            px = spl[0::2]
            py = spl[1::2]
            dqflag += threedhst.utils.region_mask(fi[1].data.shape,px,py)
    
    ##### Set DQ bit
    dqflag[np.where(dqflag > 0)] = addval
    fi[3].data+=dqflag
    ##### Write back to flt_file
    fi.flush()