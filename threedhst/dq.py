"""
3DHST.dq

Check DQ (satellite trails, scattered light)
of direct/grism exposures.


"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import os

import tkinter as tk

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import numpy as np

import threedhst

try:
    import pysao
    threedhst.options['PYSAO_INSTALLED'] = True
except ImportError:
    try:
        import pyds9 as pysao
    except:
        print('No pysao installation found.')
        threedhst.options['PYSAO_INSTALLED'] = False
    
class myDS9(pysao.ds9):
    """
myDS9(pysao.ds9)
    
    Add shortcut methods to extend `pysao.ds9`
    """
    def v(self, img, vmin=0, vmax=3):
        self.view(img)
        self.set('scale limits %f %f' %(vmin,vmax))
        
    def scale(self, min, max):
        """
scale(self, min, max)
    
    ds9.set('scale limits min max')
        """
        self.set('scale limits %f %f' %(min,max))
    #
    def match_all(self, mode='wcs'):
        """
match_all(self, alignType='wcs')

    ds9.set('match frames `alignType`')
    ds9.set('match colorbars')
    ds9.set('match scales')
        """
        self.set('match frames %s' %(mode))
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
    
    def cds_query(self, radius=1):
        """
        Search the CDS database around the center position of the DS9 window
        """
        import os
        coords = self.get('pan fk5 sexagesimal')
        if coords == ' 0 0 \n':
            print('No WCS information.  Load the image with a header.')
            return False
            
        format='+'.join(coords.replace('+','%2B').replace('-','%2D').split())
        link="http://vizier.u-strasbg.fr/viz-bin/VizieR?-c=%s&-c.rs=%.1f" %(format, radius)
        print(coords)
        os.system('open \"%s\"' %(link))
        
class checkDQ:
    """
checkDQ(asn_direct_file='ib3704050_asn.fits',
        asn_grism_file='ib3704060_asn.fits',
        path_to_flt='../RAW/'):

    `master` is an instance of Tkinter.tk()
    
    Widget application.
    
    Just display all images in ds9 and mask files if they exist:
    
    ds9 & 
    files=`grep flt files.info |awk '{print $1}'`
    for f in $files; do  newframe ../RAW/${f}; reg=`echo $f |sed "s/gz/mask.reg/"`; ds9_reg $reg; done
    # make polygon region
    xpaset -p ds9 regions save xxx_flt.fits.mask.reg
    
    """
    
    def __init__(self, asn_direct_file='ib3704050_asn.fits',
           asn_grism_file='ib3704060_asn.fits',
           path_to_flt='../RAW/', SIMPLE_DS9=True, wait_time=30, size=600):
           
        master = tk.Tk()
        self.master = master
        master.title('3D-HST')
        master.geometry('200x90-350+40')
        frame = tk.Frame(master)
        frame.pack()
        self.frame = frame
        self.path_to_flt = path_to_flt
        
        # maybe want to put ds9 outside of the method
        # to preserve state between subsequent ASN files
        self.ds9 = myDS9(wait_time=wait_time)
        self.ds9.set('width %d' %(size))
        self.ds9.set('height %d' %(size+200))
        
        if SIMPLE_DS9:
            self.ds9.set('view panner no')
            self.ds9.set('view magnifier no')
            
        self.asn_grism = threedhst.utils.ASNFile(asn_grism_file)
        self.asn_direct = threedhst.utils.ASNFile(asn_direct_file)
        print('=== Grism ===')
        self.asn_grism.showContents()
        print('=== Direct ===')
        self.asn_direct.showContents()
        
        self.nexp = len(self.asn_grism.exposures)
        self.idx = 0
        self.current = None  # set if editing grism or direct
        self.showExposures()
        self.busy = False
        
        #### Buttons
        self.button_quit = tk.Button(frame, 
            text="QUIT", fg="red", command=self.finish)
        self.button_quit.pack(side=tk.LEFT)
        
        self.button_next = tk.Button(frame,
            text="NEXT", bg="green", command=self.goNext)
        self.button_next.pack(side=tk.LEFT)

        self.button_back = tk.Button(frame,
            text="BACK", bg="yellow", command=self.goBack)
        self.button_back.pack(side=tk.LEFT)
        
        self.button_kill = tk.Button(frame,
            text="Kill exposure", bg="red", command=self.kill_exposure)
        self.button_kill.pack(side=tk.LEFT)
        
        self.button_dqdirect = tk.Button(frame,
            text="Flag DIRECT", bg="white", fg="blue", command=self.dqdirect)
        self.button_dqdirect.pack(side=tk.LEFT)
        
        self.button_dqgrism = tk.Button(frame,
            text="Flag GRISM", bg="blue", fg="white", command=self.dqgrism)
        self.button_dqgrism.pack(side=tk.LEFT)
        
        self.button_next.grid(row=0, column=0)
        self.button_back.grid(row=1, column=0)
        self.button_kill.grid(row=2, column=0)
        self.button_dqdirect.grid(row=0, column=1)
        self.button_dqgrism.grid(row=1, column=1)
        self.button_quit.grid(row=2, column=1)
        
        self.master.mainloop()
        
    def goNext(self):
        """
goNext(self)
    
    Go to next image in ASN file.
        """
        self.idx+=1
        if (self.idx == self.nexp):
            self.finish()
        else:
            self.showExposures()
    
    def goBack(self):
        """
goBack(self)

    Go back one image in ASN file.
        """
        if self.idx == 0:
            return None
        
        self.idx-=1
        self.showExposures()
    
    def kill_exposure(self):
        """
kill_exposure(self)

    Remove current exposure from the ASN lists.
        """
        out = self.asn_direct.exposures.pop(self.idx)
        out = self.asn_grism.exposures.pop(self.idx)
        self.idx-=1
        self.nexp-=1
        self.goNext()
        
    def dqdirect(self):
        """
dqdirect(self)
        """
        if not self.busy:
            self.current = 'DIRECT'
        else:
            ### If in the middle of defining polygons
            ### and other button clicked, do nothing
            if self.current == 'GRISM':
                return None
                
        self.readMask()
        # do things ... run polygon mask routine for direct image
        
    def dqgrism(self):
        """
dqgrism(self)
        """
        if not self.busy:
            self.current = 'GRISM'
        else:
            ### If in the middle of defining polygons
            ### and other button clicked, do nothing
            if self.current == 'DIRECT':
                return None
                
        self.readMask()
        # do things ...  run polygon mask routine for grism image
    
    def readMask(self):
        """
readMask(self)
    
    Get polygon regions from current DS9 and save a mask region file.
        """
        ### First click on region button
        if not self.busy:
            self.busy = True
            self.ds9.set('tile no')
            if self.current == 'DIRECT':
                self.ds9.frame(1)
                self.polyfile_root = self.asn_direct.exposures[self.idx]
            else:
                self.ds9.frame(3)
                self.polyfile_root = self.asn_grism.exposures[self.idx]
            
            self.ds9.set('zoom to fit')
            self.ds9.set('regions shape polygon')
            print('Define region polygon(s) in DS9 and click again ')
            
        else:
            self.busy = False
            self.ds9.set('regions system image')
            regions = self.ds9.get('regions source')
            
            #### If polygon regions defined, write them to the mask file
            if len(regions.split('polygon')) > 1:
                fp = open(self.polyfile_root+'_flt.fits.mask.reg','w')
                fp.write(regions)
                fp.close()
                print('Wrote mask, %s_flt.fits.mask.reg' %self.polyfile_root)
            
            self.ds9.set('tile yes')
            self.ds9.frame(3)
            self.ds9.set('zoom to fit')
            self.ds9.set('match frames image')
            
    def showExposures(self):
        """
showExposure(self)

    Display Direct and Grism exposures in ds9.
        """
        # open 4 frames and show grism/direct, SCI/DQ extensions
        # Use the `self.idx`th object from the ASN list
        flt_grism = self.asn_grism.exposures[self.idx]
        flt_direct = self.asn_direct.exposures[self.idx]
        
        fits_grism = threedhst.utils.find_fits_gz(self.path_to_flt+
                                                  flt_grism+'_flt.fits')
        fits_direct = threedhst.utils.find_fits_gz(self.path_to_flt+
                                                   flt_direct+'_flt.fits')
        
        fi_grism = pyfits.open(fits_grism)
        fi_direct = pyfits.open(fits_direct)

        ### Display SCI extension [1]
        self.ds9.frame(1)
        self.ds9.view(fi_direct[1].data)
        self.ds9.scale(-0.1,2)
        if os.path.exists(flt_direct+'_flt.fits.mask.reg'):
            self.ds9.set('regions load %s_flt.fits.mask.reg' %flt_direct)
            
        ### Display DQ extension [3]
        self.ds9.frame(2)
        self.ds9.view(fi_direct[3].data)
        self.ds9.scale(0,100)

        ### Display SCI extension [1]
        self.ds9.frame(3)
        self.ds9.view(fi_grism[1].data-np.median(fi_grism[1].data))
        self.ds9.scale(-0.1,0.6)
        if os.path.exists(flt_grism+'_flt.fits.mask.reg'):
            self.ds9.set('regions load %s_flt.fits.mask.reg' %flt_grism)
        
        ### Display DQ extension [3]
        self.ds9.frame(4)
        self.ds9.view(fi_grism[3].data)
        self.ds9.scale(0,100)

        ### DS9 Tile and return to Frame 1
        self.ds9.set('tile yes')
        self.ds9.set('tile grid')

        self.ds9.frame(3)
        self.ds9.set('zoom to fit')
        self.ds9.set('match frames image')
    
    def finish(self):
        """
finish(self)
        
    Write out ASN files and close widget.
        """
        ### write asn files
        self.asn_grism.write('self')
        self.asn_direct.write('self')
        
        ### close widget
        self.frame.quit()
        self.master.destroy()
        del(self.ds9)

def checkAllDQ(clobber=False, path_to_flt='../RAW', copy_new_files=True):
    """
checkAllDQ(clobber=False)
    
    Copy ASN files from ../RAW to ./ [DATA].  If clobber==False, 
    don't overwrite ASN files in DATA.
    
    Figure out which ASN are grism, which are direct and run checkDQ.
    
    This should be run in the `DATA` directory.
    
    Assumes that grism/direct pairs will be next to eachother 
    in the list of *asn.fits.
    """
    import glob
    import shutil
    import os
    
    old_pwd = os.getcwd()+''
    
    if copy_new_files:
        os.chdir(path_to_flt)
        asn_files = glob.glob('i*asn.fits')
        for file in asn_files:
            if (clobber is True) | (not os.path.exists('../DATA/'+file)):
                shutil.copy(file,old_pwd)
                print('Copy %s to %s' %(file, old_pwd))
            
        os.chdir(old_pwd)
    
    asn_files = glob.glob('i*asn.fits')
    
    direct_files = []
    grism_files = []
    for file in asn_files:
        print(file)
        asn = threedhst.utils.ASNFile(file)
        fits = threedhst.utils.find_fits_gz(path_to_flt+asn.exposures[0]+
                                            '_flt.fits')
        fi = pyfits.open(fits)
        head = fi[0].header
        if head['FILTER'] == 'F140W':
            direct_files.append(file)
        if head['FILTER'] == 'G141':
            grism_files.append(file)
    
    if len(direct_files) != len(grism_files):
        print("""
    checkAllDQ: Number of direct ASN files is not the same as the
                number of grism ASN files!'
            """)
        return None
    
    Npairs = len(direct_files)
    for i in range(Npairs):
        checkDQ(asn_grism_file=grism_files[i],
                asn_direct_file=direct_files[i],
                path_to_flt=path_to_flt)
                
