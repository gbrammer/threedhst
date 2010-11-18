"""
3DHST.prep_flt_files

Process RAW flt files to 

    1) Subtract background
    
    2) Align to reference (e.g. ACS-z)
    
"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import numpy as np
import pyfits
import scipy.linalg
import matplotlib.pyplot as plt

from pyraf import iraf
from iraf import stsdas,dither,slitless,axe 

INDEF = iraf.INDEF
no = iraf.no
yes = iraf.yes

import threedhst

class fit_2D_background():
    
    def __init__(self, ORDER=2, grism=False, x0=None, y0=None, DQMAX=10,
        DIRECT_SKY='/research/HST/GRISM/CONF/F140W_GOODSN_median.fits',
        GRISM_SKY='/research/HST/GRISM/CONF/G141_sky_cleaned.fits'):
        """
__init__(self, ORDER=2, grism=False, x0=None, y0=None, DQMAX=10)
    
    ORDER: Polynomial order of the fit, e.g.
        ORDER=2 - 0th order, x, y, x**2, y**2, x*y
        ORDER=3 - 0th order, x, y, x**2, y**2, x*y, x**2 y, x y**2, x**3, y**3
           
    if `grism` is True:
        0th order is median G141 image
    else:
        0th order is median F140W image
        
    x0, y0: reference coordinate for polynomical fit.  Defaults to image center
    
    DQMAX: Maximum value in FLT[DQ] extension considered OK
    
        """
        self.ORDER = ORDER
        self.x0 = x0
        self.y0 = y0
        self.DQMAX = DQMAX
        self.setup_matrices(DIRECT_SKY=DIRECT_SKY, GRISM_SKY=GRISM_SKY)
        if grism:
            self.A = self.B
            
    def setup_matrices(self,
        DIRECT_SKY='/research/HST/GRISM/CONF/F140W_GOODSN_median.fits',
        GRISM_SKY='/research/HST/GRISM/CONF/G141_sky_cleaned.fits'):
        """
setup_matrices()
    
    Setup self.A (direct) and self.B (grism) matrices for polynomial fit.
        """
        NX = 1014
        NY = 1014
        
        #### Image matrix indices
        xi,yi = np.indices((NX,NY))
        
        #### Default reference position is image center
        if self.x0 is None:
            self.x0 = NX/2.
        if self.y0 is None:
            self.y0 = NY/2.
                
        xi = (xi-self.x0*1.)/NX
        yi = (yi-self.y0*1.)/NY
        
        #### 0th order, `grism` is False, is median F140W image
        medF = pyfits.open(DIRECT_SKY)
        med = medF[0].data
        medF.close()

        NPARAM  = np.sum(np.arange(self.ORDER+2))
        self.A = np.zeros((NPARAM,NX,NY))
        self.A[0,:,: ] = med  ## zeroth order is background median

        ##### Set up polynomial grid for fit, including cross terms
        count = 1
        for pow in range(1,self.ORDER+1):
            pi = pow-1

            #### Cross terms
            while (pi > pow/2.):
                self.A[count,:,:] = xi**pi*yi**(pow-pi)
                print 'A[%d,:,:] = xi**%d*yi**%d' %(count,pi,pow-pi)
                count+=1
                self.A[count,:,:] = xi**(pow-pi)*yi**pi
                print 'A[%d,:,:] = xi**%d*yi**%d' %(count,pow-pi,pi)
                count+=1
                pi-=1
            
            #### x**pow/2 * y**pow/2 term
            if (pow/2. == np.int(pow/2.)):
                self.A[count,:,:] = xi**(pow/2)*yi**(pow/2)
                print 'A[%d,:,:] = xi**%d*yi**%d' %(count,pow/2,pow/2)
                count+=1
            
            #### x**pow, y**pow terms
            print 'A[%d,:,:] = xi**%d' %(count,pow)
            self.A[count,:,:] = xi**pow
            count+=1
            print 'A[%d,:,:] = yi**%d' %(count,pow)
            self.A[count,:,:] = yi**pow
            count+=1
        
        #### Oth order for `grism` is True is G141 median image
        #medF = pyfits.open('../CONF/WFC3.IR.G141.sky.V1.0.fits') # from aXe web
        # cleaned of value=0 pixels
        medF = pyfits.open(GRISM_SKY)
        med_g141 = medF[0].data
        medF.close()
        
        self.B = self.A*1.
        self.B[0,:,:] = med_g141*1.
        self.NPARAM = NPARAM
    
    def fit_image(self, root, A=None, overwrite=False, show=True):
        """
fit_image(self, root, A=None, overwrite=False, show=True)
    
    Fit and optionally subtract the background from a FLT image.
    
    `root` is like "ib3728d2q"
    
    `A` is a matrix computed by `setup_matrices` or `__init__`.
    
    if `overwrite` is True:
        Write the background-subtracted image to the original FLT file
        
    if `show` is True:
        Show a plot of the original, background, and bg-subtracted images.
        
        """
        import os
        
        if A is None:
            A = self.A
        
        #### Read the FLT file and get dimensions 
        #### (should always be 1014x1014 for WFC3)
        fi = pyfits.open(root+'_flt.fits',mode='update')
        IMG = fi[1].data
        DQ = fi[3].data
        NX,NY = IMG.shape
        
        #### array indices for each pixel position
        xi,yi = np.indices((NX,NY))
        xi = (xi-NX)/2./NX
        yi = (yi-NY)/2./NY
        
        #### If a segmentation image is available, read it
        #### for use as an object mask.
        if not os.path.exists(root+'_flt.seg.fits'):
            seg = IMG*0.
        else:
            fis = pyfits.open(root+'_flt.seg.fits')
            seg = fis[0].data
        
        #### Apply segmentation mask, also mask out extreme IMG values 
        #### and any pixel with DQ flag > self.DQMAX
        
        q = np.where((seg == 0) & (IMG > -1) & (IMG < 4) & (DQ < self.DQMAX)) 
        qb = np.where((seg > 0) | (IMG < -1) | (IMG > 4) | (DQ >= self.DQMAX))
        IMGb = IMG*1.
        IMGb[qb] = np.nan
        
        #### Apply mask to IMG and fit matrices
        Aq = np.transpose(A[:,q[0],q[1]])
        IMGq = IMG[q[0],q[1]]
        
        #### Get fit parameters with least-sq. fit.
        p, resid, rank, s = scipy.linalg.lstsq(Aq,IMGq)

        print p
        
        #### Create the bg fit image from the fit parameters
        IMGout = IMG*0.
        for i in range(self.NPARAM):
            IMGout += A[i,:,:]*p[i]
        print 'Done'
        
        #### Save fit parameters to an ASCII file
        fp = open(root+'_flt.polybg','w')
        for pi in p:
            fp.write('%13.5e\n' %pi)
        fp.close()
        
        #### Show the results, note that subsequent 
        #### plots aren't cleared from memory, so memory 
        #### fills up quickly with repeated calls with show=True.
        if show:
            dshow = 0.3
            plt.figure()
            plt.subplot(221)
            plt.imshow(IMGb,vmin=np.median(IMGb)-dshow,
                vmax=np.median(IMGb)+dshow)
            plt.subplot(222)
            plt.imshow(IMG-IMGout,vmin=0-dshow,vmax=dshow)
            plt.subplot(223)
            plt.imshow(IMGout,vmin=np.median(IMGb)-dshow,
                vmax=np.median(IMGb)+dshow)

            plt.subplot(224)
            plt.imshow(DQ,vmin=0,vmax=10)
        
        #### Subtract fitted background, write
        #### bg-subtracted image to original FLT file if `overwrite` is True
        FIX = IMG-IMGout
        if overwrite:
            print 'Overwrite: '+root
            fi[1].data = FIX
            fi.flush()
        
        #### Save images to self.
        self.FIX = FIX
        self.IMG = IMG
        self.MODEL = IMGout
        fi.close()
    
def make_grism_shiftfiles(direct_files='ib*050_asn.fits', 
                          grism_files='ib*060_asn.fits'):
    """
make_grism_shiftfiles(direct_files='ib*050_asn.fits', 
                      grism_files='ib*060_asn.fits')
    
    Find all of the shiftfiles determined for the direct images
    and make corresponding versions for the G141 associations.
    """
    import threedhst
    import glob
    asn_direct_files = glob.glob(direct_files)
    asn_grism_files = glob.glob(grism_files)
    for i, asn_direct_file in enumerate(asn_direct_files):
        asn_grism_file=asn_grism_files[i]
        # asn_grism_file = asn_direct_file.split('50_asn')[0]+'60_asn.fits'
        threedhst.shifts.make_grism_shiftfile(asn_direct_file, asn_grism_file)
    
def prep_all(asn_files='ib*050_asn.fits', get_shift=True, bg_skip=False,
             redo_background=True, bg_only=False, grism=False,
             ALIGN_IMAGE='../ACS/h_nz_sect*img.fits',
             skip_drz=False,final_scale=0.06, pixfrac=0.8,
             DIRECT_SKY='/research/HST/GRISM/CONF/F140W_GOODSN_median.fits',
             GRISM_SKY='/research/HST/GRISM/CONF/G141_sky_cleaned.fits',
             align_geometry='rxyscale,shift', clean=True):
    """
prep_all(asn_files='ib*050_asn.fits', get_shift=True, 
         redo_background=True, bg_only=False, grism=False)
    
    asn_files = glob.glob(asn_files)
    
    Run prep_flt on all direct or grism associations in the current directory.
    See `prep_flt` for parameter descriptions.
    
    """
    import glob
    # asn_files = glob.glob('ib*050_asn.fits')
    # if grism:
    #     asn_files = glob.glob('ib*060_asn.fits')
    asn_files = glob.glob(asn_files)
       
    for file in asn_files:
        prep_flt(asn_file=file, get_shift=get_shift, bg_skip=bg_skip,
                    redo_background=redo_background,
                    bg_only=bg_only, grism=grism, ALIGN_IMAGE=ALIGN_IMAGE,
                    skip_drz=skip_drz,final_scale=final_scale, pixfrac=pixfrac,
                    DIRECT_SKY=DIRECT_SKY, GRISM_SKY=GRISM_SKY, 
                    align_geometry=align_geometry, clean=clean)

def prep_flt(asn_file=None, get_shift=True, bg_only=False, bg_skip=False,
                redo_background=True, grism=False,
                ALIGN_IMAGE='../ACS/h_nz_sect*img.fits',
                skip_drz=False, final_scale=0.06, pixfrac=0.8,
                DIRECT_SKY='/research/HST/GRISM/CONF/F140W_GOODSN_median.fits',
                GRISM_SKY='/research/HST/GRISM/CONF/G141_sky_cleaned.fits',
                align_geometry='rxyscale,shift', clean=True):
    """
prep_flt(asn_file=None, get_shift=True, bg_only=False,
            redo_background=True, grism=False)

    
    Subtract background and align WCS of direct/grism FLT files.
    
    1) Apply the DQ masks defined in the *mask.reg files, as created
       by threedhst.dq
    
    2) First pass on background subtraction 
        
        o 0th order is median background image, either for 
          F140W or G141 [if grism is True]
        
        o 1-nth order is polynomial fit with x-y cross terms.
    
        o [if `bg_only` is True then return]
        
    3) Run tweakshifts  [if `get_shift` is True & `grism` is False]
    
    4) Run Multidrizzle with first guess tweakshifts
    
    5) Get object mask for FLT files:
    
        o Blot DRZ output to all the individual FLT frames
        
        o Run SExtractor on blot images for object (segmentation)
          mask  
    
    6) Use threedhst routines to align the F140W direct image to
       ACS reference.   [if `get_shift` is True & `grism` is False]
       
    7) Redo background subtraction with improved object mask defined
       by the segmentation images [if redo_background is True]
       
    8) Run multidrizzle again [if redo_background is True]
    
    
    """
    #import fit_2d_poly
    import threedhst
    #import threedhst.dq    
    
    if asn_file is None:
        asn_file = 'ib3728050_asn.fits'
    
    if grism:
        get_shift=False
    
    if bg_skip:
        bg_only=False
        redo_background=False
        
    ROOT_DIRECT = asn_file.split('_asn.fits')[0]
    # ALIGN_IMAGE = '../ACS/h_nz_sect*img.fits'
    
    asn = threedhst.utils.ASNFile(asn_file)
    
    #### Set up matrix for fitting
    fit = fit_2D_background(ORDER=0, grism=grism, DIRECT_SKY=DIRECT_SKY,
                            GRISM_SKY=GRISM_SKY)#, x0=507, y0=507)
    
    #### First pass background subtraction
    if not bg_skip:
        for exp in asn.exposures:
            threedhst.regions.apply_dq_mask(exp+'_flt.fits')
            fit.fit_image(exp, A=fit.A, show=False, overwrite=True)
    
    #### Stop here if only want background subtraction
    if bg_only:
        return
        
    #### First guess at shifts
    if get_shift:
        threedhst.shifts.run_tweakshifts(asn_file, verbose=True)
        threedhst.shifts.checkShiftfile(asn_file)
        
    if not skip_drz:
        startMultidrizzle(asn_file, use_shiftfile=True, 
            skysub=bg_skip, final_scale=final_scale, pixfrac=pixfrac)
        
    #### Blot combined images back to reference frame and make a 
    #### segmentation mask
    run = MultidrizzleRun((asn_file.split('_asn.fits')[0]).upper())
    
    #### ACS has entries in run file for each of two WFC chips
    flt = pyfits.open(asn.exposures[0]+'_flt.fits')
    inst = flt[0].header.get('INSTRUME').strip()
    if inst is 'ACS':
        skip=2
    else:
        skip=1
        
    for i,exp in enumerate(asn.exposures):
        run.blot_back(ii=i*skip, copy_new=(i is 0))
        make_segmap(run.flt[i])
        
    if get_shift:
        #### If shift routine gets confused, run the following instead
        #for geom in ['shift','rxyscale','shift']:
        # for geom in ['rxyscale','shift']:
        for geom in align_geometry.split(','):
            refine_shifts(ROOT_DIRECT=ROOT_DIRECT,
                          ALIGN_IMAGE=ALIGN_IMAGE,
                          fitgeometry=geom.strip(), clean=clean)
            startMultidrizzle(asn_file, use_shiftfile=True, skysub=False,
                final_scale=final_scale, pixfrac=pixfrac, driz_cr=False,
                updatewcs=False)
                
        
    #### Run BG subtraction with improved mask and run multidrizzle again
    if redo_background:
        fit = fit_2D_background(ORDER=1, grism=grism, DIRECT_SKY=DIRECT_SKY,
                                GRISM_SKY=GRISM_SKY)#, x0=507, y0=507)
        for exp in asn.exposures:
            fit.fit_image(exp, A=fit.A, show=False, overwrite=True)
        
        startMultidrizzle(asn_file, use_shiftfile=True, skysub=False,
            final_scale=final_scale, pixfrac=pixfrac, driz_cr=False,
            updatewcs=False)
    
#
def refine_shifts(ROOT_DIRECT='f160w',
                  ALIGN_IMAGE='../../ACS/h_sz*drz_img.fits',
                  fitgeometry='shift', clean=True):
    """
refine_shifts(ROOT_DIRECT='f160w',
              ALIGN_IMAGE='../../ACS/h_sz*drz_img.fits',
              fitgeometry='shift')
                
    Refine shifts by catalog matching an input multidrizzle image, 
    ROOT_DIRECT+'_drz.fits' to one or more alignment images
    """
    run = MultidrizzleRun(ROOT_DIRECT.upper())
    
    xshift, yshift, rot, scale = threedhst.shifts.align_to_reference(
                        ROOT_DIRECT,
                        ALIGN_IMAGE,
                        fitgeometry=fitgeometry, clean=clean)

    #### shifts measured in DRZ frame.  Translate to FLT frame
    drz = pyfits.open(ROOT_DIRECT+'_drz.fits')
    alpha = (180.-drz[1].header['PA_APER'])/360.*2*np.pi
    
    xsh = (xshift*np.cos(alpha)-yshift*np.sin(alpha))*np.float(run.scl)
    ysh = (xshift*np.sin(alpha)+yshift*np.cos(alpha))*np.float(run.scl)

    print 'Final shift:', xsh, ysh, drz[1].header['PA_APER']
    fp = open(ROOT_DIRECT+'_align.info','w')
    fp.write('%s %8.3f %8.3f %8.3f\n' %(ALIGN_IMAGE, xsh, ysh, rot)) 
    fp.close()
    
    #### Read the shiftfile       
    shiftF = threedhst.shifts.ShiftFile(ROOT_DIRECT+'_shifts.txt')
    
    #### Apply the alignment shifts to the shiftfile
    shiftF.xshift = list(np.array(shiftF.xshift)-xsh)
    shiftF.yshift = list(np.array(shiftF.yshift)-ysh)
    shiftF.rotate = list((np.array(shiftF.rotate)+rot) % 360)
    shiftF.scale = list(np.array(shiftF.scale)*scale)
    
    shiftF.print_shiftfile(ROOT_DIRECT+'_shifts.txt')

def startMultidrizzle(root='ib3727050_asn.fits', use_shiftfile = True,
    skysub=True, updatewcs=True, driz_cr=True, median=True,
    final_scale=0.06, pixfrac=0.8, 
    final_outnx='', final_outny='', final_rot=0., ra='', dec=''):
    """
startMultidrizzle(root='ib3727050_asn.fits', use_shiftfile = True,
                  skysub=True, final_scale=0.06, updatewcs=True, driz_cr=True,
                  median=True, final_scale=0.06, pixfrac=0.8, 
                  final_outnx='', final_outny='', final_rot=0., ra='', dec='')
    
    Run multidrizzle on an input asn table.
    
    if `use_shiftfile` is True:
        use a root+'_shifts.txt' shiftfile.
    else: 
        no shiftfile
    
    if skysub is True:
        Run multidrizzle WITH sky subtraction
    else:
        "        "       WITHOUT   "         "
        
    final_scale: Final pixel scale of output image (arcsec)
    
    """
    asn_direct_file = root #'ib3727050_asn.fits'
    
    asn_direct = threedhst.utils.ASNFile(file=asn_direct_file)
    ROOT_DIRECT = asn_direct_file.split('_asn')[0]
    
    if use_shiftfile:
        shiftfile=ROOT_DIRECT+'_shifts.txt'
    else:
        shiftfile=''
    
    #### If fewer than 4 exposures in the asn list, use
    #### a larger `pixfrac`.
    if len(asn_direct.exposures) < 4:
        pixfrac = 1.0
    
    if skysub:
        skysub=yes
    else:
        skysub=no
    
    if updatewcs:
        updatewcs=yes
    else:
        updatewcs=no
    
    if driz_cr:
        driz_cr=yes
    else:
        driz_cr=no

    if median:
        median=yes
    else:
        median=no
    
    #### Run Multidrizzle
    iraf.multidrizzle(input=asn_direct_file, \
       shiftfile=shiftfile, \
       output = '', skysub = skysub, updatewcs = updatewcs, driz_cr=driz_cr,
       final_scale = final_scale, final_pixfrac = pixfrac, median=median, 
       blot=median,
       final_outnx=final_outnx, final_outny=final_outny, 
       final_rot=final_rot, ra=ra, dec=dec)
    
    #### Delete created files    
    threedhst.process_grism.cleanMultidrizzleOutput()
    
# def go_blot_all():
#     import glob
#     
#     files=glob.glob('IB*run')
#     for file in files:
#         run = MultidrizzleRun(root=file.split('.run')[0])
#         for i in range(4):
#             try:
#                 run.blot_back(ii=i)
#             except:
#                 pass
# 
# def blot_g141(asn_file=None):
#     import os
#     if not asn_file:
#         asn_file = 'ib3704060_asn.fits'
#     
#     iraf.multidrizzle.static=no
#     iraf.multidrizzle.skysub=no
#     iraf.multidrizzle.driz_separate=no
#     iraf.multidrizzle.median=no
#     iraf.multidrizzle.median_newmasks=no
#     iraf.multidrizzle.blot=no
#     iraf.multidrizzle.driz_cr=no
#     iraf.multidrizzle.driz_combine=yes
#     startMultidrizzle(root=asn_file, use_shiftfile=False)
#     asn_root = (asn_file.split('_asn.fits')[0]).upper()
#     os.remove(asn_root+'_drz.fits')
#     
#     asn = threedhst.utils.ASNFile(asn_file)
#     #direct = (asn_file.split('060_')[0]+'050').upper()
#     run = MultidrizzleRun(root=asn_root)
#     
#     run.root = '../DATA/'+asn_file.split('_asn')[0]+'CONT'
#     # for i,exp in enumerate(asn.exposures):
#     #     run.flt[i] = exp+'_flt'
#     for ii,flt in enumerate(run.flt):
#         run.blot_back(ii=ii)
    
class MultidrizzleRun():
    """
MultidrizzleRun(root='IB3728050')
    
    Read a .run file output from MultiDrizzle.
    
    Get list of flt files and their shifts as used by multidrizzle.
    """
    def __init__(self, root='IB3728050'):
        import numpy as np
        
        runfile = root+'.run'
        self.root = root
        
        self.flt = []
        self.xsh = []
        self.ysh = []
        self.rot = []
        self.scl = 1.
        
        for line in open(runfile,'r'):
            if line.startswith('drizzle.scale'):
                self.scl = line.split()[2]
                
            if line.startswith('drizzle '):
                spl = line.split()
                self.flt.append(spl[1].split('.fits')[0])
                for tag in spl:
                    if tag.startswith('xsh'):
                        self.xsh.append(np.float(tag.split('=')[1]))
                    if tag.startswith('ysh'):
                        self.ysh.append(np.float(tag.split('=')[1]))
                    if tag.startswith('rot'):
                        self.rot.append(np.float(tag.split('=')[1]))
                        
    def blot_back(self, ii=0, SCI=True, WHT=True, copy_new=True):
        """
blot_back(self, ii=0, SCI=True, WHT=True, copy_new=True)
    
    Blot the output DRZ file back to exposure #ii pixels.
    
    if SCI is True:
        blot science extension to FLT+'.BLOT.SCI.fits'

    if WHT is True:
        blot weight extension to FLT+'.BLOT.WHT.fits'
    
    if copy_new is True:
        imcopy SCI and WHT extensions of DRZ image to separate files.
        
        """
        #flt_orig = pyfits.open('../RAW/'+self.flt[ii]+'.fits.gz')
        threedhst.process_grism.flprMulti()
        
        flt_orig = pyfits.open(self.flt[ii]+'.fits')
        exptime = flt_orig[0].header.get('EXPTIME')
        flt_orig.close()
        
        inNX = flt_orig[1].header.get('NAXIS1')
        inNY = flt_orig[1].header.get('NAXIS2')
        
        iraf.delete(self.flt[ii]+'.BLOT.*.fits')
        if copy_new:
            iraf.delete('drz_*.fits')
            # iraf.imcopy(self.root+'_drz.fits[1]','drz_sci.fits')
            # iraf.imcopy(self.root+'_drz.fits[2]','drz_wht.fits')
            
            ### NEED TO STRIP FITS HEADER
            im_drz = pyfits.open(self.root+'_drz.fits')
            sci = im_drz[1].data            
            s_hdu = pyfits.PrimaryHDU(sci)
            s_list = pyfits.HDUList([s_hdu])
            copy_keys = ['CTYPE1','CTYPE2','CRVAL1','CRVAL2','CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2','LTM1_1','LTM2_2']
            s_list[0].header.update('EXPTIME',im_drz[0].header.get('EXPTIME'))
            s_list[0].header.update('CDELT1',im_drz[1].header.get('CD1_1'))
            s_list[0].header.update('CDELT2',im_drz[1].header.get('CD2_2'))
            for key in copy_keys:
                s_list[0].header.update(key, im_drz[1].header.get(key))
            s_list.writeto('drz_sci.fits', clobber=True)
            
            wht = im_drz[2].data
            w_hdu = pyfits.PrimaryHDU(sci)
            w_list = pyfits.HDUList([w_hdu])
            copy_keys = ['CTYPE1','CTYPE2','CRVAL1','CRVAL2','CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2','LTM1_1','LTM2_2']
            w_list[0].header.update('EXPTIME',im_drz[0].header.get('EXPTIME'))
            w_list[0].header.update('CDELT1',im_drz[1].header.get('CD1_1'))
            w_list[0].header.update('CDELT2',im_drz[1].header.get('CD2_2'))
            for key in copy_keys:
                w_list[0].header.update(key, im_drz[1].header.get(key))
            w_list.writeto('drz_wht.fits', clobber=True)
            
        if SCI:
            iraf.blot(data='drz_sci.fits',
                outdata=self.flt[ii]+'.BLOT.SCI.fits', scale=self.scl,
                coeffs=self.flt[ii]+'_coeffs1.dat', xsh=self.xsh[ii], 
                ysh=self.ysh[ii], 
                rot=self.rot[ii], outnx=inNX, outny=inNY, align='center', 
                shft_un='input', shft_fr='input', in_un='cps', out_un='cps', 
                interpol='poly5', sinscl='1.0', expout=exptime, 
                expkey='EXPTIME',fillval=0.0)
        
        if WHT:
            iraf.blot(data='drz_wht.fits',
                outdata=self.flt[ii]+'.BLOT.WHT.fits', scale=self.scl,
                coeffs=self.flt[ii]+'_coeffs1.dat', xsh=self.xsh[ii], 
                ysh=self.ysh[ii], 
                rot=self.rot[ii], outnx=inNX, outny=inNY, align='center', 
                shft_un='input', shft_fr='input', in_un='cps', out_un='cps', 
                interpol='poly5', sinscl='1.0', expout=exptime, 
                expkey='EXPTIME',fillval=0.0)
                
def jitter_info():
    """
jitter_info()
    
    Get LIMBANG values from jitter files and also get 
    image stats from FLT images.  Useful for flagging 
    exposures affected by earthglow.
    """
    import glob
    import os
    import numpy as np
    
    fp = open('jitter_info.dat','w')
    
    jit_files = glob.glob('../JITTER/*jit.fits')
    jit_files = glob.glob('../RAW/ib*0_asn.fits')
    
    for file in jit_files:
        #im = pyfits.open(file)
        im = pyfits.open('../JITTER/'+ 
                         os.path.basename(file).split('_')[0]+
                         '_jit.fits')
        nExten = len(im)-1
        for ext in range(1,nExten+1):
            head = im[ext].header
            dat = im[ext].data
            flt = pyfits.open('../RAW/'+head.get('EXPNAME')[:-1]+
                              'q_flt.fits.gz')
            
            med = np.median(flt[1].data[180:300,180:300])
            
            str = ("%s %10s %8.1f %5.1f %2d %5.1f %2d  %8.2f"
               %(head.get('EXPNAME')[:-1]+'q',
               flt[0].header['FILTER'], flt[0].header['EXPTIME'],
               dat.field('LimbAng')[0], dat.field('BrightLimb')[0],
               dat.field('LimbAng')[-1], dat.field('BrightLimb')[-1],med))
            
            fp.write(str+'\n')
            print str
    
    fp.close()
    
# def go_make_segmap():
#     
#     import glob
#     
#     files = glob.glob('*BLOT.SCI.fits')
#     for file in files:
#         make_segmap(root=file.split('.BLOT')[0])
        
def make_segmap(root='ib3701ryq_flt', sigma=1):
    """
make_segmap(root='ib3701ryq_flt', sigma=1)
    
    Get a segmentation image for a flt file after creating its 
    BLOT SCI and WHT images.
    
    DETECT_THRESH = ANALYSIS_THRESH = sigma
    """
    import threedhst
    
    se = threedhst.sex.SExtractor()
    ## Set the output parameters required for aXe 
    ## (stored in [threedhst source]/data/aXe.param) 
    se.aXeParams()
    ## XXX add test for user-defined .conv file
    se.copyConvFile()
    
    se.overwrite = True
    se.options['CATALOG_NAME']    = root+'.BLOT.SCI.cat'
    se.options['CHECKIMAGE_NAME'] = root+'.seg.fits, bg.fits'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION, BACKGROUND'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = root+'.BLOT.WHT.fits'
    se.options['FILTER']    = 'Y'

    se.options['BACK_TYPE']     = 'AUTO'
    se.options['BACK_FILTERSIZE']     = '2'
    
    #### Detect thresholds (default = 1.5)
    se.options['DETECT_THRESH']    = '%f' %sigma
    se.options['ANALYSIS_THRESH']  = '%f' %sigma
    se.options['MAG_ZEROPOINT'] = '26.46'
    status = se.sextractImage(root+'.BLOT.SCI.fits')
    
