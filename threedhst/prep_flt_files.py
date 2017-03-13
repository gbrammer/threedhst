"""
3DHST.prep_flt_files

Process RAW flt files to 

    1) Subtract background
    
    2) Align to reference (e.g. ACS-z)
    
    3) Combine full direct mosaics and mosaics of grism pointings with 
       like orientation angles.
    
"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import os
import glob
import shutil

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import numpy as np

import scipy.linalg
# from scipy import polyfit, polyval
import matplotlib.pyplot as plt
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

USE_PLOT_GUI = False

import threedhst
import threedhst.grism_sky

class fit_2D_background():
    
    def __init__(self, ORDER=-1, x0=None, y0=None, DQMAX=10,
        IMAGES=['/research/HST/GRISM/CONF/G141_sky_cleaned.fits'], NX=1014, NY=1014):
        """
__init__(self, ORDER=-1, x0=None, y0=None, DQMAX=10,
         IMAGES=['/research/HST/GRISM/CONF/G141_sky_cleaned.fits'])
    
    ORDER: Polynomial order of the fit, e.g.
        ORDER=2 - 0th order, x, y, x**2, y**2, x*y
        ORDER=3 - 0th order, x, y, x**2, y**2, x*y, x**2 y, x y**2, x**3, y**3
           
    x0, y0: reference coordinate for polynomical fit.  Defaults to image center
    
    DQMAX: Maximum value in FLT[DQ] extension considered OK
    
    IMAGES: Include actual images in the model that is fit.  This can be multiple files, such as the various background images for different sky levels.
    
        """
        self.ORDER = ORDER
        self.x0 = x0
        self.y0 = y0
        self.NX = NX
        self.NY = NY
        self.DQMAX = DQMAX
        self.IMAGES = IMAGES
        self.setup_matrices()
            
    def setup_matrices(self):
        """
setup_matrices()
    
    Setup self.A matrix for polynomial fit.
        """
        NX = self.NX #1014
        NY = self.NY #1014
        
        #### Image matrix indices
        yi,xi = np.indices((NY,NX))
        
        #### Default reference position is image center
        if self.x0 is None:
            self.x0 = NX/2.
        if self.y0 is None:
            self.y0 = NY/2.
                
        xi = (xi-self.x0*1.)/NX
        yi = (yi-self.y0*1.)/NY
        
        NPARAM  = np.sum(np.arange(self.ORDER+2)) #+1 #+1
        NPARAM += len(self.IMAGES)
        self.NPARAM = NPARAM
        
        self.A = np.zeros((NPARAM,NY,NX))
        
        #### Read images to add to the "model"
        count=0
        for img in self.IMAGES:
            hdu = pyfits.open(img)
            self.A[count,:,: ] = hdu[0].data
            hdu.close()
            count += 1
        
        #### zeroth order, flat background
        if self.ORDER >= 0:
            self.A[count,:,:] += 1
            count+=1
        
        for pow in range(1,self.ORDER+1):
            pi = pow-1

            #### Cross terms
            while (pi > pow/2.):
                self.A[count,:,:] = xi**pi*yi**(pow-pi)
                print('A[%d,:,:] = xi**%d*yi**%d' %(count,pi,pow-pi))
                count+=1
                self.A[count,:,:] = xi**(pow-pi)*yi**pi
                print('A[%d,:,:] = xi**%d*yi**%d' %(count,pow-pi,pi))
                count+=1
                pi-=1
            
            #### x**pow/2 * y**pow/2 term
            if (pow/2. == np.int(pow/2.)):
                self.A[count,:,:] = xi**(pow/2)*yi**(pow/2)
                print('A[%d,:,:] = xi**%d*yi**%d' %(count,pow/2,pow/2))
                count+=1
            
            #### x**pow, y**pow terms
            print('A[%d,:,:] = xi**%d' %(count,pow))
            self.A[count,:,:] = xi**pow
            count+=1
            print('A[%d,:,:] = yi**%d' %(count,pow))
            self.A[count,:,:] = yi**pow
            count+=1
        
        # #### Oth order for `grism` is True is G141 median image
        # #medF = pyfits.open('../CONF/WFC3.IR.G141.sky.V1.0.fits') # from aXe web
        # # cleaned of value=0 pixels
        # medF = pyfits.open(GRISM_SKY)
        # med_g141 = medF[0].data
        # medF.close()
        # 
        # self.B = self.A*1.
        # self.B[0,:,:] = med_g141*1.
        # self.NPARAM = NPARAM
    
    def fit_image(self, root, A=None, overwrite=False, show=True,
                  save_fit=False, root_filename=True):
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
        import glob
        import scipy.ndimage as nd
        
        if A is None:
            A = self.A
        
        #### Read the FLT file and get dimensions 
        #### (should always be 1014x1014 for WFC3)
        if root_filename:
            fi = pyfits.open(root+'_flt.fits',mode='update')
        else:
            fi = pyfits.open(root, mode='update')
            
        try:
            IMG = fi[1].data
            DQ = fi[3].data
        except:
            IMG = fi[0].data
            DQ = np.cast[int](IMG*0)
            
        NX, NY = IMG.shape
        
        #### array indices for each pixel position
        xi,yi = np.indices((NX,NY))
        xi = (xi-NX)/2./NX
        yi = (yi-NY)/2./NY
        
        #### If a segmentation image is available, read it
        #### for use as an object mask.
        segfile = glob.glob(root+'_flt.seg.fits*')
        if len(segfile) == 0:
            seg = IMG*0.
        else:
            print('Segmentation image: %s' %segfile[0])
            fis = pyfits.open(segfile[0])
            seg = fis[0].data
        
        seg_grow = nd.maximum_filter(seg, size=3)
        
        #### Apply segmentation mask, also mask out extreme IMG values 
        #### and any pixel with DQ flag > self.DQMAX
        
        
        q = np.where((seg_grow == 0) & (IMG > -1) & (IMG < 4) & (DQ < self.DQMAX)) 
        qb = np.where((seg_grow > 0) | (IMG < -1) | (IMG > 4) | (DQ >= self.DQMAX))
        IMGb = IMG*1.
        IMGb[qb] = np.nan
        
        #### Apply mask to IMG and fit matrices
        Aq = np.transpose(A[:,q[0],q[1]])
        IMGq = IMG[q[0],q[1]]
        
        #### Get fit parameters with least-sq. fit.
        p, resid, rank, s = scipy.linalg.lstsq(Aq,IMGq)

        print(p)
        
        #### Create the bg fit image from the fit parameters
        IMGout = IMG*0.
        for i in range(self.NPARAM):
            IMGout += A[i,:,:]*p[i]
        print('Done')
        
        if 'SKY0' in list(fi[0].header.keys()):
            BG_LEVEL = fi[0].header['SKY0']
        else:
            BG_LEVEL = 0
            
        BG_LEVEL += IMGout[NY/2, NY/2]*fi[0].header['EXPTIME']
        
        fi[0].header.update('SKY0', BG_LEVEL)
        
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
        
        if save_fit:
            hdu0 = pyfits.PrimaryHDU(header=fi[0].header)
            hdu1 = pyfits.ImageHDU(data=IMG, header=fi[1].header)
            hdu2 = pyfits.ImageHDU(data=IMGb, header=fi[1].header)
            hdu3 = pyfits.ImageHDU(data=IMGout, header=fi[1].header)
            save_im = pyfits.HDUList([hdu0,hdu1,hdu2,hdu3])
            #save_im[0].data = IMGout
            id=len(glob.glob(root+'_flt.BG*fits'))
            save_im.writeto(root+'_flt.BG_%02d.fits' %(id+1), clobber=True)
        
        #### Subtract fitted background, write
        #### bg-subtracted image to original FLT file if `overwrite` is True
        FIX = IMG-IMGout
        if overwrite:
            print('Overwrite: '+root)
            fi[1].data = FIX
            fi.flush()
        
        #### Save images to self.
        self.FIX = FIX
        self.IMG = IMG
        self.MODEL = IMGout
        fi.close()

def asn_grism_background_subtract(asn_file='ibhj42040_asn.fits', nbin=8, path='./', verbose=True, savefig=True):
    """
    Run the 1-D background subtraction routine for all FLT files
    defined in an ASN file.
    """
    
    asn = threedhst.utils.ASNFile(asn_file)
    if verbose:
        print('Background:')
    for flt in asn.exposures:
        if verbose:
            print('   %s' %flt)
        test = oned_grism_background_subtract(flt, nbin=nbin, verbose=verbose, savefig=savefig)
            
    
def oned_grism_background_subtract(flt_root, nbin=8, path='./', savefig=True, force=False, verbose=True):
    """
    Collapse a WFC FLT image along the y-axis to get a 
    model of the overall background.  The structure of 
    the sky background seems to be primarily along this axis.
    
    Note that this assumes that the grism sky background has already 
    been subtracted from the FLT file, and that an object mask 
    exists with filename `flt_root` + '_flt.seg.fits[.gz]'.
    
    Input is a rootname, like 'ibhm46i3q'.
    
    `nbin` is the number of pixels to combine along X to 
    smooth out the profile.
    """
    from threedhst.utils import find_fits_gz
    
    #### Find fits or fits.gz files
    flt_file = find_fits_gz(flt_root+'_flt.fits', hard_break=True)
    seg_file = find_fits_gz(flt_root+'_flt.seg.fits', hard_break=True)
    
    #### Open fits files
    flt = pyfits.open(flt_file,'update')
    seg = pyfits.open(seg_file)
    
    #### Don't proceed if not a grism exposure
    keys = list(flt[0].header.keys())
    IS_GRISM = False
    FILTER_STRING = ""
    for key in keys:
        if key.startswith('FILTER'):
            FILTER_STRING += flt[0].header[key]+" "
            if flt[0].header[key].startswith('G'):
                IS_GRISM = True
                
    if not IS_GRISM:
        if verbose:
            print('%s is not a grism exposure (%s)' %(flt_root, FILTER_STRING))
        return False
        
    #### Don't proceed if already been done
    if 'GRIS-BG' in list(flt[1].header.keys()):
        if (flt[1].header['GRIS-BG'] == 1) | (force is False):
            if verbose:
                print('Background already subtracted from %s.' %(flt_root))
            return False
    
    #### Arrays     
    xi = np.arange(1014/nbin)*nbin+nbin/2.
    yi = xi*1.
    si = xi*1.
     
    #### Set up output plot  
    if savefig:
        if USE_PLOT_GUI:
            fig = plt.figure(figsize=[5,3],dpi=100)
        else:
            fig = Figure(figsize=[5,3], dpi=100)
        
        fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.17,
                            bottom=0.17,right=0.97,top=0.97)
        ax = fig.add_subplot(111)
    
    #### Iterate on bg subtraction
    NITER = 4
    for it in range(NITER):
        for i in np.arange(1014/nbin):
            #### Masks, object and DQ
            seg_stripe = seg[0].data[:,i*nbin:(i+1)*nbin]
            dq_stripe = flt[3].data[:,i*nbin:(i+1)*nbin]
            OK_PIXELS = (seg_stripe == 0) & ((dq_stripe & (dq_stripe*0+4096)) == 0)
            #### Data columns
            data = flt[1].data[:,i*nbin:(i+1)*nbin]
            
            for it2 in range(1):
                #### Biweight mean and sigma
                stats = threedhst.utils.biweight(data[OK_PIXELS], both=True)
                yi[i], si[i] = stats
                OK_PIXELS = OK_PIXELS & (np.abs(data-0*stats[0]) < 3*stats[1])
                        
            # ypix, xpix = np.indices(data.shape)
            # xx = (ypix-507)/2.
            # poly = polyfit(xx[OK_PIXELS], data[OK_PIXELS], 4)
            # flt[1].data[:,i*nbin:(i+1)*nbin] -= polyval(poly, xx)
            # plt.plot(ypix[OK_PIXELS],data[OK_PIXELS],marker='.', linestyle='None', alpha=0.3)
            # plt.plot(ypix[OK_PIXELS],polyval(poly,xx[OK_PIXELS]),marker='.', alpha=0.8, linestyle='None')
            
        if savefig:
            ax.plot(xi, yi, color='red', alpha=1-it*0.8/(NITER-1))
        
        #### Interpolate smoothed back to individual pixels
        xpix = np.arange(1014)
        ypix = np.interp(xpix, xi, yi)
        for i in np.arange(1014):
            flt[1].data[:,i] -= ypix[i]
    
    #### Output figure
    if savefig:
        ax.plot(xi, yi*0, linestyle='--', alpha=0.6, color='black')
        ax.set_xlabel('x pixel')
        ax.set_ylabel('e-/s')
        ax.set_xlim(-1,1015)
        if USE_PLOT_GUI:
            plt.savefig(flt_root+'_flt.residual.png')
            plt.close()
        else:
            canvas = FigureCanvasAgg(fig)
            canvas.print_figure(flt_root+'_flt.residual.png', dpi=100, transparent=False)
        
    
    #### Add a 'GRIS-BG' header keyword to the FLT[DATA] extension.
    flt[1].header.update('GRIS-BG',1)
    
    #### Dump to FITS file
    flt.flush()
    
    return True

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
             redo_background=True, bg_only=False,
             ALIGN_IMAGE='../ACS/h_nz_sect*img.fits', ALIGN_EXT = 0,
             skip_drz=False,final_scale=0.06, pixfrac=0.8,
             IMAGES=['/research/HST/GRISM/CONF/G141_sky_cleaned.fits'],
             align_geometry='shift', 
             initial_order=-1,
             clean=True,
             save_fit=False):
    """
prep_all(asn_files='ib*050_asn.fits', get_shift=True, 
         redo_background=True, bg_only=False)
    
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
                    bg_only=bg_only, ALIGN_IMAGE=ALIGN_IMAGE,
                    ALIGN_EXT=ALIGN_EXT,
                    skip_drz=skip_drz,final_scale=final_scale,
                    pixfrac=pixfrac,
                    IMAGES=IMAGES,
                    initial_order=initial_order,
                    align_geometry=align_geometry, clean=clean,
                    save_fit=save_fit)

def make_3dhst_persistence_mask(asn_direct_file='GOODS-N-45-F140W_asn.fits'):
    """
    The 3D-HST visit sequence is D-G / D-G / D-G / D-G
    """
    from threedhst.prep_flt_files import MultidrizzleRun
    
    STARS_FLAGGED = False
        
    asn = threedhst.utils.ASNFile(asn_direct_file)
    
    flt = pyfits.open(asn.exposures[0]+'_flt.fits')
    instrume = flt[0].header.get('INSTRUME').strip()
    if instrume is 'ACS':
        skip=2
    else:
        skip=1
    
    #### Blot the mosaic back to the FLT frame
    run = MultidrizzleRun((asn_direct_file.split('_asn.fits')[0]).upper())
    
    for i,exp in enumerate(asn.exposures):
        run.blot_back(ii=i*skip, copy_new=(i is 0))
    
    ### Run SExtractor on a direct image, which tells where the persistence will be
    ### in the *next* direct image
    
    ### only very bright stars in the 2nd exposure will persist on the 3rd 
    maglim = [18,18,15,18] 
    
    for iflt in range(1,4):
        root = asn.exposures[iflt]+'_flt'
    
        se = threedhst.sex.SExtractor()
        se.aXeParams()
        ## XXX add test for user-defined .conv file
        se.copyConvFile(grism=False)
        se.overwrite = True
        se.options['CATALOG_NAME']    = root+'.BLOT.SCI.cat'
        se.options['CHECKIMAGE_NAME'] = root+'.seg.fits, bg.fits'
        se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION, BACKGROUND'
        se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'        
        se.options['WEIGHT_IMAGE']    = root+'.BLOT.WHT.fits'
        se.options['FILTER']    = 'Y'
        se.options['BACK_TYPE']     = 'AUTO'
        se.options['BACK_FILTERSIZE']     = '2'
        se.options['FILTER_NAME'] = 'default.conv'
        sigma=10
        se.options['DETECT_THRESH']    = '%f' %sigma
        se.options['ANALYSIS_THRESH']  = '%f' %sigma
        se.options['MAG_ZEROPOINT'] = '26.46'
        status = se.sextractImage(root+'.BLOT.SCI.fits')
        
        cat = threedhst.sex.mySexCat(root+'.BLOT.SCI.cat')
        xpix, ypix, mag, radius = np.cast[float](cat.X_IMAGE), np.cast[float](cat.Y_IMAGE), np.cast[float](cat.MAG_AUTO), np.cast[float](cat.FLUX_RADIUS)
        
        try:
            os.remove(root+'.BLOT.SCI.fits')
            os.remove(root+'.BLOT.WHT.fits')
        except:
            pass
            
        stars = (mag < maglim[iflt]) & (radius < 2)
        idx = np.arange(len(stars))[stars]

        ### first-order polygon obtained from a reference image
        xmask = np.array([1,0,1,0,1,0,1,0])
        ymask = np.array([0,1,0,1,0,1,0,1])

        first_mask = np.cast[float]('583.52156,585.79943,722.12037,585.60725,721.63812,580.78472,585.16049,580.30247'.split(','))
        first_mask -= 542.81*xmask + 579.843*ymask + 1*ymask
                
        zeroth_mask = np.array([-1,-1,1,-1,1,1,-1,1])*3
        
        if len(idx) > 0:
            STARS_FLAGGED = True
            
            regfile = asn.exposures[iflt]+'_flt.fits.mask.reg'
            if os.path.exists(regfile):
                fp = open(regfile,'a')
            else:
                fp = open(regfile,'w')
                fp.write('image\n')
            
            for i in idx:
                #fp.write('circle(%f,%f,5)\n' %(xpix[i], ypix[i]))
                
                #### Get parameters from aXe calibration of G141
                dldp_0 = field_dependent(xpix[i], ypix[i],'8.95431E+03   9.35925E-02   0.0')
                dldp_1 = field_dependent(xpix[i], ypix[i], '4.51423E+01   3.17239E-04   2.17055E-03  -7.42504E-07   3.48639E-07   3.09213E-07')
                
                dydx_0 = field_dependent(xpix[i], ypix[i],'1.96882E+00  9.09159E-05 -1.93260E-03')
                dydx_1 = field_dependent(xpix[i], ypix[i],'1.04275E-02 -7.96978E-06 -2.49607E-06  1.45963E-09  1.39757E-08  4.84940E-10')
                
                y0 = dydx_0 + dydx_1*(xpix[i])-5
                dy = 3
                xleft = (1.05e4 - dldp_0)/dldp_1
                xright = (1.75e4 - dldp_0)/dldp_1
                first_mask = np.array([xleft,y0-dy,xright,y0-dy,xright,y0+dy,xleft,y0+dy])
                
                pstr = 'polygon(%f' %(first_mask[0]+xmask[0]*xpix[i]+ymask[0]*ypix[i])
                for j in range(1,8):
                    pstr += ',%f' %(first_mask[j]+xmask[j]*xpix[i]+ymask[j]*ypix[i])
                
                fp.write(pstr+')\n')
                
                #### zeroth order position is field dependent, parameters from
                #### the aXe calibration for G141
                dldp_0 = 459047.749023
                dldp_1 = 2324.048828
                xoff_0 = field_dependent(xpix[i], ypix[i],'-0.2400520   -0.0023144    0.0111089')-1
                x0th = (1.4e4 - (dldp_0-dldp_1*xoff_0))/dldp_1
                
                pstr = 'polygon(%f' %(zeroth_mask[0]+xmask[0]*x0th+xmask[0]*xpix[i]+ymask[0]*ypix[i])
                for j in range(1,8):
                    pstr += ',%f' %(zeroth_mask[j]+xmask[j]*x0th+xmask[j]*xpix[i]+ymask[j]*ypix[i])
                
                fp.write(pstr+')\n')
                
            fp.close()
        
    return STARS_FLAGGED

def field_dependent(xi, yi, str_coeffs):
    """ 
    Calculate field-dependent parameter for aXe conventions.
    
    a = a_0 + a_1 * xi + a_2 * yi + a_3 * xi**2 + a_4 * xi * yi + a_5 * yi**2
    
    """
    coeffs = np.cast[float](str_coeffs.split())
    xy = np.array([1,xi,yi,xi**2,xi*yi,yi**2, xi**3, xi**2*yi, xi*yi**2, yi**3])
    a = np.sum(coeffs*xy[0:len(coeffs)])
    return a

def prep_flt(asn_file=None, get_shift=True, bg_only=False, bg_skip=False,
                first_run=True, redo_background=True,
                ALIGN_IMAGE='../ACS/h_nz_sect*img.fits', ALIGN_EXT = 0, 
                skip_drz=False, final_scale=0.06, pixfrac=0.8,
                IMAGES=['/research/HST/GRISM/CONF/G141_sky_cleaned.fits'],
                align_geometry='shift', clean=True,
                initial_order=-1, save_fit=False,
                TWEAKSHIFTS_ONLY=False,
                oned_background=True, make_persistence_mask=False,
                redo_segmentation=True, 
                clean_drz=False):
    """
prep_flt(asn_file=None, get_shift=True, bg_only=False,
            redo_background=True)

    
    Subtract background and align WCS of direct/grism FLT files.
    
    1) Apply the DQ masks defined in the *mask.reg files, as created
       by threedhst.dq
    
    2) First pass on background subtraction 
        
        o 0th order is median background image, e.g. G141 sky 
        
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
       
    8) Subtract the collapsed background to fix residuals from the grism sky fit
    
    9) Run multidrizzle again [if redo_background is True]
    
    
    """
    #import fit_2d_poly
    import threedhst
    import os
    import glob
    
    #import threedhst.dq    
    
    if asn_file is None:
        asn_file = 'ib3728050_asn.fits'
    
    if bg_skip:
        bg_only=False
        redo_background=False
        
    ROOT_DIRECT = asn_file.split('_asn.fits')[0]
    # ALIGN_IMAGE = '../ACS/h_nz_sect*img.fits'
    
    asn = threedhst.utils.ASNFile(asn_file)
    flt = pyfits.open(asn.exposures[0]+'_flt.fits')

    if flt[0].header['DETECTOR'] == 'UVIS':
        bg_skip = True
        redo_background = False
        
    #### First guess at shifts
    if get_shift:

        threedhst.shifts.run_tweakshifts(asn_file, verbose=True)
        threedhst.shifts.checkShiftfile(asn_file)
        threedhst.shifts.default_rotation(asn_file, path_to_flt='./')
            
        if (not TWEAKSHIFTS_ONLY):
            for ig, geom in enumerate(align_geometry.split(',')):
                first = ig == 0
                startMultidrizzle(asn_file, use_shiftfile=True, skysub=True,
                    final_scale=final_scale, pixfrac=pixfrac, driz_cr=first,
                    updatewcs=first, clean=True, median=first)
                
                threedhst.shifts.refine_shifts(ROOT_DIRECT=ROOT_DIRECT,
                          ALIGN_IMAGE=ALIGN_IMAGE, 
                          ALIGN_EXTENSION = ALIGN_EXT,
                          fitgeometry=geom.strip(), clean=clean)
                #
                #### Put persistence mask here, still not quite working
                #if make_persistence_mask:
                #    make_3dhst_persistence_mask(asn_file)
                
            
        ### Need fresh FLT files now
        threedhst.process_grism.fresh_flt_files(asn_file)
             
    #### First pass background subtraction
    if not bg_skip:
        #### Set up matrix for fitting
        fit = fit_2D_background(ORDER=initial_order,
                                IMAGES=IMAGES)#, x0=507, y0=507)
        
        for exp in asn.exposures:
            #threedhst.regions.apply_dq_mask(exp+'_flt.fits')
            fit.fit_image(exp, A=fit.A, show=False, overwrite=True,
                          save_fit=save_fit)
    
    #### Stop here if only want background subtraction
    if bg_only:
        return
            
    if not skip_drz:
        if first_run:
            #### First MDRZ run needs native pixels to avoid flagging stars as CRs
            startMultidrizzle(asn_file, use_shiftfile=True, 
                skysub=True, final_scale=0.128254, pixfrac=1.0,
                driz_cr=True, median=True, updatewcs=True)
        
        if (final_scale != 0.128254) | (pixfrac != 1.0) | (not first_run):
            startMultidrizzle(asn_file, use_shiftfile=True, 
                skysub=bg_skip, final_scale=final_scale, pixfrac=pixfrac,
                driz_cr=False, median=False, updatewcs=False, skyuser='SKY0')
     
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
    
    if flt[0].header['DETECTOR'] == 'UVIS':
        skip=2
    
    if redo_background:
        copy_new = True
        #### Clean DRZ image 
        #asn_file='J000208-152314-F160W_asn.fits'
        if clean_drz:
            print('Clean_DRZ: clean bad pixels from DRZ mosaic')
            try:
                import scipy.ndimage as nd
                drz = pyfits.open(asn_file.replace('asn','drz'), mode='update')
                med = nd.median_filter(drz[1].data+0.1, size=4)
                ratio = (drz[1].data+0.1)/med
                bad = ratio > 5
                drz[1].data[bad] = 0.
                drz[2].data[bad] = 0
                drz.flush()
            except:
                print('Clean_DRZ failed.')
                pass
        
        for i,exp in enumerate(asn.exposures):
            if (not os.path.exists(run.flt[i]+'.seg.fits')) | redo_segmentation:
                run.blot_back(ii=i*skip, copy_new=copy_new)
                make_segmap(run.flt[i])
                copy_new=False
        
        ### Flag bright stars in segmentation map
        asn_mask = asn_file+'.mask.reg'
        if os.path.exists(asn_mask):
            threedhst.showMessage("Apply ASN mask: %s" %(asn_mask))
            for flt in asn.exposures:
                seg = flt+'_flt.seg.fits'
                threedhst.regions.apply_dq_mask(seg, extension=0, 
                    mask_file = asn_mask)
                            
    
    #### Run BG subtraction with improved mask and run multidrizzle again
    if redo_background:
        fit = fit_2D_background(ORDER=initial_order, IMAGES=IMAGES)
        for exp in asn.exposures:
            #### 2-D background, fit
            fit.fit_image(exp, A=fit.A, show=False, overwrite=True, 
                          save_fit=save_fit)
            #
            #### 1-D background, measured
            if oned_background:
                print('\n Extract 1D background \n')
                test = oned_grism_background_subtract(exp, nbin=26, savefig=True, verbose=False)
        
        startMultidrizzle(asn_file, use_shiftfile=True, skysub=True,
            final_scale=final_scale, pixfrac=pixfrac, driz_cr=False,
            updatewcs=False, median=False, clean=clean, skyuser='SKY0')
    
    #### Take out huge bad pixels from the FLT images
    if flt[0].header['DETECTOR'] == 'IR':
        for exp in asn.exposures:
            flt = pyfits.open('%s_flt.fits' %(exp), mode='update')
            bad = (flt['SCI'].data > 1.e5) | ((flt['DQ'].data & 4096) > 0)
            if bad.sum() > 0:
                flt['SCI'].data[bad] = 0
                flt['ERR'].data[bad] = 1.e6
                flt.flush()
            
            flt.close()
            
    if clean:
        files=glob.glob('*BLOT*')
        files.extend(glob.glob('drz_???.fits'))
        for file in files: 
            #print 'rm '+file
            os.remove(file)
    
def flag_data_quality():
    """ 
    Flag asteroid trails.
    """
    import threedhst.dq
    
    ####********************************************####
    ####                 AEGIS
    ####********************************************####
    os.chdir('/research/HST/GRISM/3DHST/AEGIS/PREP_FLT')
    
    threedhst.dq.checkDQ('ibhj42030_asn.fits','ibhj42040_asn.fits', 
                         path_to_flt='./')
    
    
    ####********************************************####
    ####                 COSMOS
    ####********************************************####
    os.chdir('/research/HST/GRISM/3DHST/COSMOS/PREP_FLT')
    
    threedhst.dq.checkDQ('ibhm31030_asn.fits','ibhm31040_asn.fits')
    threedhst.dq.checkDQ('ibhm44030_asn.fits','ibhm44040_asn.fits')
    threedhst.dq.checkDQ('ibhm51030_asn.fits','ibhm51040_asn.fits')
    threedhst.dq.checkDQ('ibhm53030_asn.fits','ibhm53040_asn.fits')
    
    ####********************************************####
    ####                 SN-MARSHALL
    ####********************************************####
    os.chdir('/research/HST/GRISM/3DHST/SN-MARSHALL/PREP_FLT')
    
    threedhst.dq.checkDQ('ibfuw1070_asn.fits','ibfuw1070_asn.fits', 
                         path_to_flt='./')
    
    
    ####********************************************####
    ####                 GOODS-N
    ####********************************************####
    
    threedhst.dq.checkDQ('ib3749050_asn.fits','ib3749050_asn.fits',                                               path_to_flt='./')
    threedhst.dq.checkDQ('ib3703050_asn.fits','ib3703060_asn.fits',                                               path_to_flt='./')
    
    threedhst.dq.checkDQ('GOODS-N-23-G141_asn.fits','GOODS-N-23-G141_asn.fits',                                               path_to_flt='./')
    
def process_all():
    """
    Initial processing of all 3D-HST frames
    """
    from pyraf import iraf
    from iraf import stsdas,dither,slitless,axe 

    INDEF = iraf.INDEF
    no = iraf.no
    yes = iraf.yes
    
    from threedhst.prep_flt_files import process_3dhst_pair as pair
    import threedhst.prep_flt_files
    
    ####********************************************####
    ####                 GOODS-N
    ####********************************************####
    os.chdir('/research/HST/GRISM/3DHST/GOODS-N/PREP_FLT')
    ALIGN = '../ACS/h_nz*drz_img.fits'
    
    direct = glob.glob('ib37[012]*050_asn.fits')
    grism = glob.glob('ib37[012]*060_asn.fits')
    
    for i  in range(28):
        pair(direct[i], grism[i], ALIGN_IMAGE = ALIGN, SKIP_GRISM=False,
             GET_SHIFT=True, DIRECT_HIGHER_ORDER=2)
        
    # Make mosaic
    direct_files = glob.glob('GOODS-N*[0-9]*-F140W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='GOODS-N-F140W',
                       path_to_FLT='./', run_multidrizzle=False)
    
    SCALE = 0.06
    NX = np.round(0.25/SCALE*2326)
    NY = np.round(0.25/SCALE*3919)
    threedhst.prep_flt_files.startMultidrizzle('GOODS-N-F140W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=SCALE, pixfrac=0.8, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=189.22805, dec=62.236116,
             final_outnx=NX, final_outny=NY, final_rot=45)
    
    #### Make direct image for each pointing that also include 
    #### neighboring pointings
    files = glob.glob('GOODS-N-[0-9]*-F140W_asn.fits')
    for file in files:
        pointing = file.split('_asn.fits')[0]
        threedhst.prep_flt_files.mosaic_to_pointing(mosaic_list='GOODS-N*[0-9]*F140W',
                                    pointing=pointing,
                                    run_multidrizzle=True)
                                                                    
    ####********************************************####
    ####                      AEGIS
    ####********************************************####
    os.chdir('/research/HST/GRISM/3DHST/AEGIS/PREP_FLT')
    
    ALIGN = '../NMBS/AEGIS-N2_K_sci.fits'
    pair('ibhj42030_asn.fits','ibhj42040_asn.fits', ALIGN_IMAGE = ALIGN)
    pair('ibhj43030_asn.fits','ibhj43040_asn.fits', ALIGN_IMAGE = ALIGN)
    
    ALIGN = '/research/HST/CANDELS/EGS/WIRDS/WIRDS_Ks_141927+524056_T0002.SUB2.fits'
    ALIGN = '../WIRDS/WIRDS_Ks_141927+524056_T0002.fits'

    pair('ibhj49030_asn.fits','ibhj49040_asn.fits', ALIGN_IMAGE = ALIGN)
    pair('ibhj40030_asn.fits','ibhj40040_asn.fits', ALIGN_IMAGE = ALIGN)
    pair('ibhj39030_asn.fits','ibhj39040_asn.fits', ALIGN_IMAGE = ALIGN)
    
    ### Make mosaic
    direct_files = glob.glob('AEGIS*-F140W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='AEGIS-F140W',
                       path_to_FLT='./', run_multidrizzle=False)
    threedhst.prep_flt_files.startMultidrizzle('AEGIS-F140W_asn.fits',
             use_shiftfile=True, skysub=True,
             final_scale=0.06, pixfrac=0.8, driz_cr=False,
             updatewcs=False, clean=True, median=False)
    #
    files = glob.glob('AEGIS-*-F140W_asn.fits')
    for file in files:
        pointing = file.split('_asn.fits')[0]
        threedhst.prep_flt_files.mosaic_to_pointing(mosaic_list='AEGIS-*-F140W',
                                    pointing=pointing,
                                    run_multidrizzle=True, grow=400)
    
    # threedhst.prep_flt_files.make_grism_subsets(root='AEGIS')
    
    ####********************************************####
    ####                     COSMOS
    ####********************************************####
    os.chdir('/research/HST/GRISM/3DHST/COSMOS/PREP_FLT')
    ALIGN = '../NMBS/COSMOS-1.V4.K_nosky.fits'
            
    SKIP = False
    pair('ibhm51030_asn.fits','ibhm51040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    pair('ibhm31030_asn.fits','ibhm31040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    pair('ibhm46030_asn.fits','ibhm46040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    pair('ibhm29030_asn.fits','ibhm29040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)

    pair('ibhm30030_asn.fits','ibhm30040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    pair('ibhm32030_asn.fits','ibhm32040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    pair('ibhm33030_asn.fits','ibhm33040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)

    pair('ibhm34030_asn.fits','ibhm34040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    pair('ibhm35030_asn.fits','ibhm35040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    pair('ibhm36030_asn.fits','ibhm36040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    pair('ibhm37030_asn.fits','ibhm37040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=False, SKIP_GRISM=True, GET_SHIFT=True)
    
    pair('ibhm43030_asn.fits','ibhm43040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    pair('ibhm44030_asn.fits','ibhm44040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    pair('ibhm45030_asn.fits','ibhm45040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    pair('ibhm48030_asn.fits','ibhm48040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    pair('ibhm52030_asn.fits','ibhm52040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    pair('ibhm53030_asn.fits','ibhm53040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP, GRISM_HIGHER_ORDER=-1)
    pair('ibhm54030_asn.fits','ibhm54040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    pair('ibhm55030_asn.fits','ibhm55040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    pair('ibhm56030_asn.fits','ibhm56040_asn.fits', ALIGN_IMAGE = ALIGN, IMAGES=['G141_fixed_sky.fits'], SKIP_DIRECT=SKIP)
    
    threedhst.prep_flt_files.startMultidrizzle('COSMOS-23-F140W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=0.06, pixfrac=0.8, driz_cr=False,
             updatewcs=False, clean=True, median=False)
    
    # LAEs
    threedhst.prep_flt_files.startMultidrizzle('COSMOS-1-F140W_asn.fits',
            use_shiftfile=True, skysub=False,
            final_scale=0.06, pixfrac=0.8, driz_cr=False,
            ra=150.06581, dec=2.4083496, final_outnx=420, final_outny=185,
            updatewcs=False, clean=True, median=False)

    threedhst.prep_flt_files.startMultidrizzle('COSMOS-23-G141_asn.fits',
            use_shiftfile=True, skysub=False,
            final_scale=0.06, pixfrac=0.8, driz_cr=False,
            ra=150.13029, dec=2.3016245, final_outnx=530, final_outny=300,
            updatewcs=False, clean=True, median=False)
    
    ### Make mosaic
    direct_files = glob.glob('COSMOS-*-F140W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='COSMOS-F140W',
                       path_to_FLT='./', run_multidrizzle=False)
    threedhst.prep_flt_files.startMultidrizzle('COSMOS-F140W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=0.06, pixfrac=0.8, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=150.12356, dec=2.3608425,
             final_outnx = 9355, final_outny=11501)
             
    threedhst.prep_flt_files.make_grism_subsets(root='COSMOS')
    
    #### Make direct image for each pointing that also include 
    #### neighboring pointings
    files = glob.glob('COSMOS-[0-9]*-F140W_asn.fits')
    for file in files:
        pointing = file.split('_asn.fits')[0]
        threedhst.prep_flt_files.mosaic_to_pointing(mosaic_list='COSMOS-*[0-9]*F140W',
                                    pointing=pointing,
                                    run_multidrizzle=True, grow=400)
    
    ####********************************************####
    ####                   GOODS-S
    ####********************************************####
    os.chdir('/research/HST/GRISM/3DHST/GOODS-S/PREP_FLT')
    ALIGN = '../ECDFS/MUSYC_ECDFS_BVR.fits'
    
    pair('ibhj06030_asn.fits','ibhj06040_asn.fits', ALIGN_IMAGE = ALIGN)
    pair('ibhj27030_asn.fits','ibhj27040_asn.fits', ALIGN_IMAGE = ALIGN)
    pair('ibhj28030_asn.fits','ibhj28040_asn.fits', ALIGN_IMAGE = ALIGN)

    ALIGN = '../ACS/h_sz*drz_img.fits'
    pair('ibhj24030_asn.fits','ibhj24040_asn.fits', ALIGN_IMAGE = ALIGN)
    
    ### Make mosaic
    direct_files = glob.glob('GOODS-S*-F140W_asn.fits')
    threedhst.utils.combine_asn_shifts(direct_files, out_root='GOODS-S-F140W',
                       path_to_FLT='./', run_multidrizzle=False)
    threedhst.prep_flt_files.startMultidrizzle('GOODS-S-F140W_asn.fits',
             use_shiftfile=True, skysub=False,
             final_scale=0.06, pixfrac=0.8, driz_cr=False,
             updatewcs=False, clean=True, median=False,
             ra=53.24009, dec=-27.84430, 
             final_outnx=7440, final_outny=4650)
    
    threedhst.prep_flt_files.make_grism_subsets(root='GOODS-S')
    
    ####********************************************####
    ####              SN-MARSHALL (UDS)
    ####********************************************####
    ## Shifts and direct images determined separately
    os.chdir('/research/HST/GRISM/3DHST/SN-MARSHALL/PREP_FLT')
    #### First orientation
    for i in range(1,2,2):
        ist = '%0d' %(i)
        shutil.copy('ibfuw'+ist+'070_shifts.txt',
                    'MARSHALL-'+ist+'-G141_shifts.txt')
        pair(None,'ibfuw'+ist+'070_asn.fits', ALIGN_IMAGE = None, SKIP_DIRECT=True, SKIP_GRISM=False, GRISM_HIGHER_ORDER=2)
    
    asn_list = glob.glob('MARSHALL-[135]-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(asn_list, out_root='MARSHALL-225-G141',
                    path_to_FLT='./', run_multidrizzle=False)
    threedhst.prep_flt_files.startMultidrizzle('MARSHALL-225-G141_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    
    #### Second orientation
    for i in range(2,7,2):
        ist = '%0d' %(i)
        shutil.copy('ibfuw'+ist+'070_shifts.txt',
                    'MARSHALL-'+ist+'-G141_shifts.txt')
        pair(None,'ibfuw'+ist+'070_asn.fits', ALIGN_IMAGE = None, SKIP_DIRECT=True, SKIP_GRISM=False, GRISM_HIGHER_ORDER=2)
    
    asn_list = glob.glob('MARSHALL-[246]-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(asn_list, out_root='MARSHALL-245-G141',
                    path_to_FLT='./', run_multidrizzle=False)
    threedhst.prep_flt_files.startMultidrizzle('MARSHALL-245-G141_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
        
    
    ####********************************************####
    ####         SN-PRIMO (GOODS-S, UDF)
    ####********************************************####
    ## Shifts determined separately
    os.chdir('/research/HST/GRISM/3DHST/SN-PRIMO/PREP_FLT')
    shutil.copy('G141_1026_shifts.txt','PRIMO-1026-G141_shifts.txt')
    pair(None,'G141_1026_asn.fits', ALIGN_IMAGE = None, SKIP_DIRECT=True, SKIP_GRISM=False, IMAGES=['G141_fixed_sky.fits'], GRISM_HIGHER_ORDER=2)
    
    shutil.copy('G141_1101_shifts.txt','PRIMO-1101-G141_shifts.txt')
    pair(None,'G141_1101_asn.fits', ALIGN_IMAGE = None, SKIP_DIRECT=True, SKIP_GRISM=False, IMAGES=['G141_fixed_sky.fits'], GRISM_HIGHER_ORDER=2)
    
    
    #### Direct image
    CANDELS='/research/HST/CANDELS/GOODS-S/PREP_FLT/'
    f125w_list = threedhst.shifts.find_align_images_that_overlap('PRIMO-1101-G141_drz.fits', CANDELS+'*F125W*pointing.reg', is_region=True)
    threedhst.utils.combine_asn_shifts(f125w_list, out_root='PRIMO_F125W',
                    path_to_FLT=CANDELS, run_multidrizzle=False)
    threedhst.process_grism.fresh_flt_files('PRIMO_F125W_asn.fits', from_path=CANDELS)
    threedhst.prep_flt_files.startMultidrizzle('PRIMO_F125W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False,
                 refimage='PRIMO-1101-G141_drz.fits[1]')
    # Shifts aren't quite right, due to differen scales
    threedhst.shifts.refine_shifts(ROOT_DIRECT='PRIMO_F125W',
                      ALIGN_IMAGE='F125W_drz.fits',
                      fitgeometry='shift', clean=True,
                      ALIGN_EXTENSION=1)
    #
    threedhst.prep_flt_files.startMultidrizzle('PRIMO_F125W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.7, driz_cr=False,
                 updatewcs=False, clean=True, median=False,
                 ra=53.159487, dec=-27.776371,
                 final_outnx=2980, final_outny=2380)
    
    threedhst.process_grism.clean_flt_files('PRIMO_F125W_asn.fits')
    
    ##### F160W
    f160w_list = threedhst.shifts.find_align_images_that_overlap('PRIMO-1101-G141_drz.fits', CANDELS+'*F160W*pointing.reg', is_region=True)
    threedhst.utils.combine_asn_shifts(f160w_list, out_root='PRIMO_F160W',
                    path_to_FLT=CANDELS, run_multidrizzle=False)
    threedhst.process_grism.fresh_flt_files('PRIMO_F160W_asn.fits', from_path=CANDELS)
    threedhst.prep_flt_files.startMultidrizzle('PRIMO_F160W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False,
                 refimage='PRIMO-1101-G141_drz.fits[1]')
    # Shifts aren't quite right, due to differen scales
    threedhst.shifts.refine_shifts(ROOT_DIRECT='PRIMO_F160W',
                      ALIGN_IMAGE='F125W_drz.fits',
                      fitgeometry='shift', clean=True,
                      ALIGN_EXTENSION=1)
    #
    threedhst.prep_flt_files.startMultidrizzle('PRIMO_F160W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.7, driz_cr=False,
                 updatewcs=False, clean=True, median=False,
                 ra=53.159487, dec=-27.776371,
                 final_outnx=2980, final_outny=2380)
                 
    threedhst.process_grism.clean_flt_files('PRIMO_F160W_asn.fits')
    
    ##### ACS F850LP
    acs_list = threedhst.shifts.find_align_images_that_overlap('PRIMO_F160W_drz.fits', '../ACS/h_sz*drz_img.fits')

    threedhst.shifts.matchImagePixels(input=acs_list,
                     matchImage='PRIMO_F160W_drz.fits',
                     output='PRIMO_F850LP.fits', match_extension = 1,
                     input_extension=0)
    
    files = glob.glob('[rgb].fits')
    for file in files: os.remove(file)
    s = 10**(-0.4*(26.25-25.96))
    iraf.imcalc('PRIMO_F125W_drz.fits[1]','g.fits','im1*%f' %s)
    s = 10**(-0.4*(25.84-25.96))*3.5
    iraf.imcalc('PRIMO_F850LP.fits[0]','b.fits','im1*%f' %s)
    iraf.imcalc('PRIMO_F160W_drz.fits[1]','r.fits','im1')
    
    ####********************************************####
    ####            SN-GEORGE (GOODS-S)
    ####********************************************####
    ## Shifts and direct images determined separately
    os.chdir('/research/HST/GRISM/3DHST/SN-GEORGE/PREP_FLT')
    shutil.copy('ibfug1040_shifts.txt','GEORGE-1-G141_shifts.txt')
    pair(None,'ibfug1040_asn.fits', ALIGN_IMAGE = None, SKIP_DIRECT=True, SKIP_GRISM=False,  GRISM_HIGHER_ORDER=2)

    shutil.copy('ibfug2040_shifts.txt','GEORGE-2-G141_shifts.txt')
    pair(None,'ibfug2040_asn.fits', ALIGN_IMAGE = None, SKIP_DIRECT=True, SKIP_GRISM=False,  GRISM_HIGHER_ORDER=2)
    
    asn_files = glob.glob('GEORGE-?-G141_asn.fits')
    threedhst.utils.combine_asn_shifts(asn_files, out_root='GEORGE-G141',
                    path_to_FLT='./', run_multidrizzle=False)
    
    threedhst.prep_flt_files.startMultidrizzle('GEORGE-G141_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False)
    
    ###### Direct image
    # Shifts aren't quite right, due to differen scales, use direct as reference
    threedhst.utils.combine_asn_shifts(['ibfug1030_asn.fits'],
                    out_root='GEORGE_DIRECT',
                    path_to_FLT='../RAW/', run_multidrizzle=False)
    #
    threedhst.process_grism.fresh_flt_files('GEORGE_DIRECT_asn.fits', from_path='../RAW/')
    threedhst.prep_flt_files.startMultidrizzle('GEORGE_DIRECT_asn.fits',
                 use_shiftfile=True, skysub=True,
                 final_scale=0.06, pixfrac=0.8, driz_cr=True,
                 updatewcs=True, clean=True, median=True,
                 ra=53.084992, dec=-27.835811,
                 final_outnx=3450, final_outny=3100)
                 
    threedhst.process_grism.clean_flt_files('GEORGE_DIRECT_asn.fits')
    
    
    CANDELS='/research/HST/CANDELS/GOODS-S/PREP_FLT/'
    ##### Full F125W
    
    f125w_list = threedhst.shifts.find_align_images_that_overlap('GEORGE-G141_drz.fits', CANDELS+'*F125W*pointing.reg', is_region=True)
    threedhst.utils.combine_asn_shifts(f125w_list, out_root='GEORGE_F125W',
                    path_to_FLT=CANDELS, run_multidrizzle=False)
    threedhst.process_grism.fresh_flt_files('GEORGE_F125W_asn.fits', from_path=CANDELS)
    threedhst.prep_flt_files.startMultidrizzle('GEORGE_F125W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False,
                 ra=53.084992, dec=-27.835811,
                 final_outnx=3450, final_outny=3100)
                      
    threedhst.shifts.refine_shifts(ROOT_DIRECT='GEORGE_F125W',
                      ALIGN_IMAGE='GEORGE_DIRECT_drz.fits',
                      fitgeometry='shift', clean=False,
                      ALIGN_EXTENSION=1)
    #
    threedhst.prep_flt_files.startMultidrizzle('GEORGE_F125W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.7, driz_cr=False,
                 updatewcs=False, clean=True, median=False,
                 ra=53.084992, dec=-27.835811,
                 final_outnx=3450, final_outny=3100)
                     
    threedhst.process_grism.clean_flt_files('GEORGE_F125W_asn.fits')
    
    ##### Full F160W
    
    f160w_list = threedhst.shifts.find_align_images_that_overlap('GEORGE-G141_drz.fits', CANDELS+'*F160W*pointing.reg', is_region=True)
    threedhst.utils.combine_asn_shifts(f160w_list, out_root='GEORGE_F160W',
                    path_to_FLT=CANDELS, run_multidrizzle=False)
    threedhst.process_grism.fresh_flt_files('GEORGE_F160W_asn.fits', from_path=CANDELS)
    threedhst.prep_flt_files.startMultidrizzle('GEORGE_F160W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False,
                 ra=53.084992, dec=-27.835811,
                 final_outnx=3450, final_outny=3100)
                      
    threedhst.shifts.refine_shifts(ROOT_DIRECT='GEORGE_F160W',
                      ALIGN_IMAGE='GEORGE_DIRECT_drz.fits',
                      fitgeometry='shift', clean=False,
                      ALIGN_EXTENSION=1)
    #
    threedhst.prep_flt_files.startMultidrizzle('GEORGE_F160W_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.7, driz_cr=False,
                 updatewcs=False, clean=True, median=False,
                 ra=53.084992, dec=-27.835811,
                 final_outnx=3450, final_outny=3100)
                     
    threedhst.process_grism.clean_flt_files('GEORGE_F160W_asn.fits')
    
    ##### ACS F850LP
    acs_list = threedhst.shifts.find_align_images_that_overlap('GEORGE_DIRECT_drz.fits', '../ACS/h_sz*drz_img.fits')

    threedhst.shifts.matchImagePixels(input=acs_list,
                     matchImage='GEORGE_DIRECT_drz.fits',
                     output='GEORGE_F850LP.fits', match_extension = 1,
                     input_extension=0)
    
    files = glob.glob('[rgb].fits')
    for file in files: os.remove(file)
    s = 10**(-0.4*(26.25-25.96))
    iraf.imcalc('GEORGE_F125W_drz.fits[1]','g.fits','im1*%f' %s)
    s = 10**(-0.4*(25.84-25.96))*3.5
    iraf.imcalc('GEORGE_F850LP.fits[0]','b.fits','im1*%f' %s)
    iraf.imcalc('GEORGE_F160W_drz.fits[1]','r.fits','im1')
    
    ####********************************************####
    ####              Stanford clusters in NDWFS
    ####********************************************####
    os.chdir('/research/HST/GRISM/3DHST/STANFORD/PREP_FLT')
    ALIGN = '../NMBS/AEGIS-N2_K_sci.fits'
    
    pair('ib5l06020_asn.fits','ib5l06030_asn.fits', ALIGN_IMAGE = None, GET_SHIFT=True, TWEAKSHIFTS_ONLY=True)
    pair('ib5l09020_asn.fits','ib5l09030_asn.fits', ALIGN_IMAGE = None, GET_SHIFT=True, TWEAKSHIFTS_ONLY=True)
    pair('ib5l13020_asn.fits','ib5l13030_asn.fits', ALIGN_IMAGE = None, GET_SHIFT=True, TWEAKSHIFTS_ONLY=True)
    pair('ib5l17020_asn.fits','ib5l17030_asn.fits', ALIGN_IMAGE = None, GET_SHIFT=True, TWEAKSHIFTS_ONLY=True)
    pair('ib5l01020_asn.fits','ib5l01030_asn.fits', ALIGN_IMAGE = None, GET_SHIFT=True, TWEAKSHIFTS_ONLY=True)
    pair('ib5l04020_asn.fits','ib5l04030_asn.fits', ALIGN_IMAGE = None, GET_SHIFT=True, TWEAKSHIFTS_ONLY=True)
    pair('ib5l18020_asn.fits','ib5l18030_asn.fits', ALIGN_IMAGE = None, GET_SHIFT=True, TWEAKSHIFTS_ONLY=True)
    
    
def mosaic_to_pointing(mosaic_list='GOODS-N-*F140W',
                       pointing='GOODS-N-43-F140W',
                       run_multidrizzle=True, grow=0):
    """ 
    Given an input list of ASN tables that could be combined to form a mosaic,
    find only those that overlap with the input pointing and make a new ASN
    table for that pointing that includes the other overlapping pointings.
    
    This would be identical to using the full mosaic ASN and then using the 
    individual pointing as a reference.  The only difference is that this 
    procedure only includes images that overlap with the output file to 
    save computation time.
    """
    from threedhst.shifts import find_align_images_that_overlap as overlap
     
    list = overlap(pointing+'_drz.fits', mosaic_list+'*drz.fits',
                   ALIGN_EXTENSION=1)
    
    asn_files = []
    for item in list:
        asn_files.append(item.replace('drz','asn'))
    
    threedhst.utils.combine_asn_shifts(asn_files, out_root='mostmp',
                    path_to_FLT='./', run_multidrizzle=False)
    
    #### Make the output area a bit larger than the original to allow for 
    #### extracting grism objects near the edges
    header = pyfits.getheader(pointing+'_drz.fits',1)
    NX = header.get('NAXIS1')
    NY = header.get('NAXIS2')
    RA_CENTER = header.get('CRVAL1')
    DEC_CENTER = header.get('CRVAL2')
    SCALE = np.sqrt((header.get('CD1_1')*3600.)**2+
                    (header.get('CD1_2')*3600.)**2)
    
    ANGLE = np.arctan2(header.get('CD2_1'),
                       header.get('CD2_2'))/2/np.pi*360.
    
    # print 'NX NY RA DEC SCALE'
    # print '%d %d %13.6f %13.6f %6.3f' %(NX, NY, RA_CENTER, DEC_CENTER, SCALE)
    
    threedhst.prep_flt_files.startMultidrizzle('mostmp_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=SCALE, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False,
                 final_outnx=NX+grow, final_outny=NY+grow, 
                 final_rot=0., ra=RA_CENTER, dec=DEC_CENTER) #refimage=pointing+'_drz.fits[1]')
    
    os.remove('mostmp_shifts.txt')
    os.remove('mostmp_asn.fits')
    #### Copy the outputs but not the ASN or shift files
    files=glob.glob('mostmp*')
    for file in files:
        out=file.replace('mostmp',pointing)
        shutil.move(file, out)
        
    
def make_grism_subsets(root='GOODS-S', run_multidrizzle=True, single=None):
    """
    Group grism exposures with the same orientation angle together.  
    
    Use the full direct mosaic (root+'-D_drz.fits') as a MDRZ reference
    image if it exists.
    
    """
    import threedhst.catIO as catIO
    
    info = catIO.Readfile('files.info')
    
    if root=='GOODS-S':
        for i in range(info.N):
            info.targname[i] = info.targname[i].replace('SOUTH','S')
    
    if root=='GOODS-N':
        for i in range(info.N):
            info.targname[i] = info.targname[i].replace('GNGRISM','GOODS-N-')
        
    info.pa_v3 = np.cast[int](np.round(info.pa_v3))
    angles = np.unique(np.round(info.pa_v3))
    
    #### Just redo one angle
    if single is not None:
        angles = [angles[single]]
    
    for angle in angles:
        targets = np.unique(info.targname[info.pa_v3 == angle])
        list = []
        for targ in targets:
            list.append("%s-G141_asn.fits" %(targ))
        #
        out_root = root+'-%03d' %(angle)
        print(out_root)
        
        threedhst.utils.combine_asn_shifts(list, out_root=out_root,
                   path_to_FLT='./', run_multidrizzle=False)
        #
        if run_multidrizzle:
            direct_ref = root+'-F140W_drz.fits'
            if not os.path.exists(direct_ref):
                direct_ref=''
            else:
                direct_ref+='[1]'
            
            threedhst.prep_flt_files.startMultidrizzle(out_root+'_asn.fits',
                 use_shiftfile=True, skysub=False,
                 final_scale=0.06, pixfrac=0.8, driz_cr=False,
                 updatewcs=False, clean=True, median=False,
                 refimage=direct_ref)
                 
def make_targname_asn(asn_file, newfile=True, use_filtname=True, path_to_flt='../RAW/',field='ANY', ext='flt'):
    """
    Take an ASN file like 'ibhm51030_asn.fits' and turn it into 
    'COSMOS-3-F140W_asn.fits'
    """
    asn = threedhst.utils.ASNFile(asn_file)
    
    flt_file = threedhst.utils.find_fits_gz(path_to_flt+'/'+asn.exposures[0]+'_%s.fits' %(ext))
    im = pyfits.open(flt_file)
    
    instrum = im[0].header['INSTRUME']
    if instrum == 'ACS':
        filter = ''
        for fkey in ['FILTER1','FILTER2']:
            if 'CLEAR' not in im[0].header[fkey]:
                filter=im[0].header[fkey]
    else:
        filter=im[0].header['FILTER']
        
    if filter.startswith('F'):
        type='D'
    else:
        type='G'
    
    if use_filtname:
        type=filter
        
    target = im[0].header['TARGNAME']
    #target = target.replace('SOUTH','S')
    #target = target.replace('GNGRISM','GOODS-N-')
    target = target.replace('GEORGE','GEORGE-')
    target = target.replace('ANY',field)
    
    #### 3D-HST translations
    translate = {'AEGIS-':'aegis-', 'COSMOS-':'cosmos-', 'GNGRISM':'goodsn-', 'GOODS-SOUTH-':'goodss-', 'UDS-':'uds-'}
    for key in list(translate.keys()):
        target = target.replace(key, translate[key])
    
    ## pad i < 10 with zero
    for key in list(translate.keys()):
        if translate[key] in target:
            spl = target.split('-')
            if int(spl[-1]) < 10:
                spl[-1] = '%02d' %(int(spl[-1]))
                target = '-'.join(spl)
                
    
    if target == 'MARSHALL':
        #### Add the pointing number, 1-6
        ID = asn.exposures[0][5]
        target+='-'+ID

    if target == 'PRIMO':
        #### Add the date, like "1026"
        date = ''.join(im[0].header['DATE-OBS'].split('-')[1:])
        target+='-'+date
    
    if target.startswith('GEORGE'):
        #### Add the date, like "1026"
        hour = np.int(im[0].header['TIME-OBS'].split(':')[0])
        if hour > 12:
            target='GEORGE-2'
        else:
            target='GEORGE-1'
            
    if instrum == 'ACS':
        ID = asn.exposures[0][4:6]
        target+='-'+ID

    # 
    # if target.startswith('GOODS-N'):
    #     #### Some Visits were redone
    #     date = im[0].header['DATE-OBS']
    #     if date > '2011-04-01':
    #         target = target.replace('GOODS-N','GOODS-N2')
    
    product = target+'-'+type
    print(product)

    asn.product = product
    if newfile:
        asn.write(product+'_asn.fits', clobber=True)
    return product+'_asn.fits'
    
def process_3dhst_pair(asn_direct_file='ib3706050_asn.fits',
                       asn_grism_file='ib3706060_asn.fits',
                       adjust_targname=True,
                       ALIGN_IMAGE='../ACS/h_nz_sect*img.fits',
                       ALIGN_EXTENSION=0,
                       align_geometry='shift',
                       PATH_TO_RAW='../RAW',
                       IMAGES = [os.getenv('THREEDHST')+'/CONF/G141_sky_cleaned.fits',
                                 os.getenv('THREEDHST')+'/CONF/G141wLO_fixed_sky.fits', 
                                 os.getenv('THREEDHST')+'/CONF/G141wHI_fixed_sky.fits'],
                       SKIP_GRISM=False,
                       SKIP_DIRECT=False,
                       GET_SHIFT=True,
                       TWEAKSHIFTS_ONLY=False,
                       DIRECT_HIGHER_ORDER=2,
                       GRISM_HIGHER_ORDER=1,
                       save_fit=False, 
                       second_pass=True, overall=True,
                       sky_images=['sky.G141.set001.fits', 'sky.G141.set002.fits','sky.G141.set003.fits','sky.G141.set004.fits','sky.G141.set005.fits','sky.G141.set025.fits','sky.G141.set120.fits'],
                       final_scale=0.06, clean_drz=False):
        
    #### Old sky_images=['sky_cosmos.fits', 'sky_goodsn_lo.fits', 'sky_goodsn_hi.fits', 'sky_goodsn_vhi.fits']
    import threedhst
    import threedhst.prep_flt_files
    from threedhst.prep_flt_files import make_targname_asn
    
    if (asn_direct_file is not None) & adjust_targname:
        asn_direct_file = make_targname_asn(asn_direct_file)
    
    if (asn_grism_file is not None) & adjust_targname:
        asn_grism_file = make_targname_asn(asn_grism_file)
    
    print('DIRECT: %s, GRISM: %s\n' %(asn_direct_file, asn_grism_file))
        
    ##### Direct images
    if not SKIP_DIRECT:
        
        threedhst.process_grism.fresh_flt_files(asn_direct_file,
                      from_path=PATH_TO_RAW)

        ##### Make region files for the pointing
        if not os.path.exists(asn_direct_file.replace('fits','pointing.reg')):
            threedhst.regions.asn_region(asn_direct_file)
        
        threedhst.prep_flt_files.prep_flt(asn_file=asn_direct_file,
                        get_shift=GET_SHIFT, first_run=True,
                        bg_only=False, bg_skip=False, redo_background=True,
                        ALIGN_IMAGE=ALIGN_IMAGE, 
                        ALIGN_EXT=ALIGN_EXTENSION,
                        skip_drz=False, final_scale=final_scale, pixfrac=0.8,
                        IMAGES=[],
                        align_geometry=align_geometry, clean=True,
                        initial_order=0, save_fit=save_fit,
                        TWEAKSHIFTS_ONLY=TWEAKSHIFTS_ONLY, make_persistence_mask=True, clean_drz=clean_drz)
        
        if DIRECT_HIGHER_ORDER > 0:
            threedhst.prep_flt_files.prep_flt(asn_file=asn_direct_file,
                        get_shift=False, first_run=False,
                        bg_only=False, bg_skip=False, redo_background=False,
                        skip_drz=False, final_scale=final_scale, pixfrac=0.8,
                        IMAGES=[], clean=True,
                        initial_order=DIRECT_HIGHER_ORDER, save_fit=save_fit, clean_drz=clean_drz)
        
    #### Grism images
    if not SKIP_GRISM:
        if asn_direct_file:
            threedhst.shifts.make_grism_shiftfile(asn_direct_file,
                                                  asn_grism_file)
        
        #### Have to account for the second epoch GOODS-N images 
        #### taken to fix the high background
        if asn_direct_file.startswith('GOODS-N'):
            sfd = threedhst.shifts.ShiftFile(asn_direct_file.replace('asn.fits','shifts.txt'))
            sfg = threedhst.shifts.ShiftFile(asn_grism_file.replace('asn.fits','shifts.txt'))
            x0, y0 = sfd.xshift[0], sfd.yshift[0]
            x1, y1 = x0, y0
            for j in range(len(sfd.images)):
                if sfd.images[j].startswith('ib374'):
                    x1 = sfd.xshift[j]
                    y1 = sfd.yshift[j]
                    #
                    for k in range(len(sfg.images)):
                        if sfg.images[k].startswith('ib374'):
                            sfg.xshift[k], sfg.yshift[k] = x1, y1
                        else:
                            sfg.xshift[k], sfg.yshift[k] = x0, y0
                    #
                    sfg.write(sfg.filename)
                    #    
                    break
                
        #
        asn = threedhst.utils.ASNFile(asn_grism_file)
        test = True
        for exp in asn.exposures:
            test = test & (os.path.exists(exp+'_flt.seg.fits'))
        
        if not test:
            threedhst.process_grism.fresh_flt_files(asn_grism_file,
                          from_path=PATH_TO_RAW)

            if not os.path.exists(asn_grism_file.replace('fits','pointing.reg')):
                threedhst.regions.asn_region(asn_grism_file)

            threedhst.prep_flt_files.prep_flt(asn_file=asn_grism_file,
                            get_shift=False, 
                            bg_only=False, bg_skip=False, redo_background=True,
                            skip_drz=False, final_scale=final_scale, pixfrac=0.8,
                            IMAGES=IMAGES, clean=True,
                            initial_order=-1, save_fit=save_fit, clean_drz=clean_drz)
        
        #### Now that we have the segmentation masks for the grism, do the
        #### division by the flat + subtracting the sky image
        skies = " ".join(sky_images)
        # for img in sky_images:
        #     skies += img+" "
        message = """Divide by the flat and subtract the sky images.
Imaging flat: %s
Sky images: %s""" %(threedhst.grism_sky.flat_direct.replace('//','/'), skies)

        threedhst.showMessage(message)
        
        ## Copy new FLT files
        threedhst.process_grism.fresh_flt_files(asn_grism_file, 
                     from_path=PATH_TO_RAW, preserve_dq=False)
        
        ## Run the sky background division             
        asn_grism = threedhst.utils.ASNFile(asn_grism_file)
        for exp in asn_grism.exposures:
            threedhst.grism_sky.remove_grism_sky(flt=exp+'_flt.fits', list=sky_images, path_to_sky=os.getenv('THREEDHST')+'/CONF/', verbose=True, second_pass=second_pass, overall=overall)
        
        ## Run Multidrizzle twice, the first time to flag CRs + hot pixels
        startMultidrizzle(asn_grism_file, use_shiftfile=True, skysub=False,
                final_scale=0.128, pixfrac=1.0, driz_cr=True,
                updatewcs=True, median=True, clean=True)
                       
        startMultidrizzle(asn_grism_file, use_shiftfile=True, skysub=False,
                final_scale=0.06, pixfrac=0.8, driz_cr=False,
                updatewcs=False, median=False, clean=True)
        
        # if GRISM_HIGHER_ORDER > 0:
        #     threedhst.prep_flt_files.prep_flt(asn_file=asn_grism_file,
        #                 get_shift=False, first_run=False,
        #                 bg_only=False, bg_skip=False, redo_background=False,
        #                 skip_drz=False, final_scale=0.06, pixfrac=0.8,
        #                 IMAGES=[], clean=True,
        #                 initial_order=GRISM_HIGHER_ORDER, save_fit=save_fit)
    
    threedhst.showMessage("""FLT prep done.  Run gzip *flt.fits to save disk space.""")
    
def startAstrodrizzle(root='ib3727050_asn.fits', use_shiftfile = True,
        skysub=True, updatewcs=True, driz_cr=True, median=True, final_driz=True, 
        final_scale=0.06, pixfrac=0.8, clean=True,
        final_outnx='', final_outny='', final_rot=0., ra='', dec='', 
        refimage='', unlearn=True, use_mdz_defaults=True, ivar_weights=True, rms_weights=False, build_drz=True, generate_run=False, skyuser=''):
    """
    Use astrodrizzle instead, test....
    """
    import drizzlepac
    from drizzlepac import astrodrizzle
    pass
    
    
def startMultidrizzle(root='ib3727050_asn.fits', use_shiftfile = True,
    skysub=True, updatewcs=True, driz_cr=True, median=True, final_driz=True, 
    final_scale=0.06, pixfrac=0.8, clean=True,
    final_outnx='', final_outny='', final_rot=0., ra='', dec='', 
    refimage='', unlearn=True, use_mdz_defaults=True, ivar_weights=True, rms_weights=False, build_drz=True, generate_run=False, skyuser=''):
    """
startMultidrizzle(root='ib3727050_asn.fits', use_shiftfile = True,
                  skysub=True, final_scale=0.06, updatewcs=True, driz_cr=True,
                  median=True, final_scale=0.06, pixfrac=0.8, 
                  final_outnx='', final_outny='', final_rot=0., ra='', dec='',
                  refimage='', unlearn=True)
    
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
    from pyraf import iraf
    from iraf import stsdas,dither,slitless,axe 

    INDEF = iraf.INDEF
    no = iraf.no
    yes = iraf.yes
    
    asn_direct_file = root #'ib3727050_asn.fits'
    
    asn_direct = threedhst.utils.ASNFile(file=asn_direct_file)
    ROOT_DIRECT = asn_direct_file.split('_asn')[0]
    
    if use_shiftfile:
        shiftfile=ROOT_DIRECT+'_shifts.txt'
    else:
        shiftfile=''
    
    #### If fewer than 4 exposures in the asn list, use
    #### a larger `pixfrac`.
    # if len(asn_direct.exposures) < 4:
    #     pixfrac = 1.0
    
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
    #
    if final_driz:
        driz_combine=yes
    else:
        driz_combine=no
    
    if unlearn:
        iraf.unlearn('multidrizzle')
    
    #### Set default parameters from pipeline mdz file
    if use_mdz_defaults:
        #### Read the first FLT image and read its MDRIZTAB file
        flt = pyfits.open(asn_direct.exposures[0]+'_flt.fits')
        
        #### Get the filter string, WFC3 or ACS
        if flt[0].header['INSTRUME'] == 'WFC3':
            filter=flt[0].header['FILTER']
            REF = 'iref'
        else:
            filter=(flt[0].header['FILTER1']+','+flt[0].header['FILTER2']).strip()
            REF = 'jref'
        
        mdz = pyfits.open(flt[0].header['MDRIZTAB'].replace(REF+'$',os.getenv(REF)+'/'))[1].data
        
        #### Force direct filter because parameters are a bit strange for grisms
        if filter.startswith('G1'):
            filter='F140W'
        
        if filter.startswith('G8'):
            filter='F814W'
        
        #### find 
        idx = np.where(mdz.field('filter') == filter)[0]
        if len(idx) == 0:
            filter='ANY'
            idx = np.where(mdz.field('filter') == filter)[0]
        
        #### Find right column for given "numimages" = len(exposures)  
        use = idx[0]
        for i in idx[1:]:
            if len(asn_direct.exposures) >= mdz.field('numimages')[i]:
                use = i
        
        #### Now set all of the parameters
        for param in mdz.names:
            try:
                value = mdz.field(param)[use]
                if (not np.isfinite(value)) | (value < -1000):
                    value = iraf.INDEF
                #
                iraf.dither.multidrizzle.setParam(param, value)
            except:
                #### some columns in the MDZ file aren't actually parameters, skip
                pass
        
        #### Don't use these default values from the mdrz file
        iraf.dither.multidrizzle.setParam('crbit','')
        iraf.dither.multidrizzle.setParam('combine_type','minmed')
        #iraf.dither.multidrizzle.setParam('combine_type','median')
        iraf.dither.multidrizzle.setParam('mdriztab',iraf.no)
        iraf.dither.multidrizzle.setParam('context',iraf.no)
        iraf.dither.multidrizzle.setParam('clean',iraf.no)
        iraf.dither.multidrizzle.setParam('ra','')
        iraf.dither.multidrizzle.setParam('dec','')
        iraf.dither.multidrizzle.setParam('runfile','')
        # iraf.dither.multidrizzle.setParam('driz_cr_snr','3.5 3.0')
        
    ### Set CR SNR parameter following candels
    if flt[0].header['INSTRUME'] == 'WFC3':
        #iraf.dither.multidrizzle.driz_cr_snr = '6.0 3.0'
        iraf.dither.multidrizzle.driz_cr_snr = '3.5 3.0'
        #iraf.dither.multidrizzle.driz_cr_scale = '1.6 0.7'
        ### More conservative to avoid rejecting central pixels of stars
        iraf.dither.multidrizzle.driz_cr_scale = '2.5 0.7'
        
    #
    if rms_weights:
        #### Generate inverse variance weight map, will need 
        #### flat + dark images in the iref or jref directories
        iraf.dither.multidrizzle.setParam('final_wht_type','IVM')
    
    if ivar_weights:
        #### Generate inverse variance weight map, will need 
        #### flat + dark images in the iref or jref directories
        iraf.dither.multidrizzle.setParam('final_wht_type','IVM')
    
    if build_drz:
        build=iraf.yes
    else:
        build=iraf.no
    
    if generate_run:        
        import multidrizzle
        from pydrizzle import pydrizzle, process_input
        
        md = multidrizzle.Multidrizzle(asn_direct_file, output='',
           shiftfile=asn_direct_file.replace('_asn.fits','_shifts.txt'), editpars=iraf.no, 
           skysub = 0, updatewcs = 0, driz_cr=0,
           driz_final_scale = final_scale, driz_final_pixfrac = pixfrac, median=0, 
           blot=0, driz_separate=1, static=0, driz_combine=0,
           driz_sep_outnx = final_outnx, driz_sep_outny = final_outny, 
           driz_final_outnx=final_outnx, driz_final_outny=final_outny, 
           driz_final_rot=final_rot, ra=ra, dec=dec, refimage=refimage, build=0)
           
        md.build()
        md.image_manager._setOutputFrame(md.driz_final_pars)

        assoc = pydrizzle._PyDrizzle(md.asndict, output=md.output,
                                    idckey=md.coeffs,
                                    section=md.driz_sep_pars['group'],
                                    bits_single=md.driz_sep_bits,
                                    bits_final=md.final_bits,
                                    )
        #
        drizpars = md.driz_final_pars
        _final_shape = (drizpars['outnx'],drizpars['outny'])
        _new_field = pydrizzle.SkyField(shape=_final_shape)
        _new_field.set(psize=drizpars['scale'], orient=drizpars['rot'],
                        ra=drizpars['ra'], dec=drizpars['dec'])
        #
        assoc.resetPars(field=_new_field, 
                    pixfrac=drizpars['pixfrac'], 
                    kernel=drizpars['kernel'], units=drizpars['units'] ) 
        
        runfile = asn_direct_file.replace('_asn.fits','.run')
        runlog = open(runfile,'w')

        runlog.write("drizzle.outnx = "+str(assoc.parlist[0]['outnx'])+"\n")
        runlog.write("drizzle.outny = "+str(assoc.parlist[0]['outny'])+"\n")
        runlog.write("drizzle.scale = "+str(assoc.parlist[0]['scale'])+"\n")
        runlog.write("drizzle.pixfrac = "+str(assoc.parlist[0]['pixfrac'])+"\n")
        runlog.write("drizzle.shft_fr = 'output'\n")
        runlog.write("drizzle.shft_un = 'output'\n")
        runlog.write("drizzle.in_un = "+str(assoc.parlist[0]['in_units'])+"\n")
        runlog.write("drizzle.out_un = '"+assoc.parlist[0]['units']+"'\n")
        runlog.write("drizzle.align = 'center'\n")
        runlog.write("drizzle.expkey = 'EXPTIME'\n")
        runlog.write("drizzle.fillval = "+str(assoc.parlist[0]['fillval'])+"\n")
        runlog.write("drizzle.outcont = '"+assoc.parlist[0]['outcontext']+"'\n")
        runlog.write("drizzle.kernel = '"+assoc.parlist[0]['kernel']+"'\n")
        runlog.write("\n")

        for p in assoc.parlist:
            xsh_str = "%.4f"  % p['xsh']
            ysh_str = "%.4f"  % p['ysh']
            rot_str = "%.5f"  % p['rot']

            print(("\ndrizzle "+p['data']+" "+p['outdata']+
                  " in_mask="+p['driz_mask']+" outweig="+p['outweight']+
                  " xsh="+xsh_str+" ysh="+ysh_str+" rot="+rot_str+
                  " coeffs='"+p['coeffs']+"' wt_scl='"+str(p['wt_scl'])+"'"+
                  " xgeoim='"+p['xgeoim']+"' ygeoim='"+p['ygeoim']+"'\n"))
            
            runlog.write("drizzle "+p['data']+" "+p['outdata']+
                         " in_mask="+p['driz_mask']+" outweig="+p['outweight']+
                         " xsh="+xsh_str+" ysh="+ysh_str+" rot="+rot_str+
                         " coeffs='"+p['coeffs']+"' wt_scl='"+str(p['wt_scl'])+"'"+
                         " xgeoim='"+p['xgeoim']+"' ygeoim='"+p['ygeoim']+"'\n")

        # Close the "runfile" log
        if runlog != None:
            runlog.close()
        
        ### Need to copy back the flt files
        threedhst.utils.replace_OrIg()
        threedhst.process_grism.cleanMultidrizzleOutput()
        
        return md
    
    #iraf.dither.multidrizzle.setParam('final_bits',0)
    
    #### Run Multidrizzle
    iraf.multidrizzle(input=asn_direct_file, \
       shiftfile=shiftfile, \
       output = '', skysub = skysub, updatewcs = updatewcs, driz_cr=driz_cr,
       final_scale = final_scale, final_pixfrac = pixfrac, median=median, 
       blot=median, driz_separate=median, static=median, driz_combine=driz_combine,
       driz_sep_outnx = final_outnx, driz_sep_outny = final_outny, 
       final_outnx=final_outnx, final_outny=final_outny, 
       final_rot=final_rot, ra=ra, dec=dec, refimage=refimage, build=build, 
       skyuser=skyuser)
    
    #### Delete created files    
    if clean is True:
        threedhst.process_grism.cleanMultidrizzleOutput()
            
class MultidrizzleRun():
    """
MultidrizzleRun(root='IB3728050')
    
    Read a .run file output from MultiDrizzle.
    
    Get list of flt files and their shifts as used by multidrizzle.
    """
    def __init__(self, root='IB3728050'):
        
        runfile = root+'.run'
        self.root = root
        
        self.flt = []
        self.xsh = []
        self.ysh = []
        self.rot = []
        self.scl = 1.
        self.exptime = []
        self.outnx = 1.
        self.outny = 1.
        
        for line in open(runfile,'r'):
            if line.startswith('drizzle.outnx'):
                self.outnx = int(line.split()[2])
            if line.startswith('drizzle.outny'):
                self.outny = int(line.split()[2])
            if line.startswith('drizzle.scale'):
                self.scl = line.split()[2]
            if line.startswith('drizzle '):
                spl = line.split()
                self.flt.append(spl[1].split('.fits')[0])
                self.exptime.append(-1)
                for tag in spl:
                    if tag.startswith('xsh'):
                        self.xsh.append(np.float(tag.split('=')[1]))
                    if tag.startswith('ysh'):
                        self.ysh.append(np.float(tag.split('=')[1]))
                    if tag.startswith('rot'):
                        self.rot.append(np.float(tag.split('=')[1]))
        
        self.count = len(self.flt)
        
    def blot_back(self, ii=0, SCI=True, WHT=True, copy_new=True, shape = None, ACS_CHIP=None):
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
        from pyraf import iraf
        from iraf import stsdas,dither,slitless,axe 

        INDEF = iraf.INDEF
        no = iraf.no
        yes = iraf.yes
        
        #flt_orig = pyfits.open('../RAW/'+self.flt[ii]+'.fits.gz')
        threedhst.process_grism.flprMulti()
        
        if ACS_CHIP is None:
            ACS = False
            coeffs_ext = '_coeffs1.dat'
            EXT = 1
        else:
            ACS = True
            coeffs_ext = '_coeffs%0d.dat' %(3-ACS_CHIP)
            if ACS_CHIP == 1:
                EXT = 1
            else:
                EXT = 4
                
        ## Copy the coeffs1.dat file if necessary
        if not os.path.exists(self.flt[ii]+coeffs_ext):
            coeffs = threedhst.utils.get_package_data('wfc3_coeffs1.dat')
            fp = open(self.flt[ii]+'_coeffs1.dat','w')
            fp.writelines(coeffs)
            fp.close()
            
        if self.exptime[ii] < 0:
            try:
                flt_orig = pyfits.open(self.flt[ii]+'.fits')
                exptime = flt_orig[0].header.get('EXPTIME')
                if not ACS:
                    filter = flt_orig[0].header.get('FILTER').strip()
                else:
                    filter = (flt_orig[0].header.get('FILTER1').strip(),flt_orig[0].header.get('FILTER2').strip())
                
                flt_orig.close()
            except:
                exptime = 1.
                filter='INDEF'
        else:
            exptime = self.exptime[ii]
        
        if shape is None:  
            try:
                inNX = flt_orig[EXT].header.get('NAXIS1')
                inNY = flt_orig[EXT].header.get('NAXIS2')
            except:
                inNX = 1014
                inNY = 1014
            
            shape = (inNX, inNY)
        else:
            inNX, inNY = shape
        
        #### Need to update reference position of coeffs file
        #### for an output shape different than 1014, 1014
        coeffs = self.flt[ii]+coeffs_ext #'_coeffs1.dat'
        fp = open(coeffs)
        coeffs_lines = fp.readlines()
        fp.close()
        
        if shape != (1014, 1014):
            for i, line in enumerate(coeffs_lines):
                if line.strip().startswith('refpix'):
                    ### Default to center pixel
                    coeffs_lines[i] = 'refpix %9.3f %9.3f\n' %(inNX*1./2, inNY*1./2)
        
        fp = open('tmp'+coeffs_ext,'w')
        fp.writelines(coeffs_lines)
        fp.close()
        
                
        #iraf.delete(self.flt[ii]+'.BLOT.*.fits')
        files = glob.glob(self.flt[ii]+'.BLOT.*.fits')
        for file in files:
            os.remove(file)
            
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
            w_hdu = pyfits.PrimaryHDU(1./wht)
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
                coeffs='tmp'+coeffs_ext, xsh=self.xsh[ii], 
                ysh=self.ysh[ii], 
                rot=self.rot[ii], outnx=inNX, outny=inNY, align='center', 
                shft_un='input', shft_fr='input', in_un='cps', out_un='cps', 
                interpol='poly5', sinscl='1.0', expout=exptime, 
                expkey='EXPTIME',fillval=0.0)
        
        if WHT:
            iraf.blot(data='drz_wht.fits',
                outdata=self.flt[ii]+'.BLOT.WHT.fits', scale=self.scl,
                coeffs='tmp'+coeffs_ext, xsh=self.xsh[ii], 
                ysh=self.ysh[ii], 
                rot=self.rot[ii], outnx=inNX, outny=inNY, align='center', 
                shft_un='input', shft_fr='input', in_un='cps', out_un='cps', 
                interpol='poly5', sinscl='1.0', expout=exptime, 
                expkey='EXPTIME',fillval=0.0)
            #
            im = pyfits.open(self.flt[ii]+'.BLOT.WHT.fits', mode='update')
            bad = im[0].data <= 0
            im[0].data = np.sqrt(im[0].data) #/im[0].header['EXPTIME']
            im[0].data[bad] = 1.e7
            im.flush()
            
        #iraf.delete('drz_*.fits')
        
class DRZFile(MultidrizzleRun):
    """
    Get the information from a drz file directly, rather than from a .run 
    file
    """
    def __init__(self, fitsfile='ib6o23010_drz.fits'):
        
        self.root = fitsfile.split('_drz.fits')[0]
        
        drz = pyfits.open(fitsfile)
        hdrz = drz[0].header
        self.count = 0
        for key in list(hdrz.keys()):
            if key.startswith('D') & key.endswith('XSH'):
                self.count+=1
        
        self.flt = []
        self.xsh = []
        self.ysh = []
        self.rot = []
        self.scl = 0.
        self.exptime = []
        
        for i in range(self.count):
            self.flt.append(hdrz.get('D%03dDATA' %(i+1)).split('.fits')[0])
            self.xsh.append(hdrz.get('D%03dXSH' %(i+1)))
            self.ysh.append(hdrz.get('D%03dYSH' %(i+1)))
            self.rot.append(hdrz.get('D%03dROT' %(i+1)))
            self.exptime.append(hdrz.get('D%03dDEXP' %(i+1)))
            self.scl += hdrz.get('D%03dSCAL' %(i+1))
        
        self.scl /= self.count
        
def jitter_info():
    """
jitter_info()
    
    Get LIMBANG values from jitter files and also get 
    image stats from FLT images.  Useful for flagging 
    exposures affected by earthglow.
    """
    import glob
    import os
    
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
            print(str)
    
    fp.close()


# def go_make_segmap():
#     
#     import glob
#     
#     files = glob.glob('*BLOT.SCI.fits')
#     for file in files:
#         make_segmap(root=file.split('.BLOT')[0])
        
def make_segmap(root='ib3701ryq_flt', sigma=1.1, IS_GRISM=None, grow_size=5):
    """
make_segmap(root='ib3701ryq_flt', sigma=1)
    
    Get a segmentation image for a flt file after creating its 
    BLOT SCI and WHT images.
    
    DETECT_THRESH = ANALYSIS_THRESH = sigma
    """
    import threedhst
    import scipy.ndimage as nd
    
    ## Find if image is for grism or direct image
    if IS_GRISM is None:
        flt = pyfits.open(root+'.fits')
        IS_GRISM = flt[0].header.get('FILTER').startswith('G')
        flt.close()
    
    se = threedhst.sex.SExtractor()
    ## Set the output parameters required for aXe 
    ## (stored in [threedhst source]/data/aXe.param) 
    se.aXeParams()
    ## XXX add test for user-defined .conv file
    se.copyConvFile(grism=IS_GRISM)
    
    se.overwrite = True
    se.options['CATALOG_NAME']    = root+'.BLOT.SCI.cat'
    se.options['CHECKIMAGE_NAME'] = root+'.seg.fits, bg.fits'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION, BACKGROUND'
    se.options['WEIGHT_TYPE']     = 'MAP_RMS'
    if IS_GRISM:
        se.options['WEIGHT_TYPE'] = 'NONE'
        
    se.options['WEIGHT_IMAGE']    = root+'.BLOT.WHT.fits'
    se.options['FILTER']    = 'Y'

    #se.options['BACK_TYPE']     = 'MANUAL'
    #se.options['BACK_TYPE']     = 'AUTO'
    #se.options['BACK_FILTERSIZE']     = '2'
    
    if IS_GRISM:
        se.options['FILTER_NAME'] = 'grism.conv'
    else:
        se.options['FILTER_NAME'] = threedhst.sex.USE_CONVFILE #'default.conv'
    
    #### Detect thresholds (default = 1.5)
    se.options['DETECT_THRESH']    = '%f' %sigma
    se.options['ANALYSIS_THRESH']  = '%f' %sigma
    se.options['MAG_ZEROPOINT'] = '26.46'
    status = se.sextractImage(root+'.BLOT.SCI.fits')
    
    seg = pyfits.open(root+'.seg.fits', mode='update')
    seg[0].data = nd.maximum_filter(seg[0].data, size=grow_size)
    seg.flush()
    
    if os.path.exists(root+'.seg.fits.mask.reg'):
        threedhst.regions.apply_dq_mask(root+'.seg.fits', extension=0,
           addval=100)
           
def apply_best_flat(fits_file, verbose=False, use_cosmos_flat=True, use_candels_flat=True, apply_BPM=True):
    """
    Check that the flat used in the pipeline calibration is the 
    best available.  If not, multiply by the flat used and divide
    by the better flat.
    
    Input fits_file can either be an ASN list or an individual FLT file
     """
    fits_list = [fits_file]
    
    if fits_file.find('_asn.fits') > 0:
        asn = threedhst.utils.ASNFile(fits_file)
        fits_list = []
        for exp in asn.exposures:
            fits_list.append(exp+'_flt.fits')
    
    for file in fits_list:
        im = pyfits.open(file, 'update')
        if im[0].header['INSTRUME'] == 'ACS':
            return 'ACS'
        
        USED_PFL = im[0].header['PFLTFILE'].split('$')[1]
        BEST_PFL = find_best_flat(file, verbose=False)
        
        if (use_cosmos_flat) & (im[0].header['DATE'] > '2010-08-01') & (im[0].header['FILTER'] == 'F140W'):
            #### Updated F140W flat from COSMOS
            #BEST_PFL = 'cosmos_f140w_flat.fits'
            ### Time dependent flats
            BEST_PFL = 'flat_3DHST_F140W_t1_v0.1.fits'
            
            if im[0].header['EXPSTART'] > 55793.:
                BEST_PFL = 'flat_3DHST_F140W_t2_v0.1.fits'
            
            if im[0].header['EXPSTART'] > 55911.:
                BEST_PFL = 'flat_UDF_F140W_v0.fits'
            
        if (use_candels_flat) & (im[0].header['FILTER'] == 'F125W'):
            #BEST_PFL = 'flat.F125W.fits'
            BEST_PFL = 'flat_F125W_t1_v0.3.fits'
            if im[0].header['EXPSTART'] > 55780.:
                BEST_PFL = 'flat_F125W_t2_v0.3.fits'
                
        if (use_candels_flat) & (im[0].header['FILTER'] == 'F160W'):
            #BEST_PFL = 'flat.F160W.fits'
            BEST_PFL = 'flat_F160W_t1_v0.3.fits'
            if im[0].header['EXPSTART'] > 55780.:
                BEST_PFL = 'flat_F160W_t2_v0.3.fits'
           
        IREF = os.environ["iref"]+"/"
        
        MSG = ''
        if apply_BPM:
            my_bpm = pyfits.open('%s/flat_BPM_v0.1.fits' %(os.environ['iref']))[0].data
            im[3].data[my_bpm > 0] |= (100+4096)
            im[1].data += my_bpm*1000000 ## make new bad pixels obviously bad
            ##im[3].data[im[1].data > 1000000] |= 4096
                            
            MSG = 'extra BPM: ${iref}/flat_BPM_v0.1.fits\n'
            
            if im[0].header['DATE-OBS'] > '2012-01-01':
                more_bpm = pyfits.open('%s/badpix_spars200_Nov9.fits' %(os.environ['iref']))[0].data
                im[3].data[more_bpm > 0] |= (100+4096)
                im[1].data += my_bpm*1000000 ## make new bad pixels obviously bad
                MSG += '           ${iref}/badpix_spars200_Nov9.fits\n'
            
            # print 'BP flag 1', im[3].data[im[1].data > 1.e4].min()
            
        MSG += 'PFLAT, %s: Used= %s, Best= %s' %(file, USED_PFL, BEST_PFL)
        
        if BEST_PFL is None:
            threedhst.showMessage("No PFL file found! (NEED %s)" %(USED_PFL), warn=True)
            return
                    
        BEST_PFL = os.path.basename(BEST_PFL)
                
        if USED_PFL != BEST_PFL:
            MSG += ' *'
            used = pyfits.open(IREF+USED_PFL)
            best = pyfits.open(IREF+BEST_PFL)
            
            im[1].data *= (used[1].data/best[1].data)[5:-5,5:-5]
            im[0].header.update('PFLTFILE', 'iref$'+BEST_PFL)
            #print 'BP flag 2', im[3].data[im[1].data > 1.e4].min()
            im.flush()
                    
        if verbose:
            print(MSG)

def apply_persistence_mask(flt_file, limit_sigma=0.6, filter_size=3, persistence_path='/3DHST/Spectra/Work/PERSISTENCE/All', verbose=False):
    """
    Mask WFC3-IR persistence mask if a "persist" file exists for the 
    specified FLT file.
    """
    import scipy.ndimage as nd
    persist_file = os.path.join(persistence_path, flt_file.replace('flt','persist'))
    
    if not os.path.exists(persist_file):
        if verbose:
            print('%s not found, ignoring persistence.' %(persist_file))
        #
        return 
        
    flt = pyfits.open(flt_file, mode='update')
    persist = pyfits.open(persist_file)
    
    msg = 'P'
    if threedhst.options['FLT_PERSISTENCE_SUBTRACT']:
        msg = 'Subtracted P'
        flt[1].data -= persist[1].data
        #print 'Subtracted persistence mask: %s' %(persist_file)
        
    mask = (persist[1].data > (limit_sigma*flt['ERR'].data))*1
    
    #### Unflagged cosmic rays
    cr = (nd.convolve(mask, np.ones((3,3))) == 1) & mask & ((flt['DQ'].data & 4096) == 0)
    mask[cr > 0] = 0
    mask_grow = nd.maximum_filter(mask, filter_size)
    flt['DQ'].data[mask_grow > 0] |= (100+4096)
    flt.flush()
    
    if verbose:
        print(msg+'ersistence mask: %s, S/N > %.1f [%d masked pixels]' %(persist_file, limit_sigma, mask_grow.sum()))
    
def find_best_flat(flt_fits, verbose=True): #, IREF='/research/HST/GRISM/IREF/'):
    """
    Find the most recent PFL file in $IREF for the filter used for the 
    provided FLT image.  Doesn't do any special check on USEAFTER date, just
    looks for the most-recently modified file. 
    """
    import glob
    import os.path
    import time
    
    IREF = os.environ["iref"]+"/"
    
    the_filter = pyfits.getheader(flt_fits,0).get('FILTER')
    
    pfls = glob.glob(IREF+'/*pfl.fits')
    latest = 0
    best_pfl = None
    
    for pfl in pfls:
        head = pyfits.getheader(pfl)
        if head.get('FILTER') != the_filter:
            continue    
        
        this_created = os.path.getmtime(pfl)
        if this_created > latest:
            best_pfl = pfl
            latest = this_created
            
        if verbose:
            print('%s %s %s' %(pfl, the_filter, time.ctime(latest)))
    
    return best_pfl #, the_filter, time.ctime(latest)

def prep_acs(root='jbhj01',force=False):
    """
    Apply destripe and CTE correction to ACS images.
    
    Run this in a 'RAW' directory and results will be put in ../FIXED.
    
    *** requires irafx, i.e. type 'irafx' before starting pyraf ***
    
    """
    from pyraf import iraf
    from iraf import stsdas,dither,slitless,axe 

    INDEF = iraf.INDEF
    no = iraf.no
    yes = iraf.yes
    
    from acstools import acs_destripe
    from acstools import PixCteCorr
    
    if not (os.getcwd().endswith('RAW') or force):
        threedhst.showMessage("CWD is not 'RAW'.  Run with force=True if this is OK.", warn=True)
        return
    
    #### Make a ../FIXED directory if it doesn't exist
    if not os.path.exists('../FIXED'):
        print("Making output directory ../FIXED")
        os.mkdir('../FIXED')

    files = glob.glob('%s*_flt.fits' %(root))
    print(files)
    for file in files:
        #epar acs_destripe # *flt.fits dstrp clobber+
        flt = pyfits.open(file)
        if (not os.path.exists(file.replace('flt','flt_dstrp'))) & ('PCTEFILE' not in list(flt[0].header.keys())):
            acs_destripe.clean(file,'dstrp',clobber=True)
        
    files=glob.glob('%s*dstrp*' %(root))
    print(files)
    for file in files:
        iraf.hedit(images=file+'[0]', fields='PCTEFILE',
            value='jref$pctefile_101109.fits', add=iraf.yes, verify=iraf.no, 
            update=iraf.yes, show=iraf.no)
    
    ## epar PixCteCorr # *dstrp*fits ''
    PixCteCorr.CteCorr('%s*dstrp*fits' %(root))
    files=glob.glob('%s*cte.fits' %(root))
    for file in files:
        out=file.replace('cte','flt')
        shutil.move(file,'../FIXED/'+out)
    
    ### Clean up dstrp files
    files=glob.glob('%s*dstrp*fits' %(root))
    for file in files:
        print(file)
        os.remove(file)

