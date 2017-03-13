"""
3DHST.shifts

Utilities for computing shifts and processing shiftfiles.
"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import os

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import numpy as np

import matplotlib.pyplot as plt

#### Set to False if you don't have GUI access (Mac) on the machine, i.e. running over SSH.
USE_PLOT_GUI=False
from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

import threedhst
import threedhst.prep_flt_files

#### DS9 reference catalog for use with refine_shifts
REFERENCE_CATALOG = 'sdss8'

def run_tweakshifts(asn_direct, verbose=False, clean=True):
    """
run_tweakshifts(asn_direct)
    
    asn_direct - filename of ASN table of direct images [...]_asn.fits
    
    This routine only uses dither.tweakshifts to compute the relative shifts of 
    the direct images
    """
    from pyraf import iraf
    from iraf import stsdas,dither

    no = iraf.no
    yes = iraf.yes
    INDEF = iraf.INDEF
    
    root = asn_direct.split('_asn.fits')[0]#.lower()

    try:
        os.remove(root+'_tweak.fits')
    except:
        pass        

    iraf.flpr()
    iraf.flpr()
    iraf.flpr()
    if clean:
        clean=iraf.yes
    else:
        clean=iraf.no
    
    iraf.unlearn('tweakshifts')
    
    status = iraf.tweakshifts(input=asn_direct, shiftfile='',
                     reference=root+'_tweak.fits',
                     output = root+'_shifts.txt', findmode = 'catalog',
                     gencatalog = 'daofind', sextractpars = '', 
                     undistort = yes, computesig = yes, idckey = 'idctab', \
       clean = clean, verbose = no, catfile = '', xcol = 1, ycol = 2, \
       fluxcol = 3, fluxmax = INDEF, fluxmin = INDEF, fluxunits = 'counts', \
       nbright = INDEF, refcat = '', refxcol = 1, refycol = 2, rfluxcol = 3, \
       rfluxmax = INDEF, rfluxmin = INDEF, rfluxunits = 'counts', \
       refnbright = INDEF, minobj = 15, nmatch = 30, matching = 'tolerance', \
       xyxin = INDEF, xyyin = INDEF, tolerance = 4.0, fwhmpsf = 1.5, \
       sigma = 0.0, datamin = INDEF, datamax = INDEF, threshold = 4.0, \
       nsigma = 1.5, fitgeometry = 'shift', function = 'polynomial', \
       maxiter = 3, reject = 3.0, crossref = '', margin = 50, tapersz = 50, \
       pad = no, fwhm = 7.0, ellip = 0.05, pa = 45.0, fitbox = 7, \
    Stdout=1)
    
    if verbose:
        for line in status:
            print(line)
    
    return status
    
def find_align_images_that_overlap(DIRECT_MOSAIC, ALIGN_IMAGE, ALIGN_EXTENSION=0, is_region=False, show=False):
    """
align_img_list = find_align_images_that_overlap()
    
    Look through the images defined in threedhst.options['ALIGN_IMAGE'] and 
    return a list of only those that overlap.  
    
    This was written for the GOODS fields, where the multiple GOODS ACS tiles
    can be used to align the F140W images.  ALIGN_IMAGE will be something like 
    h_nz_sect*, but you don't want to waste time swarping large images that 
    don't overlap with the target image.
    """
    import glob
    import threedhst.regions
    #ROOT_DIRECT = threedhst.options['ROOT_DIRECT']
    align_images = glob.glob(ALIGN_IMAGE)
    
    #### Get polygon of the direct mosaic edges
    px, py = threedhst.regions.wcs_polygon(DIRECT_MOSAIC, extension=1)
    if show:
        plt.plot(px, py, color='blue', alpha=0.8, linewidth=3)
        
    #### Loop through align_images and check if they overlap with the 
    #### direct mosaic
    align_img_list = []
    for align_image in align_images:
        if is_region:
            #### Read a polygon defined in a region file
            fp = open(align_image)
            lines = fp.readlines()
            fp.close()
            for line in lines:
                if line.startswith('polygon'):
                    spl = line.split('(')[1].split(')')[0].split(',')
                    qx = np.cast[float](spl[::2])
                    qy = np.cast[float](spl[1::2])
                    break
            align_image = align_image.replace('pointing.reg','fits')
        else:
            qx, qy = threedhst.regions.wcs_polygon(align_image,
                extension=ALIGN_EXTENSION)
        #    
        if threedhst.regions.polygons_intersect(px, py, qx, qy):
            if show:
                plt.plot(qx, qy, color='green', alpha=0.8, linewidth=3)
            align_img_list.append(align_image)
        else:
            if show:
                plt.plot(qx, qy, color='red', alpha=0.8, linewidth=3)
            
    return align_img_list

def refine_shifts(ROOT_DIRECT='f160w',
                  ALIGN_IMAGE='../../ACS/h_sz*drz_img.fits',
                  fitgeometry='shift', clean=True,
                  ALIGN_EXTENSION=0, shift_params=None,
                  toler=3, maxtoler=5, align_sdss_ds9=False,
                  verbose=False):
    """
refine_shifts(ROOT_DIRECT='f160w',
              ALIGN_IMAGE='../../ACS/h_sz*drz_img.fits',
              fitgeometry='shift', clean=True)
                
    Refine shifts by catalog matching an input multidrizzle image, 
    ROOT_DIRECT+'_drz.fits' to one or more alignment images
    """
        
    run = threedhst.prep_flt_files.MultidrizzleRun(ROOT_DIRECT.upper())
    
    ## radius for match is 2**toler.  Make it larger if fit comes out bad
    #toler, maxtoler = 3, 5  
    iter, MAXIT = 0, 5
    xrms, yrms = 100, 100
    if shift_params is not None:
        xshift, yshift, rot, scale = shift_params
        threedhst.showMessage('Using specified DRZ-frame shifts: %f %f %f %f' %(xshift, yshift, rot, scale))
    else:
        threedhst.showMessage('Aligning WCS to %s (%s)' %(threedhst.options['ALIGN_IMAGE'], fitgeometry))
        while ((xrms > 1) | (yrms > 1)) & (toler <= maxtoler) & (iter < MAXIT):
            iter = iter + 1
            xshift, yshift, rot, scale, xrms, yrms = threedhst.shifts.align_to_reference(
                        ROOT_DIRECT,
                        ALIGN_IMAGE,
                        fitgeometry=fitgeometry, clean=clean,
                        ALIGN_EXTENSION=ALIGN_EXTENSION,
                        toler=toler, skip_swarp=(toler > 3),
                        align_sdss_ds9=align_sdss_ds9, verbose=verbose)
            toler+=1

    #### shifts measured in DRZ frame.  Translate to FLT frame
    drz = pyfits.open(ROOT_DIRECT+'_drz.fits')
    #alpha = (180.-drz[1].header['PA_APER'])/360.*2*np.pi
    #### Get reference angle from first image in the ASN file
    asn = threedhst.utils.ASNFile(ROOT_DIRECT+'_asn.fits')
    alpha = (180.-pyfits.getheader(asn.exposures[0]+'_flt.fits',1)['PA_APER'])/360.*2*np.pi
    
    xsh = (xshift*np.cos(alpha)-yshift*np.sin(alpha))*np.float(run.scl)
    ysh = (xshift*np.sin(alpha)+yshift*np.cos(alpha))*np.float(run.scl)

    print('Final shift:', xsh, ysh, drz[1].header['PA_APER'])
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
    
    shiftF.write(ROOT_DIRECT+'_shifts.txt')
    
def plot_shifts(ROOT_DIRECT, ALIGN_IMAGE, clean=True, verbose=True, ALIGN_EXTENSION=0, toler=3, skip_swarp=False, threshold=7, force=False, drz=True, WEIGHT_IMAGE = None):
    """
    Run SExtractor on two images and match the objects to plot the shifts between them.
    
    ALIGN_IMAGE is a string that may contain wildcards, and the function will use 
    `align_img_list` to find ALIGN_IMAGEs
    """
    import glob
    
    from pyraf import iraf
    from iraf import stsdas,dither

    no = iraf.no
    yes = iraf.yes
    INDEF = iraf.INDEF
    
    if os.path.exists(ROOT_DIRECT+'_align.fits') & (not force):
        if verbose:
            print('Image %s_align.fits exists.  Skipping SWarp.' %(ROOT_DIRECT))
        skip_swarp = True
        
    if not skip_swarp:
        if drz:
            align_img_list = find_align_images_that_overlap(ROOT_DIRECT+'_drz.fits', ALIGN_IMAGE, ALIGN_EXTENSION=ALIGN_EXTENSION)
        else:
            align_img_list = glob.glob(ALIGN_IMAGE)
            
        if not align_img_list:
            print('threedhst.shifts.align_to_reference: no alignment images overlap.')
            return 0,0
        #
        try:
            os.remove(ROOT_DIRECT+'_align.fits')
        except:
            pass
        
        if drz:
            matchImagePixels(input=align_img_list, matchImage=ROOT_DIRECT+'_drz.fits', output=ROOT_DIRECT+'_align.fits', match_extension = 1, input_extension=ALIGN_EXTENSION)
            ALIGN_FITS = ROOT_DIRECT+'_align.fits'
        else:
            ALIGN_FITS = os.path.basename(ROOT_DIRECT.split('.fits')[0])+'_align.fits'
            matchImagePixels(input=align_img_list, matchImage=ROOT_DIRECT, output=ALIGN_FITS, match_extension = 0, input_extension=ALIGN_EXTENSION)
            
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CHECKIMAGE_TYPE'] = 'NONE'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    
    if drz:
        se.options['WEIGHT_IMAGE']    = ROOT_DIRECT+'_drz.fits[1]'
    else:
        if WEIGHT_IMAGE:
            se.options['WEIGHT_IMAGE']    = WEIGHT_IMAGE
        else:
            se.options['WEIGHT_TYPE']    = 'NONE'
            se.options['WEIGHT_IMAGE']    = 'NONE'
            
    se.options['FILTER']    = 'Y'
    ## Detect thresholds (default = 1.5)
    se.options['DETECT_THRESH']    = '%f' %(threshold) 
    se.options['ANALYSIS_THRESH']  = '%f' %(threshold)
    se.options['MAG_ZEROPOINT'] = str(threedhst.options['MAG_ZEROPOINT'])

    #### Run SExtractor on direct and alignment images
    ## direct image
    se.options['CATALOG_NAME']    = 'direct.cat'
    if drz:
        status = se.sextractImage(ROOT_DIRECT+'_drz.fits[0]')
        INPUT_IMAGE = ROOT_DIRECT+'_drz.fits'
    else:
        status = se.sextractImage(ROOT_DIRECT)
        INPUT_IMAGE = ROOT_DIRECT
        
    ## alignment image
    se.options['CATALOG_NAME']    = 'align.cat'
    se.options['WEIGHT_TYPE']     = 'NONE'
    status = se.sextractImage(ALIGN_FITS)

    ## Read the catalogs
    directCat = threedhst.sex.mySexCat('direct.cat')
    alignCat = threedhst.sex.mySexCat('align.cat')
    
    xshift = 0
    yshift = 0
    rot = 0
    scale = 1.
    
    xrms = 2
    yrms = 2
    
    NITER = 5
    IT = 0
    while (IT < NITER):
        IT = IT+1
        
        #### Get x,y coordinates of detected objects
        ## direct image
        fp = open('direct.xy','w')
        for i in range(len(directCat.X_IMAGE)):
            fp.write('%s  %s\n' %(directCat.X_IMAGE[i],directCat.Y_IMAGE[i]))
        fp.close()

        ## alignment image
        fp = open('align.xy','w')
        for i in range(len(alignCat.X_IMAGE)):
            fp.write('%s  %s\n' %(np.float(alignCat.X_IMAGE[i])+xshift,
                       np.float(alignCat.Y_IMAGE[i])+yshift))
        fp.close()

        iraf.flpr(); iraf.flpr(); iraf.flpr()
        
        #### iraf.xyxymatch to find matches between the two catalogs
        pow = toler*1.
        try:
            os.remove('align.match')
        except:
            pass
            
        status1 = iraf.xyxymatch(input="direct.xy", reference="align.xy", output="align.match", tolerance=2**pow, separation=0, verbose=yes, Stdout=1)
        
        while status1[-1].startswith('0'):
            pow+=1
            os.remove('align.match')
            status1 = iraf.xyxymatch(input="direct.xy", reference="align.xy", output="align.match", tolerance=2**pow, separation=0, verbose=yes, Stdout=1)
    
    #### Images are aligned, plot the offsets
    
    dx, dy, ax, ay, di, ai = np.loadtxt('align.match', unpack=True)
    
    ddx,ddy = dx-ax, dy-ay
    keep = (np.abs(ddx) < 15) & (np.abs(ddy) < 15)
    for i in range(5):
        sx, sy = threedhst.utils.biweight(ddx[keep], both=True), threedhst.utils.biweight(ddy[keep], both=True)
        keep = keep & (np.abs(ddx-sx[0]) < 5*sx[1]) & (np.abs(ddy-sy[0]) < 5*sy[1])
        
    if USE_PLOT_GUI:
        fig = plt.figure(figsize=[8,4],dpi=100)
    else:
        fig = Figure(figsize=[8,4], dpi=100)
    
    fig.subplots_adjust(wspace=0.28, hspace=0.0, left=0.08, bottom=0.14, right=0.98, top=0.98)
    
    ax = fig.add_subplot(121)
    
    ax.plot(ddx, ddy, marker='o', linestyle='None', color='black', alpha=0.4, ms=2, zorder=-1)
    ax.errorbar(sx[0],sy[0],sx[1],sy[1], marker='o', ms=0.1, color='white', linewidth=3, zorder=100)
    ax.errorbar(sx[0],sy[0],sx[1],sy[1], marker='o', ms=0.1, color='red', linewidth=2, zorder=500)
    
    ax.grid(alpha=0.5)
    dwin = 5*np.max([sx[1],sy[1]])
    ax.set_xlim(sx[0]-dwin,sx[0]+dwin)
    ax.set_ylim(sy[0]-dwin,sy[0]+dwin)
    ax.set_xlabel(r'$\Delta x$ [pix]')
    ax.set_ylabel(r'$\Delta y$ [pix]')
    ax.text(0.5,0.95,os.path.basename(INPUT_IMAGE), fontsize=9, horizontalalignment='center', transform=ax.transAxes)
    ax.text(0.5,0.9,os.path.basename(ALIGN_IMAGE), fontsize=9, horizontalalignment='center', transform=ax.transAxes)
    
    ax.text(0.5,0.1,r'$\Delta x, \Delta y = %.2f \pm %.2f, %.2f \pm %.2f)$' %(sx[0],sx[1],sy[0],sy[1]), fontsize=11, horizontalalignment='center', transform=ax.transAxes)
    
    ax = fig.add_subplot(122)
    
    ax.plot(dx[keep], dy[keep], marker='o', ms=1, linestyle='None', color='black', alpha=0.1)
    ax.quiver(dx[keep], dy[keep], ddx[keep], ddy[keep], alpha=0.5, angles='xy', headlength=1, headwidth=1, scale=30./(dx.max()-dx.min()), units='x', minlength=1)
    aa = np.array([1,1])
    ax.quiver(dx[keep].mean()*aa, dy[keep].max()*0.95*aa, 1*aa, 0*aa, alpha=0.9, angles='xy', headlength=0, headwidth=1, scale=30./(dx.max()-dx.min()), units='x', color='red')
    
    ax.set_xlabel(r'$x$ [pix]')
    ax.set_ylabel(r'$y$ [pix]')
    
    if USE_PLOT_GUI:
        fig.savefig(os.path.basename(ROOT_DIRECT.split('.fits')[0])+'_align.pdf',dpi=100,transparent=False)
    else:
        canvas = FigureCanvasAgg(fig)
        canvas.print_figure(os.path.basename(ROOT_DIRECT.split('.fits')[0])+'_align.pdf', dpi=100, transparent=False)
    
    if clean:
        rmfiles = ['SCI.fits','WHT.fits','align.cat',
               'align.map','align.match','align.reg','align.xy',
               'direct.cat','direct.reg','direct.xy',
               'drz_sci.fits','drz_wht.fits','bg.fits']
        
        for file in rmfiles:
            try:
                os.remove(file)
            except:
                pass
       
def align_to_reference(ROOT_DIRECT, ALIGN_IMAGE, fitgeometry="shift",
    clean=True, verbose=False, ALIGN_EXTENSION=0, toler=3, skip_swarp=False,
    align_sdss_ds9=False, catalog=None):
    """
xshift, yshift, rot, scale, xrms, yrms = align_to_reference()
    """        
    import os
    import glob
    import shutil
    
    from pyraf import iraf
    from iraf import stsdas,dither
    
    import threedhst
    from threedhst import catIO
    
    no = iraf.no
    yes = iraf.yes
    INDEF = iraf.INDEF
    
    #### Clean slate    
    rmfiles = ['SCI.fits','WHT.fits','align.cat','direct.cat'
               'align.map','align.match','align.reg','align.xy', 
               'direct.reg','direct.xy','ds9_align.tsv']
    
    for file in rmfiles:
        try:
            os.remove(file)
        except:
            pass
    
    if catalog is not None: 
        align_sdss_ds9 = True
                    
    #### Get only images that overlap from the ALIGN_IMAGE list    
    if not align_sdss_ds9:
        align_img_list = find_align_images_that_overlap(ROOT_DIRECT+'_drz.fits', ALIGN_IMAGE, ALIGN_EXTENSION=ALIGN_EXTENSION)
        if not align_img_list:
            print('threedhst.shifts.align_to_reference: no alignment images overlap.')
            return 0,0
    
    #### Use swarp to combine the alignment images to the same image 
    #### dimensions as the direct mosaic
    if (not skip_swarp) & (not align_sdss_ds9):
        try:
            os.remove(ROOT_DIRECT+'_align.fits')
        except:
            pass
        matchImagePixels(input=align_img_list,
                     matchImage=ROOT_DIRECT+'_drz.fits',
                     output=ROOT_DIRECT+'_align.fits', match_extension = 1,
                     input_extension=ALIGN_EXTENSION)
                     
    #### Run SExtractor on the direct image, with the WHT 
    #### extension as a weight image
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CHECKIMAGE_TYPE'] = 'NONE'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = 'WHT.fits'
    se.options['FILTER']    = 'Y'
    ## Detect thresholds (default = 1.5)
    THRESH = 10
    if align_sdss_ds9:
        if 'Vizier' not in REFERENCE_CATALOG:
            THRESH = 20
            
    se.options['DETECT_THRESH']    = '%d' %(THRESH)
    se.options['ANALYSIS_THRESH']  = '%d' %(THRESH)
    se.options['MAG_ZEROPOINT'] = str(threedhst.options['MAG_ZEROPOINT'])

    #### Run SExtractor on direct and alignment images
    ## direct image
    se.options['CATALOG_NAME']    = 'direct.cat'
    iraf.imcopy(ROOT_DIRECT+'_drz.fits[1]',"SCI.fits")
    iraf.imcopy(ROOT_DIRECT+'_drz.fits[2]',"WHT.fits")
    status = se.sextractImage('SCI.fits')

    ## Read the catalog
    directCat = threedhst.sex.mySexCat('direct.cat')

    if align_sdss_ds9:
        ### Use ds9 SDSS catalog to refine alignment
        import threedhst.dq
        import pywcs
        import threedhst.catIO as catIO
        
        wcs = pywcs.WCS(pyfits.getheader('SCI.fits', 0))
        #wcs = pywcs.WCS(pyfits.getheader('Q0821+3107-F140W_drz.fits', 1))
        
        if 'Vizier' in REFERENCE_CATALOG:
            #### Use (unstable) astroquery Vizier search
            #### CFHTLS-Deep: 'Vizier.II/317'
            VIZIER_CAT = REFERENCE_CATALOG.split('Vizier.')[1]
            print('Align to Vizier catalog: http://vizier.u-strasbg.fr/viz-bin/VizieR?-source=%s' %(VIZIER_CAT))
            
            import astroquery
            if astroquery.__version__ < '0.0.dev1078':
                from astroquery import vizier

                query = {}
                query["-source"] = VIZIER_CAT
                #query["-out"] = ["_r", "CFHTLS", "rmag"]
                query["-out"] = ["_RAJ2000", "_DEJ2000"]  ### Just RA/Dec.

                #### Center position and query radius
                r0, d0 = wcs.wcs_pix2sky([[wcs.naxis1/2., wcs.naxis2/2.]], 1)[0]
                rll, dll = wcs.wcs_pix2sky([[0, 0]], 1)[0]
                corner_radius = np.sqrt((r0-rll)**2*np.cos(d0/360.*2*np.pi)**2+(d0-dll)**2)*60.*1.5
                h = query["-c"] = "%.6f %.6f" %(r0, d0)
                query["-c.rm"] = "%.3f" %(corner_radius)  ### xxx check image size

                #### Run the query
                vt = vizier.vizquery(query)
            else:
                #### Newer astroquery
                from astroquery.vizier import Vizier
                import astropy.coordinates as coord
                import astropy.units as u
                
                Vizier.ROW_LIMIT = -1
                
                r0, d0 = wcs.wcs_pix2sky([[wcs.naxis1/2., wcs.naxis2/2.]], 1)[0]
                rll, dll = wcs.wcs_pix2sky([[0, 0]], 1)[0]
                corner_radius = np.sqrt((r0-rll)**2*np.cos(d0/360.*2*np.pi)**2+(d0-dll)**2)*60.*1.5
                #
                c = coord.ICRSCoordinates(ra=r0, dec=d0, unit=(u.deg, u.deg))
                #### something with astropy.coordinates
                # c.icrs.ra.degree = c.icrs.ra.degrees
                # c.icrs.dec.degree = c.icrs.dec.degrees
                #
                vt = Vizier.query_region(c, width=u.Quantity(corner_radius, u.arcminute), catalog=[VIZIER_CAT])[0]
                
            #### Make a region file
            ra_list, dec_list = vt['RAJ2000'], vt['DEJ2000']
            print('Vizier, found %d objects.' %(len(ra_list)))
            fp = open('%s.vizier.reg' %(ROOT_DIRECT),'w')
            fp.write('# %s, r=%.1f\'\nfk5\n' %(VIZIER_CAT, corner_radius))
            for ra, dec in zip(ra_list, dec_list):
                fp.write('circle(%.6f, %.6f, 0.5")\n' %(ra, dec))
            #
            fp.close()
        else:
            #### Use DS9 catalog
            ds9 = threedhst.dq.myDS9()
            ds9.set('file SCI.fits')
            #ds9.set('file Q0821+3107-F140W_drz.fits')
            ds9.set('catalog %s' %(REFERENCE_CATALOG))
            ### Can't find XPA access point for "copy to regions"
            ds9.set('catalog export tsv ds9_align.tsv')
            lines = open('ds9_align.tsv').readlines()
            ra_list, dec_list = [], []
            for line in lines[1:]:
                spl = line.split()
                ra, dec = float(spl[0]), float(spl[1])
                ra_list.append(ra)
                dec_list.append(dec)
            #
            del(ds9)
            
        x_image, y_image = [], []
        for ra, dec in zip(ra_list, dec_list):
            x, y = wcs.wcs_sky2pix([[ra, dec]], 1)[0]
            x_image.append(x)
            y_image.append(y)
        
        alignCat = catIO.EmptyCat()
        alignCat['X_IMAGE'] = np.array(x_image)
        alignCat['Y_IMAGE'] = np.array(y_image)
        
    else:
        ## alignment image
        se.options['CATALOG_NAME']    = 'align.cat'
        status = se.sextractImage(ROOT_DIRECT+'_align.fits')
        alignCat = threedhst.sex.mySexCat('align.cat')
    
    xshift = 0
    yshift = 0
    rot = 0
    scale = 1.
    
    xrms = 2
    yrms = 2
    
    NITER = 5
    IT = 0
    while (IT < NITER):
        IT = IT+1
        
        #### Get x,y coordinates of detected objects
        ## direct image
        fp = open('direct.xy','w')
        for i in range(len(directCat.X_IMAGE)):
            fp.write('%s  %s\n' %(directCat.X_IMAGE[i],directCat.Y_IMAGE[i]))
        fp.close()

        ## alignment image
        fp = open('align.xy','w')
        for i in range(len(alignCat.X_IMAGE)):
            fp.write('%s  %s\n' %(np.float(alignCat.X_IMAGE[i])+xshift,
                       np.float(alignCat.Y_IMAGE[i])+yshift))
        fp.close()

        iraf.flpr()
        iraf.flpr()
        iraf.flpr()
        #### iraf.xyxymatch to find matches between the two catalogs
        pow = toler*1.
        try:
            os.remove('align.match')
        except:
            pass
            
        status1 = iraf.xyxymatch(input="direct.xy", reference="align.xy",
                       output="align.match",
                       tolerance=2**pow, separation=0, verbose=yes, Stdout=1)
        
        nmatch = 0
        while status1[-1].startswith('0') | (nmatch < 10):
            pow+=1
            os.remove('align.match')
            status1 = iraf.xyxymatch(input="direct.xy", reference="align.xy",
                           output="align.match",
                           tolerance=2**pow, separation=0, verbose=yes, Stdout=1)
            #
            nmatch = 0
            for line in open('align.match'): nmatch += 1
            
        if verbose:
            for line in status1:
                print(line)
        
                
        #### Compute shifts with iraf.geomap
        iraf.flpr()
        iraf.flpr()
        iraf.flpr()
        try:
            os.remove("align.map")
        except:
            pass
            
        status2 = iraf.geomap(input="align.match", database="align.map",
                    fitgeometry=fitgeometry, interactive=no, 
                    xmin=INDEF, xmax=INDEF, ymin=INDEF, ymax=INDEF,
                    maxiter = 10, reject = 2.0, Stdout=1)
        if verbose:
            for line in status2:
                print(line)
        
        #fp = open(root+'.iraf.log','a')
        #fp.writelines(status1)
        #fp.writelines(status2)
        #fp.close()
                
        #### Parse geomap.output 
        fp = open("align.map","r")
        for line in fp.readlines():
            spl = line.split()
            if spl[0].startswith('xshift'):
                xshift += float(spl[1])    
            if spl[0].startswith('yshift'):
                yshift += float(spl[1])    
            if spl[0].startswith('xrotation'):
                rot = float(spl[1])    
            if spl[0].startswith('xmag'):
                scale = float(spl[1])    
            if spl[0].startswith('xrms'):
                xrms = float(spl[1])    
            if spl[0].startswith('yrms'):
                yrms = float(spl[1])    
            
        fp.close()
        
        #os.system('wc align.match')
        print('Shift iteration #%d, xshift=%f, yshift=%f, rot=%f, scl=%f (rms: %5.2f,%5.2f)' %(IT, xshift, yshift, rot, scale, xrms, yrms))
    
    im = pyfits.open('SCI.fits')
        
    shutil.copy('align.map',ROOT_DIRECT+'_align.map')
    shutil.copy('align.match',ROOT_DIRECT+'_align.match')
    
    #### Cleanup
    if clean:
        rmfiles = ['SCI.fits','WHT.fits','align.cat',
               'align.map','align.match','align.reg','align.xy',
               'direct.cat','direct.reg','direct.xy',
               'drz_sci.fits','drz_wht.fits','bg.fits']
        
        for file in rmfiles:
            try:
                os.remove(file)
            except:
                pass
        
    return xshift, yshift, rot, scale, xrms, yrms

def match_diagnostic_plot(root='JKCS041-2r-168-F160W'):
    """
    Make delta_x delta_y scatter plot and vector diagram for outputs from the 
    alignment script. 
    """    
    import pywcs
    drz = pyfits.getheader(root+'_drz.fits', 'SCI')
    wcs = pywcs.WCS(drz)
     
    x_ref, y_ref, x_in, y_in, i_ref, i_in = np.loadtxt(root+'_align.match', unpack=True)
    
    dx = x_in - x_ref
    dy = y_in - y_ref
    
    plt.scatter(dx, dy, alpha=0.1)
    plt.xlim(-3,3); plt.ylim(-3,3)
    
    s = 200
    
    plt.quiver(x_in, y_in, dx*s, dy*s, scale=1, units='xy', alpha=0.5)
    plt.quiver(0.05*x_in.max(), 0.05*y_in.max(), s, 0, scale=1, units='xy', alpha=0.8, label='1 pix', color='red')
    plt.legend()
    
    
def matchImagePixels(input=None,matchImage=None,output=None,
                     input_extension=0, match_extension=0):
    """
matchImagePixels(input=None,matchImage=None,output=None,
                 match_extension=0)
    
    SWarp input image to same size/scale as matchIMage.
    Output default is input_root+'.match.fits'
    
    """
    from threedhst.sex import SWarp
    
    #input = '/research/HST/GRISM/WEINER/ACS/h_nz_sect33_v2.0_drz_img.fits'
    #matchImage = 'IB3721050_SCI.fits'
    
    if not input:
        return False
    
    if not output:
        output = os.path.basename(input).split('.fits')[0]+'.match.fits'
    
    #### initialize SWarp
    sw = SWarp()
    sw._aXeDefaults()
    sw.overwrite = True
    
    #### Get first guess coordinates
    sw.swarpMatchImage(matchImage,extension=match_extension)
    status = sw.swarpImage(matchImage+'[%d]' %match_extension,mode='wait')
    os.remove('coadd.fits')
    os.remove('coadd.weight.fits')
    
    #### Recenter
    sw.swarpRecenter()
    
    #### Make the final output image
    sw.options['IMAGEOUT_NAME'] = output
    #sw.options['WEIGHTOUT_NAME'] = base+'.match.weight.fits'
    
    for i,im in enumerate(input):
        input[i] += '[%d]' %input_extension
    status = sw.swarpImage(input,mode='direct')
    os.remove('coadd.weight.fits')

def make_blank_shiftfile(asn_file='ib3721050_asn.fits',
    xshift=0, yshift=0, rot=0., scale=1.0):
    """
make_blank_shiftfile(asn_file='ib3721050_asn.fits',
        xshift=0, yshift=0, rot=0., scale=1.0)
    
    Make a shiftfile with empty shifts and rotations for all exposures defined
    in an ASN file.
    """
    import threedhst
    
    asn = threedhst.utils.ASNFile(asn_file)
    root = asn_file.split('_asn.fits')[0]
    out = root+'_shifts.txt'
    
    fp = open(out,'w')
    fp.write("""# frame: output
# refimage: %s_flt.fits[1]
# form: delta 
# units: pixels\n""" %(asn.exposures[0]))
    
    for exp in asn.exposures:
        fp.write("%s_flt.fits %10.4f %10.4f %10.4f %10.4f\n" %(exp, xshift, yshift, rot, scale))
    
    fp.close()
    
def make_grism_shiftfile(asn_direct, asn_grism):
    """
make_grism_shiftfile(asn_direct, grism_direct)
    
    Make a shiftfile for grism exposures to match
    corresponding direct images
    """
    from threedhst.utils import ASNFile
    ROOT_DIRECT = asn_direct.split('_asn.fits')[0]#.lower()
    ROOT_GRISM  =  asn_grism.split('_asn.fits')[0]#.lower()
    #### Read shiftfile and ASN table
    sf = ShiftFile(ROOT_DIRECT+'_shifts.txt')
    asn = ASNFile(asn_grism)
    
    if sf.nrows == len(asn.exposures):
        #### Assume one direct image for each grism images, so just
        #### change the image names in the shiftfile to the grism exposures
        for i,exp in enumerate(asn.exposures):
            sf.images[i] = exp+'_flt.fits'
    else:
        #### Have different number of grism and direct images.  Just use the 
        #### shifts/rotations for the first direct image
        xs = sf.xshift[0]
        ys = sf.yshift[0]
        rot = sf.rotate[0]
        scl = sf.scale[0]
        sf.images = []
        sf.xshift = []
        sf.yshift = []
        sf.rotate = []
        sf.scale = []
        for i,exp in enumerate(asn.exposures):
            sf.images.append(exp+'_flt.fits')
            sf.xshift.append(xs)
            sf.yshift.append(ys)
            sf.rotate.append(rot)
            sf.scale.append(scl)
        
        sf.nrows = len(asn.exposures)
        
    #### Write the new shiftfile
    sf.write(ROOT_GRISM+'_shifts.txt')
    
    #print "\n3DHST.shifts.make_grism_shiftfile: %s_shifts.txt\n" %ROOT_GRISM
    threedhst.showMessage('Making grism shiftfile, %s_shifts.txt' %ROOT_GRISM)
    

def checkShiftfile(asn_direct):
    """
checkShiftfile(asn_direct)
    
    Make sure that there is a line in the shiftfile for each exposure 
    in the ASN table.  Also check that no scales are zero.
    """
    from threedhst.utils import ASNFile
    asn = ASNFile(asn_direct)
    
    sf_file = asn_direct.split('_asn.fits')[0]+'_shifts.txt'
    sf = ShiftFile(sf_file)
    flag=False
    for exp in asn.exposures:
        if exp+'_flt.fits' not in sf.images:
            flag=True
            print('Exposure, %s, not in %s' %(exp,sf_file))
            #print sf.nrows
            sf.append(exp+'_flt.fits')
            #print sf.nrows
    
    #### Check if scales are zero in the shiftfile
    if 0.0 in sf.scale:
        flag = True
        print('Found scale=0 in the shiftfile, setting to default no shift/rotation\n')
        for i in range(len(sf.scale)):
            sf.xshift[i], sf.yshift[i], sf.rotate[i], sf.scale[i] = 0.0, 0.0, 0.0, 1.0
        
    if flag:
        sf.write(sf_file)
    else:       
        #print "\n3DHST.shifts.checkShiftfile: %s looks OK.\n" %sf_file
        threedhst.showMessage('Shiftfile, %s, looks OK' %sf_file)

def default_rotation(asn_direct, path_to_flt='./'):
    """
    Force rotation to be 0.09 degrees to account for shift in the current 
    IDCTAB file (uab1537ci_idc.fits).
    
    http://www.stsci.edu/hst/wfc3/documents/newsletters/STAN_04_05_2011#section2
    
    """
    sf_file = asn_direct.split('_asn.fits')[0]+'_shifts.txt'
    sf = ShiftFile(sf_file)
    
    im = pyfits.open(threedhst.utils.find_fits_gz(sf.images[0]))
    if im[0].header['IDCTAB'].startswith('iref$uab1537'):    
        for i in range(sf.nrows):
            sf.rotate[i] = 0.09
        sf.write(sf_file)
    
class ShiftFile():
    """
ShiftFile(infile)
    
    Read and manipulate a shiftfile produced by, e.g., tweakshifts
    
    Based on aXe2html.sextractcat
    """
    def __init__(self, filename):
        linelist = self.opencat(filename)
        self.headerlines = self.extractheader(linelist)
        rowlines    = self.extractrows(linelist)
        self.processrows(rowlines)
        self.nrows = len(rowlines)
        self.filename = filename
        
    def opencat(self, filename):
        """
        Input:
            filename - the name of the sextractor ascii catalog
            
        Return:
                linelist - a list wit all lines of the sextractor
                       catalog
                       
        Description:
            Opens the file and stores all the lines in a list.

        """
        listfile = open(filename,'r')
        linelist = listfile.readlines()
        listfile.close()
        
        return linelist
    
    
    def extractheader(self, linelist):
        """
    Input:
        linelist - all lines of a sextractor catalog

    Return:
        headerlines - the lines which contain header
                      information

    Description:
        Extracts the header lines from the list of lines
        Header lines have a '#' as the first digits and are
        longer than two characters. This allows for small
        cosmetical changes to the original Sextractor tables
        """
        headerlines = []
        for index in range(len(linelist)):
            oneline = linelist[index]
            if (oneline[0] == "#") and (len(oneline)) > 2:
                headerlines.append(oneline)
        return headerlines
    
    def extractrows(self, linelist):
        """
    Input:
        linelist - all lines of a sextractor catalog
        
    Return:
        rowlines - the content lines of a sextractor
                   catalog

    Description:
        Extracts the content lines from the list of lines.
        Content lines have a all lines except the ones which start with a 
        '#' and are longer than one character. This allows for
        e.g. whitespaces as the last 'rows', which often remain 
        after editing an ascii-file
        
        """
        rowlines = []
        for index in range(len(linelist)):
            oneline = linelist[index]
            if oneline[0] != "#" and len(oneline) > 1:
                rowlines.append(oneline)        
        return rowlines
    
    
    def processrows(self, rowlines):
        """
        processrows(self, rowlines)
        
        Read: image xshift yshift rotate scale from shiftfile.
            
        """
        nlines = len(rowlines)
        self.images = []
        self.xshift = []
        self.yshift = []
        self.rotate = []
        self.scale = []
        for i in range(nlines):
            line = rowlines[i].split()
            self.images.append(line[0])
            self.xshift.append(float(line[1]))
            self.yshift.append(float(line[2]))
            if len(line) > 3:
                self.rotate.append(float(line[3]))
            else:
                self.rotate.append(0.)
                self.scale.append(1.0)
            if len(line) > 4:
                self.scale.append(float(line[4]))
    
    
    def write(self, outfile):
        """write(outfile)"""
        fp = open(outfile,'w')
        fp.writelines(self.headerlines)
        for i in range(self.nrows):
            line = '%-20s %8.4f %8.4f %8.3f %8.3f\n' %(self.images[i],
                self.xshift[i], self.yshift[i], self.rotate[i], self.scale[i])
            fp.write(line)
        fp.close()
    
    def pop(self,idx):
        """
pop(self,idx)
        """
        out = self.images.pop(idx) 
        out = self.xshift.pop(idx) 
        out = self.yshift.pop(idx) 
        out = self.rotate.pop(idx) 
        out = self.scale.pop(idx) 
        self.nrows -= 1
    
    def append(self,image, xshift=0., yshift=0.,
                    rotate=0.0, scale=1.0):
        """
append(self,image, xshift=0., yxhift=0.,
                   rotate=0.0, scale=1.0)
        """
        self.images.extend([image])
        self.xshift.extend([xshift])
        self.yshift.extend([yshift])
        self.rotate.extend([rotate])
        self.scale.extend([scale])
        self.nrows += 1