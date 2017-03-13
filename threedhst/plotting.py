"""
3DHST.plotting

Utilities for plotting grism spectra.

"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import numpy as np
errs = np.seterr(divide='ignore', invalid='ignore')
errs = np.seterr(divide='ignore', invalid='ignore')

import matplotlib
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import pylab

USE_PLOT_GUI=False

from matplotlib.figure import Figure
from matplotlib.backends.backend_agg import FigureCanvasAgg

import threedhst

def defaultPlotParameters():
    """
    defaultPlotParameters()
    """
    #plt.rcParams['font.family'] = 'serif'
    #plt.rcParams['font.serif'] = ['Times']
    plt.rcParams['ps.useafm'] = True
    plt.rcParams['patch.linewidth'] = 0.
    plt.rcParams['patch.edgecolor'] = 'black'
    plt.rcParams['text.usetex'] = False
    plt.rcParams['image.cmap'] = 'gray'
    #plt.rcParams['text.usetex'] = True
    #plt.rcParams['text.latex.preamble'] = ''

def redo_1dall():
    
    import glob
    
    for dir in ['COSMOS','GOODS-N','AEGIS','SN-PRIMO','SN-GEORGE']:
        os.chdir('/research/HST/GRISM/3DHST/'+dir+'/DATA/')
        files=glob.glob('*G141_asn.fits')
        for file in files:
            threedhst.plotting.redo_1dspec(file.split('_asn')[0])
        
def redo_1dspec(ROOT_GRISM):
    """
    Run the 1D grism plots again after the plotting routine changed
    """
    import os
    sexCat = threedhst.sex.mySexCat(ROOT_GRISM+'_drz.cat')
    
    #### Read the SPC file containing the 1D extractions
    SPC = threedhst.plotting.SPCFile(ROOT_GRISM+'_2_opt.SPC.fits',
                    axe_drizzle_dir=os.environ['AXE_DRIZZLE_PATH'])
    
    threedhst.plotting.asciiSpec(SPC,root=ROOT_GRISM,path='../HTML/ascii')
    print('')
    threedhst.plotting.makeSpec1dImages(SPC, path='../HTML/images/')

def make_data_products(ROOT_DIRECT, ROOT_GRISM):
    """
make_data_products()
    
    Make the ascii spectra, thumbnails and HTML files.
    """
    import os
    import shutil
    import tarfile
    import glob
    
    threedhst.showMessage('Making output data products (1D/2D thumbnails,\n'+
        'webpages, etc.)')
    
    #### Check directory structure
    if os.path.exists('../HTML/scripts') is False:
        os.mkdir('../HTML/scripts')
        print("""
        WARNING: ../HTML/scripts/ not found.  Download from
                 http://code.google.com/p/threedhst/
              """)
    
    if os.path.exists('../HTML/images') is False:
        os.mkdir('../HTML/images')
        
    #### Read the SExtractor catalog
    sexCat = threedhst.sex.mySexCat(ROOT_GRISM+'_drz.cat')
    
    #### Read the SPC file containing the 1D extractions
    SPC = threedhst.plotting.SPCFile(ROOT_GRISM+'_2_opt.SPC.fits',
                    axe_drizzle_dir=os.environ['AXE_DRIZZLE_PATH'])
    threedhst.currentRun['SPC'] = SPC
    
    #############################################
    #### Direct image thumbnails
    #############################################
    print('\nTHREEDHST.plotting.makeThumbs: Creating direct image ' + \
          'thumbnails...\n\n')
    threedhst.plotting.makeThumbs(SPC, sexCat, path='../HTML/images/')
    
    # 
    # fptar = tarfile.open('../HTML/images/'+ROOT_GRISM+'_thumbs.tar.gz','w|gz')
    # oldwd = os.getcwd()
    # os.chdir('../HTML/images/')
    # files = glob.glob(ROOT_GRISM+'*thumb.fits.gz')
    # for file in files:
    #     fptar.add(file)
    # fptar.close()
    # os.chdir(oldwd)
        
    #############################################
    #### 2D spectra images showing the Model, contamination, and G141 spectra
    #############################################
    print('\nmakeSpec1dImages: Creating 2D spectra '+ \
          'thumbnails...\n\n')
    threedhst.plotting.makeSpec2dImages(SPC, path='../HTML/images/')

    #### Make ASCII spectra from the SPC file
    print('\n Making ASCII spectra in ../HTML/ascii/\n')
    threedhst.plotting.asciiSpec(SPC,root=ROOT_GRISM,path='../HTML/ascii')

    #############################################
    #### 1D spectra images
    #############################################
    print('\nTHREEDHST.plotting.makeSpec1dImages: Creating 1D spectra '+ \
          'thumbnails...\n\n')
    threedhst.plotting.makeSpec1dImages(SPC, path='../HTML/images/')
    
    # fptar = tarfile.open('../HTML/images/'+ROOT_GRISM+'_2D.tar.gz','w|gz')
    # oldwd = os.getcwd()
    # os.chdir('../HTML/images/')
    # files = glob.glob(ROOT_GRISM+'*2D.fits.gz')
    # for file in files:
    #     fptar.add(file)
    # fptar.close()
    # os.chdir(oldwd)
    
    #### Done with thumbnails
    threedhst.currentRun['step'] = 'MAKE_THUMBNAILS'
    
    #############################################
    #### Make tiles for Google map layout   
    #############################################
    print('\nMaking GMap tiles in ./HTML\n\n')
    
    try:
        os.mkdir('../HTML/tiles')
    except:
        pass
    
    mapParams = threedhst.gmap.makeAllTiles(ROOT_DIRECT, ROOT_GRISM)
                
    #### Done making the map tiles
    threedhst.currentRun['step'] = 'MAKE_GMAP_TILES'
    
    #############################################
    #### Copy catalog to HTML directory
    #############################################
    status = os.system('cp '+ROOT_GRISM+'_drz.cat ../HTML/')
    
    mapParams['ZOOMLEVEL']=15
    #### Make the full HTML file
    out_web = '../HTML/'+ROOT_GRISM+'_index.html'
    print('\nthreedhst.plotting.makeHTML: making webpage: %s\n' %out_web)
    threedhst.plotting.makeHTML(SPC, sexCat, mapParams, output=out_web)
    threedhst.plotting.makeCSS()
    threedhst.plotting.makeJavascript()
    
    threedhst.currentRun['step'] = 'MAKE_HTML'
        
    
def plotThumb(object_number, mySexCat, in_image = None, size = 20, scale=0.128, 
              outfile='/tmp/thumb.png', close_window=False):
    """
    plotThumb(object_number, mySexCat, in_image = None, size = 20, scale=0.128,
              outfile='/tmp/thumb.png', close_window=False)
    """
    
    if in_image is not None:
        im = pyfits.open(mySexCat.filename.split('.cat')[0]+'.fits')
        in_image = im[1].data
    
    obj_str = str(object_number)
    idx = -1
    for i, number in enumerate(mySexCat.NUMBER):
        if number == obj_str:
            idx = i
            break
    
    if idx < 0:
        print('Object \'%s\' not found in SExtractor catalog, %s.\n' %(obj_str,
                             mySexCat.filename))
        return False
    
    ##### Figure out why X,Y are swapped in mySexCat.
    ##### Maybe the orientation of in_image is rotated w.r.t catalog?
    x0 = np.round(np.float(mySexCat.X_IMAGE[idx]))
    y0 = np.round(np.float(mySexCat.Y_IMAGE[idx]))
    
    sub = in_image[y0-size:y0+size, x0-size:x0+size]
    
    if sub.shape != (size, size):
        sub = np.zeros([size, size])
       
    interp = 'nearest'
    asp = 'auto'
    sub_max = np.max(sub)
    if sub_max > 0:
        #### Minmax scaling
        vmin = -1.1*sub_max
        vmax = 0.1*sub_max
    else:
        vmin = -0.5*0.8
        vmax = 0.1*0.8
    
    flux = np.round(np.float(mySexCat.FLUX_AUTO[idx]))
    vmin = -0.03*flux
    vmax = 0.003*flux
    
    defaultPlotParameters()
    
    
    fig = plt.figure(figsize=[3,3],dpi=100)
    fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.02,
                        bottom=0.02,right=0.98,top=0.98)
    ax = fig.add_subplot(111)
    ax.imshow(0-sub, interpolation=interp,aspect=asp,vmin=vmin,vmax=vmax)    
    
    asec_pix = 1./scale
    nasec = int(size/asec_pix)
    
    ax.set_yticklabels([])
    xtick = ax.set_xticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size)
    ax.set_xticklabels([])
    ytick = ax.set_yticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size)
    
    ### Save to PNG
    if outfile:
        fig.savefig(outfile,dpi=100,transparent=False)
    
    if close_window:
        status = plt.close()
    
    #plt.show()

def plotThumbNew(object_number, mySexCat, SPCFile,
                 outfile='/tmp/thumb.png', close_window=False):
    """
    plotThumb(object_number, mySexCat, in_image = None, size = 20, scale=0.128,
              outfile='/tmp/thumb.png', close_window=False)
    """
    import os
    import glob
    
    from pyraf import iraf
    from iraf import dither
    
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    
    obj_str = str(object_number)
    idx = -1
    for i, number in enumerate(mySexCat.NUMBER):
        if number == obj_str:
            idx = i
            break
    
    if idx < 0:
        print('Object \'%s\' not found in SExtractor catalog, %s.\n' %(obj_str,
                             mySexCat.filename))
        return False
    
    #### Match thumbnail size to size of 2D spectrum
    mef = pyfits.open('../'+threedhst.options['DRIZZLE_PATH']+
                      '/'+root+'_mef_ID'+str(object_number)+'.fits')
    head = mef['SCI'].header
    size = head['NAXIS2']
    
    #drz_image = mySexCat.filename.split('.cat')[0]+'.fits'
    drz_image = threedhst.options['DIRECT_MOSAIC']
    
    drz = pyfits.open(drz_image)
    drz_header = drz['SCI'].header
    drz_y, drz_x = drz['SCI'].data.shape
    
    orient = drz_header['PA_APER']
    
    #### Need to get orientation from *GRISM* images if a 
    #### "prefab" direct image was supplied.
    if threedhst.options['PREFAB_DIRECT_IMAGE'] is not None:
        drz_grism = threedhst.options['ROOT_GRISM']+'_drz.fits'
        grism_header = pyfits.getheader(drz_grism,'SCI')
        orient = grism_header['PA_APER']
    
    pixel_scale = np.abs(drz_header['CD1_1']*3600.)
    ra_ref = mySexCat.X_WORLD[idx]
    dec_ref = mySexCat.Y_WORLD[idx]

    xpix = np.round(np.float(mySexCat.X_IMAGE[idx]))
    ypix = np.round(np.float(mySexCat.Y_IMAGE[idx]))
    
    data_img = drz_image+'[1]'
    mask_img = drz_image+'[2]'
    
    #### Imcopy a subimage, necessary for very large input DRZ files
    data_img = 'subSCI.fits'
    mask_img = 'subWHT.fits'
    
    old_files = glob.glob('sub[SW]*.fits')
    for old_file in old_files:
        #print old_file
        try:
            os.remove(old_file)
        except:
            pass
            
    # print xpix, ypix, 3*size, '\n\n'
    
    #return
    
    # print 'pix\n\n',drz_x, drz_y, xpix+3*size, ypix+3*size, 3*size
    # print (np.max([xpix-3*size,1]),
    #               np.min([xpix+3*size, drz_x]),
    #               np.max([ypix-3*size,1]),
    #               np.min([ypix+3*size, drz_y]))
    # print '\n\ntest'
    
    threedhst.process_grism.flprMulti()
       
    iraf.imcopy(drz_image+'[SCI][%d:%d,%d:%d]' %(np.max([xpix-3*size,1]),
                  np.min([xpix+3*size, drz_x-5]),
                  np.max([ypix-3*size,1]),
                  np.min([ypix+3*size, drz_y-5])), 'subSCIa.fits',
                  verbose=iraf.no)
    # iraf.imcopy(drz_image+'[WHT][%d:%d,%d:%d]' %(xpix-3*size, xpix+3*size, 
    #               ypix-3*size, ypix+3*size), 'subWHT.fits',
    #               verbose=iraf.no)
    
    #### Imcopy misses a few WCS keywords (that are zero)
    iraf.hedit('subSCIa.fits','CD1_2',0, add=iraf.yes, update=iraf.yes, verify=iraf.no, show=iraf.no)
    iraf.hedit('subSCIa.fits','CD2_1',0, add=iraf.yes, update=iraf.yes, verify=iraf.no, show=iraf.no)
    # iraf.hedit('subWHT.fits','CD1_2',0, add=iraf.yes, update=iraf.yes, verify=iraf.no, show=iraf.no)
    # iraf.hedit('subWHT.fits','CD2_1',0, add=iraf.yes, update=iraf.yes, verify=iraf.no, show=iraf.no)
    # 
    # iraf.wcscopy(images='tmpWHT.fits',refimages='subSCI.fits',
    #      verbose=iraf.no)
    
    sci = pyfits.open('subSCIa.fits')
    sci[0].data = sci[0].data*0+1
    sci.writeto('subWHTa.fits', clobber=True)

    fitsfile = (os.path.dirname(outfile)+'/'+
                os.path.basename(outfile).split('.png')[0]+'.fits')
    
    ### NEED TO STRIP FITS HEADER
    im = pyfits.open('subSCIa.fits')
    sci = im[0].data
                
    s_hdu = pyfits.PrimaryHDU(sci)
    s_list = pyfits.HDUList([s_hdu])
    s_list[0].header.update('EXPTIME',im[0].header.get('EXPTIME'))
    s_list[0].header.update('CDELT1',im[0].header.get('CD1_1'))
    s_list[0].header.update('CDELT2',im[0].header.get('CD2_2'))
    copy_keys = ['CTYPE1','CTYPE2','CRVAL1','CRVAL2','CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2','LTM1_1','LTM2_2']
    for key in copy_keys:
        try:
            s_list[0].header.update(key, im[0].header.get(key))
        except:
            s_list[0].header.update(key, 0)
            
    s_list.writeto('subSCI.fits', clobber=True)
    
    im = pyfits.open('subWHTa.fits')
    sci = im[0].data
    w_hdu = pyfits.PrimaryHDU(sci)
    w_list = pyfits.HDUList([w_hdu])
    copy_keys = ['CTYPE1','CTYPE2','CRVAL1','CRVAL2','CRPIX1','CRPIX2','CD1_1','CD1_2','CD2_1','CD2_2','LTM1_1','LTM2_2']
    w_list[0].header.update('EXPTIME',im[0].header.get('EXPTIME'))
    w_list[0].header.update('CDELT1',im[0].header.get('CD1_1'))
    w_list[0].header.update('CDELT2',im[0].header.get('CD2_2'))
    for key in copy_keys:
        try:
            w_list[0].header.update(key, im[0].header.get(key))
        except:
            s_list[0].header.update(key, 0)
    
    w_list.writeto('subWHT.fits', clobber=True)
    
    old_files = glob.glob(fitsfile+'*')
    for old_file in old_files:
        #print old_file
        os.remove(old_file)
    
    # print fitsfile, data_img, mask_img, size, ra_ref, dec_ref, pixel_scale, orient
    
    threedhst.process_grism.flprMulti()
    
    wdrizzle_fail=False
    try:
      status = iraf.wdrizzle(data = data_img, outdata = fitsfile, \
       outweig = "", outcont = "", in_mask = mask_img, 
       wt_scl = 'exptime', \
       outnx = size, outny = size, geomode = 'wcs', kernel = 'square', \
       pixfrac = 1.0, coeffs = "", lamb = 1392., xgeoim = "", ygeoim = "", \
       align = 'center', scale = 1.0, xsh = 0.0, ysh = 0.0, rot = 0.0, \
       shft_un = 'input', shft_fr = 'input', outscl = pixel_scale, \
       raref = ra_ref, decref = dec_ref, xrefpix = size/2+0.5, yrefpix = size/2+0.5, \
       orient = orient, dr2gpar = "", expkey = 'exptime', in_un = 'cps', \
       out_un = 'cps', fillval = '0', mode = 'al', Stdout=1)
      
      wdrizzle_fail = not status[-1].startswith('-Writing output')
      
    except:
        wdrizzle_fail=True
        
    #### Sometimes imcopy/wdrizzle breaks.  Run from larger DRZ file
    #### directly if it does
    if wdrizzle_fail:
        print('Redo WDRIZZLE\n\n')
        data_img = drz_image+'[1]'
        mask_img = drz_image+'[2]'
        try:
          status = iraf.wdrizzle(data = data_img, outdata = fitsfile, \
           outweig = "", outcont = "", in_mask = mask_img, 
           wt_scl = 'exptime', \
           outnx = size, outny = size, geomode = 'wcs', kernel = 'square', \
           pixfrac = 1.0, coeffs = "", lamb = 1392., xgeoim = "", ygeoim = "", \
           align = 'center', scale = 1.0, xsh = 0.0, ysh = 0.0, rot = 0.0, \
           shft_un = 'input', shft_fr = 'input', outscl = pixel_scale, \
           raref = ra_ref, decref = dec_ref, xrefpix = size/2+0.5, yrefpix = size/2+0.5, \
           orient = orient, dr2gpar = "", expkey = 'exptime', in_un = 'cps', \
           out_un = 'cps', fillval = '0', mode = 'al', Stdout=1)
        except:
            return 
            
        if not status[-1].startswith('-Writing output'):
            return 
            
    #print status[-1]+'\n'
    
    ####  !!!!!!!!!!! Should do the same here for the segmentation image, but
    ####  rotating the segmentation image doesn't work very well for fractional
    ####  pixels.  Can do kernel='point', but then get holes in the image where
    ####  pixels don't line up.
    sub_im = pyfits.open(fitsfile)
    sub = sub_im[0].data
    
    #### Look for a better platform-independent way to gzip
    if os.path.exists(fitsfile+'.gz'):
        os.remove(fitsfile+'.gz')
    os.system('gzip '+fitsfile)
    
    interp = 'nearest'
    asp = 'auto'
    sub_max = np.max(sub)
    if sub_max > 0:
        #### Minmax scaling
        vmin = -1.1*sub_max
        vmax = 0.1*sub_max
        vmin = -0.8*sub_max
        vmax = 0.08*sub_max
    else:
        vmin = -0.5*0.8
        vmax = 0.1*0.8
    
    # vmin = -0.5
    # vmax = 0.1
    
    #flux = np.round(np.float(mySexCat.FLUX_AUTO[idx]))
    #vmin = -0.03*flux
    #vmax = 0.003*flux
    
    defaultPlotParameters()
    
    if USE_PLOT_GUI:
        fig = plt.figure(figsize=[3,3],dpi=100)
    else:
        fig = Figure(figsize=[3,3], dpi=100)
        
    fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.02,
                        bottom=0.02,right=0.98,top=0.98)
    
    ax = fig.add_subplot(111)
        
    ax.imshow(0-sub, interpolation=interp,aspect=asp,vmin=vmin,vmax=vmax)    
    
    asec_pix = 1./pixel_scale
    nasec = int(size/asec_pix/2)
    
    #print nasec
    
    ax.set_yticklabels([])
    xtick = ax.set_xticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size/2)
    ax.set_xticklabels([])
    ytick = ax.set_yticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size/2)
    
    ### Save to PNG
    if outfile:
        if USE_PLOT_GUI:
            fig.savefig(outfile,dpi=100,transparent=False)
        else:
            canvas = FigureCanvasAgg(fig)
            canvas.print_figure(outfile, dpi=100, transparent=False)
    
    if close_window & (USE_PLOT_GUI):
        status = plt.close()
    
    #plt.show()
    
def makeThumbs(SPCFile, mySexCat, path='./HTML/'):
    """
    makeThumbs(SPCFile, mySexCat,path='./HTML')
    
    Run plotThumbs for each object in SPCFile
    """
    import os
    import tarfile
    
    # im = pyfits.open(mySexCat.filename.split('.cat')[0]+'.fits')
    im = pyfits.open(threedhst.options['DIRECT_MOSAIC'])
    
    dat = im[1].data
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    ids = SPCFile._ext_map+0
    ids.sort()
    for id in ids:
        idstr = '%05d' %id
        print(threedhst.noNewLine+'plotting.makeThumbs: %s_%s_thumb.png' %(root, idstr))
        plotThumbNew(id, mySexCat, SPCFile,
                  outfile=path+'/'+root+'_'+idstr+'_thumb.png',
                  close_window=True)
        # plotThumb(id, mySexCat, SPCFile,
        #           outfile=path+'/'+root+'_'+idstr+'_thumb.png',
        #           close_window=True)
        
def plot2Dspec(SPCFile, object_number, outfile='/tmp/spec2D.png',
               close_window=False, clean=True):
    """
plot2Dspec(SPCFile, object_number, outfile='/tmp/spec2D.png', 
    close_window=False, clean=True)
    
    Make plot of 2D spectrum comparing observed, model and contamination 
    spectra.  If clean=True, then plot the cleaned-spectrum rather than the 
    contamination in the bottom panel.
    """
    import os
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    
    mef = pyfits.open('../'+threedhst.options['DRIZZLE_PATH']+
                      '/'+root+'_mef_ID'+str(object_number)+'.fits')
    
    head = mef['SCI'].header
    lmin = 10800
    lmax = 16800
    
    if threedhst.options['ACS_G800L']:
        lmin = 5500
        lmax = 1.0e4
    
    if threedhst.options['GRISM_NAME'] == 'G102':
        lmin = 8000
        lmax = 1.12e4
    
    xmin = (lmin-head['CRVAL1'])/head['CDELT1']+head['CRPIX1']
    xmax = (lmax-head['CRVAL1'])/head['CDELT1']+head['CRPIX1']
    
    defaultPlotParameters()
    
    if USE_PLOT_GUI:
        fig = plt.figure(figsize=[6,4],dpi=100)
    else:
        fig = Figure(figsize=[6,4], dpi=100)
    
    fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.06,
                        bottom=0.12,right=0.99,top=0.99)
        
    interp = 'nearest'
    asp = 'auto'
    vmin = -0.6 
    vmax = 0.1
        
    mod_max = np.max(mef['MOD'].data)
    if mod_max > 0:
        vmin = -1*mod_max
        vmax = 0.1*mod_max
    else:
        vmin = -0.5*0.8
        vmax = 0.1*0.8
    
    #### Need to convert to e- / s for ACS.  aXe bug for ACS?????
    if threedhst.options['GRISM_NAME'] == 'G800L':
        drz = pyfits.getheader(threedhst.options['ROOT_GRISM']+'_drz.fits',0)
        mef['SCI'].data /= (drz.get('EXPTIME')/drz.get('NDRIZIM')*2)
        
    #vmin *= (np.float(threedhst.options['DRZSCALE'])/0.128254)**1
    #vmax *= (np.float(threedhst.options['DRZSCALE'])/0.128254)**1
    
    ax = fig.add_subplot(311)
    ax.imshow(0-mef['SCI'].data, interpolation=interp,aspect=asp,
              vmin=vmin,vmax=vmax)    
    ax.set_xlim(xmin,xmax)
    ax.set_yticklabels([])
    ax.set_xticks((np.arange(np.ceil(lmin/1000.)*1000,
       np.ceil(lmax/1000.)*1000,1000)-head['CRVAL1'])/head['CDELT1']
       +head['CRPIX1'])
    ax.set_xticklabels([])
    ax.set_ylabel(threedhst.options['GRISM_NAME'])
    
    ax = fig.add_subplot(312)
    ax.imshow(0-mef['MOD'].data, interpolation=interp, aspect=asp,
              vmin=vmin,vmax=vmax)
    ax.set_xlim(xmin,xmax)
    ax.set_yticklabels([])
    ax.set_xticks((np.arange(np.ceil(lmin/1000.)*1000,
        np.ceil(lmax/1000.)*1000,1000)-head['CRVAL1'])/head['CDELT1']
        +head['CRPIX1'])
    ax.set_xticklabels([])
    ax.set_ylabel('Model')

    ax = fig.add_subplot(313)
    if clean:
        ax.imshow(0-(mef['SCI'].data-mef['CON'].data), interpolation=interp,    
              aspect=asp, vmin=vmin,vmax=vmax)
    else:
        ax.imshow(0-mef['CON'].data, interpolation=interp, aspect=asp,
               vmin=vmin,vmax=vmax)

    ax.set_xlim(xmin,xmax)
    ax.set_yticklabels([])
    ax.set_xticks((np.arange(np.ceil(lmin/1000.)*1000,
        np.ceil(lmax/1000.)*1000,1000)-head['CRVAL1'])/head['CDELT1']
        +head['CRPIX1'])
    ax.set_xticklabels(np.arange(np.ceil(lmin/1000.)*1000,
        np.ceil(lmax/1000.)*1000,1000)/1.e4)
    
    if clean:
        ax.set_ylabel('Cleaned')
    else:
        ax.set_ylabel('Contam.')
    
    ax.set_xlabel(r'$\lambda$ [$\mu$m]')
    
    #ax.set_xticklabels([])
    #plt.show()
    
    ### Save to PNG
    if outfile:
        if USE_PLOT_GUI:
            fig.savefig(outfile,dpi=100,transparent=False)
        else:
            canvas = FigureCanvasAgg(fig)
            canvas.print_figure(outfile, dpi=100, transparent=False)
            
    if close_window & USE_PLOT_GUI:
        status = plt.close()

def makeSpec2dImages(SPCFile, path='./HTML/', add_FITS=True):
    """
    makeSpec2dImages(SPCFile, path='./HTML')
    
    Run plotObject for each object in SPCFile
    """
    import os
    import shutil
    
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    ids = SPCFile._ext_map+0
    ids.sort()
    for id in ids:
        idstr = '%05d' %id
        print(threedhst.noNewLine+'plotting.makeSpecImages: %s_%s_2D.png' %(root, idstr))
        plot2Dspec(SPCFile, id, outfile=path+'/'+root+'_'+idstr+'_2D.png',
                   close_window=True)
        
        if add_FITS:
            mef_file = '../'+threedhst.options['DRIZZLE_PATH'] + \
                       '/'+root+'_mef_ID'+str(id)+'.fits'
            out_file = path+'/'+root+'_'+idstr+'_2D.fits'
            shutil.copy(mef_file,out_file)
            ### Gzip the result
            if os.path.exists(out_file+'.gz'):
                os.remove(out_file+'.gz')
            os.system('gzip '+out_file)
        
def plot1Dspec(SPCFile, object_number, outfile='/tmp/spec.png',
               close_window=False, show_test_lines=False, own_extraction=True):
    """
    plot1Dspec(SPCFile, object_number, outfile='/tmp/spec.png', 
               close_window=False, show_test_lines=False)
    """
    import os
    import scipy.optimize
    #from scipy.interpolate import interp1d
    
    #import threedhst.plotting as pl
    #reload pl
    #self = pl.SPCFile('ib3721050_2_opt.SPC.fits')
    
    defaultPlotParameters()
    #object_number = 42
    spec = SPCFile.getSpec(object_number)
    lam  = spec.field('LAMBDA')
    flux = spec.field('FLUX')
    ferr = spec.field('FERROR') 
    sflux = np.sum(flux[(lam > 1.3e4) & (lam < 1.5e4)])
    if sflux == 0:
        sflux = 1
    
    contam = spec.field('CONTAM')
    ferr_axe = ferr*1.
    error_color=(0,1,0)
    
    #### Own extraction
    if own_extraction:
        sp1d = threedhst.spec1d.extract1D(object_number, root=SPCFile.filename.split('_2_opt')[0], path='../HTML', show=False, out2d=False)
        lam = sp1d['lam']
        flux = sp1d['flux']
        
        #### Scale it to roughly match the normalization of the aXe extraction
        scl = np.sum(flux[(lam > 1.3e4) & (lam < 1.5e4)])/sflux
        if scl == 0:
            scl = 1
        
        flux /= scl
        ferr = sp1d['error']/scl
        contam = sp1d['contam']/scl
        error_color=(0,0,1)
        
    #### Need to convert to e- / s for ACS.  aXe bug for ACS?????
    if threedhst.options['GRISM_NAME'] == 'G800L':
        drz = pyfits.getheader(threedhst.options['ROOT_GRISM']+'_drz.fits',0)
        flux /= (drz.get('EXPTIME')/drz.get('NDRIZIM')*2)
        ferr /= (drz.get('EXPTIME')/drz.get('NDRIZIM')*2)
    
    ###
    xmin = 10800
    xmax = 16800
    
    if threedhst.options['ACS_G800L']:
        xmin = 5500
        xmax = 1.0e4
    
    if threedhst.options['GRISM_NAME'] == 'G102':
        xmin = 8000
        xmax = 1.12e4
        
    sub = np.where((lam > xmin) & (lam < xmax))[0]
    ymax = np.max((flux-0*contam)[sub])
    
    ### Initialize plot
    if USE_PLOT_GUI:
        fig = plt.figure(figsize=[5,3.4],dpi=100)
    else:
        fig = Figure(figsize=[5,3.4], dpi=100)
    
    fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.10,
                        bottom=0.15,right=0.99,top=0.90)

    ### plot window
    ax = fig.add_subplot(111)
        
    ### Plot Flux and Contamination
    ax.plot(lam, flux-contam, linewidth=1.0, color='blue',alpha=0.8)
    ax.plot(lam, contam, linewidth=1.0, color='red',alpha=1)
    ax.plot(lam, flux, color='red',linewidth=1.0,alpha=0.2)
    # ax.errorbar(lam, flux-contam,
    #            yerr= ferr,ecolor='blue', ealpha=0.3,
    #            color='blue',fmt='.',alpha=0.2) #, linestyle='none')
    
    lok = (lam > xmin) & (lam < xmax)
    ax.fill_between(lam[lok], (flux-contam+ferr)[lok], (flux-contam-ferr)[lok], color=error_color, alpha=0.3)
    #print ferr_axe.shape, ferr.shape
    
    # ax.fill_between(lam[lok], (flux-contam+ferr_axe)[lok], (flux-contam-ferr_axe)[lok], color=(0,0,1), alpha=0.1)
        
    #### Search for lines
    lines = threedhst.spec1d.findLines(SPCFile, idx=object_number)
    out_lines = []
    if lines:
        sdss_lines = threedhst.spec1d.readLinelist()
        #sdss_lines.gal_weight = sdss_lines.gal_weight+1
        sdss_lines.gal_weight /= np.max(sdss_lines.gal_weight)
        sdss_lines.qso_weight /= np.max(sdss_lines.qso_weight)
        colors = ['green','orange','blue']
        weight = sdss_lines.gal_weight
        #weight = sdss_lines.qso_weight
        sdss_use = np.where(np.array(weight) > -10)[0]
                
        iline=0
        for il, line in enumerate(lines):
            if (line.flag == 'ok') & (line.type.startswith('em')):

                ax.plot(line.wave*np.array([1,1]),np.array([-1,1]),
                        color='black',linewidth=2,alpha=0.2)
                ax.plot(line.wave*np.array([1,1]),np.array([-1,1]),'--',
                        color='orange',linewidth=2,alpha=0.7)
                
                ### Fit a gaussian to the line position
                near = np.where(np.abs(lam-line.wave) < 800)[0]                                         
                if 1 == 1:                    
                #if len(near) > 2:                    
                    p0 = np.array([np.max(flux[near]), line.wave,
                         50, np.median(flux[sub]), 0.])
                    #fi = interp1d(lam[near], (flux-contam)[near], kind='cubic')       
                    xnew = np.linspace(lam[near[0]], lam[near[-1]], len(near)*3)
                    fi = np.interp(xnew, lam[near], (flux-contam)[near])
                    ok = np.isfinite(fi)
                    if len(xnew[ok]) > 10:
                        pout = threedhst.plotting.gaussfit(xnew[ok], 
                                                           fi[ok], p0)
                    else:
                        pout = -1
                        
                    #print p0, pout[2][0]
                    #print pout[2][4]
                    
                    #### If gaussian fit successful, update line parameters
                    #### with fit results
                    if pout is not -1:
                        if pout[2][4] <= 4:
                            line.wave = pout[2][0][1]
                            line.sigma = pout[2][0][2]
                            line.ew = pout[1]
                            
                            # S/N at line peak
                            npeak = np.where(np.abs(lam-line.wave) ==
                                         np.min(np.abs(lam-line.wave)))[0]
                        
                            line.sn = pout[2][0][0]/ferr[npeak]
                            line.type = 'emgauss'
                        
                            #lines[il] = line

                            ax.plot(xnew,pout[0], color='green', linewidth=4, alpha=0.3)
                    
                ### show assuming line is OII, OIII, Ha
                compare_lines = np.array([3727., 5007, 6563.])
                if show_test_lines: 
                  for iz, z_test in enumerate(line.wave/compare_lines-1):
                    print('z: %5.2f' %z_test)
                    ygauss = lam*0.
                    gsigma = 60 ## Ang
                    for l in sdss_use:
                        ygauss += np.exp(
                            -1.*( (lam - sdss_lines.wave[l]*(1+z_test) ) 
                            /gsigma)**2)*ymax*0.1*weight[l]
                    
                    scl = 0.1
                    ax.plot(lam,ygauss+scl*ymax*iz,alpha=0.7,color='white', 
                            linewidth=4)
                    ax.plot(lam,ygauss+scl*ymax*iz,alpha=1,color=colors[iline])
                
                iline += 1
                out_lines.append(line)
                
            if (line.flag == 'contam') & (line.type=='em'):
                ax.plot(line.wave*np.array([1,1]),np.array([-1,1]),'--',
                        color='orange',linewidth=2,alpha=0.1)
            
            if (line.flag == 'ok') & (line.type=='abs'):
                ax.plot(line.wave*np.array([1,1]),np.array([-1,1]),
                        color='black',linewidth=2,alpha=0.2)
                ax.plot(line.wave*np.array([1,1]),np.array([-1,1]),'--',
                        color='green',linewidth=2,alpha=0.7)
                #out_lines.append(line)
                
            if (line.flag == 'artifact') & (line.type=='abs'):
                ax.plot(line.wave*np.array([1,1]),np.array([-1,1]),'--',
                        color='green',linewidth=2,alpha=0.1)
                
    ### Show emission line wavelengths
    # zb = np.linspace(0.5,2.5,5)
    # zi = np.linspace(0.5,2.5,20)
    # yz = np.linspace(0,ymax*1.1,5)
    # lines = [3727,4862,5007,6563]
    # for line in lines:
    #     ax.plot(line*(1+zb),yz,color='black',alpha=0.1)
    #     ax.scatter(line*(1+zb),yz,marker='o',color='black',alpha=0.1)
    
    ### Axes
    #plt.semilogx(subsx=[11000,12500,15000])
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(-0.05*ymax,1.1*ymax)
    
    ### Labels
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    if plt.rcParams['text.usetex']:
        ax.set_title(r'%s: \#%d' %(root.replace('_','\_'),object_number))
        ax.set_xlabel(r'$\lambda~\left[$\AA$\right]')
        ax.set_ylabel(r'$\mathit{f}_{\lambda}$')
    else:
        ax.set_title(r'%s: #%d' %(root.replace('_','\_'),object_number))
        ax.set_xlabel(r'$\lambda$ [$\AA$]')
        ax.set_ylabel(r'$f_{\lambda}$')
        
    ### Save to PNG
    if outfile:
        if USE_PLOT_GUI:
            fig.savefig(outfile,dpi=80,transparent=False)
        else:
            canvas = FigureCanvasAgg(fig)
            canvas.print_figure(outfile, dpi=80, transparent=False)
        
    if close_window & USE_PLOT_GUI:
        plt.close()
    
    return out_lines

def test_gaussfit():
    
    fitfunc = lambda p, x: p[0]*np.exp(-(x-p[1])**2/(2.0*p[2]**2))+p[3]+p[4]*(x-p[1])/1.e4
    continuum = lambda p, x: p[3]+p[4]*(x-p[1])/1.e4
    
    x = np.arange(-10,10,0.1)
    pp = np.array([3, 0, 0.5, 1., 0.])
    yin = fitfunc(pp, x)+np.random.normal(size=len(x))*0.4
    
    pfit = threedhst.plotting.gaussfit(x,yin,pp*0.1)
    
    plt.plot(x,yin,color='red',marker='o')
    plt.plot(x,fitfunc(pfit[2][0], x),color='blue')
    plt.plot(x,continuum(pfit[2][0], x),color='green')
    
    
def gaussfit(x,y,p0):
    import scipy
    import scipy.optimize
    
    # define a gaussian fitting function where
    # p[0] = amplitude
    # p[1] = mean
    # p[2] = sigma
    fitfunc = lambda p, x: p[0]*np.exp(-(x-p[1])**2/(2.0*p[2]**2))+p[3]+p[4]*(x-p[1])/1.e4
    #fitfunc = lambda p, x: p[0]*x+p[1]
    continuum = lambda p, x: p[3]+p[4]*(x-p[1])/1.e4
    
    errfunc = lambda p, x, y: (fitfunc(p,x)-y)
    
    use = (np.isfinite(x)) & (np.isfinite(y))
    if len(x[use]) <= 3:
        return -1
        
    output = scipy.optimize.leastsq(errfunc, p0.copy(),
                     args=(x[use],y[use]),
                     full_output=True)
        
    yfit = fitfunc(output[0], x)
    cont = continuum(output[0], x)
    eqwidth = np.trapz(1-yfit/cont, x)
    
    #print output
    return [yfit, eqwidth, output]
    
def makeSpec1dImages(SPCFile, path='./HTML/'):
    """
    makeSpec1dImages(SPCFile, path='./HTML')
    
    Run plotObject for each object in SPCFile
    """
    import os
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    ids = SPCFile._ext_map+0
    ids.sort()
    
    fp = open(path+'/'+root+'_1D_lines.dat','w')
    fp.write('# id lambda sigma eqw snpeak\n# 4 parameters for each detected em. line\n')
    for id in ids:
        idstr = '%05d' %id
        print(threedhst.noNewLine+'plotting.makeSpec1dImages: %s_%s_1D.png' %(root, idstr))
        lines = plot1Dspec(SPCFile, id,
                   outfile=path+'/'+root+'_'+idstr+'_1D.png',
                   close_window=True)
        
        str = '%5d' %id
        for line in lines:
            if line.type == 'emgauss':
                str+='   %8.1f %8.1f %9.1e %7.1f' %(line.wave, line.sigma, line.ew, line.sn)
            if line.type == 'abs':
                str+='   %8.1f' %(-1*line.wave)
                
        fp.write(str+'\n')
    
    fp.close()
    
def makeHTML(SPCFile, mySexCat, mapParams,
             output='./HTML/index.html', title=None):
    """
    makeHTML(SPCFile, mySexCat, mapParams,
             output='./HTML/index.html', title=None)
    """
    import os
    from socket import gethostname as hostname
    
    if hostname().startswith('uni'):
        GMAPS_KEY = 'ABQIAAAAzSrfHr_4F2D2YfSuYQD2ZBRWJmdnNuPtXxsK1b3ql6FJMkf8bxT-OiTDeKiGIrffKkTBi-in1FVtrw'
    else:
        # localhost
        GMAPS_KEY = 'ABQIAAAA1XbMiDxx_BTCY2_FkPh06RR20YmIEbERyaW5EQEiVNF0mpNGfBSRb_rzgcy5bqzSaTV8cyi2Bgsx3g'
    
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    if not title:
        title = root
    
    #### Header
    lines = ["""
<html>
    <head>
    
    <link rel="stylesheet" href="scripts/style.css" type="text/css" id="" media="print, projection, screen" /> 

    <script type="text/javascript" src="scripts/jquery-1.4.2.min.js"></script> 

    <script type="text/javascript" src="scripts/jquery.sprintf.js"></script> 

    <script type="text/javascript" src="scripts/threedhst.js"></script> 

    <!--  www.astro.yale.edu/ --> 
    <!-- 
    <script src="http://maps.google.com/maps?file=api&amp;v=3&amp;key=ABQIAAAAzSrfHr_4F2D2YfSuYQD2ZBTGxfEdj5ixTzExLHeue1TdBmBBTxSo_kKvXxnIoUhPTW743ryzjWdouQ"
     type="text/javascript"></script> 
    --> 
    
    <!-- localhost -->
    <script src="http://maps.google.com/maps?file=api&amp;v=3&amp;key=%s" type="text/javascript"></script> 
    """ %(GMAPS_KEY)]
        
    #### Script for the Google map
    llSW = mapParams['LLSW']
    llNE = mapParams['LLNE']
    center = mapParams['LLCENTER']
    lng_offset = mapParams['LNG_OFFSET']
    
    lines.append("""
    <script type="text/javascript"> 
    
    //////////// Global variables
    var map = 0;
    var centerLat = %f;
    var centerLng = %f;
    // var offset = %f;
    var offset = 0.0;
    var zoomLevel = %f;
    var root = '%s';
    
    var myIcon = new GIcon();
   	myIcon.iconSize = new GSize(30,25);
   	myIcon.iconAnchor = new GPoint(15,14);
   		
    function initialize() {        
        if (GBrowserIsCompatible()) {
            map = new GMap2(document.getElementById("map"));
            // map.addControl(new GScaleControl());
            var topLeft = new GControlPosition(G_ANCHOR_TOP_LEFT,
                                                 new GSize(10,40));
            map.addControl(new GSmallMapControl(), topLeft);

		    ///// Add the layer tiles
		    add_mapLayer_tiles()
            
            latLng2raDec();
            
            GEvent.addListener(map, "moveend", function() {
                latLng2raDec();
                show_centerbox();
            });
            
            ///// Add the green circles around the catalog objects
            plotXmlObjects();
        }
        initialize_SED_column();
    }
    
    </script>
        """ %(center[0],center[1],lng_offset,mapParams['ZOOMLEVEL'],
              threedhst.options['ROOT_GRISM']))
    
    #### HTML Body   
    lines.append("""

</head>
<body onload="initialize()" onunload="GUnload()">
    
    <div id="map"></div>
    
    <div id="title">
        %s
        <a href="./%s_drz.cat" class="dl"> cat </a>
        <a href="./images/%s_thumbs.tar.gz" class="dl"> thumbs </a>
        <a href="./ascii/%s_spec.tar.gz" class="dl"> 1Dspec </a>
        <a href="./images/%s_2D.tar.gz" class="dl"> 2Dspec </a>
        <a href="./images/%s_1D_lines.dat" class="dl"> lines </a>
        <a href="%s.threedhst.param" class="dl"> info </a>
    </div>
    
    <img src="scripts/3dHST.png" id="logo">
    
    <div onclick="javascript:switch_layout()" id="switchbox">
        ||
    </div>
    
    <div onclick="javascript:toggle_markers()" id="markerbox"> </div>
    
    <div id="centerbox"></div>
    
    <div id="coords">
    <form onsubmit="return false" style="display:inline;">
        <input type="text" value="#id" class="cinput" id="idInput" maxlength="4" onchange="centerOnID()"/>
        <input type="text" value="00:00:00.00" class="cinput" id="raInput" maxlength="11" onchange="centerOnInput()"/>
        <input type="text" value="+00:00:00.0" class="cinput" id="decInput" maxlength="11" onchange="centerOnInput()"/>
        </form>
        <a id="vizierLink" href="http://vizier.u-strasbg.fr/viz-bin/VizieR?-c=12:36:36.85+%%2B62:06:58.7&-c.rs=1" target="_blank"><img src="scripts/glass.png" style="width:12px; margin:0px; padding:0px"></a>
    </div>
    
    """ %(title, title, title, title, title, title, title))
    
    MAG_COLUMN = 'MAG_F1392W'
    for col in mySexCat.column_names:
        if col.startswith('MAG_F'):
            MAG_COLUMN = col
    MAG_COL_ID = mySexCat.searchcol(MAG_COLUMN)
    
    lines.append("""
    
    <div id="content" style="width:840px;height:100%%;overflow:auto">
    
    <table id="myTable" cellspacing="1" class="tablesorter"> 
    <thead> 
    <tr id="header_row"> 
    <th width=50px>ID 
		<span style="cursor:n-resize">
			<a onclick="javascript:pageUp();"> <b>+</b> </a>
		</span>
		<span style="cursor:s-resize">
			<a onclick="javascript:pageDown();"> <b>-</b> </a>
		</span>
    </th>
        <th width=70px>R.A.</th> 
        <th width=70px>Dec.</th> 
        <th width=100px>mag %s</th> 
        <th width=140px>Thumb</th> 
        <th width=210px>1D Spec</th> 
        <th width=210px>2D Spec</th> 
    </tr> 
    </thead> 
    <tbody> 
    """ %(MAG_COLUMN.split('MAG_')[1]))
        
    # for id in SPCFile._ext_map:
    for id in [SPCFile._ext_map[0]]:
        for idx,num in enumerate(mySexCat.NUMBER):
            if num == str(id):
                break
        
        ra  = mySexCat.X_WORLD[idx]
        dec = mySexCat.Y_WORLD[idx]
        mag = mySexCat.getcol(MAG_COL_ID)[idx]
        
        img = '%s_%05d' %(root,id)
        lines.append("""
        <tr id=""> 
            <td id="tab_id" onclick="javascript:recenter(0,0)">
                  %d
            </td>
            <td id="tab_ra">%13.6f</td> 
            <td id="tab_dec">%13.6f</td> 
            <td id="tab_mag">%6.2f</td> 
            <td id="tab_thumb">
                <a id="a_thumb" href='images/%s_thumb.fits.gz'><img id="img_thumb" src='images/%s_thumb.png' width=133px></a></td> 
            <td id="tab_1D">
                <a id="a_1D" href='ascii/%s.dat'><img id="img_1D" src='images/%s_1D.png' width=200px title='ascii'></a></td> 
            <td id="tab_2D">
                <a id="a_2D" href='images/%s_2D.fits.gz'><img id="img_2D" src='images/%s_2D.png' width=200px></a></td> 
        </tr> 
        """ %(id,
              np.float(ra),np.float(dec),np.float(mag),
              img,img,img,img,img,img))
    
    lines.append("""
    </tbody>
    </div> <!-- content -->
    </body>
    </html>""")
    
    if not output:
        output='HTML/index.html'
    
    fp = open(output,'w')
    fp.writelines(lines)
    fp.close()
    
def makeCSS(path="../HTML/scripts", title_size=18):
    """
makeCSS(path="../HTML/scripts")
    """
    
    fp = open(path+"/style.css","w")
    fp.write("""
    #title {
        border: solid 1px black;
        font-family:Verdana;
        font-size: %0dpt;
        font-weight:bold;
        position:absolute;
        top:5px;
        left:5px;
        width:1087px;
        padding:9px 0px 0px 49px;
        height:45px;
        vertical-align: middle;
        background: #f6f6f6;
    }

    #logo {
        position:absolute;
        height:40px;
        top:10px;
        left:10px;
    }

    a.dl {
        font-family:Verdana;
        font-size:6pt;
        color:black;
        text-decoration: none;
        padding-left:10px;
    }

    a.dl:hover {
        color:red;
        text-decoration: none;
    }

    #map {
        position: absolute;
        top:55px;
        left:5px;
        width: 300px;
        height: 300px;
        border: solid 1px black;
        vertical-align: middle;
        background: white;
        text-align: center;
    }

    #switchbox {
        background: white;
        color: black;
        font-family:Verdana;
        font-size: 6pt;
        border: 1px solid black;
        opacity: 0.8;
        position:absolute;
        top:63px;
        left:28px;
        width: 10px;
        padding-left:2px;
        height: 12px;
    }
    
    #markerbox {
        border: 2px solid #00FF03;
        opacity: 0.8;
        position:absolute;
        top:80px;
        left:28px;
        width: 8px;
        padding-left:2px;
        height: 10px;
    }
    
    #coords {
        background: white;
        color: black;
        font-family:Verdana;
        font-size: 8pt;
        text-align: center;
        border: 1px solid black;
        opacity: 0.8;

        position: absolute;
        top:41px;
        left:145px;
        width: 214px;
        height:16px;
        padding-bottom:2px;

        -moz-border-radius-topleft: 8px;
        -webkit-border-top-left-radius: 8px;
        -moz-border-radius-topright: 8px; 
        -webkit-border-top-right-radius: 8px;
    }

    #centerbox {
        width:20px;
        height:20px;

        position: absolute;
        top:100px;
        left:100px;
        opacity:0.0;
        border: 2px solid yellow;
    }

    .cinput {
        border: none;
        font-size:8pt;
        font-family:Verdana;
        width:75px;
    }

    #idInput {
        width:32px;
        text-align:right;
    }

    #content {
        width:840px;
        position:absolute;
        top:55px;
        left:301px;
        height:100%%;
        overflow:auto;
        border: solid 1px black;
    }


    /* tables */
    table.tablesorter {
    	font-family:Verdana;
    	background-color: #CDCDCD;
    	font-size: 8pt;
    /*  width: 70%%;
    */	text-align: left;
    }
    table.tablesorter thead tr th, table.tablesorter tfoot tr th {
    	background-color: #e6EEEE;
    	border: 1px solid #FFF;
    	font-size: 8pt;
    	padding: 4px;
    }
    table.tablesorter thead tr .header {
    	background-image: url(bg.gif);
    	background-repeat: no-repeat;
    	background-position: center right;
    	cursor: pointer;
    }
    table.tablesorter tbody td {
    	color: #3D3D3D;
    	padding: 4px;
    	background-color: #FFF;
    	vertical-align: top;
    }
    table.tablesorter tbody tr.odd td {
    	background-color:#F0F0F6;
    }
    table.tablesorter thead tr .headerSortUp {
    	background-image: url(asc.gif);
    }
    table.tablesorter thead tr .headerSortDown {
    	background-image: url(desc.gif);
    }
    table.tablesorter thead tr .headerSortDown, table.tablesorter thead tr .headerSortUp {
    background-color: #8dbdd8;
    }
    """ %(title_size))
    fp.close()
    
    
def asciiSpec(SPCFile, root="spec", path="../HTML/ascii"):
    """
asciiSpec(SPCFile, root="spec", path="../HTML/ascii")
    
    Put ASCII spectra in HTML/ascii.
    """
    import os
    import gzip
    import tarfile
    
    try:
        os.mkdir(path)
    except:
        pass
    
    ids = SPCFile._ext_map+0
    ids.sort()
    
    #### hack: errors seem too large when you use 
    #### drizzle with updated DRZSCALE parameter.
    #### scale by DRZSCALE/0.128254
    #ERROR_SCALE = np.float(threedhst.options['DRZSCALE'])/0.128254
    ########    ---- this was a bug fixed by M. Kuemmel 12/16/10  ----
    
    for id in ids:
        spec = SPCFile.getSpec(id)
        lam  = spec.field('LAMBDA')
        flux = spec.field('FLUX')
        ferr = spec.field('FERROR')
        contam = spec.field('CONTAM')
        
        sp1d = threedhst.spec1d.extract1D(id, root=SPCFile.filename.split('_2_opt')[0], path='../HTML', show=False, out2d=False)
        lam2 = sp1d['lam']
        flux2 = sp1d['flux']
        ferr2 = sp1d['error']
        contam2 = sp1d['contam']
        
        out = root+'_%05d.dat' %id
        print(threedhst.noNewLine+out)
        
        # fp = gzip.open(path+'/'+out+'.gz','wb')
        fp = open(path+'/'+out,'w')
        fp.write('# lam flux error contam flux2 error2 contam2\n')
        
        NL = len(lam)
        for i in range(NL):
            fp.write('%11.5e %10.3e %10.3e %10.3e %10.3e %10.3e %10.3e\n' %(lam[i],flux[i],ferr[i],contam[i], flux2[i], ferr2[i], contam2[i]))
        fp.close()
        
    ### make a tarfile
    oldwd = os.getcwd()
    os.chdir(path)
    fptar = tarfile.open(root+'_spec.tar.gz','w|gz')
    for id in ids:
        out = root+'_%05d.dat' %id
        fptar.add(out)
    fptar.close()
    os.chdir(oldwd)
    
class SPCFile(object):
    """
    SPCFile
    
    Class for reading and plotting spectra from an aXe
    
    """
    def _getPath(self):
        """
        _getPath()
        
        Figure out if we're in the root directory or in DATA
        """
        import os
        if os.getcwd().split('/')[-1] == 'DATA':
            self.path = '../'+self.axe_drizzle_dir+'/'
        elif os.getcwd().split('/')[-1] == self.axe_drizzle_dir:
            self.path = './'
        else:
            self.path = self.axe_drizzle_dir+'/'
    
    def _mapFitsExtensions(self):
        """
_mapFitsExtensions()
        
    Figure out which object corresponds to which extension in the SPC.fits file
        """
        for i in range(self.N_ext):
            self._ext_map[i] = np.int(
               self.fits[i+1].header['EXTNAME'].split('BEAM_')[1][:-1])
            
    def __init__(self, filename='ib3721050_2_opt.SPC.fits',
                 axe_drizzle_dir='DRIZZLE_G141'):
        """
__init__(filename='ib3721050_2_opt.SPC.fits',
         axe_drizzle_dir='DRIZZLE_G141')
        """
        self.filename = filename
        self.axe_drizzle_dir = axe_drizzle_dir
        self._getPath()
        self.fits = pyfits.open(self.path+filename)
        
        self.N_ext = len(self.fits)-1
        self._ext_map = np.arange(self.N_ext)*0
        self._mapFitsExtensions()
        
    def getSpec(self, object_number):
        """
        getSpec(self, object_number)
        """
        if object_number in self._ext_map:
            idx = np.where(self._ext_map == object_number)[0][0]
            return self.fits[idx+1].data
        else:
            print("Object #%d not found in %s." %(object_number, self.filename))
            return False

def makeJavascript(path="../HTML/scripts"):
    """
make_Javascript(path="../HTML/scripts")
    """

    fp = open(path+"/threedhst.js","w")
    fp.write("""
    /////////////////// Ready function
    $(document).ready(function() {
    	switch_layout();		
    });

    var layout = -1; // Start with horizontal layout
    function switch_layout() {
        if (layout == 0) {
            layout = 1;
            vertical_layout();
            $("#switchbox").text("||");
            $("#switchbox").css("cursor","s-resize");  
            map.checkResize();                    
        } else {
            layout = 0;
            horizontal_layout();
            $("#switchbox").text("=");
            $("#switchbox").css("cursor","e-resize");
            map.checkResize();          
        }
    }

    ///// Turn on/off green object markers
    var markers_on=1;
    function toggle_markers() {
    	if (markers_on == 0) {
    	    $("#markerbox").css("border","2px solid #00FF03");  
            markers_on = 1;
            for (var i = 0; i < marker_list.length; i++) {
    			marker_list[i].show();
    	    }
    	} else {
    	    $("#markerbox").css("border","2px solid #AAAAAA");  
    	    markers_on=0
            for (var i = 0; i < marker_list.length; i++) {
    			marker_list[i].hide();
    	    }
    	}
    }

    ////// Map at left, spectral images at right
    function vertical_layout() {

    	//$("#title").css("width",1087);
    	$("#title").css("width",$(window).width()-12-
    	    parseFloat($("#title").css("padding-left"))+
    	    parseFloat($("#title").css("padding-right")));

    	$("#content").css("height",$(window).height()-60-5);
    	//$("#content").css("width","840px");
    	$("#content").css("width",$(window).width()-295+2);
    	$("#content").css("top","60px");
    	$("#content").css("left","301px");

    	$("#map").css("height",$(window).height()-60);	
    	$("#map").css("width",300-5);	

        $("#coords").css("left",
    	    parseInt($("#map").css("width"))/2.-
    	    parseInt($("#coords").css("width"))/2.+10);

    	c = $("#centerbox");
        c.css("left",parseFloat($("#map").css("left"))+
                     parseFloat($("#map").css("width"))/2.-
                     parseFloat(c.css("width"))/2.-0);

        c.css("top",parseFloat($("#map").css("top"))+
                     parseFloat($("#map").css("height"))/2.-
                     parseFloat(c.css("height"))/2.-0);

    	addRowSet();

    }

    ////// Map on top, spectra info on bottom
    function horizontal_layout() {

    	$("#title").css("width",$(window).width()-12-
    	    parseFloat($("#title").css("padding-left"))+
    	    parseFloat($("#title").css("padding-right")));

    	$("#content").css("height",170);
    	$("#content").css("width",$("#title").width()+
    	    parseFloat($("#title").css("padding-left"))+
    	    parseFloat($("#title").css("padding-right")));
    	//alert(parseFloat($("#title").css("padding-left"))+20);
    	$("#content").css("top",$(window).height()-170);
        $("#content").css("left",$("#map").css("left"));

    	$("#map").css("height",$(window).height()-60-170-11);
    	$("#map").css("width",$("#title").width()+
    	    parseFloat($("#title").css("padding-left"))+
    	    parseFloat($("#title").css("padding-right")));    

    	$("#coords").css("left",
    	    parseInt($("#map").css("width"))/2.-
    	    parseInt($("#coords").css("width"))/2.);    

    	c = $("#centerbox");
        c.css("left",parseFloat($("#map").css("left"))+
                     parseFloat($("#map").css("width"))/2.-
                     parseFloat(c.css("width"))/2.-0);

        c.css("top",parseFloat($("#map").css("top"))+
                     parseFloat($("#map").css("height"))/2.-
                     parseFloat(c.css("height"))/2.-0);

    	clearRows();

    }

    /////////////////////
    /////  Map utilities
    /////////////////////
    
    ///// Add the layer tiles
    function add_mapLayer_tiles() {
        
        var copyright = new GCopyright(1,
             new GLatLngBounds(new GLatLng(-2,-2),
                               new GLatLng(2,2)),
                               zoomLevel, "3D-HST");
        
        var copyrightCollection = new GCopyrightCollection('Map Data:');
        copyrightCollection.addCopyright(copyright);
        
        // Direct image tiles
        CustomGetDirectTileUrl=function(a,b){
            return "tiles/"+root+"_d_"+a.x+"_"+a.y+"_"+b+".png"
        }
        var tilelayersDirect = [new GTileLayer(copyrightCollection,
                                      zoomLevel-2,zoomLevel+1)];
        tilelayersDirect[0].getTileUrl = CustomGetDirectTileUrl;
        var custommapDirect = new GMapType(tilelayersDirect, 
               new GMercatorProjection(zoomLevel+2), "Direct");
        map.addMapType(custommapDirect);

        // Grism image tiles
        CustomGetGrismTileUrl=function(a,b){
            return "tiles/"+root+"_g_"+a.x+"_"+a.y+"_"+b+".png"
        }
        var tilelayersGrism = [new GTileLayer(copyrightCollection,
                                      zoomLevel-2,zoomLevel+1)];
        tilelayersGrism[0].getTileUrl = CustomGetGrismTileUrl;
        var custommapGrism = new GMapType(tilelayersGrism, 
               new GMercatorProjection(zoomLevel+2), "Grism");
        map.addMapType(custommapGrism);
        
        // Model image tiles
        CustomGetModelTileUrl=function(a,b){
            return "tiles/"+root+"_m_"+a.x+"_"+a.y+"_"+b+".png"
        }
        var tilelayersModel = [new GTileLayer(copyrightCollection,
                                      zoomLevel-2,zoomLevel+1)];
        tilelayersModel[0].getTileUrl = CustomGetModelTileUrl;
        var custommapModel = new GMapType(tilelayersModel, 
               new GMercatorProjection(zoomLevel+2), "Model");
        map.addMapType(custommapModel);    
        
        // Can't remove all three for some reason
        map.removeMapType(G_NORMAL_MAP);
        map.removeMapType(G_HYBRID_MAP);
        map.removeMapType(G_SATELLITE_MAP);
        
        map.addControl(new GMapTypeControl());
        
        // Set map center
        map.setCenter(new GLatLng(0.0, offset), zoomLevel,
         custommapDirect);
		
    }
    
    ///// Flash a little box to show where the center of the map is, 
    ///// which corresponds to the listed coordinates
    function show_centerbox() {
        c = $("#centerbox");
        document.getElementById("centerbox").style.display = "block";
        c.animate({ opacity: 1.0}, 300, function() { });        
        c.animate({ opacity: 0.0}, 300, function() {
            document.getElementById("centerbox").style.display = "none";
        });
    }

    ///// Convert RA/Dec to Map Lat/Lng coordinates, which are centered around 0,0
    function latLng2raDec() {
        var mapcenter = map.getCenter();
        var dec = mapcenter.lat()+centerLat;                       
        var dsign = "+";
        var hex = "%2B"
        if (dec < 0) {
            dsign = "-";
            hex = "%2D";
        }
        dec = Math.abs(dec);
        var ded = parseInt(dec);
        var dem = parseInt((dec-ded)*60);
        var des = parseInt(((dec-ded)*60-dem)*60);
        var dess = parseInt((((dec-ded)*60-dem)*60-des)*10);
        if (ded < 10) {ded = "0"+ded;} 
        if (dem < 10) {dem = "0"+dem;} 
        if (des < 10) {des = "0"+des;} 
        var decstr = ded+":"+dem+":"+des+"."+dess;
        document.getElementById("decInput").value = dsign + decstr;

        var ra = ((360-mapcenter.lng()/Math.cos(centerLat/360.*2*3.14159)+offset-centerLng)/360.*24);
        var rah = parseInt(ra);
        var ram = parseInt((ra-rah)*60);
        var ras = parseInt(((ra-rah)*60-ram)*60);
        var rass = parseInt((((ra-rah)*60-ram)*60-ras)*100);
        if (rah < 10) {rah = "0"+rah;} 
        if (ram < 10) {ram = "0"+ram;} 
        if (ras < 10) {ras = "0"+ras;} 
        if (rass < 10) {rass = "0"+rass;} 
        var rastr = rah+":"+ram+":"+ras+"."+rass;
        document.getElementById("raInput").value = rastr;   

        document.getElementById("vizierLink").href = "http://vizier.u-strasbg.fr/viz-bin/VizieR?-c=" + rastr + "+" + hex + decstr +  "&-c.rs=1";        

    }

    ////// Center the map after changing the coordinates in the input box
    function centerOnInput() {
        var rastr = document.getElementById("raInput").value;
        var rsplit = rastr.split(":");
        if (rsplit.length != 3) rsplit = rastr.split(" ");
        var ra = parseFloat(rastr);

        if (rsplit.length == 3) {
            ra = (parseInt(rsplit[0])+
                parseInt(rsplit[1])/60.+
                parseFloat(rsplit[2])/3600.)/24.*360;
        } 

        var decstr = document.getElementById("decInput").value;
        var dsplit = decstr.split(":");
        if (dsplit.length != 3) dsplit = decstr.split(" ");
        var dec = parseFloat(decstr);
        if (rsplit.length == 3) {
            dec = Math.abs(parseInt(dsplit[0])+
                Math.abs(parseInt(dsplit[1])/60.)+
                Math.abs(parseFloat(dsplit[2])/3600.));

            /// Don't know why, but need to do this twice
            dec = Math.abs(parseInt(dsplit[0]))+
            parseInt(dsplit[1])/60.+
            parseFloat(dsplit[2])/3600.;

            if (parseFloat(dsplit[0]) < 0) {
                dec *= -1;
            }
        }

        recenter(ra,dec);
        latLng2raDec();
    }

    //// Pan to object ID entered in the input box, show spectra + info
    function centerOnID() {
        var id = document.getElementById("idInput").value;
        for (j=0; j < ids.length; j++) {
            if (ids[j] == id) {
                ROWSTART = j;
                setFirstRow();
                recenter(0,0);
                break;
            }
        }
    }

    ////// Globals
    var ra_list = [];
    var de_list = [];
    var mag_list = [];
    var lats = [];
    var lngs = [];
    var ids = [];
    var nObject = 0;
    var marker_list = [];
    
    var SED_COLUMN = 1;
    var ROWSTART = 0;
	var NSHOW = 25;
	
	
    ///// Read objects from XML file and plot regions
    function plotXmlObjects() {
        GDownloadUrl(root+".xml", function(data) {
            var xml = GXml.parse(data);
            var markers = xml.documentElement.getElementsByTagName("marker");
            nObject = markers.length;

            for (var i = 0; i < markers.length; i++) {
                // Read from XML
                var id = markers[i].getAttribute("id");
                var ra = markers[i].getAttribute("ra");
                var dec = markers[i].getAttribute("dec");
                var mag = markers[i].getAttribute("mag");
                var lat = dec-centerLat;
                var lng = ((360-ra)-centerLng+offset)*Math.cos(centerLat/360.*2*3.14159);

    			ra_list.push(ra);
    			de_list.push(dec);
    			mag_list.push(mag);
                lats.push(lat);
                lngs.push(lng);
                ids.push(id);

                // The marker
                myIcon.image = "scripts/circle.png"; // ".php?id="+id;
                markerOptions = { icon:myIcon, title:id};
                var point = new GLatLng(lat,lng);
                var marker = new GMarker(point, markerOptions);

                // The only thing passed to the listener is LatLng(), 
                // so need to find which object is closest to the clicked marker.
                GEvent.addListener(marker, "click", function(self) {
                    //alert(self.lat()+' '+nObject+' '+lats[0]);
                    var matchID = 0;
                    var lat = self.lat();
                    var lng = self.lng();
                    for (var j=0;j<nObject;j++) {
                        var dj = (lats[j]-lat)*(lats[j]-lat)+
                                 (lngs[j]-lng)*(lngs[j]-lng) ;
                        if (dj == 0) {
                            matchID = ids[j];
                            break;
                        }
                    }
                    //alert(matchID);
                    //window.location = '#i'+matchID;

    				ROWSTART = j;
    				setFirstRow();

    				if (layout == 1) {
    					clearRows();
    					addRowSet();
    				}
                });

                marker_list.push(marker);
                map.addOverlay(marker);            
            }
        });
    }

    function recenter(ra,dec) {
    	if ((ra+dec) == 0) {
    		ra = document.getElementById("tab_ra").innerText;
    		dec = document.getElementById("tab_dec").innerText;
    	}

    	var lat = dec-centerLat;
        var lng = ((360-ra)-centerLng+offset)*Math.cos(centerLat/360.*2*3.14159)
        map.setCenter(new GLatLng(lat,lng)); //, zoomLevel);
    }

    /////////////////////
    /////  Dynamic control of the table showing the spectra/thumbnails
    /////////////////////
    function setFirstRow() {
    	if (ROWSTART < 0) ROWSTART=0;
    	j = ROWSTART;

    	var matchID = 0+$.sprintf('%d', ids[j]);
    	if (ids[j] < 1000) matchID = 0+matchID;
    	if (ids[j] < 100) matchID = 0+matchID;
    	if (ids[j] < 10) matchID = 0+matchID;

    	document.getElementById("tab_id").innerText = ids[j];
    	document.getElementById("tab_ra").innerText = $.sprintf('%14.6f', ra_list[j]);
    	document.getElementById("tab_dec").innerText = $.sprintf('%14.6f', de_list[j]);
    	document.getElementById("tab_mag").innerText = $.sprintf('%6.2f', mag_list[j]);

    	document.getElementById("img_thumb").src = 'images/'+root+'_' + matchID + '_thumb.png';
    	document.getElementById("a_thumb").href = 'images/'+root+'_' + matchID + '_thumb.fits.gz';

    	document.getElementById("img_1D").src = 'images/'+root+'_' + matchID + '_1D.png';
    	document.getElementById("a_1D").href = 'ascii/'+root+'_' + matchID + '.dat';

    	document.getElementById("img_2D").src = 'images/'+root+'_' + matchID + '_2D.png';
    	document.getElementById("a_2D").href = 'images/'+root+'_' + matchID + '_2D.fits.gz';

        if (SED_COLUMN == 1) {
            document.getElementById("img_SED").src = 'SED/'+root+'_' + matchID + '_SED.png';
        	document.getElementById("a_SED").href = 'SED/'+root+'_' + matchID + '_SED.png';

        }

    }

    //////// Clear all but the first row for horizonal layout and redrawing
    function clearRows() {
    	var theTable =  document.getElementById("myTable");
    	var NROW = theTable.rows.length;
    	if (NROW > 2) {
    		for (var i = NROW-1; i > 1; i--) {
    			theTable.deleteRow(i);
    		}
    	}
    }

    /////// Add NSHOW rows, starting with index 'ROWSTART'
    function addRowSet() {
    	if (ROWSTART < 0) {
    		ROWSTART=0;
    	}
    	stop = ROWSTART + NSHOW;
    	if (stop >= ids.length) stop = ids.length-1;
    	if (layout == 1) {
    		for (var i=ROWSTART+1; i<=stop; i++) {
    			addRow(i);
    		}
    	}
    }

    ////// Add a row to the display table
    function addRow(idx) {
    	var theTable =  document.getElementById("myTable");
    	var NROW = theTable.rows.length;
    	var newRow = theTable.insertRow(NROW);

    	thisID = 0+ids[idx]
    	if (ids[idx] < 1000) thisID = 0+thisID;
    	if (ids[idx] < 100) thisID = 0+thisID;
    	if (ids[idx] < 10) thisID = 0+thisID;

    	var tid = newRow.insertCell(0);
    	tid.innerHTML = "<a onclick='javascript:recenter("+$.sprintf('%14.6f', ra_list[idx])+","+$.sprintf('%14.6f', de_list[idx])+");'>"+ids[idx]+"</a>";

    	var tra = newRow.insertCell(1);
    	tra.innerText = $.sprintf('%14.6f', ra_list[idx]);

    	var tdec = newRow.insertCell(2);
    	tdec.innerText = $.sprintf('%14.6f', de_list[idx]);

    	var tmag = newRow.insertCell(3);
    	tmag.innerText = $.sprintf('%6.2f', mag_list[idx]);

    	var tthumb = newRow.insertCell(4);
    	tthumb.innerHTML = "<a href='images/"+root+"_"+thisID+ "_thumb.fits.gz'> <img src='images/"+root+"_"+thisID+"_thumb.png'  width=133px></a>";

    	var t1D = newRow.insertCell(5);
    	t1D.innerHTML = "<td><a href='ascii/"+root+"_"+thisID+".dat'><img src='images/"+root+"_"+thisID+"_1D.png' width=200px title='ascii'></a></td>"

    	var t2D = newRow.insertCell(6);
    	t2D.innerHTML = "<td><a href='ascii/"+root+"_"+thisID+"_2D.fits.gz'><img src='images/"+root+"_"+thisID+"_2D.png' width=200px title='2D'></a></td>"

        if (SED_COLUMN == 1) {
        	var tSED = newRow.insertCell(7);
        	tSED.innerHTML = "<td><a href='SED/"+root+"_"+thisID+"_SED.png'><img src='SED/"+root+"_"+thisID+"_SED.png' width=200px title='2D'></a></td>"
        }

    }

    ///// Page down through objects on vertical view
    function pageDown() {
        if (layout == 1) {
            ROWSTART+=NSHOW;
        } else {
            ROWSTART+=1;
        }
        setFirstRow();
        clearRows();
        addRowSet();
    }

    ///// Page up on vertical view
    function pageUp() {
        if (layout == 1) {
            ROWSTART-=NSHOW;
        } else {
            ROWSTART-=1;
        }
        setFirstRow();
        clearRows();
        addRowSet();
    }
    
    ///// Add an additional column for plots with BB SEDs compared
    ///// to the spectra.
    function initialize_SED_column() {
        if (SED_COLUMN == 1) {
        	var theTab = document.getElementById("myTable");
        	
        	//// <th> column, called "SED"
        	var tabHead = theTab.tHead;
        	thObj = document.createElement("th");
        	thObj.style.width="210px";
        	thObj.innerHTML = "SED";
        	tabHead.rows[0].appendChild(thObj);
        	
        	//// Data column, with a link and an image
        	var newCell = theTab.tBodies[0].rows[0].insertCell(7);
        	newCell.id = "tab_sed";
        	var newLink = document.createElement("a");
        	newLink.id="a_SED";
        	var newImg = document.createElement("img");
        	newImg.id="img_SED";
        	newImg.style.width="200px";
        	newImg.src = "SED/"+root+"_0001.png";
        	newLink.appendChild(newImg);
        	newCell.appendChild(newLink);

            //// newCell.innerHTML = "<a id=\"a_sed\"><img id=\"img_sed\" src=\"SED/"+root+"_0001.png\" width=200px></a>";
        }
    }
    
    """)
    fp.close()
        