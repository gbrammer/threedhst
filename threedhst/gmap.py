"""
Convert Multidrizzle drz.fits (North up) to Google Maps tiles
"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import os
import glob

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import numpy as np

# Specifies the size of the map (in pixels).
TILE_SIZE = 256
MAP_SIZE = [TILE_SIZE,TILE_SIZE]    
# This is the Maps API key for running on localhost:8080
MAP_KEY = 'ABQIAAAA1XbMiDxx_BTCY2_FkPh06RR20YmIEbERyaW5EQEiVNF0mpNGfBSRb' \
    '_rzgcy5bqzSaTV8cyi2Bgsx3g'

def makeGMapTiles(fitsfile=None,outPath=None,tileroot='direct', extension=1,
                  zmin=-0.1, zmax=1, rgb_params=(5, 3, -0.05), verbose=False, 
                  rgb_clip=True):
    """
makeGMapTiles(fitsfile=None,outPath=None,tileroot='direct', extension=1,
              zmin=-0.1, zmax=1)
    
    Make Google map tiles for an input FITS image, which is assumed to be
    North-up, East-left, like normal Multidrizzle output.
    
    "fitsfile" can be a single FITS image or a comma-separated list of three
    images (R,G,B) that will be used to generate a 3-color image.  
    """
    import pywcs
    #import fitsimage
    
    if fitsfile is None:
        fitsfile = 'ib3721050_drz.fits'
    
    if outPath is None:
        outPath = '/tmp/'
    
    RGB = False
    
    if ',' in fitsfile:
        fitsfile = fitsfile.split(',')
        RGB = True
        if len(fitsfile) != 3:
            threedhst.showMessage("`fitsfile` must be either a single image filename or a comma-separated list of 3 images for RGB", warn=True)
            raise ValueError
    else:
        fitsfile=[fitsfile]
        
    # print fitsfile, outPath
    
    ### Read the FITS file
    fi = pyfits.open(fitsfile[0])
    head = fi[extension].header
    data = [fi[extension].data]
    if RGB:
        for file in fitsfile[1:]:
            fi = pyfits.open(file)
            data.append(fi[extension].data)
            
    xsize, ysize = data[0].shape #[1]
    
    ### Image corners in Lat/Lon
    wcs = pywcs.WCS(head)
    llSE = radec2latlon(wcs.wcs_pix2sky([[wcs.naxis1,1]],1)[0])
    llNE = radec2latlon(wcs.wcs_pix2sky([[wcs.naxis1,wcs.naxis2]],1)[0])
    llSW = radec2latlon(wcs.wcs_pix2sky([[1,1]],1)[0])
    llNW = radec2latlon(wcs.wcs_pix2sky([[1,wcs.naxis2]],1)[0])
    llCenter = (llSW+llNE)/2.
    # print llNE,llSW
    # print llCenter
    
    lng_offset = 90
    lng_offset = 0.
    
    params = {}
    params['LNG_OFFSET'] = lng_offset
    params['LLNE'] = llNE*1.
    params['LLSW'] = llSW*1.
    params['LLCENTER'] = llCenter*1.
    
    #makeMapHTML(llSW,llNE,lng_offset=lng_offset)
    
    llSW[1] += lng_offset-llCenter[1]
    llSE[1] += lng_offset-llCenter[1]
    llNW[1] += lng_offset-llCenter[1]
    llNE[1] += lng_offset-llCenter[1]
    llCenter[1] += lng_offset-llCenter[1]
    
    cdelt = np.cos(llCenter[0]/360.*2*np.pi)
    llSW[1] *= cdelt
    llSE[1] *= cdelt
    llNW[1] *= cdelt
    llNE[1] *= cdelt
    
    llSW[0] -= llCenter[0]
    llSE[0] -= llCenter[0]
    llNW[0] -= llCenter[0]
    llNE[0] -= llCenter[0]
    llCenter[0] -= llCenter[0]
    
    ### Get Google Map pixel/tile coordinates
    m = MercatorProjection()
    bounds = [llSW,llNE]
    view_size = [wcs.naxis1,wcs.naxis2]
    zoomLevel = m.CalculateBoundsZoomLevel(bounds, view_size)
    pixScale = 1./np.array(m.pixels_per_lon_degree)-np.abs(head['CD2_2'])
    zoomLevel = (np.where(np.abs(pixScale) == np.min(np.abs(pixScale))))[0][0]
    params['ZOOMLEVEL'] = zoomLevel
    
    pixSW = m.FromLatLngToPixel(llSW,zoomLevel)
    pixSE = m.FromLatLngToPixel(llSE,zoomLevel)
    pixNW = m.FromLatLngToPixel(llNW,zoomLevel)
    pixNE = m.FromLatLngToPixel(llNE,zoomLevel)
    pixCenter = m.FromLatLngToPixel(llCenter,zoomLevel)
    
    ### Padding to make the output image size
    ### multiples of TILE_SIZE
    padL = (pixSW.tilex-np.floor(pixSW.tilex))*TILE_SIZE
    padR = (np.ceil(pixNE.tilex)-pixNE.tilex)*TILE_SIZE
    # print padL,padR
    
    #padR = (pixNE.tilex-np.floor(pixNE.tilex))*TILE_SIZE
    #padL = (np.ceil(pixSW.tilex)-pixSW.tilex)*TILE_SIZE
    #print padL,padR
    
    ### Need to shift the image padding for the
    ### coords to come out right.  The expression 
    ### below works empirically for my test case, but not sure why.
    dd = (padL-padR)/2; padR+=dd; padL-=dd
    
    padT = (pixNE.tiley-np.floor(pixNE.tiley))*TILE_SIZE
    padB = (np.ceil(pixSW.tiley)-pixSW.tiley)*TILE_SIZE
    
    dx = pixNE.x-pixSW.x
    dy = pixSW.y-pixNE.y
    
    pixPerDeg = ysize/(llNE[1]-llSW[1])
    pixRatio = m.pixels_per_lon_degree[zoomLevel]/pixPerDeg
    
    #data_copy = congrid(data,(dy,dx))
    #data_copy = data
    #print xsize-dy
    
    fix_LR = (xsize-dy)/2.
    if fix_LR == np.round(fix_LR):
        padL+=fix_LR
        padR+=fix_LR
    else:
        padL+=fix_LR

    fix_TB = (ysize-dx)/2.
    if fix_TB == np.round(fix_TB):
        padT+=fix_TB
        padB+=fix_TB
    else:
        padT+=fix_TB
        
    # padT+=(ysize-dx)/2.
    # padB+=(ysize-dx)/2.
    dx, dy = ysize, xsize
    
    #print pixNE.tilex-pixSW.tilex, pixNE.tiley-pixSW.tiley, padL, padR, padT, padB, xsize, ysize, dx, dy
    
    #data_copy.resize((int(ysize*pixRatio),int(xsize*pixRatio)))
    #data_copy.resize((dy,dx))
    fullx = padL+padR+dx
    fully = padT+padB+dy
    full_images = []
    for ch in range(len(data)):
        full_image = np.zeros((fully,fullx))    
        if (np.int(padT) == 0) | (np.int(padR) == 0):
            if (np.int(padT)+np.int(padR)) == 0:
                full_image[padB:, padR:] = data[ch]*1
            if (np.int(padT) == 0) & (np.int(padR) != 0):
                full_image[padB:, padL:-padR] = data[ch]*1
            if (np.int(padT) != 0) & (np.int(padR) == 0):
                full_image[padB:-padT, padL:] = data[ch]*1
        else:
            full_image[padB:-padT, padL:-padR] = data[ch]*1
        
        full_images.append(full_image)
        
    # print pixRatio, dx/xsize, fullx/256., fully/256.
    
    NX = (padL+padR+dx)*1./TILE_SIZE
    NY = (padT+padB+dy)*1./TILE_SIZE
    
    tileX0 = int(pixSW.tilex)
    tileY0 = int(pixNE.tiley)
        
    for i in np.arange(NX):
        for j in np.arange(NY):
            #i,j = 0,0
            subs = []
            for ch in range(len(data)):
                sub = full_images[ch][fully-(j+1)*TILE_SIZE:fully-j*TILE_SIZE,
                             i*TILE_SIZE:(i+1)*TILE_SIZE]
                subs.append(sub)
                
            if (sub.shape[0] == 0) | (sub.shape[1] == 0):
                continue
                
            outfile = outPath+'%s_%d_%d_%d.png' %(tileroot,
                            tileX0+i,tileY0+j,zoomLevel)
            
            if len(subs) == 1:
                subim = data2image(sub, zmin=zmin, zmax=zmax)
                subim.save(outfile)
            else:
                if rgb_params is not None:
                    luptonRGB(subs[0]*zmax[0], subs[1]*zmax[1], subs[2]*zmax[2], Q=rgb_params[0], alpha=rgb_params[1], m0=rgb_params[2], m1=1, shape=None, filename=outfile, ds9=None, verbose=False, rgb_clip=rgb_clip)
                else:
                    linearRGB(subs[0], subs[1], subs[2], shape=None, filename=outfile)
                    
            if verbose: 
                print('threedhst.gmap: %s' %(outfile))
            #print outfile
            
    return params

def makeOtherTiles(reference_image='ib6o23020_drz.fits', 
                   reference_ext = 1,
      other_image='../../ECDFS_DR1.FITS', 
      other_ext = 0, zmin=-0.1, zmax=1,
      outPath='../HTML/tiles/', 
      tileroot='VLA', 
      swarpmode='wait',
      rescale=True):
    """
makeOtherTiles
    
    Xray, radio etc.
    """
    import threedhst
    
    #### Swarp the `other_image` to the `reference_image` pixels
    sw = threedhst.sex.SWarp()
    sw._aXeDefaults()
    sw.swarpMatchImage(reference_image, extension=reference_ext)
    sw.swarpImage(reference_image+'[%d]' %(reference_ext),mode='wait')    
    ## Refine center position from SWarp's own output
    sw.swarpRecenter()                  
    sw.swarpImage(other_image+'[%d]' %(other_ext), mode=swarpmode)
    
    co = pyfits.open('coadd.fits')
    im = pyfits.open(reference_image)
    im[1].data = co[0].data
    im.writeto('other_match.fits', clobber=True)
    
    #### Size of reference image in arcsec
    Xasec = sw.NX*np.float(sw.options['PIXEL_SCALE'])
    Yasec = sw.NY*np.float(sw.options['PIXEL_SCALE'])
    
    m = threedhst.gmap.MercatorProjection()
    aperpix = 1./np.array(m.pixels_per_lon_degree)*3600
    
    ########### Prepare map tiles for different zoom levels:
    ########### 0.07 x [1,2,4,8] arcsec/pix
    for aper in range(13,17):
        ### base image
        
        threedhst.showMessage("%s\nMap tiles, zoom level: %10.6f arcsec/pix"
                              %(other_image, aperpix[aper]))
                
        sw.options['IMAGE_SIZE']='0'
        sw.options['PIXELSCALE_TYPE']='MANUAL'
        sw.options['PIXEL_SCALE']='%10.6f' %aperpix[aper]
        
        #### Direct
        sw.swarpImage('other_match.fits[1]', mode=swarpmode)
        im = pyfits.open('coadd.fits')
        #im[0].data *= 4**(aper-15)
        if rescale:
            if aper <= 14:
                im[0].data /= 3**np.abs(aper-15)
            if aper == 16:
                im[0].data *= 4
            
        im.writeto('scale.fits', clobber=True)
        mapParamsD = threedhst.gmap.makeGMapTiles(fitsfile='scale.fits',
                                                 outPath=outPath,
                                                 tileroot=tileroot,
                                                 extension=0,
                                                 zmin=zmin, zmax=zmax)
        os.remove('scale.fits')
        
def makeAllTiles(ROOT_DIRECT, ROOT_GRISM, zmin=-0.1, zmax=1, verbose=False,
                 PARAM_ONLY=False):
    """
mapParams = makeAllTiles(ROOT_DIRECT, ROOT_GRISM, zmin=-0.1, zmax=1, 
                         PARAM_ONLY=False)
    """
    import threedhst
    
    threedhst.showMessage("""
Make all of the google map tiles.  This involves running SWarp
multiple times for each zoom level and image combination.
""")

    #### Make XML file of the catalog, coordinates and ID number
    threedhst.gmap.makeCatXML(catFile=ROOT_GRISM+'_drz.cat',
                              xmlFile='../HTML/'+ROOT_GRISM+'.xml',
                              SPCFile = threedhst.currentRun['SPC'])
    
    threedhst.gmap.makeCirclePNG(outfile='../HTML/scripts/circle.php')            
    
    #### Need to swarp the drz images to the correct pixel scale for 
    #### the pixel size of the google map tiles
    m = threedhst.gmap.MercatorProjection()
    aperpix = 1./np.array(m.pixels_per_lon_degree)*3600
    sw = threedhst.sex.SWarp()
    sw._aXeDefaults()
    # sw.swarpMatchImage(ROOT_GRISM.lower()+'_drz.fits')
    sw.swarpMatchImage(threedhst.options['DIRECT_MOSAIC'])
    
    ########### Prepare map tiles for different zoom levels:
    ########### 0.07 x [1,2,4,8] arcsec/pix
    aper_list = list(range(13,17))
    if PARAM_ONLY:
        aper_list = [16]
        
    for aper in aper_list:
        ### base image
        
        threedhst.showMessage("Map tiles, zoom level: %10.6f arcsec/pix"
                              %aperpix[aper])
        
        sw.options['IMAGE_SIZE']='0'
        sw.options['PIXELSCALE_TYPE']='MANUAL'
        sw.options['PIXEL_SCALE']='%10.6f' %aperpix[aper]
        
        zmi = zmin
        zma = zmax
                
        #### Direct
        # sw.swarpImage(ROOT_DIRECT.lower()+'_drz.fits[1]', mode='wait')
        sw.swarpImage(threedhst.options['DIRECT_MOSAIC']+'[1]', mode='wait')
        im = pyfits.open('coadd.fits')
        #im[0].data *= 4**(aper-15)
        if aper <= 14:
            im[0].data /= 2
        im.writeto('scale.fits', clobber=True)
        mapParamsD = threedhst.gmap.makeGMapTiles(fitsfile='scale.fits',
                                                 outPath='../HTML/tiles/',
                                                 tileroot=ROOT_GRISM+'_d',
                                                 extension=0,
                                                 zmin=zmi, zmax=zma,
                                                 verbose=verbose)
        
        #### Get map parameters from high-resolution image
        if (aper == 16):
            mapParams=mapParamsD.copy()
            if PARAM_ONLY:
                return mapParams
                
        if aper == 16:
            zmi*=0.2
            zma*=0.2
            
        #### Grism
        sw.swarpImage(ROOT_GRISM+'_drz.fits[1]', mode='wait')
        im = pyfits.open('coadd.fits')
        #im[0].data *= 4**(aper-15)
        if aper <= 14:
            im[0].data /= 2
        im.writeto('scale.fits', clobber=True)
        mapParamsG = threedhst.gmap.makeGMapTiles(fitsfile='scale.fits',
                                                 outPath='../HTML/tiles/',
                                                 tileroot=ROOT_GRISM+'_g',
                                                 extension=0,
                                                 zmin=zmi, zmax=zma,
                                                 verbose=verbose)
        
        #### Model
        sw.swarpImage(ROOT_GRISM+'CONT_drz.fits[1]', mode='wait')
        im = pyfits.open('coadd.fits')
        #im[0].data *= 4**(aper-15)
        if aper <= 14:
            im[0].data /= 2
        im.writeto('scale.fits', clobber=True)
        mapParamsM = threedhst.gmap.makeGMapTiles(fitsfile='scale.fits',
                                                 outPath='../HTML/tiles/',
                                                 tileroot=ROOT_GRISM+'_m',
                                                 extension=0,
                                                 zmin=zmi, zmax=zma,
                                                 verbose=verbose)
        
    # #### direct tiles
    # mapParamsD = threedhst.gmap.makeGMapTiles(fitsfile=
    #                                          ROOT_GRISM.lower()+'_drz.fits',
    #                                          outPath='../HTML/tiles/',
    #                                          tileroot=ROOT_GRISM+'_d')
    # #### grism tiles
    # mapParamsG = threedhst.gmap.makeGMapTiles(fitsfile=
    #                                          ROOT_GRISM.lower()+'_drz.fits',
    #                                          outPath='../HTML/tiles/',
    #                                          tileroot=ROOT_GRISM+'_g')
    # #### model tiles                           
    # mapParamsM = threedhst.gmap.makeGMapTiles(fitsfile=
    #                                          ROOT_GRISM.lower()+'CONT_drz.fits',
    #                                          outPath='../HTML/tiles/',
    #                                          tileroot=ROOT_GRISM+'_m')
    
    #### Done making the map tiles
    threedhst.currentRun['step'] = 'MAKE_GMAP_TILES'
    try:
        os.remove('scale.fits')
        os.remove('coadd.fits')
        os.remove('coadd.weight.fits')
    except:
        pass
        
    return mapParams
        
def makeCirclePNG(outfile='circle.php'):
    """
makeCirclePNG(outfile=None)
    
    Simple icon for overlay on Google Map.  Optionally includes
    an object label on top of a circle.
    
    Example:
    >>> makeCirclePNG(outfile='~/Sites/circle.php')
    
    Then point to 'http://localhost/~[user]/circle.php?id=100' in a web browser 
    (with web sharing and PHP enabled, if viewing locally)
    """
    
    PHPstring = """<?php
    $string = $_GET['id'];
    // Create a blank image.
    $image = imagecreatetruecolor(30, 25);
    // Make the background transparent
    $black = imagecolorallocate($im, 0, 0, 0);
    imagecolortransparent($image, $black);
    // Choose a color for the ellipse.
    $green = imagecolorallocate($image, 0, 255, 0);
    // Draw the ellipse.
    imageellipse($image, 14.5, 14.5, 10, 10, $green);
    // Add the ID number
    $px     = (imagesx($image) - 4.5 * strlen($string)) / 2;
    imagestring($image, 0, $px, 1.5, $string, $green);
    // Output the image.
    header("Content-type: image/png");
    imagepng($image);
    ?>
    """
    
    if outfile:
        fp = open(outfile,'w')
        fp.write(PHPstring)
        fp.close()
    else:
        print(PHPstring)


def data2image(data,zmin=-0.1,zmax=0.5):
    """
data2image(data,zmin=-0.1,zmax=0.5)
    
    Take a 2D data array and send it to a PIL Image 
    object after (linear) z-scaling.
    
    Parts taken from the `fitsimage` class in wcs2kml.
    
    """ 
    from PIL import Image
    # array sizes
    xsize = data.shape[1]
    ysize = data.shape[0]
    # copy of data
    fits_data = data*1.
    fits_data = np.where(fits_data > zmin, fits_data, zmin)
    fits_data = np.where(fits_data < zmax, fits_data, zmax)
    scaled_data = (fits_data - zmin) * (255.0 / (zmax - zmin)) + 0.5
    # convert to 8 bit unsigned int
    scaled_data = scaled_data.astype("B")
    # create the image
    image = Image.frombuffer("L", (xsize, ysize), scaled_data, "raw", "L", 0, 0)
    return image

def linearRGB(R, G, B, filename='junk.png', shape=None):
    """
    Use linear scaling rather than more complicated lupton scaling.  
    Input images are assumed to already be scaled 0..1
    """
    import Image
    
    R = np.clip(R, 0, 1)
    G = np.clip(G, 0, 1)
    B = np.clip(B, 0, 1)
    
    im = Image.merge('RGB', (Image.fromarray((R[::-1,:]*255).astype('int8'), mode='L'), Image.fromarray((G[::-1,:]*255).astype('int8'), mode='L'), Image.fromarray((B[::-1,:]*255).astype('int8'), mode='L')))
    
    if shape is not None:
        im = im.resize(shape)
    
    im.save(filename)
    
def luptonRGB(imr, img, imb, Q=5, alpha=3, m0=-0.05, m1=1, shape=None, filename='junk.png', ds9=None, verbose=False, rgb_clip=True):
    """
    Make a 3 color image scaled with the color clipping and 
    asinh scaling from Lupton et al. (2004)
    """   
    import Image
    
    I = (imr+img+imb-3*m0)/3.
    fI = np.arcsinh(alpha*Q*I)/Q
    M = m0 + np.sinh(Q*1.)/(alpha*Q)
    #ds9.v(fI, vmin=0, vmax=1)
    if verbose:
        print('min, max = %f, %f' %(m0, M))
        
    fI[I < m0] = 0
    R = np.maximum(imr-m0, 0)*fI/I
    G = np.maximum(img-m0, 0)*fI/I
    B = np.maximum(imb-m0, 0)*fI/I
        
    min_RGB = np.minimum(np.minimum(R,G),B)
    zero = min_RGB < 0
    zero = fI < 0
    R[zero] = 0.
    G[zero] = 0.
    B[zero] = 0.
    
    R[R < 0] = 0
    G[G < 0] = 0
    B[B < 0] = 0
    
    max_RGB = np.maximum(np.maximum(R,G),B)
    if rgb_clip:
        clip = max_RGB > 1
        R[clip] = R[clip]/max_RGB[clip]
        G[clip] = G[clip]/max_RGB[clip]
        B[clip] = B[clip]/max_RGB[clip]
    else:
        R[R > 1] = 1.
        G[G > 1] = 1.
        B[B > 1] = 1.

    if ds9 is not None:
        #ds9.set('rgb True')
        v1=1
        ds9.set('rgb lock colorbar')
        ds9.set('rgb red'); ds9.v(R, vmin=0, vmax=v1); ds9.set('scale linear')
        ds9.set('rgb green'); ds9.v(G, vmin=0, vmax=v1); ds9.set('scale linear')
        ds9.set('rgb blue'); ds9.v(B, vmin=0, vmax=v1); ds9.set('scale linear')
        return True
        
    #rgb = np.array([R,G,B]).T
    #im = Image.merge('RGB', (Image.fromarray((R[::-1,:]*255).astype('uint8')), Image.fromarray((G[::-1,:]*255).astype('uint8')), Image.fromarray((B[::-1,:]*255).astype('uint8'))))
    im = Image.merge('RGB', (Image.fromarray((R[::-1,:]*255).astype('int8'), mode='L'), Image.fromarray((G[::-1,:]*255).astype('int8'), mode='L'), Image.fromarray((B[::-1,:]*255).astype('int8'), mode='L')))
    
    if shape is not None:
        im = im.resize(shape)
    
    im.save(filename)
  
def radec2latlon(radec):
    """
radec2latlon(radec)
    
    Convert R.A./Dec to Lat./Lon.
    
    """
    #latlon = np.zeros(2.)
    latlon = np.array([radec[1],360.-radec[0]])
    #latlon = np.array([radec[1],(360.-radec[0])/np.cos(radec[1]/360.*2*np.pi)])
    #latlon = np.array([radec[1],radec[0]])
    return latlon
    
def makeCatXML(catFile=None, xmlFile=None, SPCFile=None):
    """
makeCatXML(catFile=None,xmlFile=None)
    
    Make XML file suitable for reading into a Google map from 
    a SExtractor catalog
    """
    import threedhst
    #catFile='ib3721050_drz.cat'
    if not catFile:
        print('makeCatXML: no `catFile` specified')
        return None
    ### Read the catalog and get id, ra, and dec columns
    cat = threedhst.sex.mySexCat(catFile)
    colID = cat.columns[cat.searchcol('NUMBER')].entry
    colRA = cat.columns[cat.searchcol('X_WORLD')].entry
    colDEC = cat.columns[cat.searchcol('Y_WORLD')].entry
    
    if SPCFile is None:
        spec_list = np.array(np.cast[int](cat.NUMBER))
    else:
        spec_list = SPCFile._ext_map
        
    ### Mag column, sorted
    col = 'MAG_AUTO'
    for col in cat.column_names:
        if col.startswith('MAG_F'):
            break
    colMag = cat.columns[cat.searchcol(col)].entry
    
    mag = np.array(np.cast[float](colMag))
    sort_idx = mag.argsort()
    
    ### Make XML string
    xmlString = '<markers>'
    for j in range(cat.nrows):
        i = sort_idx[j]
        if int(cat.NUMBER[i]) in spec_list:
            xmlString+='<marker id="%s" ra="%s" dec="%s" mag="%s"/>' %(colID[i],
                      colRA[i], colDEC[i], colMag[i])
                      
    xmlString+='</markers>'
    ### Print to output file
    if xmlFile:
        fp = open(xmlFile,'w')
        fp.write(xmlString)
        fp.close()
    

class Point():
    """
Stores a simple (x,y) point.  It is used for storing x/y pixels.
    
Attributes:
    x: An int for a x value.
    y: An int for a y value.
        
http://code.google.com/p/google-ajax-examples/source/browse/trunk/nonjslocalsearch/localSearch.py
    
    """
    def __init__(self, x, y):
        self.x = x
        self.y = y
        self.tilex = x*1./TILE_SIZE
        self.tiley = y*1./TILE_SIZE
        
    def ToString(self):
        return '(%s, %s)' % (self.x, self.y)
    
    def Equals(self, other):
        if other is None :
            return false
        else:
            return (other.x == self.x and other.y == self.y)

class MercatorProjection():
  """
MercatorProjection

Calculates map zoom levels based on bounds or map points.

This class contains functions that are required for calculating the zoom  
level for a point or a group of points on a static map.  Usually the Maps API 
does the zoom for you when you specify a point, but not on static maps.

Attributes:
    pixels_per_lon_degree: A list for the number of pixels per longitude 
      degree for each zoom.
    pixels_per_lon_radian: A list for the number of pixels per longitude
      radian for each zoom.
    pixel_origo: List of number of x,y pixels per zoom.
    pixel_range: The range of pixels per zoom.
    pixels: Number of pixels per zoom.
    zoom_levels: A list of numbers representing each zoom level to test.
    
http://code.google.com/p/google-ajax-examples/source/browse/trunk/nonjslocalsearch/localSearch.py
    
  """
  def __init__(self, zoom_levels=18):
    self.pixels_per_lon_degree = []
    self.pixels_per_lon_radian = []
    self.pixel_origo = []
    self.pixel_range = []
    self.pixels = TILE_SIZE
    zoom_levels = list(range(0, zoom_levels))
    for z in zoom_levels:
      origin = self.pixels / 2
      self.pixels_per_lon_degree.append(self.pixels / 360)
      self.pixels_per_lon_radian.append(self.pixels / (2 * np.pi))
      self.pixel_origo.append(Point(origin, origin))
      self.pixel_range.append(self.pixels)
      self.pixels = self.pixels * 2
    
  def CalcWrapWidth(self, zoom):
    return self.pixel_range[zoom]
    
  def FromLatLngToPixel(self, lat_lng, zoom):
    """Given lat/lng and a zoom level, returns a Point instance.

    This method takes in a lat/lng and a _test_ zoom level and based on that it 
    calculates at what pixel this lat/lng would be on the map given the zoom 
    level.  This method is used by CalculateBoundsZoomLevel to see if this 
    _test_ zoom level will allow us to fit these bounds in our given map size.

    Args:
      lat_lng: A list of a lat/lng point [lat, lng]
      zoom: A list containing the width/height in pixels of the map.

    Returns:
      A Point instance in pixels.
    """
    o = self.pixel_origo[zoom]
    x = round(o.x + lat_lng[1] * self.pixels_per_lon_degree[zoom])
    siny = Bound(np.sin(DegreesToRadians(lat_lng[0])), 
        -0.9999, 0.9999)
    y = round(o.y + 0.5 * np.log((1 + siny) / 
        (1 - siny)) * -self.pixels_per_lon_radian[zoom])
    return Point(x, y)

  def CalculateBoundsZoomLevel(self, bounds, view_size):
    """Given lat/lng bounds, returns map zoom level.

    This method is used to take in a bounding box (southwest and northeast 
    bounds of the map view we want) and a map size and it will return us a zoom 
    level for our map.  We use this because if we take the bottom left and 
    upper right on the map we want to show, and calculate what pixels they 
    would be on the map for a given zoom level, then we can see how many pixels 
    it will take to display the map at this zoom level.  If our map size is 
    within this many pixels, then we have the right zoom level.

    Args:
      bounds: A list of length 2, each holding a list of length 2. It holds
        the southwest and northeast lat/lng bounds of a map.  It should look 
        like this: [[southwestLat, southwestLat], [northeastLat, northeastLng]]
      view_size: A list containing the width/height in pixels of the map.

    Returns:
      An int zoom level.
    """
    zmax = 18
    zmin = 0
    bottom_left = bounds[0]
    top_right = bounds[1]
    backwards_range = list(range(zmin, zmax))
    backwards_range.reverse()
    for z in backwards_range:
      bottom_left_pixel = self.FromLatLngToPixel(bottom_left, z)
      top_right_pixel = self.FromLatLngToPixel(top_right, z)
      if bottom_left_pixel.x > top_right_pixel.x :
        bottom_left_pixel.x -= self.CalcWrapWidth(z)
      if abs(top_right_pixel.x - bottom_left_pixel.x) <= view_size[0] \
          and abs(top_right_pixel.y - bottom_left_pixel.y) <= view_size[1] :
        return z
    return 0

def Bound(value, opt_min, opt_max):
    """
    Returns value if in min/max, otherwise returns the min/max.
    
  Args:
    value: The value in question.
    opt_min: The minimum the value can be.
    opt_max: The maximum the value can be.
    
  Returns:
    An int that is either the value passed in or the min or the max.
    
  http://code.google.com/p/google-ajax-examples/source/browse/trunk/nonjslocalsearch/localSearch.py
    """
    if opt_min is not None:
        value = max(value, opt_min)
    if opt_max is not None:
        value = min(value, opt_max)
    return value

def DegreesToRadians(deg):
    """
    http://code.google.com/p/google-ajax-examples/source/browse/trunk/nonjslocalsearch/localSearch.py
    """
    return deg * (np.pi / 180)

def congrid(a, newdims, method='linear', centre=False, minusone=False):
    '''Arbitrary resampling of source array to new dimension sizes.
    Currently only supports maintaining the same number of dimensions.
    To use 1-D arrays, first promote them to shape (x,1).
    
    Uses the same parameters and creates the same co-ordinate lookup points
    as IDL''s congrid routine, which apparently originally came from a VAX/VMS
    routine of the same name.

    method:
    neighbour - closest value from original data
    nearest and linear - uses n x 1-D interpolations using
                         scipy.interpolate.interp1d
    (see Numerical Recipes for validity of use of n 1-D interpolations)
    spline - uses ndimage.map_coordinates

    centre:
    True - interpolation points are at the centres of the bins
    False - points are at the front edge of the bin

    minusone:
    For example- inarray.shape = (i,j) & new dimensions = (x,y)
    False - inarray is resampled by factors of (i/x) * (j/y)
    True - inarray is resampled by(i-1)/(x-1) * (j-1)/(y-1)
    This prevents extrapolation one element beyond bounds of input array.
    '''
    import scipy.interpolate
    import scipy.ndimage
    
    if not a.dtype in [n.float64, n.float32]:
        a = n.cast[float](a)

    m1 = n.cast[int](minusone)
    ofs = n.cast[int](centre) * 0.5
    old = n.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print("[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions.")
        return None
    newdims = n.asarray( newdims, dtype=float )
    dimlist = []

    if method == 'neighbour':
        for i in range( ndims ):
            base = n.indices(newdims)[i]
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        cd = n.array( dimlist ).round().astype(int)
        newa = a[list( cd )]
        return newa

    elif method in ['nearest','linear']:
        # calculate new dims
        for i in range( ndims ):
            base = n.arange( newdims[i] )
            dimlist.append( (old[i] - m1) / (newdims[i] - m1) \
                            * (base + ofs) - ofs )
        # specify old dims
        olddims = [n.arange(i, dtype = n.float) for i in list( a.shape )]

        # first interpolation - for ndims = any
        mint = scipy.interpolate.interp1d( olddims[-1], a, kind=method )
        newa = mint( dimlist[-1] )

        trorder = [ndims - 1] + list(range( ndims - 1))
        for i in range( ndims - 2, -1, -1 ):
            newa = newa.transpose( trorder )

            mint = scipy.interpolate.interp1d( olddims[i], newa, kind=method )
            newa = mint( dimlist[i] )

        if ndims > 1:
            # need one more transpose to return to original dimensions
            newa = newa.transpose( trorder )

        return newa
    elif method in ['spline']:
        oslices = [ slice(0,j) for j in old ]
        oldcoords = n.ogrid[oslices]
        nslices = [ slice(0,j) for j in list(newdims) ]
        newcoords = n.mgrid[nslices]

        newcoords_dims = list(range(n.rank(newcoords)))
        #make first index last
        newcoords_dims.append(newcoords_dims.pop(0))
        newcoords_tr = newcoords.transpose(newcoords_dims)
        # makes a view that affects newcoords

        newcoords_tr += ofs

        deltas = (n.asarray(old) - m1) / (newdims - m1)
        newcoords_tr *= deltas

        newcoords_tr -= ofs

        newa = scipy.ndimage.map_coordinates(a, newcoords)
        return newa
    else:
        print("Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported.")
        return None
        
#
def parseImageString(IMAGE_STRING="test.fits[1]*1.", default_extension=1):
    """
    The input image string is image.fits[extension]*scale.  
    
    Output is a tuple with (image.fits, extension, scale).
    """
    input_list = IMAGE_STRING.split(',')
    image, extension, scale = [], [], []
    for input in input_list:
        sp_scl = input.split('*')
        if len(sp_scl) == 2:
            scale.append(float(sp_scl[1]))
        else:
            scale.append(1.)

        sp_ext = sp_scl[0].split('[')
        if len(sp_ext) == 2:
            extension.append(int(sp_ext[1][:-1]))
        else:
            extension.append(default_extension)

        image.append(sp_ext[0])
    
    return image, extension, scale
    
def makeImageMap(FITS_IMAGES, extension=1, zmin=-0.1, zmax=1, verbose=True,
                 path=os.getenv('HOME')+'/Sites/FITS/', tileroot='tile',
                 aper_list=[15], polyregions=None, rgb_params=(5, 3, -0.05),
                 invert=False, rgb_clip=True):
    """
    Make a google map viewer for a FITS image.
    
    FITS_IMAGES is a [list] of images that are used to generate the map.  
    The first image is used to define the extent of the output map.
    
    Entries in FITS_IMAGES can be simple filenames and can also contain
    extensions and scale values, specified like
    
        FITS_IMAGES = ['image1.fits[1]*10']
    
    The list entries can also be comma-separated list of three images 
    (R,G,B) that will be used to generate a 3-color image following 
    Lupton et al. (2004):
    
        FITS_IMAGES = ['red.fits[1]*10,green.fits[0]*2,blue.fits[0]*1.']
      
    """
    import threedhst    
    
    #### Make into lists if not already
    #if FITS_IMAGES.__class__ == ''.__class__:
    if isinstance(FITS_IMAGES, str):
        FITS_IMAGES = [FITS_IMAGES]
    
    roots = tileroot
    #if tileroot.__class__ == ''.__class__:
    if isinstance(tileroot, str):
        if len(FITS_IMAGES) > 1:
            roots = []
            for i,im in enumerate(FITS_IMAGES):
                roots.append(tileroot+'%d' %(i))
        else:
            roots = [tileroot]
    tileroot = roots
    
    print(tileroot)
    
    #### Need to swarp the drz images to the correct pixel scale for 
    #### the pixel size of the google map tiles
    m = threedhst.gmap.MercatorProjection()
    aperpix = 1./np.array(m.pixels_per_lon_degree)*3600
    sw = threedhst.sex.SWarp()
    sw._aXeDefaults()
    
    im0, ext0, scale0 = parseImageString(FITS_IMAGES[0], default_extension=extension)
    sw.swarpMatchImage(im0[0], extension=ext0[0], verbose=verbose)
    
    ########### Prepare map tiles for different zoom levels:
    ########### 0.07 x [1,2,4,8] arcsec/pix
    #aper_list = range(12,18)

    NX, NY = sw.options['IMAGE_SIZE'].split(',')
    NX, NY = int(NX), int(NY)
    NATIVE_SCALE = float(sw.options['PIXEL_SCALE'])
    
    for aper in aper_list:
        ### base image
        
        threedhst.showMessage("Map tiles, zoom level: %10.6f arcsec/pix"
                              %aperpix[aper])
        
        
        sw.options['IMAGE_SIZE']='%d,%d' %(np.round(NX*NATIVE_SCALE/aperpix[aper]),
                                           np.round(NY*NATIVE_SCALE/aperpix[aper]))
        sw.options['PIXELSCALE_TYPE']='MANUAL'
        sw.options['PIXEL_SCALE']='%10.6f' %aperpix[aper]
        
        for i,im in enumerate(FITS_IMAGES):
            imi, exti, scalei = parseImageString(FITS_IMAGES[i],
                                    default_extension=extension)
            
            print(imi, exti, scalei)
                
            #### Direct
            # sw.swarpImage(ROOT_DIRECT.lower()+'_drz.fits[1]', mode='wait')
            if len(imi) == 3:
                ### RGB
                channels = ['r','g','b']
                for ch in range(3):
                    sw.swarpImage(imi[ch]+'[%0d]' %(exti[ch]), mode='wait')
                
                    im = pyfits.open('coadd.fits')
                    #if aper <= 15:
                    im[0].data /= 4.**(16-aper)
                    print('Scale: %f' %(4.**(16-aper)))
                    
                    im.writeto('ch_%s.fits' %(channels[ch]), clobber=True)
                
                fitsfile='ch_r.fits,ch_g.fits,ch_b.fits'
                zmi = None
                zma = np.array(scalei)*(0.06/aperpix[aper])**1.5 ## soften a bit
                
            else:
                sw.swarpImage(imi[0]+'[%0d]' %(exti[0]), mode='wait')
                im = pyfits.open('coadd.fits')
                if aper <= 14:
                    im[0].data /= 4
                
                if invert:
                    im[0].data = 0-im[0].data
                    zma = -zmin/scalei[0]
                    zmi = -zmax/scalei[0]
                else:
                    zmi = zmin/scalei[0]
                    zma = zmax/scalei[0]
                    
                im.writeto('scale.fits', clobber=True)
                fitsfile='scale.fits'
                
            mapParamsD = threedhst.gmap.makeGMapTiles(fitsfile=fitsfile,
                                                     outPath=path+'tiles/',
                                                     tileroot=tileroot[i],
                                                     extension=0,
                                                     zmin=zmi, zmax=zma,
                                                     rgb_params=rgb_params, rgb_clip=rgb_clip,
                                                     verbose=verbose)
        
        #### Get map parameters from high-resolution image
        if (aper == aper_list[-1]):
            mapParams=mapParamsD.copy()
    
    mapParams['ZOOM_RANGE'] = [np.min(aper_list), np.max(aper_list)]
    
    threedhst.gmap.makeMapHTML(FITS_IMAGES, mapParams, output=path+'map.html', polyregions=polyregions)
    threedhst.plotting.makeCSS(path=path+'scripts/', title_size=9)
    threedhst.gmap.makeJavascript(path=path+'scripts/', tileroot=tileroot, mapParams=mapParams)
    threedhst.gmap.makePolygons(mapParams, polyregions=polyregions, path=path+'scripts/')
    
    ####
    files = glob.glob('threedhst_auto-*swarp')
    files.extend(['scale.fits', 'coadd.fits', 'coadd.weight.fits', 'ch_r.fits', 'ch_g.fits', 'ch_b.fits'])
    for file in files:
        try:
            os.remove(file)
        except:
            pass
        
#
def makePolygons(mapParams=None, polyregions=None, path='./', color='#00aa00', alpha=0.5, linewidth=2):
    """
    Translate ds9 polygons to Javascript for the map.
    """
    if polyregions is None:
        return False
    
    center = np.cast[float](mapParams['LLCENTER'])
    
    label = """
<?php
        $string = $_GET['text'];
    	$font_size = 5;
        // Create a blank image.
        $image = imagecreatetruecolor(ImageFontWidth($font_size) * strlen($string), 16);
    	//$string = ImageFontHeight($font_size);
        // Make the background transparent
        $black = imagecolorallocate($im, 0, 0, 0);
        imagecolortransparent($image, $black);
        // Choose a color for the ellipse.
        $green = imagecolorallocate($image, 0, 255, 0);
        // Add the text
        $px     = (imagesx($image) - ImageFontWidth($font_size) * strlen($string)) / 2;
        imagestring($image, $font_size, $px, 1.5, $string, $green);
        // Output the image.
        header("Content-type: image/png");
        imagepng($image);
?>
    """
    
    fp = open(path+'/threedhst.js')
    js_lines = fp.readlines()
    js_lines.append("\nfunction addPolylines() {\n")
    fp.close()
    
    #files=polyregions.split(',')
    files = polyregions
    
    N = len(files)
    counter=1
    for i in range(N):
        fp = open(files[i])
        reg_lines = fp.readlines()
        fp.close()
        
        for line in reg_lines:
            if line.startswith('polygon'):
                poly = np.cast[float](line.split('(')[1].split(')')[0].split(','))
                NP = len(poly)/2
                polystr = "\nvar polyL%d = new GPolyline([\n" %(counter)
                list = list(range(NP))
                list.append(0)
                for j in list:
                    lat = poly[j*2+1]-center[0]
                    lng = ((360-poly[j*2])-center[1])*np.cos(poly[j*2+1]/360.*2*np.pi)
                    polystr += "  new GLatLng(%f, %f),\n" %(lat, lng)
                
                polystr = polystr[:-2]+"\n], \"%s\", %d, %f);\n" %(color, linewidth, alpha)
                polystr += "map.addOverlay(polyL%d);\nmarker_list.push(polyL%d);" %(counter, counter)
                js_lines.append(polystr)
                counter+=1
                example = """
                var polyline = new GPolyline([
                  new GLatLng(37.4419, -122.1419),
                  new GLatLng(37.4519, -122.1519)
                ], "#ff0000", 10);
                map.addOverlay(polyline);
                marker_list.push(polyline);
                """
            if line.strip('#').strip().startswith('text'):
                poly = np.cast[float](line.split('(')[1].split(')')[0].split(','))
                label = line.split('text={')[1].split('}')[0]
                
                j=0
                lat = poly[j*2+1]-center[0]
                lng = ((360-poly[j*2])-center[1])*np.cos(poly[j*2+1]/360.*2*np.pi)
                
                marker="""

var label_text = "%s";
var labelIcon = new GIcon();
labelIcon.iconSize = new GSize(9*label_text.length, 15);
labelIcon.iconAnchor = new GPoint(9*label_text.length/2, 7.5);

labelIcon.image = "scripts/label_marker.php?text="+label_text;
markerOptions = { icon:labelIcon };
var point = new GLatLng(%f,%f);
var marker = new GMarker(point, markerOptions);
// GEvent.addListener(marker, "click", function() {window.open('http://xxx/%s/yyy');});

marker_list.push(marker);
map.addOverlay(marker);            
                """ %(label, lat, lng, label)
                
                add_link = """
                GEvent.addListener(marker, "click", function() {
                    window.open("http://www.google.com","mywindow", "");
                });
                """
                js_lines.append(marker)
                
                counter+=1
                
                
    js_lines.append("\n}\n");
    
    js_fp = open(path+'/threedhst.js','w')
    js_fp.writelines(js_lines)
    js_fp.close()
    
def makeMapHTML(FITS_IMAGES, mapParams, output='./HTML/index.html', title=None, polyregions=None):
    """

    """
    import os
    from socket import gethostname as hostname
    
    if hostname().startswith('uni'):
        GMAPS_KEY = 'ABQIAAAAzSrfHr_4F2D2YfSuYQD2ZBRWJmdnNuPtXxsK1b3ql6FJMkf8bxT-OiTDeKiGIrffKkTBi-in1FVtrw'
    else:
        # localhost
        GMAPS_KEY = 'ABQIAAAA1XbMiDxx_BTCY2_FkPh06RR20YmIEbERyaW5EQEiVNF0mpNGfBSRb_rzgcy5bqzSaTV8cyi2Bgsx3g'
    
    #### New general key for all referers
    GMAPS_KEY = 'AIzaSyAMtpy9tGl609hpmt30fPKjSf0Jl34D0PM'
        
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
    <script src="http://maps.google.com/maps?file=api&amp;v=2&amp;key=%s" type="text/javascript"></script> 
    """ %(GMAPS_KEY)]
        
    #### Script for the Google map
    llSW = mapParams['LLSW']
    llNE = mapParams['LLNE']
    center = mapParams['LLCENTER']
    lng_offset = mapParams['LNG_OFFSET']
    
    if polyregions is not None:
        addpolylines="addPolylines()"
    else:
        addpolylines=""
        
    lines.append("""
    <script type="text/javascript"> 
    
    //////////// Global variables
    var map = 0;
    var centerLat = %f;
    var centerLng = %f;
    // var offset = %f;
    var offset = 0.0;
    var zoomLevel = %f;
    var root = "tile";
    
    var myIcon = new GIcon();
	myIcon.iconSize = new GSize(30, 25);
	myIcon.iconAnchor = new GPoint(14.5, 14.5);
	
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
            
            %s
            
            ///// Add the green circles around the catalog objects
            //plotXmlObjects();
        }
        //initialize_SED_column();
        //// Check for position in URL name
        var url_split = document.URL.split('#');
        if (url_split.length > 1) {
            var rd = url_split[1].split(',');
            if (rd.length == 2) {
                document.getElementById('raInput').value = rd[0];
                document.getElementById('decInput').value = rd[1];
                centerOnInput();
            }
        }
    }
    
    </script>
    """ %(center[0],center[1],lng_offset,mapParams['ZOOMLEVEL'], addpolylines))
    
    FITS_LIST = []
    for im in FITS_IMAGES:
        FITS_LIST.append(os.path.basename(im))
    
    #### HTML Body   
    lines.append("""
    
    <!-- Dummy text for search/replace insertion -->
    <!-- xxx -->
    
</head>
<body onload="initialize()" onunload="GUnload()">
    
    <div id="map"></div>
    
    <div id="title">
        %s
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
    
    """ %(', '.join(FITS_LIST)))
    
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
    
#
def makeJavascript(tileroot=['tile'], path="../HTML/scripts", mapParams=None):
    """
make_Javascript(path="../HTML/scripts")
    """
    
    if tileroot.__class__ == ''.__class__:
        tileroot = [tileroot]
    
    if mapParams == None:
        mapParams = {}
        mapParams['ZOOM_RANGE'] = [12,17]
    
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
        } else {
            layout = 0;
            horizontal_layout();
            $("#switchbox").text("=");
            $("#switchbox").css("cursor","e-resize");
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

    	$("#map").css("height",$(window).height()-60-11);
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

    	//clearRows();

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
    """)
    
    for i,root in enumerate(tileroot):
        fp.write("""
        // Direct image tiles
        CustomGetTileUrl%0d=function(a,b){
            return "tiles/%s_"+a.x+"_"+a.y+"_"+b+".png"
        }
        var tilelayers%0d = [new GTileLayer(copyrightCollection,%d,%d)];
        tilelayers%0d[0].getTileUrl = CustomGetTileUrl%0d;
        var customMap%0d = new GMapType(tilelayers%0d, 
               new GMercatorProjection(%d), "%s");
        map.addMapType(customMap%0d);
    """ %(i, root, i, mapParams['ZOOM_RANGE'][0], mapParams['ZOOM_RANGE'][1], i, i, i, i, mapParams['ZOOM_RANGE'][1]+1, root, i))
    
    fp.write("""
        // Can't remove all three for some reason
        map.removeMapType(G_NORMAL_MAP);
        map.removeMapType(G_HYBRID_MAP);
        map.removeMapType(G_SATELLITE_MAP);
        
        map.addControl(new GMapTypeControl());
        
        // Set map center
        map.setCenter(new GLatLng(0.0, offset), zoomLevel,
         customMap0);
		
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
    }""" )

    fp.write("""
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
        if (dsplit.length == 3) {
            dec = Math.abs(parseInt(dsplit[0])+
                Math.abs(parseInt(dsplit[1])/60.)+
                Math.abs(parseFloat(dsplit[2])/3600.));

            /// Don't know why, but need to do this twice
            dec = Math.abs(parseInt(dsplit[0]))+
            parseInt(dsplit[1])/60.+
            parseFloat(dsplit[2])/3600.;

            // if (parseFloat(dsplit[0]) < 0) {
            //     dec *= -1;
            // }
            if (dsplit[0][0] == "-") {
                dec *= -1;
            }            
        }

        recenter(ra,dec);
        latLng2raDec();
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
	var NSHOW = 25;""")
	
    fp.close()
