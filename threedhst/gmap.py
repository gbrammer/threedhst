"""
Convert Multidrizzle drz.fits (North up) to Google Maps tiles
"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import os
import numpy as np

# Specifies the size of the map (in pixels).
TILE_SIZE = 256
MAP_SIZE = [TILE_SIZE,TILE_SIZE]    
# This is the Maps API key for running on localhost:8080
MAP_KEY = 'ABQIAAAA1XbMiDxx_BTCY2_FkPh06RR20YmIEbERyaW5EQEiVNF0mpNGfBSRb' \
    '_rzgcy5bqzSaTV8cyi2Bgsx3g'

def makeGMapTiles(fitsfile=None,outPath=None,tileroot='direct', extension=1,
                  zmin=-0.1, zmax=1, verbose=False):
    """
makeGMapTiles(fitsfile=None,outPath=None,tileroot='direct', extension=1,
              zmin=-0.1, zmax=1)
    
    Make Google map tiles for an input FITS image, which is assumed to be
    North-up, East-left, like normal Multidrizzle output.
    """
    import pyfits
    import pywcs
    #import fitsimage
    import numpy as np
    
    if not fitsfile:
        fitsfile = 'ib3721050_drz.fits'
    if not outPath:
        outPath = '/tmp/'
    
    # print fitsfile, outPath
    
    ### Read the FITS file
    fi = pyfits.open(fitsfile)
    head = fi[extension].header
    data = fi[extension].data
    #data = np.fliplr(fi[1].data)
    xsize, ysize = data.shape #[1]
    
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
    data_copy = data
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
    full_image = np.zeros((fully,fullx))    
    full_image[padB:-padT, padL:-padR] = data_copy
    
    # print pixRatio, dx/xsize, fullx/256., fully/256.
    
    NX = (padL+padR+dx)*1./TILE_SIZE
    NY = (padT+padB+dy)*1./TILE_SIZE
    
    tileX0 = int(pixSW.tilex)
    tileY0 = int(pixNE.tiley)
        
    for i in range(NX):
        for j in range(NY):
            #i,j = 0,0
            sub = full_image[fully-(j+1)*TILE_SIZE:fully-j*TILE_SIZE,
                             i*TILE_SIZE:(i+1)*TILE_SIZE]
            subim = data2image(sub, zmin=zmin, zmax=zmax)
            outfile = outPath+'%s_%d_%d_%d.png' %(tileroot,
                            tileX0+i,tileY0+j,zoomLevel)
            subim.save(outfile)
            if verbose: 
                print 'threedhst.gmap: %s' %(outfile)
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
    import pyfits
    import numpy as np
    
    #### Swarp the `other_image` to the `reference_image` pixels
    sw = threedhst.sex.SWarp()
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
    
def makeAllTiles(ROOT_DIRECT, ROOT_GRISM, zmin=-0.1, zmax=1, verbose=False):
    """
mapParams = makeAllTiles(ROOT_DIRECT, ROOT_GRISM, zmin=-0.1, zmax=1)
    """
    import threedhst
    import pyfits
    import numpy as np
    
    threedhst.showMessage("""
Make all of the google map tiles.  This involves running SWarp
multiple times for each zoom level and image combination.
""")

    #### Make XML file of the catalog, coordinates and ID number
    threedhst.gmap.makeCatXML(catFile=ROOT_DIRECT.lower()+'_drz.cat',
                              xmlFile='../HTML/'+ROOT_DIRECT+'.xml')
    
    threedhst.gmap.makeCirclePNG(outfile='../HTML/scripts/circle.php')            
    
    #### Need to swarp the drz images to the correct pixel scale for 
    #### the pixel size of the google map tiles
    m = threedhst.gmap.MercatorProjection()
    aperpix = 1./np.array(m.pixels_per_lon_degree)*3600
    sw = threedhst.sex.SWarp()
    # sw.swarpMatchImage(ROOT_DIRECT.lower()+'_drz.fits')
    sw.swarpMatchImage(threedhst.options['DIRECT_MOSAIC'])
    
    ########### Prepare map tiles for different zoom levels:
    ########### 0.07 x [1,2,4,8] arcsec/pix
    for aper in range(13,17):
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
                                                 tileroot=ROOT_DIRECT+'_d',
                                                 extension=0,
                                                 zmin=zmi, zmax=zma,
                                                 verbose=verbose)
        
        #### Get map parameters from high-resolution image
        if (aper == 16):
            mapParams=mapParamsD.copy()
        
        if aper == 16:
            zmi*=0.2
            zma*=0.2
            
        #### Grism
        sw.swarpImage(ROOT_GRISM.lower()+'_drz.fits[1]', mode='wait')
        im = pyfits.open('coadd.fits')
        #im[0].data *= 4**(aper-15)
        if aper <= 14:
            im[0].data /= 2
        im.writeto('scale.fits', clobber=True)
        mapParamsG = threedhst.gmap.makeGMapTiles(fitsfile='scale.fits',
                                                 outPath='../HTML/tiles/',
                                                 tileroot=ROOT_DIRECT+'_g',
                                                 extension=0,
                                                 zmin=zmi, zmax=zma,
                                                 verbose=verbose)
        
        #### Model
        sw.swarpImage(ROOT_GRISM.lower()+'CONT_drz.fits[1]', mode='wait')
        im = pyfits.open('coadd.fits')
        #im[0].data *= 4**(aper-15)
        if aper <= 14:
            im[0].data /= 2
        im.writeto('scale.fits', clobber=True)
        mapParamsM = threedhst.gmap.makeGMapTiles(fitsfile='scale.fits',
                                                 outPath='../HTML/tiles/',
                                                 tileroot=ROOT_DIRECT+'_m',
                                                 extension=0,
                                                 zmin=zmi, zmax=zma,
                                                 verbose=verbose)
        
    # #### direct tiles
    # mapParamsD = threedhst.gmap.makeGMapTiles(fitsfile=
    #                                          ROOT_DIRECT.lower()+'_drz.fits',
    #                                          outPath='../HTML/tiles/',
    #                                          tileroot=ROOT_DIRECT+'_d')
    # #### grism tiles
    # mapParamsG = threedhst.gmap.makeGMapTiles(fitsfile=
    #                                          ROOT_GRISM.lower()+'_drz.fits',
    #                                          outPath='../HTML/tiles/',
    #                                          tileroot=ROOT_DIRECT+'_g')
    # #### model tiles                           
    # mapParamsM = threedhst.gmap.makeGMapTiles(fitsfile=
    #                                          ROOT_GRISM.lower()+'CONT_drz.fits',
    #                                          outPath='../HTML/tiles/',
    #                                          tileroot=ROOT_DIRECT+'_m')
    
    #### Done making the map tiles
    threedhst.currentRun['step'] = 'MAKE_GMAP_TILES'
    
    return mapParams
    
def makeMapHTML(llSW, llNE, lng_offset=90):
    """
makeHTML(llSW, llNE, lng_offset=90)

    Make webpage container for holding the Google map.
    """
    
    center = (llSW+llNE)/2.
    
    web = """<html> 
    <head> 
        <meta http-equiv="content-type" content="text/html; charset=UTF-8"/> 
        <title>Google Maps</title> 
        <script src="http://maps.google.com/maps?file=api&amp;v=3&amp;key=ABQIAAAA1XbMiDxx_BTCY2_FkPh06RR20YmIEbERyaW5EQEiVNF0mpNGfBSRb_rzgcy5bqzSaTV8cyi2Bgsx3g" type="text/javascript"></script> 

        <script type="text/javascript"> 

function initialize() {        
    if (GBrowserIsCompatible()) {
        var map = new GMap2(document.getElementById("map"));
        map.addControl(new GScaleControl());
        // map.addControl(new GSmallMapControl());
        // map.addControl(new GMapTypeControl());
        var copyright = new GCopyright(1,
             new GLatLngBounds(new GLatLng(%f,%f),
                               new GLatLng(%f,%f)),
                               14, "3D-HST");
        var copyrightCollection = new GCopyrightCollection('Map Data:');
        copyrightCollection.addCopyright(copyright);

        CustomGetTileUrl=function(a,b){
            return "direct_"+a.x+"_"+a.y+"_"+b+".jpg"
        }
        var tilelayers = [new GTileLayer(copyrightCollection,14,14)];
        tilelayers[0].getTileUrl = CustomGetTileUrl;
        var custommap = new GMapType(tilelayers, 
               new GMercatorProjection(15), "FITS");
        map.addMapType(custommap);
        map.setCenter(new GLatLng(%f,  %f), 14, custommap);
          
        plotXmlObjects(map);
    }
}

// Globals
var myIcon = new GIcon();
myIcon.iconSize = new GSize(40, 40);
myIcon.iconAnchor = new GPoint(19.5, 19.5);
var lats = [];
var lngs = [];
var ids = [];
var nObject = 0;

// Read objects from XML file and plot regions
function plotXmlObjects(map, centerLng, offset) {
    GDownloadUrl("cat.xml", function(data) {
        var xml = GXml.parse(data);
        var markers = xml.documentElement.getElementsByTagName("marker");
        nObject = markers.length;
        
        for (var i = 0; i < markers.length; i++) {
            // Read from XML
            var id = markers[i].getAttribute("id");
            var ra = markers[i].getAttribute("ra");
            var dec = markers[i].getAttribute("dec");
            var lat = dec;
            var lng = (360-ra)-centerLng+offset;
            lats.push(lat);
            lngs.push(lng);
            ids.push(id);
            
            // The marker
            myIcon.image = "circle.php?id="+id;
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
                window.location = '#i'+matchID;
            });
            map.addOverlay(marker);            
        }
    });
}
        </script>
        </head> 
      <body onload="initialize()" onunload="GUnload()"> 
        <div id="map" style="width: 300px; height: 300px"></div> 
      </body> 
    </html>
    """ %(llSW[0],llSW[1]-center[1]+lng_offset,
                  llNE[0],llNE[1]-center[1]+lng_offset,
                  center[0],lng_offset)
    
    outfile = '/Users/gbrammer/Sites/map/ASTR/map.html'
    outfile = './map.html'
    
    fp = open(outfile,'w')
    fp.write(web)
    fp.close()
    print outfile
    
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
        print PHPstring


def data2image(data,zmin=-0.1,zmax=0.5):
    """
data2image(data,zmin=-0.1,zmax=0.5)
    
    Take a 2D data array and send it to a PIL Image 
    object after (linear) z-scaling.
    
    Parts taken from the `fitsimage` class in wcs2kml.
    
    """ 
    from PIL import Image
    import numpy as np
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
    
def radec2latlon(radec):
    """
radec2latlon(radec)
    
    Convert R.A./Dec to Lat./Lon.
    
    """
    import numpy as np
    #latlon = np.zeros(2.)
    latlon = np.array([radec[1],360.-radec[0]])
    #latlon = np.array([radec[1],(360.-radec[0])/np.cos(radec[1]/360.*2*np.pi)])
    #latlon = np.array([radec[1],radec[0]])
    return latlon
    
def makeCatXML(catFile=None,xmlFile=None):
    """
makeCatXML(catFile=None,xmlFile=None)
    
    Make XML file suitable for reading into a Google map from 
    a SExtractor catalog
    """
    import threedhst
    #catFile='ib3721050_drz.cat'
    if not catFile:
        print 'makeCatXML: no `catFile` specified'
        return None
    ### Read the catalog and get id, ra, and dec columns
    cat = threedhst.sex.mySexCat(catFile)
    colID = cat.columns[cat.searchcol('NUMBER')].entry
    colRA = cat.columns[cat.searchcol('X_WORLD')].entry
    colDEC = cat.columns[cat.searchcol('Y_WORLD')].entry
    ### Make XML string
    xmlString = '<markers>'
    for i in range(cat.nrows):
        xmlString+='<marker id="%s" ra="%s" dec="%s" />' %(colID[i],
                      colRA[i],colDEC[i])
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
    zoom_levels = range(0, zoom_levels)
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
    backwards_range = range(zmin, zmax)
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
    import numpy as n
    import scipy.interpolate
    import scipy.ndimage
    
    if not a.dtype in [n.float64, n.float32]:
        a = n.cast[float](a)

    m1 = n.cast[int](minusone)
    ofs = n.cast[int](centre) * 0.5
    old = n.array( a.shape )
    ndims = len( a.shape )
    if len( newdims ) != ndims:
        print "[congrid] dimensions error. " \
              "This routine currently only support " \
              "rebinning to the same number of dimensions."
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

        trorder = [ndims - 1] + range( ndims - 1 )
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

        newcoords_dims = range(n.rank(newcoords))
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
        print "Congrid error: Unrecognized interpolation type.\n", \
              "Currently only \'neighbour\', \'nearest\',\'linear\',", \
              "and \'spline\' are supported."
        return None