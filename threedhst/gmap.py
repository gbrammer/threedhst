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
MAP_SIZE = [256,256]    
# This is the Maps API key for running on localhost:8080
MAP_KEY = 'ABQIAAAA1XbMiDxx_BTCY2_FkPh06RR20YmIEbERyaW5EQEiVNF0mpNGfBSRb' \
    '_rzgcy5bqzSaTV8cyi2Bgsx3g'

def makeGMapTiles(fitsimage):
    import pyfits
    import pywcs
    import fitsimage
    
    fitsimage = 'ib3721050_drz.fits'
    
    ### Read the FITS file
    fi = pyfits.open(fitsimage)
    head = fi[1].header
    data = fi[1].data
    
    ### Image corners in Lat/Lon
    wcs = pywcs.WCS(head)
    llSW = radec2latlon(wcs.wcs_pix2sky([[1,1]],1)[0])
    llNW = radec2latlon(wcs.wcs_pix2sky([[wcs.naxis1,1]],1)[0])
    llSE = radec2latlon(wcs.wcs_pix2sky([[1,wcs.naxis2]],1)[0])
    llNE = radec2latlon(wcs.wcs_pix2sky([[wcs.naxis1,wcs.naxis2]],1)[0])
    
    ### Get Google Map pixel/tile coordinates
    m = gmap.MercatorProjection()
    bounds = [llSW,llNE]
    view_size = [wcs.naxis1,wcs.naxis2]
    zoomLevel = m.CalculateBoundsZoomLevel(bounds, view_size)
    pixCenter = m.FromLatLngToPixel((llSW+llNE)/2.,zoomLevel)
    tilex = pixCenter.x/256
    tiley = pixCenter.y/256
    
    pixSW = m.FromLatLngToPixel(llSW,zoomLevel)
    pixSE = m.FromLatLngToPixel(llSE,zoomLevel)
    pixNW = m.FromLatLngToPixel(llNW,zoomLevel)
    pixNE = m.FromLatLngToPixel(llNE,zoomLevel)
    
    image.save('/tmp/junk.png')
    
    return None
    
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
    latlon = np.zeros(2.)
    latlon[0] = radec[1]
    latlon[1] = 360-radec[0]
    return latlon
    
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
        self.tilex = x/256.
        self.tiley = y/256.
        
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
    self.pixels = 256
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
