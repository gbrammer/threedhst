"""
3DHST.regions

Utilities for processing DS9 region files.

"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import os
import string
import time

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import numpy as np

import threedhst

def asn_region(asn_file, path_to_flt='./'):
    """
asn_region(asn_file)
    
Create a DS9 region file for the exposures defined in an ASN file.
    
    """
    ##### Output file
    output_file = asn_file.split('.fits')[0]+'.pointing.reg'
    fp = open(output_file,'w')
    fp.write('fk5\n') ### WCS coordinates
    ##### Read ASN file
    asn = threedhst.utils.ASNFile(asn_file)
    NEXP = len(asn.exposures)
    RAcenters  = np.zeros(NEXP)
    DECcenters = np.zeros(NEXP)
    ##### Loop through exposures and get footprints
    for i, exp_root in enumerate(asn.exposures):
        flt_file =threedhst.utils.find_fits_gz(path_to_flt + '/' + exp_root.lower()+'_flt.fits', hard_break = True)
        
        #head = pyfits.getheader(exp_root.lower()+'_flt.fits')
        head = pyfits.getheader(flt_file)
        if head.get('INSTRUME') == 'ACS':
            extensions=[1,4]
        else:
            extensions=[1]
        
        for ext in extensions:
            regX, regY = wcs_polygon(flt_file,extension=ext)
            line = "polygon(%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f)"  %(regX[0],regY[0],regX[1],regY[1],regX[2],regY[2],regX[3],regY[3])

            RAcenters[i] += np.mean(regX)/len(extensions)
            DECcenters[i] += np.mean(regY)/len(extensions)
            
            fp.write(line+' # color=magenta\n')
        
    ##### Text label with ASN filename
    fp.write('# text(%10.6f,%10.6f) text={%s} font="Helvetica 14 normal" color=magenta\n' \
        %(np.mean(RAcenters),np.mean(DECcenters),
          asn_file.split('_asn.fits')[0]))
    fp.close()
    #print '3D-HST / ASN_REGION: %s\n' %(output_file)
    threedhst.showMessage('Create region file, %s.' %output_file)
    
def fits_regions(fits_list, output_file='list.reg', force_extension=None):
    ##### Output file
    fp = open(output_file,'w')
    fp.write('fk5\n') ### WCS coordinates
    
    NEXP = len(fits_list)
    
    ##### Loop through exposures and get footprints
    for i, exp_root in enumerate(fits_list):
        flt_file =threedhst.utils.find_fits_gz(exp_root, hard_break = True)
        
        #head = pyfits.getheader(exp_root.lower()+'_flt.fits')
        head = pyfits.getheader(flt_file)
        if (head.get('INSTRUME') == 'ACS') | ('UVIS' in head.get('APERTURE')):
            extensions=[1,4]
        else:
            extensions=[1]
        
        if force_extension is not None:
            extensions = force_extension
            
        RAcenters, DECcenters = 0.,0.
        
        for ext in extensions:
            regX, regY = wcs_polygon(flt_file,extension=ext)
            line = "polygon(%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f,%10.6f)"  %(regX[0],regY[0],regX[1],regY[1],regX[2],regY[2],regX[3],regY[3])
            fp.write(line+' # color=magenta\n')

            RAcenters += np.mean(regX)/len(extensions)
            DECcenters += np.mean(regY)/len(extensions)

            ##### Text label with ASN filename
            fp.write('# text(%10.6f,%10.6f) text={%s} font="Helvetica 9 normal" color=magenta\n' \
                %(RAcenters,DECcenters, exp_root.split('.fits')[0]))
            
    fp.close()
    threedhst.showMessage('Create region file, %s.' %output_file)
    
def make_zeroth(sexCat, outfile='zeroth.reg'):
    """
make_zeroth(sexCat, outfile='zeroth.reg')
    
    Make a DS9 region file with the shapes defined by the parameters in the
    SExtractor catalog but with the centers shifted to the expected position of
    the zeroth spectral order.
    """
    import pywcs
    
    x_world = sexCat.X_WORLD
    y_world = sexCat.Y_WORLD
    a_col = sexCat.A_WORLD
    b_col = sexCat.B_WORLD
    theta_col = sexCat.THETA_WORLD
    asec = 3600.
    pp = '"'
    theta_sign = -1
    
    ROOT_DIRECT = os.path.basename(sexCat.filename).split('_drz.cat')[0]
    asn_direct  = threedhst.utils.ASNFile(file=ROOT_DIRECT+'_asn.fits')
    
    flt_file = threedhst.utils.find_fits_gz(asn_direct.exposures[0]+
                                            '_flt.fits')
                                            
    flt = pyfits.open(flt_file)
    wcs = pywcs.WCS(flt[1].header)
    world_coord = []
    for i in range(len(x_world)):
        world_coord.append([np.float(x_world[i]),np.float(y_world[i])])
    xy_coord = wcs.wcs_sky2pix(world_coord,0)
    ## conf: XOFF_B -192.2400520   -0.0023144    0.0111089
    for i in range(len(x_world)):
        xy_coord[i][0] += (-192.2400520 - 0.0023144*xy_coord[i][0] +
                            0.0111089*xy_coord[i][1])
    world_coord = wcs.wcs_pix2sky(xy_coord,0)
    fp = open(outfile,'w')
    fp.write('fk5\n')
    for i in range(len(x_world)):
        line = "ellipse(%s, %s, %6.2f%s, %6.2f%s, %6.2f)\n" %(world_coord[i][0],
              world_coord[i][1], 
              float(a_col[i])*asec, pp,
              float(b_col[i])*asec, pp, float(theta_col[i])*theta_sign)
        fp.write(line)
    fp.close()
    
def trim_edge_objects(sexCat):
    """
    
    """
    import pywcs
    
    conf = threedhst.process_grism.Conf(threedhst.options['CONFIG_FILE'])
    beam_x = np.cast[float](np.array(conf.params['BEAMA'].split()))
    beam_width = 10
    
    ROOT_DIRECT = os.path.basename(sexCat.filename).split('_drz.cat')[0]
    asn_direct  = threedhst.utils.ASNFile(file=ROOT_DIRECT+'_asn.fits')
    
    flt_file = threedhst.utils.find_fits_gz(asn_direct.exposures[0]+
                                            '_flt.fits')                                        
    flt = pyfits.open(flt_file)
    wcs_flt = pywcs.WCS(flt[1].header)
    
    drz = pyfits.open(ROOT_DIRECT+'_drz.fits')
    wcs_drz = pywcs.WCS(drz[1].header)
    drz_size = drz[1].data.shape
    
    NOBJ = sexCat.nrows*1.
    
    kill = np.zeros(NOBJ)
    
    for idx in range(NOBJ):
        print(threedhst.noNewLine+'Trim edge objects: %d of %d' %(idx+1,NOBJ))
        
        rd0_drz = np.cast[float](np.array([[sexCat.X_WORLD[idx],
                                            sexCat.Y_WORLD[idx]]]))

        xy_flt = wcs_flt.wcs_sky2pix(rd0_drz,0)[0]
        xy_poly_flt = [ [xy_flt[0]+beam_x[0], xy_flt[1]-beam_width/2.], 
                        [xy_flt[0]+beam_x[1], xy_flt[1]-beam_width/2.],
                        [xy_flt[0]+beam_x[1], xy_flt[1]+beam_width/2.],
                        [xy_flt[0]+beam_x[0], xy_flt[1]+beam_width/2.] ]
        
        xy_poly_flt = []
        for dx in range(beam_x[0],beam_x[1]):
            xi = xy_flt[0]+dx
            for dy in range(-1*beam_width/2., beam_width/2):
                yi = xy_flt[1]+dy
                xy_poly_flt.append([xi,yi])
                
        rd_poly_flt = wcs_flt.wcs_pix2sky(xy_poly_flt,0)
        
        xy_poly_drz = np.round(wcs_drz.wcs_sky2pix(rd_poly_flt,0))
        # str = 'image\npolygon('
        # for i in range(4):
        #     str+='%7.2f,%7.2f,' %(xy_poly_drz[i][0],xy_poly_drz[i][1])
        # fp = open('first.reg','w')
        # fp.write(str[:-2]+')')
        # fp.close()
        
        px = np.cast[int](xy_poly_drz[:,0])
        py = np.cast[int](xy_poly_drz[:,1])
        use = np.where((px > 0) & (px < drz_size[1]) &
                       (py > 0) & (py < drz_size[0]))[0]
        if len(use) == 0:
            kill[idx] = 1
            continue
        
        px = px[use]
        py = py[use]
        
        ntot = len(xy_poly_drz)
        ### pixels that fall off the edge of the drz image
        nbad = ( len(np.where(drz[1].data[py,px] == 0)[0]) +
                 len(xy_poly_drz)-len(use) )
        #print '%d %5.2f' %(sexCat.id[idx], nbad*1./ntot)
        if nbad*1./ntot > 0.5:
            kill[idx] = 1.
    
    id = sexCat.id+0
    for i in range(NOBJ):
        if kill[i] == 1:
            sexCat.popItem(id[i])
    
def wcs_polygon(fits_file, extension=1, use_pywcs=False):
    """    
X, Y = wcs_polygon(fits_file, extension=1)
    
Calculate a DS9/region polygon from WCS header keywords.  
    
Will try to use pywcs.WCS.calcFootprint if pywcs is installed.  Otherwise
will compute from header directly.
    
    """
    ##### Open the FITS file
    hdulist = pyfits.open(fits_file) 
    ##### Get the header
    try:
        sci = hdulist[extension].header
    except IndexError:
        print('ERROR 3D-HST/wcs_polygon:\n'+\
              'Extension #%d out of range in %s' %(extension, fits_file))
        raise
    
    #### Try to use pywcs if it is installed
    pywcs_exists = True
    try:
        import pywcs
    except:
        pywcs_exists = False   
    
    if pywcs_exists & use_pywcs:
        wcs = pywcs.WCS(sci)
        footprint = wcs.calcFootprint()
        regX = footprint[:,0]    
        regY = footprint[:,1]    
        return regX, regY
    
    #### Do it by hand if no pywcs    
    NAXIS = [sci['NAXIS1'],sci['NAXIS2']]
    CRPIX = [sci['CRPIX1'],sci['CRPIX2']]
    CRVAL = [sci['CRVAL1'],sci['CRVAL2']]
    
    CD1_1 = sci['CD1_1']
    CD2_2 = sci['CD2_2']
    try:
        CD1_2 = sci['CD1_2']
        CD2_1 = sci['CD2_1']
    except:
        CD1_2 = 0.
        CD2_1 = 0.
        
    cosDec = np.cos(CRVAL[1]/180*np.pi)
    ##### Make region polygon from WCS keywords
    regX = CRVAL[0] + \
            ((np.array([0,NAXIS[0],NAXIS[0],0])-CRPIX[0])*CD1_1 +                        
             (np.array([0,0,NAXIS[1],NAXIS[1]])-CRPIX[1])*CD1_2) / cosDec
    
    regY = CRVAL[1] + \
            ((np.array([0,NAXIS[0],NAXIS[0],0])-CRPIX[0])*CD2_1 +         
             (np.array([0,0,NAXIS[1],NAXIS[1]])-CRPIX[1])*CD2_2)
             
    return regX, regY
    
def region_mask(shape,px,py):
    """
mask = region_mask(image.shape,px,py)
    
Make a mask image where pixels within the polygon defined by px_i, py_i
are set to 1.  This is the same algorithm as in :ref:`point_in_polygon`
but with array orders switched around to be much more efficient.
    
Note: something like this could be used to flag grism 0th order contaminants
    """
    NX=shape[0]
    NY=shape[1]
    y,x = np.mgrid[1:NX+1,1:NY+1]
    ##### Close polygons
    NPOLY = px.shape[0]
    tmp_px = np.append(px,px[0])
    tmp_py = np.append(py,py[0])
    theta = np.zeros((NX,NY),dtype=np.float)
    for i in np.arange(NPOLY):
        ##### Dot, cross products
        X1 = tmp_px[i] - x
        Y1 = tmp_py[i] - y 
        X2 = tmp_px[i+1] - x
        Y2 = tmp_py[i+1] - y
        dp = X1*X2 + Y1*Y2
        cp = X1*Y2 - Y1*X2
        theta += np.arctan2(cp,dp) 
    ##### Set up mask
    dq = np.zeros((NX,NY),dtype=np.int)
    flag_idx = np.where(np.abs(theta) > np.pi)
    dq[flag_idx] = 1
    return dq

class Polyreg():
    """
    Class for a polygon region to store the vertices, any text if it exists
    and a method for testing whether a specified coordinate lies within
    the polygon.
    """
    def __init__(self, line=None):
        self.px = []
        self.py = []
        self.text = ''
        if line is not None:
            self.parse(line)
            
    def parse(self, line='polygon(0,0,1,0,1,1,0,1) # text={test}'):
        import re
        spl = np.cast[float](re.split('[()]', line)[1].split(','))
        self.px = spl[0::2]
        self.py = spl[1::2]
        spl = line.split('text')
        if len(spl) > 1:
            self.text = re.split('[{}]', spl[1])[1]
            
    def contains(self, x, y):
        return point_in_polygon(x,y,self.px, self.py)
        
def parse_polygons(file='goods-s_ACSz.reg'):
    regions = []
    fp = open(file)
    lines = fp.readlines()
    fp.close()
    for line in lines:
        if line.startswith('polygon'):
            regions.append(Polyreg(line))
    
    return regions

def apply_dq_mask(flt_file, extension=3, mask_file=None, addval=4096, fk5=False, verbose=False):
    """
apply_dq_mask(flt_file, addval=2048)

    Read mask polygons from `flt_file`+'.mask.reg', if available,
    and apply to the DQ extension of `flt_file`.
    
    DQnew = DQold + `addval` within the polygon.
    
    Rewrite with Pyregion
    """
    import pyregion
    
    try:
        if mask_file is None:
            mask_file = flt_file.split('.gz')[0]+'.mask.reg'
            
        fp = open(mask_file,'r')
    except:
        return None
    #
    print('Applying mask from %s.mask_reg' %(flt_file.split('.gz')[0]))
    regions = ''.join(fp.readlines())
    fp.close()
    
    if not flt_file.endswith('.gz'):
        mode = 'update'
    else:
        mode = 'readonly'
        
    fi = pyfits.open(threedhst.utils.find_fits_gz(flt_file.split('.gz')[0]),
                     mode=mode)
                         
    ##### Set DQ bit
    r = pyregion.open(mask_file).as_imagecoord(header=fi['SCI'].header)
    mask = r.get_mask(hdu=fi['SCI'])
    fi['DQ'].data |= (addval*mask)
    
    ##### Write back to flt_file
    if not flt_file.endswith('.gz'):
        fi.flush()
    else:
        fi.writeto(flt_file.split('.gz')[0], clobber=True)
        os.remove(flt_file)
        os.system('gzip %s' %(flt_file.split('.gz')[0]))

def apply_dq_mask_old(flt_file, extension=3, mask_file=None, addval=4096, fk5=False, verbose=False):
    """
apply_dq_mask(flt_file, addval=2048)

    Read mask polygons from `flt_file`+'.mask.reg', if available,
    and apply to the DQ extension of `flt_file`.
    
    DQnew = DQold + `addval` within the polygon.
    """
    if fk5:
        import pywcs
        
    try:
        if mask_file is None:
            mask_file = flt_file.split('.gz')[0]+'.mask.reg'
            
        fp = open(mask_file,'r')
    except:
        return None
    #
    print('Applying mask from %s.mask_reg' %(flt_file.split('.gz')[0]))
    regions = ''.join(fp.readlines())
    fp.close()
    
    if not flt_file.endswith('.gz'):
        mode = 'update'
    else:
        mode = 'readonly'
        
    fi = pyfits.open(threedhst.utils.find_fits_gz(flt_file.split('.gz')[0]),
                     mode=mode)
    dqflag = np.zeros(fi[extension].data.shape,dtype=np.int)
    ##### Loop through user-defined regions
    for region in regions.split('\n'):
        if region.strip().startswith('polygon'):
            #region = 'polygon(375.05333,642.2,465.18667,642.2,751.36,
            # 709.8,393.08,326.73333,210.56,461.93333,465.18667,
            # 552.06667,375.05333,552.06667,221.82667,509.25333)'
            if verbose:
                print(region)
            
            spl = np.float_(np.array(
                     region[region.find('(')+1:region.find(')')].split(',')
                     ))
            px = spl[0::2]
            py = spl[1::2]
            if fk5:
                pra = px
                pdec = py
                wcs = pywcs.WCS(fi[extension].header)
                xy = wcs.wcs_sky2pix(pra, pdec, 1)
                px = xy[0]
                py = xy[1]
                
            dqflag += threedhst.regions.region_mask(fi[extension].data.shape,px,py)
    
    ##### Set DQ bit
    dqflag[np.where(dqflag > 0)] = addval
    fi[extension].data+=dqflag
    ##### Write back to flt_file
    
    if not flt_file.endswith('.gz'):
        fi.flush()
    else:
        fi.writeto(flt_file.split('.gz')[0], clobber=True)
        os.remove(flt_file)
        os.system('gzip %s' %(flt_file.split('.gz')[0]))
        

def point_in_polygon(x,y,px,py):
    """
test = point_in_polygon(x,y,px,py)
    
Test if coordinates (x,y) are inside polygon defined by (px, py)
    
<http://www.dfanning.com/tips/point_in_polygon.html>, translated to Python
    
    """    
    N = px.shape[0]
    ##### Close polygons
    tmp_px = np.append(px,px[0])
    tmp_py = np.append(py,py[0])
    ##### Counters
    i = np.arange(N)
    ip = np.arange(N)+1
    ##### Dot, cross products
    X1 = tmp_px[i] - x
    Y1 = tmp_py[i] - y 
    X2 = tmp_px[ip] - x
    Y2 = tmp_py[ip] - y
    dp = X1*X2 + Y1*Y2
    cp = X1*Y2 - Y1*X2
    theta = np.arctan2(cp,dp)
    if np.abs(np.sum(theta)) > np.pi:
        return True
    else:
        return False
    
    # return np.abs(np.sum(theta)) > np.pi

class PointXY():
    def __init__(self, x, y):
        self.x = x
        self.y = y

def ccw(A,B,C):
    """
ccw(A, B, C)

    Test if moving from point A to B to C moves in a counter-clockwise
    direction.
    
    Points are (P.x, P.y)
    
    http://www.bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
    """
    return (C.y-A.y)*(B.x-A.x) > (B.y-A.y)*(C.x-A.x)

def intersect(A,B,C,D):
    """
intersect(A, B, C, D)
    
    Test if line segments AB and CD intersect.
    
    Points are (P.x, P.y)
    
    http://www.bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
        """
    return ccw(A,C,D) != ccw(B,C,D) and ccw(A,B,C) != ccw(A,B,D)

def polygons_intersect(px, py, qx, qy):
    """
polygons_intersect(px, py, qx, qy)
    
    Test if two polygons intersect.  First check if line segments intersect.  If
    not, check if all points of one polygon are within another.
    """
    ### check if polygons are closed
    if (px[-1] != px[0]) | (py[-1] != py[0]):
        px = np.append(px,px[0])
        py = np.append(py,py[0])
    if (qx[-1] != qx[0]) | (qy[-1] != qy[0]):
        qx = np.append(qx,qx[0])
        qy = np.append(qy,qy[0])
    
    ### Test line segments.  Return true if any segments intersect
    NP = len(px)
    NQ = len(qx)
    for ip in range(NP-1):
        for iq in range(NQ-1):
            if intersect(PointXY(px[ip], py[ip]), PointXY(px[ip+1], py[ip+1]),
                         PointXY(qx[iq], qy[iq]), PointXY(qx[iq+1], qy[iq+1])):
                return True
    
    ### Test if first vertex of one polygon is within the other.  If it is, and 
    ### the intersection tests above were false, then all vertices have to be 
    ### within the test polygon
    if point_in_polygon(px[0],py[0],qx,qy):
        return True
    
    if point_in_polygon(qx[0],qy[0],px,py):
        return True
    
    ### Regions are (almost) identical
    if np.sum(np.abs(px-qx)+np.abs(py-qy))/np.mean(np.abs(px-np.mean(px))) < 0.1:
        return True
        
    ### All tests failed, so polygons don't intersect
    return False

def test_intersect():
    import glob
    
    ### GOODS-N
    px, py = wcs_polygon('ib3701050_drz.fits')
    zfiles = glob.glob('../ACS/h_nz_sect*drz_img.fits')
    for zfile in zfiles:
        qx, qy = wcs_polygon(zfile, extension=0)
        if threedhst.regions.polygons_intersect(px, py, qx, qy):
            print(zfile)
    