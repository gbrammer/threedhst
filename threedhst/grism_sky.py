"""

New correction for the grism sky:
    1) Divide by an imaging flat-field (F140W)
    2) Divide by a master sky image
    3) Subtract the residual flux parallel to the image columns
    4) Subtract the overall biweight mean
    
"""

__version__ = "$Rev: 199 $"
# $URL: https://threedhst.googlecode.com/svn/threedhst/gmap.py $
# $Author: gbrammer $
# $Date: 2011-05-22 01:59:38 -0400 (Sun, 22 May 2011) $

import matplotlib.pyplot as plt

import numpy as np
import threedhst
import threedhst.prep_flt_files
import pyfits
import glob
import os

IREF = os.getenv('iref')

# try:
#     flat_f140 = pyfits.open(IREF+'/uc721143i_pfl.fits')
#     #print 'Make grism_sky_flat'
#     #flat_f140 = pyfits.open(IREF+'/flat.IR_avg.fits')
#     flat_f140 = pyfits.open(IREF+'cosmos_f140w_flat.fits')
#     flat_g141 = pyfits.open(IREF+'/u4m1335mi_pfl.fits')
#     flat_master_g141 = flat_g141[1].data[5:1019,5:1019] / flat_f140[1].data[5:1019, 5:1019]
#     flat_master_g141[flat_master_g141 <= 0] = 5
#     flat_master_g141[flat_master_g141 > 5] = 5
#     #
#     flat_f105 = pyfits.open(IREF+'/uc72113oi_pfl.fits')
#     #print 'Make grism_sky_flat'
#     flat_g102 = pyfits.open(IREF+'/u4m1335li_pfl.fits')
#     flat_master_g102 = flat_g102[1].data[5:1019,5:1019] / flat_f105[1].data[5:1019, 5:1019]
#     flat_master_g102[flat_master_g102 <= 0] = 5
#     flat_master_g102[flat_master_g102 > 5] = 5    
# except:
#     print '\nthreedhst.grism_sky: Flat-field files (uc721143i_pfl.fits) not found in IREF: %s\n' %(IREF)
#     flat_master_g141 = np.ones((1014,1014))
#     flat_master_g102 = np.ones((1014,1014))

xprofile = None
yprofile = None

flat_grism = None
flat_direct = IREF

def set_grism_flat(grism='G141', verbose=True):
    import threedhst.grism_sky as bg
    
    if bg.flat_grism == grism:
        return True
    
    if verbose:
        print 'Set flat for grism: %s' %(grism)
    
    if grism == 'G141':
        #flat_f140 = pyfits.open(IREF+'/uc721143i_pfl.fits')
        #flat_f140 = pyfits.open(IREF+'cosmos_f140w_flat.fits')
        flat_f140 = pyfits.open(IREF+'/flat_3DHST_F140W_t1_v0.1.fits')
        flat_g141 = pyfits.open(IREF+'/u4m1335mi_pfl.fits')
        flat = flat_g141[1].data[5:1019,5:1019] / flat_f140[1].data[5:1019, 5:1019]
        flat[flat <= 0] = 5
        flat[flat > 5] = 5
        bg.flat = flat
        bg.flat_grism = 'G141'
        bg.flat_direct = flat_f140.filename()
        
    else:
        flat_f105 = pyfits.open(IREF+'/uc72113oi_pfl.fits')
        flat_g102 = pyfits.open(IREF+'/u4m1335li_pfl.fits')
        flat = flat_g102[1].data[5:1019,5:1019] / flat_f105[1].data[5:1019, 5:1019]
        flat[flat <= 0] = 5
        flat[flat > 5] = 5
        bg.flat = flat
        bg.flat_grism = 'G102'
        bg.flat_direct = flat_f105.filename()
        
    return True
        
def remove_grism_sky(flt='ibhm46ioq_flt.fits', list=['sky_cosmos.fits', 'sky_goodsn_lo.fits', 'sky_goodsn_hi.fits', 'sky_goodsn_vhi.fits'],  path_to_sky = '../CONF/', out_path='./', verbose=False, plot=False, flat_correct=True, sky_subtract=True, second_pass=True, overall=True):
    """ 
    Process a (G141) grism exposure by dividing by the F140W imaging flat-field
    and then subtracting by a master sky image.  
    
    v1.6: list=['sky_cosmos.fits', 'sky_goodsn_lo.fits', 'sky_goodsn_hi.fits', 'sky_goodsn_vhi.fits']
    
    testing: list=['sky.G141.set001.fits','sky.G141.set002.fits','sky.G141.set003.fits','sky.G141.set004.fits','sky.G141.set005.fits','sky.G141.set025.fits','sky.G141.set120.fits']
    """
    
    import threedhst.grism_sky as bg
    #import scipy.signal as sign
    
    # flt = '../../GOODS-N/RAW/ib3708ilq_flt.fits.gz'
    im = pyfits.open(flt)
    bg.set_grism_flat(grism=im[0].header['FILTER'])
    
    segfile = os.path.basename(flt.replace('.fits','.seg.fits')).replace('.gz','')
    if os.path.exists(segfile):
        seg = pyfits.open(segfile)[0].data
        use_biweight=False
    else:
        seg = np.zeros(im[1].data.shape)
        use_biweight=True
    
    xin, yin = bg.profile(flt, extension=1, flatcorr=True, biweight=use_biweight)
    #yin /= threedhst.utils.biweight(yin[(np.abs(xin-507) < 50) & np.isfinite(yin)])
    
    if plot:
        plt.plot(xin, yin, color='black', linewidth=2)
    
    #### Loop through sky images and find the one whose column profile most
    #### closely matches the input image
    chi2 = 1.e10
    keep = None
    for sky in list:
        xsky, ysky = bg.profile(flt=path_to_sky+sky, extension=0, flatcorr=False, biweight=True)
        ysky /= np.mean(ysky[np.abs(xsky-507) < 50])
        #
        ok = np.isfinite(ysky) & np.isfinite(yin) & (yin*ysky != 0)
        a = np.sum((ysky*yin)[ok])/np.sum((ysky*ysky)[ok])
        if plot:
            plt.plot(xsky, ysky*a)
        #
        chi2_i = np.sum((ysky[ok]*a-yin[ok])**2)
        if verbose:
            print sky, chi2_i
        #
        if chi2_i < chi2:
            chi2 = chi2_i*1
            keep = sky
     
    if keep is None:
        keep = 'sky_goodsn_vhi.fits'
            
    #### The best sky image
    sk = pyfits.open(path_to_sky+keep)
    sk[0].data[sk[0].data == 0] = 1.
    sk[0].data[~np.isfinite(sk[0].data)] = 1.
    
    flat = bg.flat*1.
    
    #### Only flat correction
    dq_ok = (im[3].data & (4+32+16+512+2048+4096)) == 0
    mask = (seg == 0) & dq_ok
    if plot:
        corr = im[1].data*flat#/sk[0].data
        corr -= threedhst.utils.biweight(corr[mask], mean=True)
        ds9.frame(1)
        ds9.v(corr, vmin=-0.5,vmax=0.5)
    
    if flat_correct is False:
        flat = flat*0+1
    
    if sky_subtract is False:
        sk[0].data = sk[0].data*0+1
        
    #### Divide by the sky flat
    #corr = im[1].data*flat/sk[0].data
    # #### Show the result
    # if plot:
    #     ds9.frame(2)
    #     ds9.v(corr-threedhst.utils.biweight(corr[mask], mean=True), vmin=-0.5,vmax=0.5)

    ### Instead, subtract the sky flat
    sky_stats = threedhst.utils.biweight((im[1].data*flat/sk[0].data)[mask], both=True)
    corr = im[1].data*flat-sky_stats[0]*sk[0].data
    # #### Show the result
    # if plot:
    #     ds9.frame(3)
    #     ds9.v(corr, vmin=-0.5,vmax=0.5)
    
    #### Put the result in the FLT data extension
    im[1].data = corr*1.
     
    #### Need to write an output file to use `profile`
    im.writeto(out_path+os.path.basename(flt).replace('.gz',''), clobber=True)
    xin, yin = bg.profile(out_path+os.path.basename(flt).replace('.gz',''), extension=1, flatcorr=False, biweight=True)
        
    im = pyfits.open(out_path+os.path.basename(flt).replace('.gz',''), mode='update')

    #### Subtract the residual difference between the observed and master sky
    resid = np.dot(np.ones((1014,1)), threedhst.utils.medfilt(yin, 41).reshape(1,1014))
    if second_pass:
        im[1].data -= resid
    
    #### Subtract the overall biweight mean
    if overall:
        im[1].data -= threedhst.utils.biweight(im[1].data[mask], mean=True)
    
    #### Add a header keyword and write to the output image
    im[0].header.update('GRISMSKY',keep,comment='Image used for sky subtraction')
    im[0].header.update('SKYSCALE',sky_stats[0],comment='Scale factor of sky')
    bad = ~np.isfinite(im[1].data)
    im[1].data[bad] = 1
    im[3].data[bad] = im[3].data[bad] | 32
    im.flush()
    
    #### Show the final result, compare to the earlier version in PREP_FLT
    if plot:
        ds9.frame(3)
        ds9.v(im[1].data, vmin=-0.5,vmax=0.5)
    
        chk = pyfits.open(threedhst.utils.find_fits_gz(flt.replace('RAW','PREP_FLT').replace('.gz','')))
        ds9.frame(4)
        ds9.v(chk[1].data, vmin=-0.5,vmax=0.5)
    
def profile(flt='ibhm46ioq_flt.fits', extension=1, flatcorr=True, biweight=False):
    """
    Get a cut across the columns of a FLT image, optionally masking objects and DQ 
    pixels.  
    
    If `flatcorr` is True, then divide by the F140W flat.  
    
    If `biweight`, then the output is the biweight mean of each column.  
    Otherwise, it's just the mean.
    """
    import threedhst.grism_sky as bg
    
    im = pyfits.open(flt)
    
    if flatcorr:
        im[extension].data *= flat
    
    segfile = os.path.basename(flt.replace('.fits','.seg.fits')) #.replace('.gz','')
    if os.path.exists(segfile):
        seg = pyfits.open(segfile)[0].data
    else:
        seg = im[extension].data*0
    
    shp = im[extension].data.shape
    xpix = np.arange(shp[0])
    
    if '_flt' in flt:
        dq_ok = (im[3].data & (4+32+16+512+2048+4096)) == 0
    else:
        dq_ok = np.isfinite(im[extension].data)
        
    mask = (~dq_ok) | (seg >= 1)
    N = np.ones(shp)
    
    im[extension].data[mask] = 0
    N[mask] = 0
    
    ypix = np.sum(im[extension].data, axis=0) / np.sum(N, axis=0)
    if biweight:
        for i in range(shp[0]):
            column = im[extension].data[:,i]
            ypix[i] = threedhst.utils.biweight(column[column != 0], mean=True)
    #
    bg.xprofile, bg.yprofile = xpix, ypix
    
    return xpix, ypix #, ylo, yhi

def show_profile():
    import threedhst.grism_sky as bg
    """
    Look at images collapsed along columns to separate into groups with
    similar patterns.
    """
    
    #### COSMOS
    flt_files = glob.glob('ibhm*flt.seg.fits')
    PATH = '/3DHST/Spectra/Work/COSMOS/RAW/'
    GZ = '.gz'
    
    fp = open('COSMOS.g141.list')
    flt_files = fp.readlines()
    fp.close()
    for i in range(len(flt_files)):
        flt_files[i] = flt_files[i][:-1].replace('msk','flt')
    
    
    N = len(flt_files)
    profiles = np.zeros((N, 1014))
    for i,flt in enumerate(flt_files):
        flt = flt.replace('.seg','')
        if os.path.exists(flt+'.mask.reg'):
            continue
        #
        print flt
        xi, yi = bg.profile(flt=PATH+flt+GZ)
        profiles[i,:] = yi
    
    norm = np.zeros(N)
    test = norm > 0
    for i in range(N):
        yi = profiles[i,:]
        norm[i] = np.mean(yi[np.abs(xi-507) < 50])
        test[i] = np.median(yi[np.abs(xi-40) < 10]/norm[i]) < 1.95
        if test[i]:
            p = plt.plot(xi,yi/norm[i],color=(norm[i]/3.3,0,0), alpha=0.1)
        else:
            norm[i] = 0
            
    profiles_norm = profiles / np.dot(norm.reshape(N,1), np.ones((1,1014)))
    avg = np.mean(profiles_norm[norm != 0, :], axis=0)
    plt.plot(xi, avg, color='blue', alpha=0.5)
    
    # for i in range(N):
    #     yi = profiles[i,:]*1.
    #     if yi.sum() == 0:
    #         continue
    #     #
    #     yi-=0.
    #     nor = np.mean(yi[np.abs(xi-307) < 50])
    #     p = plt.plot(xi,yi/nor,color=(norm[i]/6,0,0), alpha=0.1)
    
    plt.ylim(0.92,1.04)
    plt.xlim(-10,1024)
    plt.savefig('COSMOS.png')
        
    #### GOODS-N
    flt_files = glob.glob('ib37*seg.fits')
    PATH = '/research/HST/GRISM/3DHST/GOODS-N/RAW/'
    Ng = len(flt_files)
    
    profiles_g = np.zeros((Ng, 1014))
    for i,flt in enumerate(flt_files):
        flt = flt.replace('.seg','')
        if os.path.exists(flt+'.mask.reg'):
            continue
        #
        print flt
        xi, yi = bg.profile(flt=PATH+flt+GZ)
        profiles_g[i,:] = yi
    
    xi = np.arange(1014)
    norm_g = np.zeros(Ng)
    test = norm_g > 0
    for i in range(Ng):
        yi = profiles_g[i,:]
        norm_g[i] = np.mean(yi[np.abs(xi-507) < 50])
        # very hi
        test[i] =  np.median(yi[np.abs(xi-200) < 20]/norm_g[i]) > 1.02
        # lo
        #test[i] =  np.median(yi[np.abs(xi-200) < 20]/norm_g[i]) < 1.01
        #test[i] =  test[i] & (np.median(yi[np.abs(xi-40) < 10]/norm_g[i]) > 0.96)
        # hi
        #test[i] =  test[i] & (np.median(yi[np.abs(xi-40) < 10]/norm_g[i]) < 0.96)
        #
        if test[i]:
            p = plt.plot(xi,yi/norm_g[i],color=(0,0,norm_g[i]/1.8), alpha=0.1)
        else:
            norm_g[i]*=0
    #
    plt.ylim(0.92,1.04)
    plt.xlim(-10,1024)
    
    profiles_norm_g = profiles_g / np.dot(norm_g.reshape(Ng,1), np.ones((1,1014)))
    avg_g = np.mean(profiles_norm_g[norm_g != 0, :], axis=0)
    plt.plot(xi, avg_g, color='green', alpha=0.5)

    plt.ylim(0.92,1.04)
    plt.xlim(-10,1024)
    
    #### AEGIS
    flt_files = glob.glob('ibhj[4]*seg.fits')
    PATH = '/research/HST/GRISM/3DHST/AEGIS/RAW/'
    Na = len(flt_files)
    
    profiles_a = np.zeros((Na, 1014))
    for i,flt in enumerate(flt_files):
        flt = flt.replace('.seg','')
        if os.path.exists(flt+'.mask.reg'):
            continue
        #
        print flt
        xi, yi = bg.profile(flt=PATH+flt+GZ)
        profiles_a[i,:] = yi
    
    xi = np.arange(1014)
    norm_a = np.zeros(Na)
    test = norm_a > 0
    for i in range(Na):
        yi = profiles_a[i,:]
        norm_a[i] = np.mean(yi[np.abs(xi-507) < 50])
        # very hi
        test[i] =  np.median(yi[np.abs(xi-200) < 20]/norm_a[i]) < 1.52
        # lo
        #test[i] =  test[i] & (np.median(yi[np.abs(xi-40) < 10]/norm_a[i]) > 0.96)
        # hi
        #test[i] =  test[i] & (np.median(yi[np.abs(xi-40) < 10]/norm_a[i]) < 0.96)
        #
        if test[i]:
            p = plt.plot(xi,yi/norm_a[i],color=(0,0,norm_a[i]/1.8), alpha=0.1)
        else:
            norm_a[i]*=0
    #
    plt.ylim(0.92,1.04)
    plt.xlim(-10,1024)
    
def make_bg(GZ='.gz'):
    """
    Make average background images with object masks
    """
    files = glob.glob('ibhm*flt.seg.fits')
    PATH = '/research/HST/GRISM/3DHST/COSMOS/RAW/'
    PATH = '/3DHST/Spectra/Work/COSMOS/RAW/'
    
    fp = open('COSMOS.g141.list')
    files = fp.readlines()
    fp.close()
    for i in range(len(flt_files)):
        files[i] = files[i][:-1].replace('msk','flt')
    
    files = glob.glob('ib37*flt.seg.fits')
    PATH = '/research/HST/GRISM/3DHST/GOODS-N/RAW/'
    
    #### Direct flat-field
    flat = pyfits.open(IREF+'/uc721143i_pfl.fits')[1].data[5:-5,5:-5]
    flat[flat <= 0] = 5
    flat[flat > 5] = 5
    
    #### Candels
    os.chdir('/Users/gbrammer/CANDELS/Flats/')
    files = np.array(glob.glob('ib*flt.seg.fits'))
    PATH = '/Users/gbrammer/CANDELS/UDS/RAW/'
    
    info = catIO.Readfile(PATH+'../PREP_FLT/files.info')
    
    files = files[info.filter == 'F125W']
    flat = pyfits.open(IREF+'/uc72113qi_pfl.fits')[1].data[5:-5,5:-5]
    
    NF = len(files)
    idx = np.arange(NF)
    X = np.zeros((NF,1014.**2))
        
    ## Otherwise get it from "show_profile" above
    test = idx > -10
    
    for j,i in enumerate(idx):
        if ~test[i]:
            continue
        #
        fi = files[i].replace('.seg','')
        if not os.path.exists(fi.replace('flt','flt.seg')+GZ):
            continue
        #    
        if os.path.exists(fi+'.mask.reg'):
            continue
        #
        print '%d %s' %(i, files[i])
        flt = pyfits.open(PATH+fi+'.gz')
        flt[1].data *= flat
        ### Segmentation mask
        masked = pyfits.open(fi.replace('flt','flt.seg')+GZ)[0].data == 0
        ### DQ mask, hot pixels and the "death star"
        dq_ok = (flt[3].data & (4+32+16)) == 0
        #
        ok = masked & np.isfinite(flt[1].data) & (dq_ok)
        flt[1].data /= np.median(flt[1].data[ok])
        flt[1].data[(ok == False)] = 0
        X[j,:] = flt[1].data.flatten()

    #### Average
    nsum = np.sum(X != 0, axis=0).reshape(1014,1014)
    avg = np.sum(X, axis=0).reshape(1014,1014)/nsum
     
    ### Fill empty pixels with no input images
    sky = avg
    x,y = np.where((np.isfinite(sky) == False) | (sky/flat > 1.15))
    NX = len(x)
    pad = 1
    for i in range(NX):
        xi = x[i]
        yi = y[i]
        sub = sky[xi-pad:xi+pad+2,yi-pad:yi+pad+2]
        if (np.sum(sub) != 0.0):
            sky[xi,yi] = np.median(sub[np.isfinite(sub)])
    
    still_bad = (np.isfinite(sky) == False) | (sky <= 0.01)
    sky[still_bad] = flat[still_bad]
    
    # bad_flat = (flat < 0.5)
    # sky[bad_flat] = flat[bad_flat]
        
    im_sky = pyfits.PrimaryHDU(data=sky)
    im_n = pyfits.ImageHDU(data=nsum)
    im = pyfits.HDUList([im_sky, im_n])
    im.writeto('sky.fits', clobber=True)
    
    #### for DIRECT flat
    flatim = pyfits.open(IREF+'/uc721143i_pfl.fits')
    flatim[1].data[5:-5,5:-5] = sky
    flatim[3].data[5:-5,5:-5] = nsum
    #flatim.writeto('/research/HST/GRISM/IREF/cosmos_f140w_flat.fits', clobber=True)
    
def regenerate_segmaps():
    import os
    import glob
    import threedhst
    
    os.chdir('/research/HST/GRISM/3DHST/COSMOS/PREP_FLT')
    os.chdir('/research/HST/GRISM/3DHST/GOODS-N/DATA')
    os.chdir('/research/HST/GRISM/3DHST/AEGIS/DATA')
    os.chdir('/research/HST/GRISM/3DHST/GOODS-S/DATA')
    
    IS_GRISM=True
    
    #### New F140W flat
    IS_GRISM=False
    
    drz = glob.glob('*G141.run')
    for d in drz:
        run = threedhst.prep_flt_files.MultidrizzleRun(d.replace('.run',''))
        for i in range(4):
            run.blot_back(ii=i, copy_new=(i is 0), WHT=False)
            threedhst.prep_flt_files.make_segmap(run.flt[i], IS_GRISM=IS_GRISM)
        files=glob.glob('*BLOT*')
        for file in files:
            os.remove(file)
    
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS/BACKGROUND_PCA')  
                     
def get_segmap():
    
    fp = open('cosmos_g141.list')
    exp_list = fp.readlines()
    fp.close()
    
    asn = threedhst.utils.ASNFile('/research/HST/GRISM/3DHST/COSMOS/RAW/ibhm31040_asn.fits')
    
    flat_f140 = pyfits.open('/research/HST/GRISM/IREF//uc721143i_pfl.fits')
    flat_g141 = pyfits.open('/research/HST/GRISM/IREF//u4m1335mi_pfl.fits')
    
    for i in range(0,len(exp_list),4):
        asn.product='junk'
        asn.exposures = []
        for j in range(4):
            asn.exposures.append(exp_list[i+j].split('_flt')[0])
        
        asn.write('junk_asn.fits')
        threedhst.shifts.make_blank_shiftfile('junk_asn.fits')
        
        threedhst.process_grism.fresh_flt_files('junk_asn.fits', from_path='/research/HST/GRISM/3DHST/COSMOS/RAW/')
        
        for exp in asn.exposures:
            im = pyfits.open(exp+'_flt.fits','update')
            im[1].data *= flat_g141[1].data[5:1019,5:1019] / flat_f140[1].data[5:1019, 5:1019]
            im.flush()
            
        threedhst.prep_flt_files.startMultidrizzle('junk_asn.fits', use_shiftfile=True, 
            skysub=False, final_scale=0.128254, pixfrac=1.0,
            driz_cr=False, median=False, updatewcs=False)
        
        run = threedhst.prep_flt_files.MultidrizzleRun('junk')
        for i,exp in enumerate(asn.exposures):
            run.blot_back(ii=i, copy_new=(i is 0))
            threedhst.prep_flt_files.make_segmap(run.flt[i])
        
        cleanup()
        
def cleanup():
    import glob
    import os
    files=glob.glob('*coeffs?.dat')
    files.extend(glob.glob('*BLOT*'))
    files.extend(glob.glob('drz*'))
    files.extend(glob.glob('bg.fits'))
    for file in files:
        os.remove(file)
        