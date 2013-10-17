"""
Align and background-subtract direct and grism FLT images using Drizzlepac/Astrodrizzle rather
than Multidrizzle

Seems to be working:
====================
runAstroDrizzle    - wrapper around AstroDrizzle (irafx, 2013-09-12)
ablot_segmentation - Use drizzlepac.ablot to blot a segmentation image to an output DRZ or FLT frame 

Future ideas:
=============
Run tweakreg, astrodrizzle, tweakreg, tweakback, astrodrizzle for alignment
Make grism segmentation mask from direct image using polygons or convolution kernel
Make direct segmentation mask 
"""
import os

import pyfits
import matplotlib.pyplot as plt
import numpy as np
import scipy.ndimage as nd

import stsci.convolve

import threedhst

try:
    #### Own modified version
    import mydrizzlepac
    from mydrizzlepac import astrodrizzle, tweakreg, tweakback
except:
    import drizzlepac
    from drizzlepac import astrodrizzle, tweakreg, tweakback
    
def get_fine_alignment(asn='ibhm04030_asn.fits'):
    """
    Cross correlate shifts for fine alignment within a visit
    """
    #os.chdir('/Users/brammer/3DHST/Spectra/Work/UDS/PREP_Astrodrizzle')
    #threedhst.process_grism.fresh_flt_files(asn)
    
    #tweakreg.TweakReg(asn, refimage=None, updatewcs=True, writecat=True, clean=False, verbose=True, updatehdr=False, shiftfile=True, outshifts=asn.replace('_asn.fits', '_shifts0.txt'), conv_width=3, threshold=10, peakmax=50)
    
    #### Generate drizzled images
    threedhst.prep_flt_astrodrizzle.runAstroDrizzle(asn, clean=False, static=False, skysub=False, driz_separate=True, median=False, blot=False, driz_cr=False, driz_combine=False)
    
    #### Make SExtractor catalog and use segmentation image as cross-correlation taper mask
    a = threedhst.utils.ASNFile(asn)
    root = a.product.lower()
    
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CATALOG_NAME']    = '%s_fine.cat' %(root)
    se.options['CHECKIMAGE_NAME'] = '%s_seg.fits' %(root)
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = '%s_single_wht.fits' %(a.exposures[0])
    se.options['FILTER']    = 'Y'
    se.options['DETECT_MINAREA']    = '10'
    se.options['DETECT_THRESH']    = '7.0'
    se.options['ANALYSIS_THRESH']  = '7.0'
    se.options['MAG_ZEROPOINT'] = '25'
    se.options['DEBLEND_NTHRESH'] = '64'
    se.options['DEBLEND_MINCONT'] = '0.00002'
    se.options['PHOT_APERTURES'] = '1'
    se.options['GAIN']  = '2.4'
    
    status = se.sextractImage('%s_single_sci.fits' %(a.exposures[0]))
    threedhst.sex.sexcatRegions('%s_fine.cat' %(root), '%s_fine.reg' %(root), format=2)
    
    ### Taper mask
    seg = pyfits.open('%s_seg.fits' %(root))[0].data
    mask = nd.gaussian_filter((seg > 0)*1., 3)
    
    ref = pyfits.open(a.exposures[0]+'_single_sci.fits')[0].data*mask
    sh = ref.shape
    yi, xi = np.indices(sh)
    r = np.sqrt((yi-sh[0]/2.)**2+(xi-sh[1]/2.)**2)
    N = len(a.exposures)
    xmax, ymax = np.zeros(N), np.zeros(N)
    
    for i in range(N):
        print 'Cross corr.: %s' %(a.exposures[i])
        im = pyfits.open(a.exposures[i]+'_single_sci.fits')[0].data*mask
        # if i == 1:
        #     im = nd.shift(pyfits.open(a.exposures[i]+'_single_sci.fits')[0].data, (0.65,0.32))*mask
        # #
        cross = stsci.convolve.correlate2d(ref, im, output=None, mode='nearest', cval=0.0, fft=1)
        window = cross*(r < 4)
        xmax[i] = (xi*window).sum()/window.sum()
        ymax[i] = (yi*window).sum()/window.sum()
        #### Make diagnostic figure
        f = plt.figure(figsize=(4,4))
        f.subplots_adjust(left=0.01, bottom=0.01, right=0.99, top=0.99)
        ax = f.add_subplot(111)
        ax.imshow(window, interpolation='Nearest')
        ax.plot(xmax[i]*np.ones(2), [0,1000], color='white', alpha=0.6)
        ax.plot([0,1000], ymax[i]*np.ones(2), color='white', alpha=0.6)
        ax.set_xlim(sh[1]/2.-5, sh[1]/2.+5)
        ax.set_ylim(sh[0]/2.-5, sh[0]/2.+5)
        ax.set_xticklabels([]); ax.set_yticklabels([])
        ax.text(0.5,0.95, a.exposures[i], ha='center', va='top', transform = ax.transAxes, color='white')
        f.savefig(a.exposures[i]+'_xcor.png')
        plt.close(f)
        
    ### Write log file  
    fp = open('%s_fine.txt' %(root), 'w')
    fp.write('# image   x_shift   y_shift\n')
    for i in range(N):
        fp.write('%s  %6.3f  %6.3f\n' %(a.exposures[i], xmax[i]-xmax[0], ymax[i]-ymax[0]))
    #
    fp.close()
    
def apply_fine_alignment(asn, min_shift=0.2, max_shift=0.7, keep=True, verbose=True):
    """
    Apply fine alignment shifts to FLT headers
    """
    from threedhst import catIO
    
    a = threedhst.utils.ASNFile(asn)
    root = a.product.lower()
    
    if (not os.path.exists('%s_fine.txt' %(root))) | (not keep):
        get_fine_alignment(asn)
    
    fine = catIO.Readfile('%s_fine.txt' %(root), save_fits=False)
    if verbose:
        print 'Applying fine shifts: '
        
    for i in range(fine.N):
        flt = pyfits.open(fine.image[i]+'_flt.fits', mode='update')
        xs, ys = np.abs(fine.x_shift[i]), np.abs(fine.y_shift[i])
        vstr = '%s  %6.3f  %6.3f' %(fine.image[i], fine.x_shift[i], fine.y_shift[i])
        if (np.minimum(xs, ys) > min_shift) & (np.maximum(xs, ys) < max_shift):
            for ext in range(5):
                flt[ext+1].header['CRPIX1'] -= fine.x_shift[i]
                flt[ext+1].header['CRPIX2'] -= fine.y_shift[i]
            #
            vstr += ' *'
        #
        if verbose:
            print vstr
            
        flt.flush()
    
def align_to_reference(asn, matchImage='/Users/brammer/3DHST/Ancillary/Mosaics/uds-f160w-astrodrizzle-v4.0_drz_sci.fits'):
    """
    Align to external reference using tweakshifts
    """
    threedhst.prep_flt_astrodrizzle.runAstroDrizzle(asn, clean=False, static=True, skysub=True, driz_separate=True, median=True, blot=True, driz_cr=True, driz_combine=True, build=True, final_rot=0.)
    
    threedhst.shifts.matchImagePixels(input=[matchImage], matchImage=asn.replace('_asn', '_drz'), output=asn.replace('_asn', '_align'), input_extension=0, match_extension=1)
    
    a = threedhst.utils.ASNFile(asn)
    root = a.product.lower()
    
    se = threedhst.sex.SExtractor()
    se.aXeParams()
    se.copyConvFile()
    se.overwrite = True
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['WEIGHT_TYPE']     = 'MAP_WEIGHT'
    se.options['WEIGHT_IMAGE']    = '%s_drz.fits[1]' %(root)
    se.options['FILTER']    = 'Y'
    se.options['DETECT_MINAREA']    = '10'
    se.options['DETECT_THRESH']    = '4.0'
    se.options['ANALYSIS_THRESH']  = '4.0'
    se.options['MAG_ZEROPOINT'] = '25'
    se.options['DEBLEND_NTHRESH'] = '64'
    se.options['DEBLEND_MINCONT'] = '0.00002'
    se.options['PHOT_APERTURES'] = '1'
    se.options['GAIN']  = '2.4'
    
    #### Detect on input image
    se.options['CATALOG_NAME']    = '%s_wcs.cat' %(root)
    se.options['CHECKIMAGE_NAME'] = '%s_wcs_seg.fits' %(root)
    status = se.sextractImage('%s_drz.fits[0]' %(root))
    threedhst.sex.sexcatRegions('%s_wcs.cat' %(root), '%s_wcs.reg' %(root), format=2)
    
    #### Detect on alignment image
    se.options['CATALOG_NAME']    = '%s_align.cat' %(root)
    se.options['CHECKIMAGE_NAME'] = '%s_align_seg.fits' %(root)
    status = se.sextractImage(asn.replace('_asn', '_align'))
    threedhst.sex.sexcatRegions('%s_align.cat' %(root), '%s_align.reg' %(root), format=2)
    
    #### 'wcs' is in (x,y) pixels and 'align' is in (ra,dec)
    for ext, col in zip(['wcs', 'align'], [(1,2), (3,4)]):
        c = np.loadtxt('%s_%s.cat' %(root, ext))
        N = c.shape[0]
        fp = open('%s_%s.xy' %(root, ext), 'w')
        for i in range(N): 
            fp.write('%.5f %.5f %.3f\n' %(c[i,col[0]], c[i,col[1]], c[i,17]))
        fp.close()
            
    fp = open('%s_wcs.list' %(root), 'w')
    fp.write('%s_drz.fits %s_wcs.xy\n' %(root, root))
    fp.close()
    
    tweakreg.TweakReg(asn.replace('_asn', '_drz'), refimage=asn.replace('_asn', '_align'), updatewcs=False, writecat=False, updatehdr=False, shiftfile=True, outshifts='%s_wcs.shifts' %(root), catfile='%s_wcs.list' %(root), xcol=1, ycol=2, fluxcol=3, fluxunits='mag', xyunits='pixels', refcat='%s_align.xy' %(root), refxcol=1, refycol=2, rfluxcol=3, refxyunits='degrees', rfluxunits='mag')
    
    #### Check wcs.shifts
    tweakreg.TweakReg(asn.replace('_asn', '_drz'), refimage=asn.replace('_asn', '_align'), updatewcs=False, writecat=False, updatehdr=True, wcsname='TWEAK2', catfile='%s_wcs.list' %(root), xcol=1, ycol=2, fluxcol=3, fluxunits='mag', xyunits='pixels', refcat='%s_align.xy' %(root), refxcol=1, refycol=2, rfluxcol=3, refxyunits='degrees', rfluxunits='mag')
    
    threedhst.shifts.matchImagePixels(input=[matchImage], matchImage=asn.replace('_asn', '_drz'), output=asn.replace('_asn', '_align'), input_extension=0, match_extension=1)

    test = threedhst.shifts.align_to_reference(root.upper(), '%s_align.fits' %(root), fitgeometry="shift",
        clean=True, verbose=False, ALIGN_EXTENSION=0, toler=3, skip_swarp=True,
        align_sdss_ds9=False)    
    #
    mydrizzlepac.updatehdr.updatewcs_with_shift('%s_drz.fits' %(a.product), '%s_align.fits' %(root), xsh=test[0], ysh=test[1], rot=test[2], force=True, verbose=True)
    
    tweakback.tweakback('%s_drz.fits' %(a.product), input=None, origwcs=None, newname=None, wcsname=None, extname='SCI', force=False, verbose=False)
    threedhst.prep_flt_astrodrizzle.runAstroDrizzle(asn, clean=False, static=False, skysub=False, driz_separate=False, median=False, blot=False, driz_cr=False, driz_combine=True, build=True, final_rot=0.)
    
def xxx():
    #os.chdir('/Users/brammer/3DHST/Spectra/Work/Perlmutter/PREP_FLT')
    #asn = 'MACSJ1720+3536-F140W_asn.fits'
    #REF_IMAGE = '../CLASH/hlsp_clash_hst_wfc3ir_macs1720_f160w_v1_drz.fits'
    
    #### First run of Astrodrizzle to flag CRs and subtract simple background
    tweakreg.TweakReg(asn, refimage=None, updatewcs=True, writecat=True, updatehdr=True)
    threedhst.prep_flt_astrodrizzle.runAstroDrizzle(asn, skysub=True, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_combine=True, build=False)
    
    tweakreg.TweakReg(asn, refimage=REF_IMAGE, updatewcs=True, writecat=True)
    tweakreg.TweakReg(asn.replace('_asn','_drz_sci'), refimage=REF_IMAGE, updatewcs=True, writecat=True, clean=True, updatehdr=True, headerlet=False, shiftfile=False)
    tweakback.tweakback(asn.replace('_asn','_drz_sci'), input=None, verbose=True)
    threedhst.prep_flt_astrodrizzle.runAstroDrizzle(asn, skysub=True, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_combine=True, build=False, wcskey='TWEAK_1')
    
    ### doesn't work
    
def runAstroDrizzle(asn_file='ibhm04030_asn.fits', final_scale=0.06, pixfrac=0.8, clean=True, use_mdrztab=True, ivar_weights=True, rms_weights=False, build_drz=True, driz_cr_snr='5.0 4.0', driz_cr_scale = '1.2 0.7', **more_params):
    """
    Wrapper to run AstroDrizzle
    
    Steps: 
        skysub=True, driz_separate=True, driz_sep_wcs=True, 
        median=True, blot=True, driz_cr=True, driz_combine=True
        
    """
    import threedhst
    asn = threedhst.utils.ASNFile(file=asn_file)
    
    #### Set default parameters from pipeline mdz file
    if use_mdrztab:
        #### Read the first FLT image and read its MDRIZTAB file
        flt = pyfits.open(asn.exposures[0]+'_flt.fits')
        
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
            if len(asn.exposures) >= mdz.field('numimages')[i]:
                use = i
        
        #### Now set all of the parameters
        param_dict = {}
        for param in mdz.names:
            try:
                value = mdz.field(param)[use]
                if (not np.isfinite(value)) | (value < -1000):
                    value = ''
                #
                param_dict[param] = value
            except:
                #### some columns in the MDZ file aren't actually parameters, skip
                pass
        #
        param_dict['crbit'] = 4096
        param_dict['combine_type'] = 'minmed'
        #param_dict['combine_type'] = 'median'
        param_dict['mdriztab'] = True
        param_dict['expkeyword'] = 'EXPTIME'
             
    #
    param_dict['driz_cr_snr'] = '5.0 4.0'
    param_dict['driz_cr_scale'] = '1.2 0.7'
    
    if flt[0].header['INSTRUME'] == 'WFC3':
        #iraf.dither.multidrizzle.driz_cr_snr = '6.0 3.0'
        driz_cr_snr = '3.5 3.0'
        #iraf.dither.multidrizzle.driz_cr_scale = '1.6 0.7'
        ### More conservative to avoid rejecting central pixels of stars
        driz_cr_scale = '2.5 0.7'
        param_dict['rnkeyword'] = 'READNSEA,READNSEB,READNSEC,READNSED'
        param_dict['gnkeyword'] = 'ATODGNA,ATODGNB,ATODGNC,ATODGND'
        
    #
    if rms_weights:
        #### Generate inverse variance weight map, will need 
        #### flat + dark images in the iref or jref directories
        param_dict['final_wht_type'] = 'RMS'
    
    if ivar_weights:
        #### Generate inverse variance weight map, will need 
        #### flat + dark images in the iref or jref directories
        param_dict['final_wht_type'] = 'IVM'
    
    bad_keys = ['subsky', 'numimages', 'ra', 'dec', 'crbitval']
    for key in bad_keys:
        if key in param_dict.keys():
            status = param_dict.pop(key)
    
    # param_dict['runfile'] = asn_file.replace('_asn.fits','_asn')
    # for key in more_params.keys():
    #     param_dict[key] = more_params[key]
    print 
    
    param_dict['final_wcs'] = True
    param_dict['final_scale'] = final_scale
    param_dict['final_pixfrac'] = pixfrac
    param_dict['clean'] = clean
    
    astrodrizzle.AstroDrizzle(input=asn_file, runfile=asn_file.replace('_asn.fits','_asn'), final_wcs=True, final_scale=final_scale, final_pixfrac=pixfrac, clean=clean, final_wht_type=param_dict['final_wht_type'], driz_cr_snr=driz_cr_snr, driz_cr_scale=driz_cr_scale, **more_params)
    
def test_grism_mask():
    """
    try convolving direct image with a kernel to make the grism
    """
    
    threedhst.prep_flt_astrodrizzle.runAstroDrizzle('ibhm04030_asn.fits', rms_weights=True, final_scale=0.128)
    threedhst.prep_flt_astrodrizzle.runAstroDrizzle('ibhm04040_asn.fits', rms_weights=True, final_scale=0.128)
    
    direct = pyfits.open('IBHM04030_drz.fits') #[1].data
    grism = pyfits.open('IBHM04040_drz.fits') #[1].data    

    sigma = 1./np.sqrt(grism[2].data)
    ok = np.isfinite(sigma)
    direct[1].data[~ok] = 0
    
    xkern = np.arange(-207,207)
    kernel = np.zeros((15, xkern.shape[0]))
    kernel[:, (xkern > 25) & (xkern < 196)] = 1.
    kernel[:, (xkern > -207) & (xkern < -177)] = 1.
    
    test = nd.convolve(direct[1].data, kernel/kernel.sum()*3)
    
    gr = grism[1].data*1.
    mask = test/np.median(sigma) > 4
    
    #### Can "ablot" the mask directly back to the FLT frame
    
def test_ablot_seg():
    threedhst.prep_flt_astrodrizzle.ablot_segmentation(segimage='UDS_F160W_seg.fits', refimage='UDS-4-F140W_drz.fits', output='blot_drz_seg.fits')
    threedhst.prep_flt_astrodrizzle.ablot_segmentation(segimage='UDS_F160W_seg.fits', refimage='ibhm04a7q_flt.fits', output='blot_flt_seg.fits')
    
def ablot_segmentation(segimage='UDS_F160W_seg.fits', refimage='UDS-4-F140W_drz.fits', output='blot_seg.fits', seg_ext=0, ref_ext=1, clobber=True):
    """
    Try blotting a segmentation image
    """
    #os.chdir("/Users/brammer/3DHST/Spectra/Work/UDS/PipelineDemo")

    # segimage='UDS_F160W_seg.fits'; refimage='UDS-4-F140W_drz.fits'; output='blot_seg.fits'; seg_ext=0; ref_ext=1; clobber=True

    import pyfits
    #import drizzlepac
    #from drizzlepac import astrodrizzle
    import stwcs
    
    #segimage = 'UDS_F160W_seg.fits'
    #segimage = 'blot_seg_full.fits'
    seg = pyfits.open(segimage)
    ref = pyfits.open(refimage)
    
    seg_wcs = stwcs.wcsutil.HSTWCS(seg, ext=seg_ext)
    ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=ref_ext)
    
    use_coeffs = refimage.split('.gz')[0].endswith('_flt.fits')
    ### Blot segmentation image
    blotted = astrodrizzle.ablot.do_blot(seg[seg_ext].data, seg_wcs, ref_wcs, 1, coeffs=use_coeffs, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)
    
    ### Blot "ones" image to get relative pixel weights
    ones = np.cast[np.int8](seg[seg_ext].data > 0)
    blotted_ones = astrodrizzle.ablot.do_blot(ones, seg_wcs, ref_wcs, 1, coeffs=use_coeffs, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)
    blotted_ones[blotted_ones == 0] = 1
    
    ### Scale by relative pixel areas to get input segmentatino values
    ratio = blotted / blotted_ones
    
    ### Force integers
    segm = np.cast[np.int32](np.round(ratio))
    
    ### Write output
    pyfits.writeto(output, data=segm, header=ref[ref_ext].header, clobber=clobber)
    
def test_blot_clean():
    """
    See if ablot can work like SWarp to blot to a rotated output frame
    """
    
    os.chdir('/Users/brammer/3DHST/Spectra/Work/AEGIS/Test_Astrodrizzle')
    asn = 'ibhj56030_asn.fits'
    threedhst.prep_flt_astrodrizzle.runAstroDrizzle(asn, rms_weights=True, final_scale=0.06)
    
    matchImage = '/Users/brammer/3DHST/Ancillary/Mosaics/aegis-f160w-astrodrizzle-v4.0_drz_sci.fits'
    matchImage = '/Users/brammer/3DHST/Ancillary/Mosaics/uds-f160w-astrodrizzle-v4.0_drz_sci.fits'
    
    #threedhst.shifts.matchImagePixels(input=[matchImage], matchImage=asn.replace('_asn', '_drz'), output=asn.replace('_asn', '_align'), input_extension=0, match_extension=1)
    
    import drizzlepac
    from drizzlepac import astrodrizzle
    import stwcs

    ref = pyfits.open(matchImage)
    next = pyfits.open(asn.replace('_asn', '_drz'))
    
    ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=0)
    flt_wcs = stwcs.wcsutil.HSTWCS(next, ext=1)
    
    root = asn.split('_asn')[0]
    a = threedhst.utils.ASNFile(asn)
    
    blotted = astrodrizzle.ablot.do_blot(ref[0].data, ref_wcs, flt_wcs, 1, coeffs=False, interp='poly5', sinscl=1.0, stepsize=10, wcsmap=None)
    pyfits.writeto('%s_align.fits' %(root), data=blotted, header=next[1].header, clobber=True)
    
    #### Alignment
    # test = threedhst.shifts.align_to_reference(root, 'align_blot.fits', fitgeometry="rotate",
    #     clean=True, verbose=False, ALIGN_EXTENSION=0, toler=3, skip_swarp=True,
    #     align_sdss_ds9=False)    
    #
    test = threedhst.shifts.align_to_reference(a.product, '%s_align.fits' %(root), fitgeometry="rxyscale",
        clean=True, verbose=False, ALIGN_EXTENSION=0, toler=3, skip_swarp=True, align_sdss_ds9=False)    
    #
    mydrizzlepac.updatehdr.updatewcs_with_shift('%s_drz.fits' %(a.product), '%s_align.fits' %(root), xsh=test[0], ysh=test[1], scale=test[3], rot=test[2], force=True, verbose=True)
    mydrizzlepac.tweakback.tweakback('%s_drz.fits' %(a.product), input=None, origwcs=None, newname=None, wcsname=None, extname='SCI', force=False, verbose=False)
    