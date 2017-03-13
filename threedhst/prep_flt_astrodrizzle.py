"""
Align and background-subtract direct and grism FLT images using 
Drizzlepac/Astrodrizzle rather than Multidrizzle

xxx still need full "pair" wrapper to put it all together, including adding direct shifts to grism exposures

"""
import os

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import numpy as np

import matplotlib.pyplot as plt
import scipy.ndimage as nd

import stsci.convolve

import threedhst
from threedhst import catIO

def clean_wcsname(flt='ibhj15wyq_flt.fits', wcsname='TWEAK', ACS=False, WFPC2=False):
    """
    Workaround for annoying TweakReg feature of not overwriting WCS solns
    """
    im = pyfits.open(flt, mode='update')
    if ACS:
        exts = [1,4]
    elif WFPC2:
        exts = [1,2,3,4]
    else:
        exts = [1]
    
    for ext in exts:
        header = im[ext].header
        for key in header:
            if key.startswith('WCSNAME'):
                if header[key] == wcsname:
                    wcs_ext = key[-1]
                    if key == 'WCSNAME':
                        header[key] = 'X' + wcsname+'X'
                        #im.flush()
        #
        for key in ['WCSNAME', 'WCSAXES', 'CRPIX1', 'CRPIX2', 'CDELT1', 'CDELT2', 'CUNIT1', 'CUNIT2', 'CTYPE1', 'CTYPE2', 'CRVAL1', 'CRVAL2', 'LONPOLE', 'LATPOLE', 'CRDER1', 'CRDER2', 'CD1_1', 'CD1_2', 'CD2_1', 'CD2_2', 'FITNAME', 'NMATCH', 'RMS_RA', 'RMS_DEC']:
            try:
                header.remove(key+wcs_ext)
            except:
                #print key
                pass
    
    im.flush()
    
def drzTweakReg(sci='goodss-34-F140W_drz_sci.fits', master_catalog='goodss_radec.dat', threshold=20, apply=True):
    import drizzlepac
    from drizzlepac import tweakback
    from stwcs import updatewcs
    
    se = threedhst.sex.SExtractor()
    se.options['WEIGHT_IMAGE'] = sci.replace('sci','wht')
    se.options['WEIGHT_TYPE'] = 'MAP_WEIGHT'
    #
    se.params['X_IMAGE'] = True; se.params['Y_IMAGE'] = True
    se.params['X_WORLD'] = True; se.params['Y_WORLD'] = True
    se.params['MAG_AUTO'] = True
    #
    se.options['CATALOG_NAME'] = sci+'.align.cat'
    se.options['FILTER'] = 'N'
    se.options['DETECT_THRESH'] = '%f' %(threshold)
    se.options['ANALYSIS_THRESH'] = '%f' %(threshold)
    #
    se.sextractImage(sci)
    threedhst.sex.sexcatRegions(sci+'.align.cat', sci+'.align.reg', format=2)
    
    c = catIO.Table(sci+'.align.cat', format='ascii.sextractor')
    #c.ra = c['X_WORLD']
    #c.dec = c['Y_WORLD']
    m = catIO.CoordinateMatcher(c, ra_column='X_WORLD', dec_column='Y_WORLD')
    r0, d0 = np.loadtxt(master_catalog, unpack=True)
    
    ### clip list to nearby objects
    rmed, dmed = np.median(c['X_WORLD']), np.median(c['Y_WORLD'])
    delta = np.sqrt((r0-rmed)**2/np.cos(dmed/180*np.pi)**2+(d0-dmed)**2)*60.
    nearby = delta < 8 # arcmin
    r0, d0 = r0[nearby], d0[nearby]
    
    dr, idx = m.match_list(r0, d0)
    
    dx = (c['X_WORLD'][idx]-r0)*np.cos(d0/180*np.pi)*3600
    dy = (c['Y_WORLD'][idx]-d0)*3600

         
    x0 = (c['X_WORLD'][idx]-np.median(c['X_WORLD']))*np.cos(d0/180*np.pi)*3600
    y0 = (c['Y_WORLD'][idx]-np.median(c['Y_WORLD']))*3600
    
    ok = dr < 1.5 
    if ok.sum() == 0:
        threedhst.showMessage('No matches found within 1.5".')
        return False
    
    # plt.scatter(x0[ok], y0[ok], color='black')
    # for i in np.arange(len(ok))[ok]:
    #     plt.plot(x0[i]+np.array([0, dx[i]*20]), y0[i]+np.array([0, dy[i]*20]), color='black')
     
    dra = (c['X_WORLD'][idx]-r0)
    dde = (c['Y_WORLD'][idx]-d0) 
    rshift, dshift = np.median(dra[ok]), np.median(dde[ok])
    
    fp = open(sci.split('.fits')[0]+'.align.dat','w')
    lines = ['# dx dy xrms yrms N\n# %s %s\n' %(sci, master_catalog), '%f %f %f %f %d\n' %(np.median(dx[ok]), np.median(dy[ok]), np.std(dx[ok]), np.std(dy[ok]), ok.sum())]
    fp.writelines(lines)
    fp.close()
    threedhst.showMessage(''.join(lines))
    
    if not apply:
        print('Not applying shift.  Re-run with apply=true to apply them.')
        return rshift, dshift
        
    for fits in [sci.replace('sci','wht'), sci]:
        print('Update WCS: %s' %(fits))
        im = pyfits.open(fits, mode='update')
        im[0].header['CRVAL1'] -= rshift
        im[0].header['CRVAL2'] -= dshift
        im.flush()
    
    im = pyfits.open(sci)
    for i in range(im[0].header['NDRIZIM']):
        flt_im = im[0].header['D%03dDATA' %(i+1)].split('[')[0]
        print('Update WCS: %s' %(flt_im))
        flt = pyfits.open(flt_im, mode='update')
        for ext in [1,2]:
            flt[ext].header['CRVAL1'] -= rshift
            flt[ext].header['CRVAL2'] -= dshift

        flt.flush()
        
    #return rshift, dshift
    
def runTweakReg(asn_file='GOODS-S-15-F140W_asn.fits', master_catalog='goodss_radec.dat', final_scale=0.06, ACS=False, threshold=5, WFPC2=False):
    """
    Wrapper around tweakreg, generating source catalogs separately from 
    `findpars`.
    """
    import glob
    import shutil
    
    import drizzlepac
    from drizzlepac import tweakreg
    from stwcs import updatewcs
    
    import threedhst.prep_flt_astrodrizzle
    
    asn = threedhst.utils.ASNFile(asn_file)
    
    if ACS:
        NCHIP=2
        sci_ext = [1,4]
        wht_ext = [2,5]
        ext = 'flc'
        dext = 'crclean'
    else:
        NCHIP=1
        sci_ext = [1]
        wht_ext = [2]
        ext = 'flt'
        dext = 'flt'
    
    ### Generate CRCLEAN images
    for exp in asn.exposures:
        updatewcs.updatewcs('%s_%s.fits' %(exp, ext))
    
    has_crclean = True
    for exp in asn.exposures:
        has_crclean &= os.path.exists('%s_crclean.fits' %(exp))
    
    threedhst.showMessage('# exposures: %d' %(len(asn.exposures)))
    
    if not has_crclean: 
        if len(asn.exposures) == 1:
            drizzlepac.astrodrizzle.AstroDrizzle(asn_file, clean=False, context=False, preserve=False, skysub=True, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_cr_corr=False, driz_combine=True)
            shutil.copy('%s_%s.fits' %(asn.exposures[0], ext), '%s_crclean.fits' %(asn.exposures[0]))
        else:
            drizzlepac.astrodrizzle.AstroDrizzle(asn_file, clean=False, context=False, preserve=False, skysub=True, driz_separate=True, driz_sep_wcs=True, median=True, blot=True, driz_cr=True, driz_cr_corr=True, driz_combine=True)
        
    #### Make SExtractor source catalogs in *each* flt
    for exp in asn.exposures:
        #updatewcs.updatewcs('%s_%s.fits' %(exp, ext))
        for i in range(NCHIP):
            se = threedhst.sex.SExtractor()
            se.options['WEIGHT_IMAGE'] = '%s_%s.fits[%d]' %(exp, dext, wht_ext[i]) #-1)
            se.options['WEIGHT_TYPE'] = 'MAP_RMS'
            #
            se.params['X_IMAGE'] = True; se.params['Y_IMAGE'] = True
            se.params['MAG_AUTO'] = True
            #
            se.options['CATALOG_NAME'] = '%s_%s_%d.cat' %(exp, ext, sci_ext[i])
            se.options['FILTER'] = 'N'
            se.options['DETECT_THRESH'] = '%f' %(threshold)
            se.options['ANALYSIS_THRESH'] = '%f' %(threshold)
            #
            se_img = '%s_%s.fits[%d]' %(exp, dext, sci_ext[i]) #-1)
            print(se_img)
            se.sextractImage(se_img)
            threedhst.sex.sexcatRegions('%s_%s_%d.cat' %(exp, ext, sci_ext[i]), '%s_%s_%d.reg' %(exp, ext, sci_ext[i]), format=1)
    
    #### TweakReg catfile
    asn_root = asn_file.split('_asn')[0]
    catfile = '%s.catfile' %(asn_root)
    fp = open(catfile,'w')
    ncat = 0
    for exp in asn.exposures:
        c = catIO.Table('%s_%s_%d.cat' %(exp, ext, sci_ext[i]))
        if len(c) < 20:
            continue
        else:
            ncat += 1
            
        line = '%s_%s.fits' %(exp, ext)
        for i in range(NCHIP):
            line += ' %s_%s_%d.cat' %(exp, ext, sci_ext[i])
        
        fp.write(line + '\n')
    
    fp.close()
    
    #### First run AstroDrizzle mosaic
    #drizzlepac.astrodrizzle.AstroDrizzle(asn_file, clean=True, context=False, preserve=False, skysub=True, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_combine=True)
    
    #### Make room for TWEAK wcsname
    for exp in asn.exposures:
        threedhst.prep_flt_astrodrizzle.clean_wcsname(flt='%s_%s.fits' %(exp, ext), wcsname='TWEAK', ACS=ACS)
    
    #### Main run of TweakReg
    if ACS:
        refimage = '%s_drc_sci.fits' %(asn_root)
    else:
        refimage = '%s_drz_sci.fits' %(asn_root)
        
    if ncat >= 1:
        tweakreg.TweakReg(asn_file, refimage=refimage, updatehdr=True, updatewcs=True, catfile=catfile, xcol=2, ycol=3, xyunits='pixels', refcat=master_catalog, refxcol=1, refycol=2, refxyunits='degrees', shiftfile=True, outshifts='%s_shifts.txt' %(asn_root), outwcs='%s_wcs.fits' %(asn_root), searchrad=5, tolerance=12, wcsname='TWEAK', interactive=False, residplot='No plot', see2dplot=False, clean=True, headerlet=True, clobber=True)
    
    #### Run AstroDrizzle again
    if ACS:
        if len(asn.exposures) == 1:
            drizzlepac.astrodrizzle.AstroDrizzle(asn_file, clean=True, final_scale=final_scale, final_pixfrac=0.8, context=False, resetbits=4096, final_bits=576, preserve=False, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_cr_corr=False)
        else:
            drizzlepac.astrodrizzle.AstroDrizzle(asn_file, clean=True, final_scale=final_scale, final_pixfrac=0.8, context=False, resetbits=4096, final_bits=576, preserve=False)
            
    else:
        if len(asn.exposures) == 1:
            drizzlepac.astrodrizzle.AstroDrizzle(asn_file, clean=True, final_scale=final_scale, final_pixfrac=0.8, context=False, resetbits=4096, final_bits=576, driz_sep_bits=576, preserve=False, driz_cr_snr='8.0 5.0', driz_cr_scale = '2.5 0.7', driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_cr_corr=False) # , final_wcs=True, final_rot=0)
        else:
            drizzlepac.astrodrizzle.AstroDrizzle(asn_file, clean=True, final_scale=final_scale, final_pixfrac=0.8, context=False, resetbits=4096, final_bits=576, driz_sep_bits=576, preserve=False, driz_cr_snr='8.0 5.0', driz_cr_scale = '2.5 0.7') # , final_wcs=True, final_rot=0)
        
    for exp in asn.exposures:
        files=glob.glob('%s*coo' %(exp))
        files.extend(glob.glob('%s*crclean.fits' %(exp)))
        for file in files:
            os.remove(file)
    
def align_drizzled(images=['MACS2129-35-F814W_drc_sci.fits', 'MACS2129-36-F814W_drc_sci.fits']):
    
    from astropy.table import Table as table
    from drizzlepac import astrodrizzle, tweakreg, tweakback
    
    for image in images:
        root = image.split('_sci.fits')[0]
        se = threedhst.sex.SExtractor()
        se.options['WEIGHT_IMAGE'] = '%s_wht.fits[0]' %(root)
        se.options['WEIGHT_TYPE'] = 'MAP_WEIGHT'
        #
        se.params['X_IMAGE'] = True; se.params['Y_IMAGE'] = True
        se.params['MAG_AUTO'] = True
        #
        se.options['CATALOG_NAME'] = '%s_sci.cat' %(root)
        se.options['FILTER'] = 'N'
        se.options['DETECT_THRESH'] = '5'
        se.options['ANALYSIS_THRESH'] = '5'
        #
        se.sextractImage('%s_sci.fits[0]' %(root))
        threedhst.sex.sexcatRegions('%s_sci.cat' %(root), '%s_sci.reg' %(root), format=1)
        #
        t = table.read('%s_sci.cat' %(root), format='ascii.sextractor')
        np.savetxt('%s_sci.xy' %(root), np.array([t['X_IMAGE'], t['Y_IMAGE']]).T, fmt='%.7f')
        fp = open('%s_sci.catfile' %(root), 'w')
        fp.write('%s_sci.fits %s_sci.xy\n' %(root, root))
        fp.close()
        
    reference = '%s_sci.xy' %(images[0].split('_sci.fits')[0])
    
    for image in images[1:]:
        root = image.split('_sci.fits')[0]
        tweakreg.TweakReg(image, refimage=images[0], updatehdr=True, updatewcs=True, catfile='%s_sci.catfile' %(root), xcol=1, ycol=2, xyunits='pixels', refcat=reference, refxcol=1, refycol=2, refxyunits='pixels', shiftfile=False, searchrad=5, tolerance=12, wcsname='TWEAK3', interactive=False, residplot='No plot', see2dplot=False, clean=True, headerlet=True, clobber=True)
        tweakback.tweakback(image)
    
    pass
    
def subtract_flt_background(root='GOODN-N1-VBA-F105W', scattered_light=False, sex_background=False, order=2):
    """
    Subtract polynomial background
    """
    import scipy.optimize
    
    import astropy.units as u
    
    from astropy.table import Table as table
    
    import stwcs
    from stwcs import updatewcs
    
    import drizzlepac
    from drizzlepac import astrodrizzle, tweakreg, tweakback
    
    import threedhst
    
    asn = threedhst.utils.ASNFile(root+'_asn.fits')
    for exp in asn.exposures:
        updatewcs.updatewcs('%s_%s.fits' %(exp, 'flt'))

    if not os.path.exists('%s_drz_sci.fits' %(root)):        
        if len(asn.exposures) == 1:
            drizzlepac.astrodrizzle.AstroDrizzle(root+'_asn.fits', clean=False, context=False, preserve=False, skysub=True, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_cr_corr=False, driz_combine=True)
        else:
            drizzlepac.astrodrizzle.AstroDrizzle(root+'_asn.fits', clean=False, context=False, preserve=False, skysub=True, driz_separate=True, driz_sep_wcs=True, median=True, blot=True, driz_cr=True, driz_cr_corr=True, driz_combine=True)
    
    se = threedhst.sex.SExtractor()
    se.options['WEIGHT_IMAGE'] = '%s_drz_wht.fits' %(root)
    se.options['WEIGHT_TYPE'] = 'MAP_WEIGHT'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION,BACKGROUND'
    se.options['CHECKIMAGE_NAME'] = '%s_drz_seg.fits,%s_drz_bkg.fits' %(root, root)
    se.options['BACK_TYPE'] = 'AUTO'
    se.options['BACK_SIZE'] = '256'
    #
    se.params['X_IMAGE'] = True; se.params['Y_IMAGE'] = True
    se.params['MAG_AUTO'] = True
    #
    se.options['CATALOG_NAME'] = '%s_drz_sci.cat' %(root)
    se.options['FILTER'] = 'Y'
    se.copyConvFile()
    se.options['FILTER_NAME'] = 'gauss_4.0_7x7.conv'
    se.options['DETECT_THRESH'] = '0.8'
    se.options['ANALYSIS_THRESH'] = '0.8'
    #
    se.options['MEMORY_OBJSTACK'] = '30000'
    se.options['MEMORY_PIXSTACK'] = '3000000'
    se.options['MEMORY_BUFSIZE'] = '2048'
    
    se.sextractImage('%s_drz_sci.fits' %(root))
    #threedhst.sex.sexcatRegions('%s_flt.cat' %(exp), '%s_flt.reg' %(exp), format=1)
    
    #### Blot segmentation map to FLT images for object mask
    asn = threedhst.utils.ASNFile('%s_asn.fits' %(root))
    
    #print 'Read files...'
    ref = pyfits.open('%s_drz_sci.fits' %(root))
    ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=0)

    seg = pyfits.open('%s_drz_seg.fits' %(root))    
    #### Fill ref[0].data with zeros for seg mask
    #seg_data = ref[0].data
    #seg_data[seg[0].data == 0] = 0
    seg_data = np.cast[np.float32](seg[0].data)
    
    bkg_data = pyfits.open('%s_drz_bkg.fits' %(root))[0].data
      
    yi, xi = np.indices((1014,1014))
    if scattered_light:        
        bg_components = np.ones((4,1014,1014))
        bg_components[1,:,:] = xi/1014.*2
        bg_components[2,:,:] = yi/1014.*2
        bg_components[3,:,:] = pyfits.open(os.getenv('THREEDHST') + '/CONF/G141_scattered_light.fits')[0].data
        #### Use flat-field itself for images affected by full-field 
        #### persistence from the tungsten lamp
        if scattered_light == 2:
            bg_components[3,:,:] = pyfits.open(os.getenv('iref') + 'flat_UDF_F140W_v0.fits')[1].data[5:-5,5:-5]
            
        NCOMP=4
    else:
        # bg_components = np.ones((3,1014,1014))
        # bg_components[1,:,:] = xi/1014.*2
        # bg_components[2,:,:] = yi/1014.*2
        # NCOMP=3
        #
        if order == 2:
            NCOMP=6
            bg_components = np.ones((NCOMP,1014,1014))
            bg_components[1,:,:] = (xi-507)/507.
            bg_components[2,:,:] = (yi-507)/507.
            bg_components[3,:,:] = ((xi-507)/507.)**2
            bg_components[4,:,:] = ((yi-507)/507.)**2
            bg_components[5,:,:] = (xi-507)*(yi-507)/507.**2
        else:
            NCOMP=3
            bg_components = np.ones((NCOMP,1014,1014))
            bg_components[1,:,:] = (xi-507)/507.
            bg_components[2,:,:] = (yi-507)/507.
            
    bg_flat = bg_components.reshape((NCOMP,1014**2))
    
    #### Loop through FLTs, blotting reference and segmentation
    models = []
    for exp in asn.exposures:
        flt = pyfits.open('%s_flt.fits' %(exp)) #, mode='update')
        flt_wcs = stwcs.wcsutil.HSTWCS(flt, ext=1)
        
        ### segmentation        
        print('Segmentation image: %s_blot.fits' %(exp))
        blotted_seg = astrodrizzle.ablot.do_blot(seg_data+0, ref_wcs, flt_wcs, 1, coeffs=True, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)
        
        blotted_bkg = 0.
        if sex_background:
            blotted_bkg = astrodrizzle.ablot.do_blot(bkg_data+0, ref_wcs, flt_wcs, 1, coeffs=True, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)
            flt[1].data -= blotted_bkg
            
        mask = (blotted_seg == 0) & (flt['DQ'].data == 0) & (flt[1].data > -1) & (xi > 10) & (yi > 10) & (xi < 1004) & (yi < 1004) 
        mask &= np.isfinite(flt[1].data) & np.isfinite(flt[2].data)
        mask &= (flt[1].data < 5*np.median(flt[1].data[mask]))
        data_range = np.percentile(flt[1].data[mask], [2.5, 97.5])
        mask &= (flt[1].data >= data_range[0]) & (flt[1].data <= data_range[1])
        data_range = np.percentile(flt[2].data[mask], [0.5, 99.5])
        mask &= (flt[2].data >= data_range[0]) & (flt[2].data <= data_range[1])
        
        ### Least-sq fit for component normalizations
        data = flt[1].data[mask].flatten()
        wht = (1./flt[2].data[mask].flatten())**2
        templates = bg_flat[:, mask.flatten()]
        p0 = np.zeros(NCOMP)
        p0[0] = np.median(data)
        obj_fun = threedhst.grism_sky.obj_lstsq
        print('XXX: %d' %(mask.sum()))
        popt = scipy.optimize.leastsq(obj_fun, p0, args=(data, templates, wht), full_output=True, ftol=1.49e-8/1000., xtol=1.49e-8/1000.)
        xcoeff = popt[0]
        model = np.dot(xcoeff, bg_flat).reshape((1014,1014))
        models.append(model)
        
        # add header keywords of the fit components
        flt = pyfits.open('%s_flt.fits' %(exp), mode='update')
        flt[1].data -= model+blotted_bkg
        for i in range(NCOMP):
            if 'BGCOMP%d' %(i+1) in flt[0].header:
                flt[0].header['BGCOMP%d' %(i+1)] += xcoeff[i]
            else:
                flt[0].header['BGCOMP%d' %(i+1)] = xcoeff[i]                
        
        flt.flush()
        coeff_str = '  '.join(['%.4f' %c for c in xcoeff])
        threedhst.showMessage('Background subtraction, %s_flt.fits:\n\n  %s' %(exp, coeff_str))
        
    return models
    
def copy_adriz_headerlets(direct_asn='GOODS-S-15-F140W_asn.fits', grism_asn='GOODS-S-15-G141_asn.fits', order=None, force=False, ACS=False):
    """
    Copy Tweaked WCS solution in direct image to the paired grism exposures.
    
    If same number of grism as direct exposures, match the WCS headers
    directly.  If not, just get the overall shift from the first direct 
    exposure and apply that to the grism exposures.
    """
    import stwcs
    from stwcs import updatewcs
    import drizzlepac
    
    direct = threedhst.utils.ASNFile(direct_asn)
    grism = threedhst.utils.ASNFile(grism_asn)
    
    Nd = len(direct.exposures)
    Ng = len(grism.exposures)
    
    if ACS:
        NCHIP=2
        sci_ext = [1,4]
        ext = 'flc'
    else:
        NCHIP=1
        sci_ext = [1]
        ext = 'flt'
        
    if Nd == Ng:
        if order is None:
            order = list(range(Nd))
            
        for i in range(Nd):
            imd = pyfits.open('%s_%s.fits' %(direct.exposures[i], ext))
            #img = pyfits.open('%s_%s.fits' %(grism.exposures[i]))
            #
            for sci in sci_ext:
                #sci_ext=1
                direct_WCS = stwcs.wcsutil.HSTWCS(imd, ext=sci)
                #
                drizzlepac.updatehdr.update_wcs('%s_%s.fits' %(grism.exposures[order[i]], ext), sci, direct_WCS, verbose=True)    
    else:
        #### Get overall shift from a shift-file and apply it to the 
        #### grism exposures
        sf = threedhst.shifts.ShiftFile(direct_asn.replace('_asn.fits', '_shifts.txt'))
        imd = pyfits.open(direct_asn.replace('asn','wcs'))
        print(imd.filename())
        direct_WCS = stwcs.wcsutil.HSTWCS(imd, ext='wcs')
        #
        for i in range(Ng):
            img = pyfits.open('%s_%s.fits' %(grism.exposures[i], ext))
            if 'WCSNAME' in img[1].header:
                if img[1].header['WCSNAME'] == 'TWEAK':
                    if force is False:
                        threedhst.showMessage('"TWEAK" WCS already found in %s_flt.fits.\nRun copy_adriz_headerlets with force=True to force update the shifts' %(grism.exposures[i]), warn=True)
                        continue
            #
            updatewcs.updatewcs('%s_%s.fits' %(grism.exposures[i], ext))
            drizzlepac.updatehdr.updatewcs_with_shift('%s_%s.fits' %(grism.exposures[i], ext), direct_WCS, rot=sf.rotate[0], scale=sf.scale[0], xsh=sf.xshift[0], ysh=sf.yshift[0], wcsname='TWEAK')
  
def subtract_acs_grism_background(asn_file='RXJ2248-08-G800L_asn.fits', final_scale=None):

    import glob
    import os
    
    import drizzlepac
    from drizzlepac import tweakreg
    from stwcs import updatewcs
    
    import threedhst.prep_flt_astrodrizzle
    
    ### First pass to flag CRs and make crclean images
    drizzlepac.astrodrizzle.AstroDrizzle(asn_file, clean=False, context=False, preserve=False, skysub=True, driz_separate=True, driz_sep_wcs=True, median=True, blot=True, driz_cr=True, driz_cr_corr=True, driz_combine=True)
    
    ### add mdrizsky back to values, fit on crclean but subtract from flc
    asn = threedhst.utils.ASNFile(asn_file)
    
    try:
        sky1 = pyfits.open(os.getenv('THREEDHST') + '/CONF/ACS.WFC.CHIP1.msky.1.smooth.fits')
        sky2 = pyfits.open(os.getenv('THREEDHST') + '/CONF/ACS.WFC.CHIP2.msky.1.smooth.fits')    
    except:
        sky1 = pyfits.open(os.getenv('THREEDHST') + '/CONF/ACS.WFC.CHIP1.msky.1.fits')
        sky2 = pyfits.open(os.getenv('THREEDHST') + '/CONF/ACS.WFC.CHIP2.msky.1.fits')    

    skies = [sky1, sky2]
    extensions = [1,4] ### SCI extensions
    
    for exp in asn.exposures:
        flt = pyfits.open(exp+'_flc.fits', mode='update')
        #crc = pyfits.open(exp+'_crclean.fits')
        ### Loop through ACS chips
        for j in [0,1]:
            ext = extensions[j]
            #
            exptime = flt[0].header['EXPTIME']
            mask = flt['dq', j+1].data == 0
            ratio = flt['sci', j+1].data/exptime/skies[j][0].data
            med = np.median(ratio[mask])
            #
            mask2 = mask & ((flt['sci',j+1].data-med*exptime)/flt['err',j+1].data < 3)
            med2 = np.median(ratio[mask2])
            #
            flt['sci', j+1].data -= med2*exptime*skies[j][0].data     
            flt['sci', j+1].header['MDRIZSKY'] = 0. #.update('MDRIZSKY', med2*exptime)
            print('%s chip%d: %.3f' %(exp, j+1, med2))
        #
        ### Write updates
        flt.flush()
    
    drizzlepac.astrodrizzle.AstroDrizzle(asn_file, clean=True, skysub=False, skyuser='MDRIZSKY', final_scale=final_scale, final_pixfrac=0.8, context=False, resetbits=4096, final_bits=576, preserve=False) # , final_wcs=True, final_rot=0)
    
def subtract_grism_background(asn_file='GDN1-G102_asn.fits', PATH_TO_RAW='../RAW/', final_scale=0.06, visit_sky=True, column_average=True, mask_grow=18, first_run=True, sky_iter=1):
    """
    Subtract master grism sky from FLTs
    """
    import os
    import scipy.ndimage as nd
    import pyregion
    
    from drizzlepac import astrodrizzle
    import drizzlepac
    
    from stwcs import updatewcs
    import stwcs
    
    import threedhst.grism_sky as bg
    
    asn = threedhst.utils.ASNFile(asn_file)
    root = asn_file.split('_asn')[0]
            
    sky_images = {'G141':['zodi_G141_clean.fits', 'excess_lo_G141_clean.fits', 'G141_scattered_light.fits'],
                  'G102':['zodi_G102_clean.fits', 'excess_G102_clean.fits']}
    #
    # sky_images = {'G141':['zodi_G141_clean.fits', 'excess_lo_G141_clean.fits', 'G141_scattered_light_v2.fits'],
    #               'G102':['zodi_G102_clean.fits', 'excess_G102_clean.fits']}
    
    # ### Don't use scattered light
    # sky_images = {'G141':['zodi_G141_clean.fits', 'excess_lo_G141_clean.fits'],
    #               'G102':['zodi_G102_clean.fits', 'excess_G102_clean.fits']}
    # 
    # ## Use aXe images
    # sky_images = {'G141':['WFC3.IR.G141.sky.V1.0.flat.fits', 'WFC3.IR.G141.sky.V1.0.flat.fits'],
    #               'G102':['zodi_G102_clean.fits', 'excess_G102_clean.fits']}
    
    if first_run:
        ### Rough background subtraction
        threedhst.process_grism.fresh_flt_files(asn_file, from_path=PATH_TO_RAW, preserve_dq=False)
        flt = pyfits.open('%s_flt.fits' %(asn.exposures[0]))
        GRISM = flt[0].header['FILTER']
        bg.set_grism_flat(grism=GRISM, verbose=True)
    
        zodi = pyfits.open(os.getenv('THREEDHST')+'/CONF/%s' %(sky_images[GRISM][0]))[0].data
    
        for exp in asn.exposures:
            updatewcs.updatewcs('%s_flt.fits' %(exp))
            flt = pyfits.open('%s_flt.fits' %(exp), mode='update')
            #flt = pyfits.open('%s_flt.fits' %(exp))
            flt[1].data *= bg.flat
            #
            mask = (flt['DQ'].data == 0)
            data_range = np.percentile(flt[1].data[mask], [20, 80])
            mask &= (flt[1].data >= data_range[0]) & (flt[1].data <= data_range[1]) & (flt[2].data != 0) & np.isfinite(flt[1].data) & np.isfinite(flt[2].data)
            ### Least-sq fit for component normalizations
            data = flt[1].data[mask].flatten()
            wht = (1./flt[2].data[mask].flatten())**2
            zodi_mask = zodi[mask].flatten()
            coeff_zodi = np.sum(data*zodi_mask*wht)/np.sum(zodi_mask**2*wht)
            flt[1].data -= zodi*coeff_zodi
            flt.flush()
            threedhst.showMessage('Rough background for %s (zodi): %0.4f' %(exp, coeff_zodi))
            #templates = bg_flat[:, mask.flatten()]
        
        ### Run astrodrizzle to make DRZ mosaic, grism-SExtractor mask
        drizzlepac.astrodrizzle.AstroDrizzle(asn_file, clean=True, context=False, preserve=False, skysub=True, driz_separate=True, driz_sep_wcs=True, median=True, blot=True, driz_cr=True, driz_combine=True, final_wcs=True, final_rot=0, resetbits=4096, final_bits=576, driz_sep_bits=576, driz_cr_snr='8.0 5.0', driz_cr_scale = '2.5 0.7')
                
    else:
        flt = pyfits.open('%s_flt.fits' %(asn.exposures[0]))
        GRISM = flt[0].header['FILTER']
        bg.set_grism_flat(grism=GRISM, verbose=True)
    
        
    se = threedhst.sex.SExtractor()
    se.options['WEIGHT_IMAGE'] = '%s_drz_wht.fits' %(root)
    se.options['WEIGHT_TYPE'] = 'MAP_WEIGHT'
    se.options['CHECKIMAGE_TYPE'] = 'SEGMENTATION'
    se.options['CHECKIMAGE_NAME'] = '%s_drz_seg.fits' %(root)
    #
    se.params['X_IMAGE'] = True; se.params['Y_IMAGE'] = True
    se.params['MAG_AUTO'] = True
    #
    se.options['CATALOG_NAME'] = '%s_drz_sci.cat' %(root)
    se.options['FILTER'] = 'Y'
    se.copyConvFile(grism=True)
    se.options['FILTER_NAME'] = 'grism.conv'
    se.options['DETECT_THRESH'] = '0.7'
    se.options['ANALYSIS_THRESH'] = '0.7'
    #
    se.sextractImage('%s_drz_sci.fits' %(root))
    
    #### Blot segmentation map to FLT images for object mask
    ref = pyfits.open('%s_drz_sci.fits' %(root))
    #ref_wcs = stwcs.wcsutil.HSTWCS(ref, ext=0)

    seg = pyfits.open('%s_drz_seg.fits' %(root))
    seg_data = np.cast[np.float32](seg[0].data)
    ref_wcs = stwcs.wcsutil.HSTWCS(seg, ext=0)
            
    #### Loop through FLTs, blotting reference and segmentation
    threedhst.showMessage('%s: Blotting grism segmentation masks.' %(root))
        
    for exp in asn.exposures:
        flt = pyfits.open('%s_flt.fits' %(exp))
        flt_wcs = stwcs.wcsutil.HSTWCS(flt, ext=1)
        ### segmentation
        #print 'Segmentation image: %s_blot.fits' %(exp)
        blotted_seg = astrodrizzle.ablot.do_blot(seg_data+0, ref_wcs, flt_wcs, 1, coeffs=True, interp='nearest', sinscl=1.0, stepsize=10, wcsmap=None)
        seg_grow = nd.maximum_filter((blotted_seg > 0)*1, size=8)
        pyfits.writeto('%s_flt.seg.fits' %(exp), header=flt[1].header, data=seg_grow, clobber=True)
        
    if first_run:
        ### Run background subtraction scripts
        threedhst.process_grism.fresh_flt_files(asn_file, from_path=PATH_TO_RAW, preserve_dq=False)
        for exp in asn.exposures:
            updatewcs.updatewcs('%s_flt.fits' %(exp))
            #threedhst.grism_sky.remove_grism_sky(flt=exp+'_flt.fits', list=sky_images[GRISM], path_to_sky=os.getenv('THREEDHST')+'/CONF/', verbose=True, second_pass=True, overall=True)
    
    if visit_sky:
        threedhst.grism_sky.remove_visit_sky(asn_file=asn_file, list=sky_images[GRISM], add_constant=False, column_average=(column_average) & (sky_iter == 1), mask_grow=mask_grow, flat_correct=first_run)
        if (sky_iter > 1) & (~first_run):
            for i in range(1, sky_iter):
                threedhst.grism_sky.remove_visit_sky(asn_file=asn_file, list=sky_images[GRISM], add_constant=False, column_average=column_average & (i == (sky_iter-1)), mask_grow=mask_grow, flat_correct=False)
    else:
        for exp in asn.exposures:
            threedhst.grism_sky.remove_grism_sky(flt='%s_flt.fits' %(exp), list=sky_images[GRISM],  path_to_sky = os.getenv('THREEDHST')+'/CONF/', out_path='./', verbose=False, plot=False, flat_correct=first_run, sky_subtract=True, second_pass=column_average, overall=True, combine_skies=False, sky_components=True, add_constant=False)
            
    ### Astrodrizzle again to reflag CRs and make cleaned mosaic
    drizzlepac.astrodrizzle.AstroDrizzle(asn_file, clean=True, skysub=False, skyuser='MDRIZSKY', final_wcs=True, final_scale=final_scale, final_pixfrac=0.8, context=False, resetbits=4096, final_bits=576, driz_sep_bits=576, preserve=False, driz_cr_snr='8.0 5.0', driz_cr_scale='2.5 0.7') # , final_wcs=True, final_rot=0)
    
    
def get_vizier_cat(image='RXJ2248-IR_sci.fits', ext=0, catalog="II/246"):
    """
    Get a list of RA/Dec coords from a Vizier catalog that can be used
    for WCS alignment.
    
    `catalog` is any catalog ID recognized by Vizier, e.g.: 
        "II/328/allwise": WISE
        "II/246": 2MASS
        "V/139": SDSS DR9
    """
    import threedhst.dq
    import astropy.wcs as pywcs
    from astropy.table import Table as table
    import astropy.io.fits as pyfits
    
    import astroquery
    from astroquery.vizier import Vizier
    import astropy.coordinates as coord
    import astropy.units as u
    
    im = pyfits.open(image)
    
    wcs = pywcs.WCS(im[ext].header)
    #wcs = pywcs.WCS(pyfits.getheader('Q0821+3107-F140W_drz.fits', 1))

    Vizier.ROW_LIMIT = -1
            
    r0, d0 = wcs.wcs_pix2world([[im[ext].header['NAXIS1']/2., im[ext].header['NAXIS2']/2.]], 1)[0]
    foot = wcs.calc_footprint()
    
    corner_radius = np.sqrt((foot[:,0]-r0)**2/np.cos(d0/360.*2*np.pi)**2 + (foot[:,1]-d0)**2).max()*60*1.1

    try:
        c = coord.ICRS(ra=r0, dec=d0, unit=(u.deg, u.deg))
    except:
        c = coord.ICRSCoordinates(ra=r0, dec=d0, unit=(u.deg, u.deg))
        
    #### something with astropy.coordinates
    # c.icrs.ra.degree = c.icrs.ra.degrees
    # c.icrs.dec.degree = c.icrs.dec.degrees
    #
    vt = Vizier.query_region(c, radius=u.Quantity(corner_radius, u.arcminute), catalog=[catalog])
    if not vt:
        threedhst.showMessage('No matches found in Vizier %s @ (%.6f, %.6f).\n\nhttp://vizier.u-strasbg.fr/viz-bin/VizieR?-c=%.6f+%.6f&-c.rs=8' %(catalog, r0, d0, r0, d0), warn=True)
        return False
    
    vt = vt[0]
            
    #### Make a region file
    ra_list, dec_list = vt['RAJ2000'], vt['DEJ2000']
    print('Vizier, found %d objects in %s.' %(len(ra_list), catalog))
    
    fp = open('%s.vizier.radec' %(image.split('.fits')[0]), 'w')
    fpr = open('%s.vizier.reg' %(image.split('.fits')[0]), 'w')
    
    fp.write('# %s, r=%.1f\'\n' %(catalog, corner_radius))
    fpr.write('# %s, r=%.1f\'\nfk5\n' %(catalog, corner_radius))
    for ra, dec in zip(ra_list, dec_list):
        fp.write('%.7f %.7f\n' %(ra, dec))
        fpr.write('circle(%.6f, %.6f, 0.5")\n' %(ra, dec))
    
    fpr.close()
    fp.close()
    
    return True
    
    
def prep_direct_grism_pair(direct_asn='goodss-34-F140W_asn.fits', grism_asn='goodss-34-G141_asn.fits', radec=None, raw_path='../RAW/', mask_grow=18, scattered_light=False, final_scale=None, skip_direct=False, ACS=False, jump=False, order=2, get_shift=True, align_threshold=20, column_average=True, sky_iter=3, run_acs_lacosmic=False):
    """
    Process both the direct and grism observations of a given visit
    """
    import threedhst.prep_flt_astrodrizzle as prep
    import drizzlepac
    from stwcs import updatewcs
    
    import time
    
    t0 = time.time()
    
    #direct_asn='goodss-34-F140W_asn.fits'; grism_asn='goodss-34-G141_asn.fits'; radec=None; raw_path='../RAW/'
    #radec = os.getenv('THREEDHST') + '/ASTRODRIZZLE_FLT/Catalog/goodss_radec.dat'
    
    ################################
    #### Direct image processing
    ################################
    
    #### xx add astroquery 2MASS/SDSS workaround for radec=None
    
    if not skip_direct:

        #### Get fresh FLTS from ../RAW/
        asn = threedhst.utils.ASNFile(direct_asn)
        if ACS:
            for exp in asn.exposures:
                print('cp %s/%s_flc.fits.gz .' %(raw_path, exp))
                os.system('cp %s/%s_flc.fits.gz .' %(raw_path, exp))
                os.system('gunzip -f %s_flc.fits.gz' %(exp))
                
                if run_acs_lacosmic:
                    try:
                        import lacosmicx
                        status = True
                    except:
                        print('import lacosmicx failed!')
                        status = False
                    
                    if status:
                        im = pyfits.open('%s_flc.fits' %(exp), mode='update')
                        for ext in [1,2]:
                            indata = im['SCI',ext].data
                            #inmask = im['DQ',ext].data > 0

                            if im['SCI',ext].header['BUNIT'] == 'ELECTRONS':
                                gain = 1
                            else:
                                gain = 1./im[0].header['EXPTIME']

                            if 'MDRIZSK0' in im['SCI',ext].header:
                                pssl = im['SCI',ext].header['MDRIZSK0']
                            else:
                                pssl = 0.
                            
                            if 'FLASHLVL' in im[0].header:
                                pssl += im[0].header['FLASHLVL']
                                sig_scale = 1.8
                            else:
                                sig_scale = 1.
                                
                            out = lacosmicx.lacosmicx(indata, inmask=None, 
                                    sigclip=3.5*sig_scale, sigfrac=0.2,
                                    objlim=7.0, gain=gain,
                                    readnoise=im[0].header['READNSEA'], 
                                    satlevel=np.inf, pssl=pssl, niter=5,
                                    sepmed=True, cleantype='meanmask',
                                    fsmode='median', psfmodel='gauss',
                                    psffwhm=2.5,psfsize=7, psfk=None,
                                    psfbeta=4.765, verbose=True)
                        
                            crmask, cleanarr  = out
                            im['DQ',ext].data |= 16*crmask
                            
                            ### Low pixels
                            if im[0].header['INSTRUME'] == 'WFC3':
                                bad = im['SCI',ext].data < -4*im['ERR',ext].data
                                im['DQ',ext].data |= 16*bad
                                
                        im.flush()
                        
        else:
            threedhst.process_grism.fresh_flt_files(direct_asn, from_path=raw_path)
        
        if (not ACS):
            #### Subtract WFC3/IR direct backgrounds
            prep.subtract_flt_background(root=direct_asn.split('_asn')[0], scattered_light=scattered_light, order=order)
            #### Flag IR CRs again within runTweakReg
        
        #### Run TweakReg
        if (radec is None) & (not ACS):
            print(len(asn.exposures))
            
            if len(asn.exposures) > 1:
                drizzlepac.astrodrizzle.AstroDrizzle(direct_asn, clean=True, final_scale=None, final_pixfrac=0.8, context=False, final_bits=576, preserve=False, driz_cr_snr='5.0 4.0', driz_cr_scale = '2.5 0.7') 
            else:
                drizzlepac.astrodrizzle.AstroDrizzle(direct_asn, clean=True, final_scale=None, final_pixfrac=1, context=False, final_bits=576, preserve=False, driz_separate=False, driz_sep_wcs=False, median=False, blot=False, driz_cr=False, driz_cr_corr=False, driz_combine=True) 
        else:
            if get_shift:
                prep.runTweakReg(asn_file=direct_asn, master_catalog=radec, final_scale=None, ACS=ACS, threshold=align_threshold)
        
        #### Subtract background of direct ACS images
        if ACS:
            for exp in asn.exposures:
                flc = pyfits.open('%s_flc.fits' %(exp), mode='update')
                if 'SUB' in flc[0].header['APERTURE']:
                    extensions = [1]
                else:
                    extensions = [1,4]
                    
                for ext in extensions:
                    threedhst.showMessage('Subtract background from %s_flc.fits[%d] : %.4f' %(exp, ext, flc[ext].header['MDRIZSKY']))
                    flc[ext].data -= flc[ext].header['MDRIZSKY']
                    flc[ext].header['MDRIZSK0'] = flc[ext].header['MDRIZSKY']
                    flc[ext].header['MDRIZSKY'] = 0.
                #
                flc.flush()
        else:
            pass
            #### Do this later, gives segfaults here???
            #prep.subtract_flt_background(root=direct_asn.split('_asn')[0], scattered_light=scattered_light)
            #### Flag CRs again on BG-subtracted image
            #drizzlepac.astrodrizzle.AstroDrizzle(direct_asn, clean=True, final_scale=None, final_pixfrac=0.8, context=False, final_bits=576, preserve=False, driz_cr_snr='5.0 4.0', driz_cr_scale = '2.5 0.7') # ,
        
    ################################
    #### Grism image processing
    ################################
    
    if grism_asn:
        asn = threedhst.utils.ASNFile(grism_asn)
        if ACS:
            for exp in asn.exposures:
                print('cp %s/%s_flc.fits.gz .' %(raw_path, exp))
                os.system('cp %s/%s_flc.fits.gz .' %(raw_path, exp))
                os.system('gunzip -f %s_flc.fits.gz' %(exp))
                updatewcs.updatewcs('%s_flc.fits' %(exp))

            prep.copy_adriz_headerlets(direct_asn=direct_asn, grism_asn=grism_asn, ACS=True)
            prep.subtract_acs_grism_background(asn_file=grism_asn, final_scale=None)
        else:
            #### Remove the sky and flag CRs
            ## with mask from rough zodi-only subtraction
            prep.subtract_grism_background(asn_file=grism_asn, PATH_TO_RAW='../RAW/', final_scale=None, visit_sky=True, column_average=False, mask_grow=mask_grow, first_run=True)
            ## Redo making mask from better combined image
            prep.subtract_grism_background(asn_file=grism_asn, PATH_TO_RAW='../RAW/', final_scale=final_scale, visit_sky=True, column_average=column_average, mask_grow=mask_grow, first_run=False, sky_iter=sky_iter)
                        
            #### Copy headers from direct images
            if radec is not None:
                prep.copy_adriz_headerlets(direct_asn=direct_asn, grism_asn=grism_asn, ACS=False)
                #### Run CR rejection with final shifts
                drizzlepac.astrodrizzle.AstroDrizzle(grism_asn, clean=True, skysub=False, final_wcs=True, final_scale=final_scale, final_pixfrac=0.8, context=False, final_bits=576, driz_sep_bits=576, preserve=False, driz_cr_snr='8.0 5.0', driz_cr_scale='2.5 0.7') # driz_cr_snr='5.0 4.0', driz_cr_scale = '2.5 0.7')
                
    if not grism_asn:
        t1 = time.time()
        threedhst.showMessage('direct: %s\n\nDone (%d s).' %(direct_asn, int(t1-t0)))
    else:
        t1 = time.time()
        threedhst.showMessage('direct: %s\ngrism: %s\n\nDone (%d s).' %(direct_asn, grism_asn, int(t1-t0)))
    
def test():
    import threedhst.prep_flt_astrodrizzle as prep
    
    asn_file = 'HOPS183B-10-140-G141_asn.fits'
    #asn_file = 'HOPS183B-10-140-F160W_asn.fits'
    
    threedhst.process_grism.fresh_flt_files(asn_file, from_path='../RAW/', preserve_dq=False)

    asn = threedhst.utils.ASNFile(asn_file)
    for exp in asn.exposures:
        prep.subtract_fixed_background(flt_file='%s_flt.fits' %(exp), path_to_raw='../RAW/')
        
def subtract_fixed_background(flt_file='ickw10upq_flt.fits', path_to_raw='../RAW/'):
    
    import astropy.io.fits as pyfits
    from threedhst import catIO
    import glob
    
    import mywfc3.zodi
    
    zodi_scl = mywfc3.zodi.flt_zodi(flt_file, verbose=False)
    
    grism_sky_images = {'G141':['zodi_G141_clean.fits', 'excess_lo_G141_clean.fits'],
                        'G102':['zodi_G102_clean.fits', 'excess_G102_clean.fits']}
    
    grism_flat = {'G141':'uc721143i_pfl.fits',
                  'G102':'uc72113oi_pfl.fits'}
    
        
    flt = pyfits.open(flt_file, mode='update')
    
    if 'FIXEDBG' in list(flt[1].header.keys()):
        if flt[1].header['FIXEDBG']:
            print('Background already subtracted!')
            return True
        
    if 'G1' in flt[0].header['FILTER']:

        bg = np.loadtxt(glob.glob('%s/%s*ramp.dat' %(path_to_raw, flt_file.split('_')[0][:-1]))[0])

        time = bg[:,0]
        dt = np.append(time[0], np.diff(time))
        
        excess_scl = np.sum(dt*(bg[:,1]-bg[:,1].min()))/time[-1]
        
        flat = pyfits.open('%s/%s' %(os.getenv('iref'), grism_flat[flt[0].header['FILTER']]))
        flt[1].data /= flat[1].data[5:-5,5:-5]
        
        zodi_file = grism_sky_images[flt[0].header['FILTER']][0]
        zodi = pyfits.open('%s/%s' %(os.getenv('iref'), zodi_file))
        flt[1].data -= zodi[0].data * zodi_scl
        flt[0].header['GSKY01'] = (zodi_scl, 'Master sky: %s' %(zodi_file))
        
        excess_file = grism_sky_images[flt[0].header['FILTER']][1]
        excess = pyfits.open('%s/%s' %(os.getenv('iref'), excess_file))
        flt[1].data -= excess[0].data * excess_scl
        flt[0].header['GSKY02'] = (excess_scl, 'Master sky: %s' %(excess_file))
        
        print('%s zodi: %.2f, excess: %.2f' %(flt[0].header['FILTER'], flt[0].header['GSKY01'], flt[0].header['GSKY02']))
        
        flt[1].header['MDRIZSKY'] = 0.
        
    else:
        flt[1].header['MDRIZSKY'] = 0

        flt[1].header['ZODISKY'] = zodi_scl
        flt[1].data -= flt[1].header['ZODISKY']
        flt.flush()
        
        print('%s zodi: %.2f' %(flt[0].header['FILTER'], flt[1].header['ZODISKY']))
    
    flt[1].header['FIXEDBG'] = True
        
    flt.flush()
    
    return True
    