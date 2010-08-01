"""
3DHST.combine_phot

Combine grism spectrum with broad-band photometry from other sources.
"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import os
import matplotlib.pyplot as pyplot
import numpy as np
import numpy.random as nprand
import threedhst

## Global        
phot_cat = None
grism_cat = None
grism_SPC = None

def matchRaDec(ra0, dec0, ralist, declist):
    """
idx, dr = matchRaDec(ra0, dec0, ralist, declist)
    
    `idx`: index of closest match
    `dr`:  distance of match, in arcsec
    
    """
    cosdec = np.cos(np.float(dec0)/360*2*np.pi)
    dr = np.sqrt(((ralist-np.float(ra0))*cosdec)**2+
                  (declist-np.float(dec0))**2)*3600.
    mat = np.where(dr == np.min(dr))[0]
    return mat[0], dr[mat[0]]
    
###############################################################################                 
####                                                                       ####
####                             CDF-South                                 ####
####                                                                       ####
###############################################################################

def init_CDFS():
    """
init_CDFS
    
    Read catalog and SPC.
    """
    drz_dir = 'DRZREJ_G141'
    
    cp = threedhst.combine_phot
    
    cp.phot_cat = threedhst.utils.ReadASCIICat('FIREWORKS/FIREWORKS_phot.cat')
    root_direct = 'ib6o23020'
    cp.grism_cat = threedhst.sex.mySexCat('DATA/'+root_direct+'_drz.cat')
    cp.grism_SPC = threedhst.plotting.SPCFile(root_direct+'_2_opt.SPC.fits',
                                                axe_drizzle_dir=drz_dir)
    
    cp.get_sed = cp.sed_CDFS
    
def sed_CDFS(idx, flam=True):
    """
lc, flux, err = cdfs_sed(idx)
    
    >>> cat = threedhst.utils.ReadASCIICat('FIREWORKS_phot.cat')
    >>> lc, flux, err = cdfs_sed(cat, 100)
    
    """
    cat = phot_cat
    abzp = 23.86
    
    lc = np.array([3663.9438, 4327.7378, 4599.7564, 5378.7950, 5957.6070, 
                   6516.0485, 7705.5403, 8658.7066,  9052.0768, 12378.455,
                   16517.513, 21681.252, 35634.260, 45110.187, 57593.123,
                   79594.189])
                   
    flux = np.array([cat.u38_colf[idx], cat.b435_colf[idx], cat.b_colf[idx],
                     cat.v_colf[idx], cat.v606_colf[idx], cat.r_colf[idx],
                     cat.i775_colf[idx], cat.i_colf[idx], cat.z850_colf[idx], 
                     cat.j_colf[idx], cat.h_colf[idx], cat.ks_colf[idx], 
                     cat.n3p6um_colf[idx], cat.n4p5um_colf[idx], 
                     cat.n5p8um_colf[idx], cat.n8p0um_colf[idx]])
    
    err = np.array([cat.u38_colfe[idx], cat.b435_colfe[idx], cat.b_colfe[idx],
                    cat.v_colfe[idx], cat.v606_colfe[idx], cat.r_colfe[idx],
                    cat.i775_colfe[idx], cat.i_colfe[idx], cat.z850_colfe[idx], 
                    cat.j_colfe[idx], cat.h_colfe[idx], cat.ks_colfe[idx], 
                    cat.n3p6um_colfe[idx], cat.n4p5um_colfe[idx], 
                    cat.n5p8um_colfe[idx], cat.n8p0um_colfe[idx]])
    
    ### zeropoint: microJy -> AB25
    # flux *= 10**(-0.4*(abzp-25))
    # err *= 10**(-0.4*(zbzp-25))
    if flam:
        fconvert = 10**(-0.4*(abzp+48.6))*3.e18/lc**2  
    else:
        fconvert = 10**(-0.4*(abzp+48.6))
        
    return lc, flux*fconvert, err*fconvert
    
def phot_grism(grism_idx=None, phot_idx=None):
    """
phot_grism(grism_idx=None, phot_idx=None)
    """
    from scipy import interpolate
    
    cp = threedhst.combine_phot
    
    #### Match object in grism and photometric catalog
    if grism_idx:
        phot_idx, dr = cp.matchRaDec(cp.grism_cat.ra[grism_idx], 
                                  cp.grism_cat.dec[grism_idx],
                                  cp.phot_cat.ra, cp.phot_cat.dec)
    else:
        if phot_idx:
            grism_idx, dr = cp.matchRaDec(cp.phot_cat.ra[phot_idx],
                                       cp.phot_cat.dec[phot_idx],
                                       cp.grism_cat.ra, cp.grism_cat.dec)
        else:
            print 'Neither `grism_idx` nor `phot_idx` are set.'
            return False
    
    if dr > 1.1:
        return False
        
    phot_z = threedhst.utils.ReadASCIICat('FIREWORKS/FIREWORKS_redshift.cat')
    
    outstr =  ' %d %d %5.2f %6.3f %6.3f\n' %(cp.grism_cat.id[grism_idx],
           cp.phot_cat.id[phot_idx],dr,
           phot_z.zsp[phot_idx], phot_z.zph_best[phot_idx])
    
    print outstr
    
    ### Get broad-band SED, grism spectrum
    lc, flux, err = cp.get_sed(phot_idx)
    
    grism = cp.grism_SPC.getSpec(cp.grism_cat.id[grism_idx])
    if grism is False:
        return False
        
    grism.corr = grism.FLUX-grism.CONTAM

    #### 
    wave_limits = [1.1e4, 1.68e4]
    mask = grism.LAMBDA*0.
    keep = np.where((grism.LAMBDA > wave_limits[0]) & 
                    (grism.LAMBDA < wave_limits[1]))[0]
    mask[keep] = 1.
    
    #### Compare detected emission line test redshifts to BB zphot/zspec
    lines = threedhst.spec1d.findLines(cp.grism_SPC, 
                        idx=cp.grism_cat.id[grism_idx], show=False)
    if lines:
        for line in lines:
            #### Mask lines for normalization
            flag = np.where(np.abs(grism.LAMBDA-line.wave) < 100)[0]
            mask[flag] = 0.
            if (line.flag == 'ok') & (line.type == 'em'):
                str = '%9.2e' %line.wave
                for mat in [3727, 5007, 6563.]:
                    str+='%6.3f' %(line.wave/mat-1)
                print str
                outstr+=str+'\n'
                
    use = np.where( (grism.LAMBDA > wave_limits[0]) & 
                    (grism.LAMBDA < wave_limits[1]) )
    
    use_mask = np.where(mask > 0)
    
    #### simple interpolation normalization.  should mask out lines
    finterp = interpolate.interp1d(lc,flux)
    einterp = interpolate.interp1d(lc,err)
    grism_to_phot = np.sum(finterp(grism.LAMBDA[use_mask])*grism.corr[use_mask]) / np.sum(finterp(grism.LAMBDA[use_mask])**2)
    
    # grism_to_phot = 10**(-0.4*( (23.86-2.5*np.log10(cp.phot_cat.h_colf[phot_idx])) - 
    #             np.float(cp.grism_cat.MAG_F1392W[grism_idx]) ))
        
    pyplot.semilogx(np.array([1,1]),np.array([1,1]))
    
    pyplot.plot(grism.LAMBDA[use], grism.corr[use]/grism_to_phot,
                    color='red')
    pyplot.errorbar(lc, flux, yerr=err, fmt='o',color='blue')
    pyplot.xlim(3000, 9.e4)
    ymax = np.max(flux)
    pyplot.ylim(-0.1*ymax,1.1*ymax)
    pyplot.text(2.e4,0.95*ymax, outstr, verticalalignment='top')
    return True
    
def check_z15():
    """
    """
    from matplotlib.backends.backend_pdf import PdfPages
    cp = threedhst.combine_phot
    
    hmag = 23.86-2.5*np.log10(cp.phot_cat.h_colf/cp.phot_cat.ks_colf* cp.phot_cat.ks_totf)
    
    jh = -2.5*np.log10(cp.phot_cat.j_colf/cp.phot_cat.h_colf)
    #jh = -2.5*np.log10(cp.phot_cat.z850_colf/cp.phot_cat.h_colf)
    
    zlimit = 1.1
    phot_z = threedhst.utils.ReadASCIICat('FIREWORKS/FIREWORKS_redshift.cat')
    q = np.where((phot_z.zph_best > zlimit) &
                  (cp.phot_cat.ra > min(cp.grism_cat.ra)) &
                  (cp.phot_cat.ra < max(cp.grism_cat.ra)) &
                  (cp.phot_cat.dec > min(cp.grism_cat.dec)) &
                  (cp.phot_cat.dec < max(cp.grism_cat.dec)) & 
                  (hmag < 25) & (jh > 0.8) )[0]
    
    ### sort on redshift
    q = q[phot_z.zph_best[q].argsort()]
    
    ### sort on H mag
    q = q[hmag[q].argsort()]
    
    pp = PdfPages('check_z15.pdf')
    
    for i,phot_idx in enumerate(q):

        if i > 200:
            continue
        
        fig = pyplot.figure(figsize=[7,4]) #,dpi=100)
        fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.08,
                            bottom=0.125,right=0.99,top=0.93)
        
        print '\n #%d \n' %i
        status = cp.phot_grism(phot_idx=phot_idx)
        if status:
            pp.savefig(fig)
        pyplot.close()
        
    pp.close()
        
    