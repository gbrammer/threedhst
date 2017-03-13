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

try:
    import astropy.io.fits as pyfits
except:
    import pyfits

import numpy as np
import numpy.random as nprand

import threedhst

## Global        
phot_cat = None
grism_cat = None
grism_SPC = None
fp = None

# G141
XMIN = 1.1e4
XMAX = 1.68e4

def make_filter_file(SPC, path='./EAZY/', NSTART=220):
    """
make_filter_file(SPC)
    
    Make an EAZY filter file for "filters" defined from grism spectrum
    wavelengths.
    """
    import scipy.integrate as integrate
    
    ## Number of points per "filter"
    NSPECFILT = 25
    ## Gaussian filter
    filtx = np.arange(NSPECFILT)/(NSPECFILT-1.)*8-4
    filty = np.exp(-0.5*filtx**2)/np.sqrt(2*np.pi)
    ## FWHM of line ~180 A
    filtx*=100
    ## Normalize 
    filty /= integrate.simps(filty,filtx)
    
    ## Make RES file
    spec = SPC.getSpec(SPC._ext_map[0])
    res = open(path+'grism.FILTER.RES','w')
    trans = open(path+'grism.translate','w')
    
    use = np.where((spec.LAMBDA >= XMIN) & (spec.LAMBDA < XMAX))[0]
    for i,l0 in enumerate(spec.LAMBDA[use]):
        res.write('%d fl%d\n' %(NSPECFILT, l0/10))
        for k in range(NSPECFILT):
            res.write('%4d %10.3e %10.3e\n' %(k+1, l0+filtx[k], filty[k]))
            
        trans.write('fl%d F%d\n' %(l0/10, NSTART+i))
        trans.write('el%d E%d\n' %(l0/10, NSTART+i))
        
    res.close()
    
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
    ROOT_DIRECT = 'ib6o23020'
    cp.grism_cat = threedhst.sex.mySexCat('DATA/'+ROOT_DIRECT+'_drz.cat')
    cp.grism_SPC = threedhst.plotting.SPCFile(ROOT_DIRECT+'_2_opt.SPC.fits',
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
    
def read_FIREWORKS_lines():
    """
fwhead, lines = read_FIREWORKS_lines()
    
    Read the header and content lines of the FIREWORKS catalog.
    """
    fp = open('FIREWORKS/FIREWORKS_phot.cat','r')
    lines = fp.readlines()
    fp.close()
    header = lines[0]
    
    skip = 0
    while lines[skip].startswith('#'):
        skip+=1
    return header, lines[skip:]
    
def check_FIREWORKS():
    """
    """
    from matplotlib.backends.backend_pdf import PdfPages
    cp = threedhst.combine_phot
    
    threedhst.plotting.defaultPlotParameters()
    pyplot.rcParams['text.usetex'] = False
    
    hmag = 23.86-2.5*np.log10(cp.phot_cat.h_colf/cp.phot_cat.ks_colf* cp.phot_cat.ks_totf)
    
    jh = -2.5*np.log10(cp.phot_cat.j_colf/cp.phot_cat.h_colf)
    jh = -2.5*np.log10(cp.phot_cat.z850_colf/cp.phot_cat.h_colf)
    
    zlimit = 0.0
    phot_z = threedhst.utils.ReadASCIICat('FIREWORKS/FIREWORKS_redshift.cat')
    q = np.where((phot_z.zsp > zlimit) &
                  (cp.phot_cat.ra >= min(cp.grism_cat.ra)) &
                  (cp.phot_cat.ra <= max(cp.grism_cat.ra)) &
                  (cp.phot_cat.dec >= min(cp.grism_cat.dec)) &
                  (cp.phot_cat.dec <= max(cp.grism_cat.dec)) & 
                  (hmag < 35) )[0]
    
    ### sort on redshift
    q = q[phot_z.zph_best[q].argsort()]
    
    ### sort on H mag
    q = q[hmag[q].argsort()]
    
    ### sort on ID
    q = q[phot_z.id[q].argsort()]
    
    pp = PdfPages('check_FIREWORKS.pdf')
    
    cp.fp = open('EAZY/grism_spec.cat','w')
    fwhead, lines = read_FIREWORKS_lines()
    fwhead = fwhead[:-1]+' z_spec z_phot'
    spec = grism_SPC.getSpec(grism_SPC._ext_map[0])
    use = np.where((spec.LAMBDA >= XMIN) & (spec.LAMBDA < XMAX))[0]
    for l0 in spec.LAMBDA[use]:
        fwhead+=' fl%d el%d' %(l0/10., l0/10.)
    cp.fp.write(fwhead+'\n')
    
    f = 1
    for i,phot_idx in enumerate(q):

        fig = pyplot.figure(figsize=[7*f,4*f]) #,dpi=100)
        fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.08,
                            bottom=0.145,right=0.99,top=0.93)
        
        print('\n #%d \n' %i)
        status = cp.phot_grism(phot_idx=phot_idx)
        if status:
            pp.savefig(fig)
            #pyplot.savefig('junk.svgz', dpi=60)
        pyplot.close()
    
    cp.fp.close() # EAZY/grism_spec.cat
    pp.close()
        
def check_bright():
    from matplotlib.backends.backend_pdf import PdfPages
    cp = threedhst.combine_phot
    
    threedhst.plotting.defaultPlotParameters()
    pyplot.rcParams['text.usetex'] = False
    
    hmag = np.cast[float](cp.grism_cat.MAG_F1392W)
    
    q = np.where(hmag < 23)[0]
    q = np.where((hmag > 23) & (hmag < 24))[0]
    q = np.where(hmag > 24)[0]
    q = q[hmag[q].argsort()]
    pp = PdfPages('check_bright.pdf')
    
    f = 1.
    for i,grism_idx in enumerate(q):
        
        fig = pyplot.figure(figsize=[7*f,4*f]) #,dpi=100)
        fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.08,
                            bottom=0.145,right=0.99,top=0.93)
        
        print('\n #%d \n' %i)
        status = cp.phot_grism(grism_idx=grism_idx)
        if status:
            pp.savefig(fig)
            #pyplot.savefig('junk.svgz', dpi=60)
        pyplot.close()
        
    pp.close()

###### General

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
            print('Neither `grism_idx` nor `phot_idx` are set.')
            return False
    
    if dr > 1.5:
        print('dr = %6.2f' %dr)
        return False
        
    phot_z = threedhst.utils.ReadASCIICat('FIREWORKS/FIREWORKS_redshift.cat')
    
    outstr =  ' %d %d %5.2f %6.3f %6.3f\n' %(cp.grism_cat.id[grism_idx],
           cp.phot_cat.id[phot_idx],dr,
           phot_z.zsp[phot_idx], phot_z.zph_best[phot_idx])
    
    print(outstr)
    
    ### Get broad-band SED, grism spectrum
    lc, flux, err = cp.get_sed(phot_idx)
    
    grism = cp.grism_SPC.getSpec(cp.grism_cat.id[grism_idx])
    if grism is False:
        return False
        
    grism.corr = grism.FLUX-grism.CONTAM

    #### 
    mask = grism.LAMBDA*0.
    keep = np.where((grism.LAMBDA > XMIN) & 
                    (grism.LAMBDA < XMAX))[0]
    mask[keep] = 1.
    
    #### Compare detected emission line test redshifts to BB zphot/zspec
    lines = threedhst.spec1d.findLines(cp.grism_SPC, 
                        idx=cp.grism_cat.id[grism_idx], show=False)
    if lines:
        for line in lines:
            #### Mask lines for normalization
            flag = np.where(np.abs(grism.LAMBDA-line.wave) < 150)[0]
            mask[flag] = 0.
            if (line.flag == 'ok') & (line.type == 'em'):
                str = '%14.2f' %(np.float(line.wave)/1.e4)
                str += r'$\mu\mathrm{m}$'
                for mat in [3727, 5007, 6563.]:
                    str+='%6.3f' %(line.wave/mat-1)
                print(str)
                outstr+=str+'\n'
                
    use = np.where( (grism.LAMBDA > XMIN) & 
                    (grism.LAMBDA < XMAX) )
    
    bad = np.where((grism.FLUX == 0) | (np.isnan(grism.FLUX)))[0]
    if len(bad) > 0:
        mask[bad] = 0.
        
    use_mask = np.where(mask > 0)
    bad_mask = np.where(mask == 0)
    ### nothing in spectrum to show
    if len(use_mask) == 0:
        return False
    
    #### simple interpolation normalization. 
    finterp = interpolate.interp1d(lc,flux)
    einterp = interpolate.interp1d(lc,err)
    grism_to_phot = np.sum(finterp(grism.LAMBDA[use_mask])*grism.corr[use_mask]) / np.sum(finterp(grism.LAMBDA[use_mask])**2)
    
    # grism_to_phot = 10**(-0.4*( (23.86-2.5*np.log10(cp.phot_cat.h_colf[phot_idx])) - 
    #             np.float(cp.grism_cat.MAG_F1392W[grism_idx]) ))
        
    pyplot.semilogx(np.array([1,1]),np.array([1,1]))
    
    pyplot.errorbar(lc, flux, yerr=err, fmt='o',color='0.95', 
                    ecolor='0.7', alpha=0.8, markersize=10)
                    
    pyplot.plot(grism.LAMBDA[use], grism.CONTAM[use]/grism_to_phot,
                    color='red', alpha=0.4, linewidth=2)
    
    pyplot.plot(grism.LAMBDA[use], grism.corr[use],
                    color='black', alpha=0.1, linewidth=2)
    
    pyplot.plot(grism.LAMBDA[use], grism.corr[use]/grism_to_phot,
                    color='blue', alpha=0.6, linewidth=2)
                    
    pyplot.xlim(3000, 9.e4)
    ymax = np.max(flux)
    pyplot.ylim(-0.1*ymax,1.1*ymax)
    pyplot.text(2.e4,0.95*ymax, outstr, verticalalignment='top')
    pyplot.ylabel(r'$f_\lambda$')
    pyplot.xlabel(r'$\lambda$')
    
    ### Make grism_spec.cat
    fwhead, lines = read_FIREWORKS_lines()
    fnu_fact = grism.LAMBDA**2/3.e18/1.e-29
    line = lines[phot_idx][:-1]+' %5.3f %5.3f' %(phot_z.zsp[phot_idx],
                      phot_z.zph_best[phot_idx])
    
    out_flux = grism.corr/grism_to_phot*fnu_fact
    out_err = grism.FERROR/grism_to_phot*fnu_fact
    out_flux[bad_mask] = -99
    out_err[bad_mask] = -99
    
    out_use = np.where((grism.LAMBDA >= XMIN) & (grism.LAMBDA < XMAX))[0]
    for i in out_use:
        line+=' %11.3e %11.3e' %(out_flux[i], out_err[i])
    
    cp.fp.write(line+'\n')
    
    return True
        