import os
import pyfits
import numpy as np
import glob

import matplotlib.pyplot as plt
from pyraf import iraf
from iraf import iraf

import threedhst
import threedhst.eazyPy as eazy
import threedhst.catIO as catIO

def go_all_SED():
    import threedhst.analysis
    threedhst.analysis.cosmos_SED_plots()
    threedhst.analysis.goods_SED_plots()

def aegis_SED_plots(grism_root='AEGIS-4-G141'):
    import threedhst.analysis #.aegis_SED_plots()

    os.chdir('/research/HST/GRISM/3DHST/AEGIS')
    
    ####### COSMOS
    MAIN_OUTPUT_FILE = 'aegis-n2.v4.6'
    OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/'
    
    cat = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'../'+MAIN_OUTPUT_FILE+'.cat')
    cat.kmag = 25.-2.5*np.log10(cat.field('ktot'))
    cat.star_flag = cat.field('star_flag')
    
    zout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
    fout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/../aegis-n2.bc03.v4.6.fout')
        
    print 'Read catalog...'
    grismCat = threedhst.sex.mySexCat('DATA/'+grism_root+'_drz.cat')
    for col in grismCat.column_names:
        if col.startswith('MAG_F'):
            grismCat.MAG = grismCat[col]
    #grismCat.MAG = grismCat.MAG_F1392W
    #grismCat.MAG = grismCat.MAG_F806W
    
    print 'Read SPC...'
    SPC = threedhst.plotting.SPCFile(grism_root+'_2_opt.SPC.fits',
                axe_drizzle_dir='DRIZZLE_'+threedhst.options['GRISM_NAME'])
    
    print 'Matched catalog'
    threedhst.analysis.match_to_NMBS(grism_root=grism_root, 
                  SPC = SPC, cat = cat,
                  grismCat = grismCat, zout = zout, fout = fout, 
                  OUTPUT = './HTML/SED/'+grism_root+'_match.cat')
    
    for id in grismCat.id:
        status = specphot(id=id, grism_root=grism_root, SPC = SPC,
            cat = cat,
            grismCat = grismCat, zout = zout, fout = fout, 
            OUT_PATH = './HTML/SED/', OUT_FILE_FORMAT=True, Verbose=False,
            MAIN_OUTPUT_FILE = MAIN_OUTPUT_FILE,
            OUTPUT_DIRECTORY = OUTPUT_DIRECTORY,
            CACHE_FILE = 'Same')
    
    
def cosmos_SED_plots(roots=['orient1','orient2']):
    
    os.chdir('/research/HST/GRISM/3DHST/COSMOS')
    import threedhst.analysis
    
    ####### COSMOS
    MAIN_OUTPUT_FILE = 'cosmos-1.v4.6'
    OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/'
    
    cat = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'../'+MAIN_OUTPUT_FILE+'.cat')
    cat.kmag = 25.-2.5*np.log10(cat.field('ktot'))
    cat.star_flag = cat.field('star_flag')
    
    zout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
    fout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/../cosmos-1.bc03.v4.6.fout')
    
    for grism_root in roots:
    
        print 'Read catalog...'
        grismCat = threedhst.sex.mySexCat('DATA/'+grism_root+'_drz.cat')
        grismCat.MAG = grismCat.MAG_F1392W
        #grismCat.MAG = grismCat.MAG_F806W
        
        print 'Read SPC...'
        SPC = threedhst.plotting.SPCFile(grism_root+'_2_opt.SPC.fits',
                    axe_drizzle_dir='DRIZZLE_'+threedhst.options['GRISM_NAME'])
        
        print 'Matched catalog'
        threedhst.analysis.match_to_NMBS(grism_root=grism_root, 
                      SPC = SPC, cat = cat,
                      grismCat = grismCat, zout = zout, fout = fout, 
                      OUTPUT = './HTML/SED/'+grism_root+'_match.cat')
        
        for id in grismCat.id:
            status = specphot(id=id, grism_root=grism_root, SPC = SPC,
                cat = cat,
                grismCat = grismCat, zout = zout, fout = fout, 
                OUT_PATH = './HTML/SED/', OUT_FILE_FORMAT=True, Verbose=False,
                MAIN_OUTPUT_FILE = MAIN_OUTPUT_FILE,
                OUTPUT_DIRECTORY = OUTPUT_DIRECTORY,
                CACHE_FILE = 'Same')
    
def goods_SED_plots(grism_root='GOODS-N-12-G141'):
    import threedhst.analysis
    
    os.chdir('/research/HST/GRISM/3DHST/GOODS-N')
    
    ######## GOODS-N
    MAIN_OUTPUT_FILE = 'photz'
    OUTPUT_DIRECTORY = '/research/HST/GRISM/3DHST/GOODS-N/MODS/EAZY/OUTPUT/'
    
    cat = catIO.Readfile('MODS/FAST/mods.cat')
    cat.addColumn('kmag',25-2.5*np.log10(cat.ks_ap))
    #cat.addColumn('kmag',cat.ks_ap)
    cat.addColumn('star_flag',cat.kmag*0.)
    
    fp = open('mods.reg','w')
    fp.write('fk5\n')
    for i in range(len(cat.ra)):
        fp.write('circle(%f,%f,1")\n' %(cat.ra[i], cat.dec[i]))
    fp.close()
    
    zout = catIO.Readfile('MODS/FAST/mods.zout')
    fout = catIO.ReadASCIICat('MODS/FAST/mods.bc03.fout')
    #fout = catIO.Readfile('MODS/FAST/mods.zout')
    #fout.addColumn('lmass',fout.z_peak*0.-1)
    
    print grism_root
    
    grismCat = threedhst.sex.mySexCat('DATA/'+grism_root+'_drz.cat')
    grismCat.MAG = grismCat.MAG_F1392W
    
    SPC = threedhst.plotting.SPCFile(grism_root+'_2_opt.SPC.fits',
                axe_drizzle_dir='DRIZZLE_G141')
    #
    threedhst.analysis.match_to_NMBS(grism_root=grism_root, 
                  SPC = SPC, cat = cat,
                  grismCat = grismCat, zout = zout, fout = fout, 
                  OUTPUT = './HTML/SED/'+grism_root+'_match.cat')
    
    for id in grismCat.id:
        specphot(id=id, grism_root=grism_root, SPC = SPC, cat = cat,
            grismCat = grismCat, zout = zout, fout = fout, 
            OUT_PATH = './HTML/SED/', OUT_FILE_FORMAT=True, Verbose=False,
            MAIN_OUTPUT_FILE = MAIN_OUTPUT_FILE,
            OUTPUT_DIRECTORY = OUTPUT_DIRECTORY,
            CACHE_FILE = 'Same')
    
def go_cdfs_SED():
    import threedhst.analysis
    # os.chdir("/research/HST/GRISM/3DHST/ERS")
    # threedhst.analysis.cdfs_SED_plots(grism_root='ib6o23020')
    # os.chdir("/research/HST/GRISM/3DHST/GOODS-S-SN")
    # threedhst.analysis.cdfs_SED_plots(grism_root='F125W_1026')
    # threedhst.analysis.cdfs_SED_plots(grism_root='F125W_1101')
    
    os.chdir("/research/HST/GRISM/3DHST/GOODS-S-SN")
    threedhst.analysis.cdfs_SED_plots(grism_root='ibhj06030')
    
    # ### G102
    # os.chdir("/research/HST/GRISM/3DHST/ERS")
    # threedhst.analysis.cdfs_SED_plots(grism_root='ib6o21020')
    
def cdfs_SED_plots(grism_root='ib6o23020'):
    import threedhst.analysis
    
    ####### FIREWORKS
    MAIN_OUTPUT_FILE = 'fireworks'
    OUTPUT_DIRECTORY = './FIREWORKS/'
    
    cat = threedhst.analysis.catIO.ReadASCIICat('FIREWORKS/FIREWORKS_phot.cat')
        
    cat.ra = cat.field('RA')
    cat.dec = cat.field('DEC')
    cat.id = cat.field('ID')
    cat.kmag = 23.86-2.5*np.log10(cat.field('Ks_totf'))
    cat.star_flag = cat.kmag*0
    
    zout = threedhst.analysis.catIO.ReadASCIICat('FIREWORKS/fireworks.zout')
    fout = threedhst.analysis.catIO.ReadASCIICat('FIREWORKS/fireworks.fout')
    ### Used 25 rather than 23.86 for FAST zeropoint.  Masses need to increase
    ### [[ fixed ]]
    #fout.lmass += np.log10(10**(0.4*(25-23.86)))
    #fout.lsfr  += np.log10(10**(0.4*(25-23.86)))
    
    print 'Read catalog...'
    grismCat = threedhst.sex.mySexCat('DATA/'+grism_root+'_drz.cat')
    for col in grismCat.column_names:
        if col.startswith('MAG_F'):
            grismCat.MAG = grismCat[col]
            
    # try:
    #     grismCat.MAG = grismCat.MAG_F1249W
    # except:
    #     grismCat.MAG = grismCat.MAG_F986W
        
    print 'Read SPC...'
    try:
        SPC = threedhst.plotting.SPCFile(grism_root+'_2_opt.SPC.fits',
                axe_drizzle_dir='DRIZZLE_G141')
    except:
        SPC = threedhst.plotting.SPCFile(grism_root+'_2_opt.SPC.fits',
                axe_drizzle_dir='DRIZZLE_G102')
        
    print 'Matched catalog'
    threedhst.analysis.match_to_NMBS(grism_root=grism_root, 
                  SPC = SPC, cat = cat,
                  grismCat = grismCat, zout = zout, fout = fout, 
                  OUTPUT = './HTML/SED/'+grism_root+'_match.cat')
    
    for id in grismCat.id:
        specphot(id=id, grism_root=grism_root, SPC = SPC, cat = cat,
            grismCat = grismCat, zout = zout, fout = fout, 
            OUT_PATH = './HTML/SED/', OUT_FILE_FORMAT=True, Verbose=False,
            MAIN_OUTPUT_FILE = MAIN_OUTPUT_FILE,
            OUTPUT_DIRECTORY = OUTPUT_DIRECTORY,
            CACHE_FILE = 'Same')

def go_UDS():
    import threedhst.analysis
    
    threedhst.analysis.uds_SED_plots(grism_root='MARSHALL-a')
    threedhst.analysis.uds_SED_plots(grism_root='MARSHALL-b')
    

def uds_SED_plots(grism_root='direct1'):
    import threedhst.analysis
    
    MAIN_OUTPUT_FILE = 'uds'
    OUTPUT_DIRECTORY = './UDSPHOT/'
    
    cat = threedhst.analysis.catIO.ReadASCIICat('UDSPHOT/uds.cat')
        
    cat.ra = cat.field('RA')
    cat.dec = cat.field('DEC')
    cat.id = cat.field('ID')
    cat.kmag = 25-2.5*np.log10(cat.field('K_totf'))
    cat.star_flag = cat.kmag*0
    
    zout = threedhst.analysis.catIO.ReadASCIICat('UDSPHOT/uds.zout')
    fout = threedhst.analysis.catIO.ReadASCIICat('UDSPHOT/uds.fout')
    
    print 'Read catalog...'
    grismCat = threedhst.sex.mySexCat('DATA/'+grism_root+'_drz.cat')
    for column in grismCat.column_names:
        if column.startswith('MAG_F'):
            grismCat.MAG = grismCat[column]
    
    print 'Read SPC...'
    SPC = threedhst.plotting.SPCFile(grism_root+'_2_opt.SPC.fits',
                axe_drizzle_dir='DRIZZLE_G141')
    
    print 'Matched catalog'
    threedhst.analysis.match_to_NMBS(grism_root=grism_root, 
                  SPC = SPC, cat = cat,
                  grismCat = grismCat, zout = zout, fout = fout, 
                  OUTPUT = './HTML/SED/'+grism_root+'_match.cat')
    
    for id in grismCat.id:
        threedhst.analysis.specphot(id=id, grism_root=grism_root, SPC = SPC, 
            cat = cat,
            grismCat = grismCat, zout = zout, fout = fout, 
            OUT_PATH = './HTML/SED/', OUT_FILE_FORMAT=True, Verbose=False,
            MAIN_OUTPUT_FILE = MAIN_OUTPUT_FILE,
            OUTPUT_DIRECTORY = OUTPUT_DIRECTORY,
            CACHE_FILE = 'Same')
    
def convolveWithThumb(id, lambdaz, temp_sed, SPC, oned=True, xint=None):
    """ 
    
    Convolve the best-fit eazy spectrum with the object shape
    
    """
    import scipy.convolve as conv
    
    thumb_file = 'HTML/images/'+SPC.filename.split('_2_opt')[0]+'_%05d' %(id) +'_thumb.fits.gz'
    try:
        thumb = pyfits.open(thumb_file)
    except:
        return -1
    
    size = thumb[0].data.shape
    
    twod_file = 'HTML/images/'+SPC.filename.split('_2_opt')[0]+'_%05d' %(id) +'_2D.fits.gz'
    twod = pyfits.open(twod_file)
    model1D = np.matrix(twod[5].data.sum(axis=1))
    model2D = np.array(np.dot(np.transpose(model1D),np.ones((1,size[0]))))
    
    if oned:
        profile = np.sum(thumb[0].data*model2D,axis=0)
        profile /= profile.sum()
    else:
        profile = thumb[0].data*model2D
        for i in len(model1D):
            profile[i,:] /= profile[i,:].sum()
    
    LSM = size[0]/2
    xmin = 3000
    xmax = 2.4e4
    
    q = np.where((lambdaz > (xmin-LSM)) & (lambdaz < (xmax+LSM)))[0]
    lambdaz = lambdaz[q]
    temp_sed = temp_sed[q]
    
    ### convolve with some factor of the object size
    #LSM = np.sqrt(LSM**2+(0.1*R[grism_idx]/0.128*46.5)**2)
    
    temp_sed_sm = temp_sed*1.
    for i in range(len(lambdaz)):
        temp_sed_sm[i] = np.trapz(1./np.sqrt(2*np.pi*50**2)*np.exp(-0.5*(lambdaz-lambdaz[i])**2/50**2)*temp_sed, lambdaz)
    
    #
    if xint is None:
        xint = np.arange(xmin-LSM*23.25,xmax+LSM*23.25,23.25)
    
    yint = np.interp(xint, lambdaz, temp_sed_sm)
    
    if oned:
        temp_sed_conv = conv(yint, profile, mode='same')
        # temp_sed_conv = yint*0.
        # for i in range(LSM,len(xint)-LSM):
        #     temp_sed_conv[i] = np.sum(yint[i-LSM:i+LSM]*profile)
    else:
        NX, NY = len(yint), len(model1D)
        temp_sed_conv = np.zeros((NY,NX))
        for i in range(LSM,len(xint)-LSM):
            temp_sed_conv[:,i] = conv(yint, profile[i,:].flatten()) #np.dot(yint[i-LSM:i+LSM],profile).reshape((NY,))
            
    
    thumb.close()
    
    return xint, temp_sed_conv
    
def specphot(id=69, grism_root='ibhm45030',
    MAIN_OUTPUT_FILE = 'cosmos-1.v4.6',
    OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/',
    CACHE_FILE = 'Same', Verbose=False,
    SPC = None, cat=None, grismCat = None,
    zout = None, fout = None, OUT_PATH='/tmp/', OUT_FILE_FORMAT=True,
    OUT_FILE='junk.png'):
    """
specphot(id)
    
    Get NMBS photometry and compare G141 spectrum of COSMOS 
    object #id
    """
    import scipy.interpolate as interpol
        
    ### 69, 54!
    
    xxx = """
    id=199
    grism_root='ibhm48'
    MAIN_OUTPUT_FILE = 'cosmos-1.v4.6'
    OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/'
    CACHE_FILE = 'Same'
    """
    
    #### Get G141 spectrum
    if Verbose:
        print 'Read SPC'
    
    if SPC is None:
        SPC = threedhst.plotting.SPCFile(grism_root+'_2_opt.SPC.fits',
                    axe_drizzle_dir='DRIZZLE_G141')
                    
    spec = SPC.getSpec(id)
    if spec is False:
        return False
        
    xmin = 3000
    xmax = 2.4e4
    
    lam = spec.field('LAMBDA')
    flux = spec.field('FLUX')
    ffix = flux-spec.field('CONTAM')
    ferr = spec.field('FERROR') #*0.06/0.128254
        
    if Verbose:
        print 'Read grism catalog'
        
    #### Read the grism catalog and get coords of desired object
    if grismCat is None:
        grismCat = threedhst.sex.mySexCat('DATA/'+grism_root+'_drz.cat')
    
    #### Source size
    R = np.sqrt(np.cast[float](grismCat.A_IMAGE)*np.cast[float](grismCat.B_IMAGE))
    grism_idx = np.where(grismCat.id == id)[0][0]
    
    Rmatch = R[grism_idx]*1.
    #print 'R=%f"' %(Rmatch)
    ra0 = grismCat.ra[grismCat.id == id][0]
    de0 = grismCat.dec[grismCat.id == id][0]
    
    #### Read EAZY outputs and get info for desired object
    if cat is None:
        cat = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'../'+MAIN_OUTPUT_FILE+'.cat')
    
    dr = np.sqrt((cat.ra-ra0)**2*np.cos(de0/360.*2*np.pi)**2+(cat.dec-de0)**2)*3600.
    
    
    nmbs_idx = np.where(dr == np.min(dr))[0][0]
    
    drMatch = dr[nmbs_idx]*1.
    #print 'dr = %7.2f\n' %(drMatch)
    #print drMatch, np.min(dr)
    
    if drMatch > 2:
        return False
        
    if Verbose:
        print 'Read zout'
    if zout is None:    
        zout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
    if fout is None:
        fout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/../cosmos-1.m05.v4.6.fout')
    
    if Verbose:
        print 'Read binaries'
        
    lambdaz, temp_sed, lci, obs_sed, fobs, efobs = \
        eazy.getEazySED(nmbs_idx, MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
                          OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
                          CACHE_FILE = CACHE_FILE)
         
    # zgrid, pz = eazy.getEazyPz(nmbs_idx, MAIN_OUTPUT_FILE=MAIN_OUTPUT_FILE, \
    #                                OUTPUT_DIRECTORY=OUTPUT_DIRECTORY, \
    #                                CACHE_FILE = CACHE_FILE)
            
    #### Smooth best-fit spectrum by gaussian with width LSM A
    # LSM = 50.
    # 
    # q = np.where((lambdaz > (xmin-LSM)) & (lambdaz < (xmax+LSM)))[0]
    # lambdaz = lambdaz[q]
    # temp_sed = temp_sed[q]
    # 
    # ### convolve with some factor of the object size
    # #LSM = np.sqrt(LSM**2+(0.1*R[grism_idx]/0.128*46.5)**2)
    # 
    # if Verbose:
    #     print 'Smooth spectrum to %6.2f A' %(LSM)
    #     
    # temp_sed_sm = temp_sed*1.
    # for i in range(len(lambdaz)):
    #     temp_sed_sm[i] = np.trapz(1./np.sqrt(2*np.pi*LSM**2)*np.exp(-0.5*(lambdaz-lambdaz[i])**2/LSM**2)*temp_sed, lambdaz)
    # 
    # #
    try:
        lambdaz, temp_sed_sm = threedhst.analysis.convolveWithThumb(id, lambdaz, temp_sed, SPC)
    except:
        temp_sed_sm = temp_sed*1.
        
    if Verbose: 
        print 'Normalize spectrum'
        
    #### Normalize G141 spectrum
    interp = interpol.interp1d(lambdaz, temp_sed_sm, kind='linear')

    q = np.where((lam > 1.08e4) & (lam < 1.68e4) & (flux > 0))[0]
    #### G102
    if lam.min() < 9000:
        q = np.where((lam > 0.8e4) & (lam < 1.13e4) & (flux > 0))[0]
    
    #### ACS G800L
    if lam.min() < 5000:
        q = np.where((lam > 0.55e4) & (lam < 1.0e4) & (flux > 0))[0]
        
    if len(q) == 0:
        return False

    yint = interp(lam[q])
        
    anorm = np.sum(yint*ffix[q])/np.sum(ffix[q]**2)
    
    if Verbose:
        print 'Start plot'
        
    #### Make the plot
    threedhst.plotting.defaultPlotParameters()
    
    xs=5.8
    ys = xs/4.8*3.2
    fig = plt.figure(figsize=[xs,ys]) #,dpi=100)
    fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.13*4.8/xs,
                        bottom=0.15*4.8/xs,right=1.-0.02*4.8/xs,top=1-0.10*4.8/xs)
    
    ymax = np.max((ffix[q])*anorm)
    
    if Verbose:
        print 'Make the plot'
        
    plt.plot(lambdaz, temp_sed_sm, color='red')
    # plt.errorbar(lam[q], ffix[q]*anorm, yerr=ferr[q]*anorm, color='blue', alpha=0.8)
    plt.plot(lam[q],ffix[q]*anorm, color='blue', alpha=0.2, linewidth=1)
    
    #### Show own extraction
    sp1d = threedhst.spec1d.extract1D(id, root=grism_root, path='./HTML', show=False, out2d=False)
    lam = sp1d['lambda']
    flux = sp1d['flux']
    ffix = sp1d['flux']-sp1d['contam']-sp1d['background']
    ferr = sp1d['error']
    anorm = np.sum(yint*ffix[q])/np.sum(ffix[q]**2)
    plt.plot(lam[q],ffix[q]*anorm, color='blue', alpha=0.6, linewidth=1)
    
    #### Show photometry + eazy template
    plt.errorbar(lci, fobs, yerr=efobs, color='orange', marker='o', markersize=10, linestyle='None', alpha=0.4)
    plt.plot(lambdaz, temp_sed_sm, color='red', alpha=0.4)

    plt.ylabel(r'$f_{\lambda}$')
    
    if plt.rcParams['text.usetex']:
        plt.xlabel(r'$\lambda$ [\AA]')
        plt.title('%s: \#%d, z=%4.1f' 
            %(SPC.filename.split('_2_opt')[0].replace('_','\_'),id,
              zout.z_peak[nmbs_idx]))
    else:
        plt.xlabel(r'$\lambda$ [$\AA$]')
        plt.title('%s: #%d, z=%4.1f' 
            %(SPC.filename.split('_2_opt')[0].replace('_','\_'),id,
              zout.z_peak[nmbs_idx]))
        
    #kmag = 25-2.5*np.log10(cat.ktot[nmbs_idx])
    kmag = cat.kmag[nmbs_idx]
    
    ##### Labels
    label = 'ID='+r'%s   K=%4.1f  $\log M$=%4.1f' %(np.int(cat.id[nmbs_idx]),
        kmag, fout.field('lmass')[nmbs_idx])
        
    plt.text(5e3,1.08*ymax, label, horizontalalignment='left',
      verticalalignment='bottom')
    
    
    label = 'R=%4.1f"' %(drMatch)
    if drMatch > 1.1:
        label_color = 'red'
    else:
        label_color = 'black'
    plt.text(2.2e4,1.08*ymax, label, horizontalalignment='right',
      color=label_color, verticalalignment='bottom')
    
    plt.xlim(xmin,xmax)
    plt.ylim(-0.1*ymax,1.2*ymax)
    
    if Verbose:
        print 'Save the plot'
    
    if OUT_FILE_FORMAT:
        out_file = '%s_%05d_SED.png' %(grism_root, id)
    else:
        out_file = OUT_FILE
    
    plt.savefig(OUT_PATH+'/'+out_file)
    
    noNewLine = '\x1b[1A\x1b[1M'
    print noNewLine+OUT_PATH+'/'+out_file
    
    if Verbose:
        print 'Close the plot window'
        
    plt.close()
    
def match_to_NMBS(grism_root='ibhm45',
        MAIN_OUTPUT_FILE = 'cosmos-1.v4.6',
        OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/',
        CACHE_FILE = 'Same', Verbose=False,
        SPC = None, cat=None, grismCat = None,
        zout = None, fout = None, OUTPUT='/tmp/match.cat'):

    
    import threedhst.analysis
    # 
    # MAIN_OUTPUT_FILE = 'cosmos-1.v4.6'
    # OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/'
    # CACHE_FILE = 'Same'
    
    if cat is None:
        cat = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'../'+MAIN_OUTPUT_FILE+'.cat')
    if zout is None:
        zout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
    if fout is None:
        fout = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/../cosmos-1.m05.v4.6.fout')
    
    rfUV = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.153-155.rf')
    
    z_peak = zout.field('z_peak')
    # kmag = 25-2.5*np.log10(cat.field('ktot'))
    kmag = cat.kmag
    
    uv = -2.5*np.log10(rfUV.field('L153')/rfUV.field('L155'))
    
    if grismCat is None:
        grismCat = threedhst.sex.mySexCat('DATA/'+grism_root+'_drz.cat')
    
    drs = grismCat.ra*1.
    ids = np.cast[int](grismCat.ra*1.)
    
    cosfact = np.cos(np.median(grismCat.dec)/360*2*np.pi)
    
    for i in range(len(drs)):
        dr = np.sqrt((cat.ra-grismCat.ra[i])**2*cosfact**2+
                     (cat.dec-grismCat.dec[i])**2)*3600
        mat = np.where(dr == np.min(dr))[0][0]
        drs[i] = dr[mat]
        ids[i] = mat
    
    fp = open(OUTPUT,'w')
    fp.write("# id_f140w  mag_f140w  id_phot  mag_Ktot  Rmatch  z_peak  UmV  logM star_flag\n# Rmatch = match distance in arcsec (should be < 1)\n# %s\n" %(MAIN_OUTPUT_FILE))
    for i in range(len(drs)):
        j = ids[i]
        line = "%5d %6.2f %8d %6.2f %5.1f %6.2f %6.2f %5.1f %d\n" %(grismCat.id[i],
                  np.float(grismCat.MAG[i]), np.int(cat.id[j]),
                  kmag[j], drs[i], 
                  z_peak[j], uv[j], fout.field('lmass')[j], 
                  np.int(cat.star_flag[j]))
        if grismCat.id[i] in SPC._ext_map:
            fp.write(line)
    fp.close()

def show_massive_galaxies(masslim=10.5, maglim=23.5, zrange=(0,5), 
    use_kmag=False):        
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS')
    
    matches = []
    
    # matches.extend(glob.glob('../COSMOS/HTML/SED/COSMOS*match.cat'))
    # matches.extend(glob.glob('../AEGIS/HTML/SED/AEGIS**match.cat'))
    matches.extend(glob.glob('../GOODS-N/HTML/SED/GOODS-N*match.cat'))
    
    # matches.extend(glob.glob('../GOODS-S/HTML/SED/*match.cat'))

    # matches.extend(glob.glob('../ERS/HTML_v1.1/SED/*match.cat'))
    # matches.extend(glob.glob('../GOODS-S-SN/HTML_v1.0/SED/*match.cat'))
    # matches.extend(glob.glob('../SN-GEORGE/HTML/SED/*match.cat'))
    # 
    # matches.extend(glob.glob('../SN-MARSHALL/HTML_v1.0/SED/*match.cat'))
        
    fplist = open('massive.coords','w')
    fplist.write('ID   ra   dec\n')
    fpreg = open('massive.reg', 'w')
    fpreg.write('fk5\n')
    
    fp = open('massive.html','w')
    fp.write("""
    <html>
    <head>
    <link rel="stylesheet" href="http://localhost/~gbrammer/COSMOS/scripts/style.css" type="text/css" id="" media="print, projection, screen" /> 
    
    <script type="text/javascript" src="http://localhost/~gbrammer/COSMOS/scripts/jquery-1.4.2.min.js"></script>
    
    <script type="text/javascript" src="http://localhost/~gbrammer/COSMOS/scripts/jquery.tablesorter.min.js"></script> 
    
    <script type="text/javascript" id="js">
    
    // Add ability to sort the table
    $(document).ready(function() {
        $.tablesorter.defaults.sortList = [[2,2]]; 
        $("table").tablesorter({
                // pass the headers argument and assing a object
                headers: {
                        // assign the secound column (we start counting zero)
                        4: {
                                sorter: false
                        },
                        5: {
                                sorter: false
                        },
                        6: {
                                sorter: false
                        },
                        7: {
                                sorter: false
                        },
                }
        });        
    });
    </script>
    
    </head>
    <body>
    <table id="myTable" cellspacing="1" class="tablesorter"> 
    <thead>
        <th> Grism id </th>
        <th> Mag_WFC3 </th>
        <th> z </th>
        <th> logM </th>
        <th> Thumb </th>
        <th> 2D </th>
        <th> 1D </th>
        <th> SED </th>
    </thead>
    <tbody>
    """)
    
    
    for match in matches:
        root_path = os.path.dirname(match).split('SED')[0]
        root = os.path.basename(match).split('_match.cat')[0]
        
        c = catIO.Readfile(match)
        xml_file = root_path+root+'.xml'
        xml = catIO.readMarkerXML(xml_file)
        #print c.keys()
        
        select_mag = c.mag_f140w
        if use_kmag:
            #print c.keys()
            select_mag = c.mag_ktot
            
        if ('ERS' in root_path) | ('GOODS-S' in root_path):
            c.logm -= 0*-0.4*(23.86-25)
            
        use = np.where((c.star_flag < 1) & (c.logm > masslim) & (c.rmatch < 1) & (select_mag < maglim) & (c.z_peak > zrange[0]) & (c.z_peak < zrange[1]))[0]
        
        use = use[np.argsort(c.logm[use])]
        
        for i in use:
            ID = c.id_f140w[i]
            
            file="%s_%05d" %(root, ID)
            
            fplist.write("%15s %14.6f %14.6f\n" %(file, xml[ID].ra, xml[ID].dec))
            fpreg.write('circle(%14.6f,%14.6f,1") # text={%s} color=magenta  \n' %(xml[ID].ra, xml[ID].dec, file))
            
            fp.write("""
        <tr>
            <td> %s<br> %13.6f %13.6f <br> %4.1f </td>
            <td> %5.1f </td>
            <td> %5.2f </td>
            <td> %5.2f </td>
            <td> <img src=%s/images/%s_thumb.png height=180px> </td>
            <td> <img src=%s/images/%s_2D.png height=180px> </td>
            <td> <img src=%s/images/%s_1D.png height=180px> </td>
            <td> <img src=%s/SED/%s_SED.png height=180px> </td>
        </tr>
        """ %(file, xml[ID].ra, xml[ID].dec, select_mag[i],
              c.mag_f140w[i], c.z_peak[i], c.logm[i],
              root_path, file,
              root_path, file,
              root_path, file,
              root_path, file))
            
    fp.write("</tbody></table></body></html>")
    fp.close()
    fplist.close()
    fpreg.close()
    
def testMatch():
    """
    Match the z=2.2 GOODS-N H-a emitters from Tadaki et al. 
    http://arxiv.org/abs/1012.4860v1
    """
    import threedhst.analysis
    
    ra_haem = [189.0141,188.9930,189.0906,189.2159,189.2582,189.2364,189.2204,189.2235,189.2152,189.2709,189.3259,189.3564]
    dec_haem = [62.1760,62.1988,62.2480,62.2513,62.2639,62.2788,62.2907,62.2901,62.3071,62.2908,62.2815,62.3194]
    
    #### LAE - HU et al.
    # ra_haem = [189.215240,189.216751,189.056106,189.254120,189.399750,189.324677,189.033004,189.456543,189.342285,189.045471,189.366013,189.320419]
    # dec_haem = [62.32683,62.36460,62.12994,62.35397,62.23944,62.29974,62.14394,62.22942,62.26277,62.17144,62.19613,62.23344]
    
    matches = []
    for i in range(len(ra_haem)):
        matches.extend(threedhst.analysis.matchObject(ra_haem[i], dec_haem[i]))
        print i, len(matches)
    
    fp = open('ha_emitters.html','w')
    fp.write("""
    <html>
    <head>
    <link rel="stylesheet" href="http://localhost/~gbrammer/COSMOS/scripts/style.css" type="text/css" id="" media="print, projection, screen" /> 
    
    <script type="text/javascript" src="http://localhost/~gbrammer/COSMOS/scripts/jquery-1.4.2.min.js"></script>
    
    <script type="text/javascript" src="http://localhost/~gbrammer/COSMOS/scripts/jquery.tablesorter.min.js"></script> 
    
    <script type="text/javascript" id="js">
    
    // Add ability to sort the table
    $(document).ready(function() {
        $.tablesorter.defaults.sortList = [[1,1]]; 
        $("table").tablesorter({
                // pass the headers argument and assing a object
                headers: {
                        // assign the secound column (we start counting zero)
                        3: {
                                sorter: false
                        },
                        4: {
                                sorter: false
                        },
                        5: {
                                sorter: false
                        },
                        6: {
                                sorter: false
                        },
                }
        });        
    });
    </script>
    
    </head>
    <body>
    <table id="myTable" cellspacing="1" class="tablesorter"> 
    <thead>
        <th> Grism id </th>
        <th> RA </th>
        <th> Dec </th>
        <th> Thumb </th>
        <th> 2D </th>
        <th> 1D </th>
        <th> SED </th>
    </thead>
    <tbody>
    """)
    
    for match in matches:
            fp.write("""
        <tr>
            <td> %s </td>
            <td> %f </td>
            <td> %f </td>
            <td> <img src=%s/images/%s_thumb.png height=180px> </td>
            <td> <img src=%s/images/%s_2D.png height=180px> </td>
            <td> <img src=%s/images/%s_1D.png height=180px> </td>
            <td> <img src=%s/SED/%s_SED.png height=180px> </td>
        </tr>
        """ %(match[1],
              match[3], match[4],
              match[0], match[1],
              match[0], match[1],
              match[0], match[1],
              match[0], match[1]))
    
    fp.write("</tbody></table></body></html>")
    fp.close()
        
    
def matchObject(ra_in, dec_in):
    
    os.chdir('/research/HST/GRISM/3DHST/ANALYSIS')
    
    catalogs = glob.glob('../GOODS-N/HTML_v1.0/*.cat')
    catalogs.extend(glob.glob('../COSMOS/HTML_v1.0/*.cat'))
    catalogs.extend(glob.glob('../ERS/HTML/*.cat'))
    catalogs.extend(glob.glob('../GOODS-S-SN/HTML/*.cat'))
    
    out = []
    for catalog in catalogs:
        
        cat = threedhst.sex.mySexCat(catalog)
        dir_root = os.path.dirname(catalog) #.split('/')[1]
        cat_root = os.path.basename(catalog).split('_drz')[0]
        
        dr = np.sqrt((cat.ra-ra_in)**2*np.cos(dec_in/360.*2*np.pi)**2+(cat.dec-dec_in)**2)*3600.
        idx = np.where(dr == np.min(dr))[0][0]
        drMatch = dr[idx]*1.
        if drMatch < 1:
            out.append([dir_root, '%s_%05d' %(cat_root, np.int(cat.NUMBER[idx])),
                        drMatch, ra_in, dec_in])
    
    return out

def test_line_histogram():
    
    import threedhst
    import matplotlib.pyplot as plt
    
    f = threedhst.spec1d.readLinesDat('HTML_v1.0/images/orient1_1D_lines.dat')
    q = (f.sn > 1.3) & (f.sn < 20) & (f.ndet < 5) & (f.sigma > 20) & (f.sigma < 100)
    
    q2 = q & (f.wave > 1.45e4) & (f.wave < 1.465e4)
    

def check_masses():
    import threedhst.catIO as catIO
    import matplotlib.pyplot as plt
    
    #### COSMOS
    OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/'
    MAIN_OUTPUT_FILE = 'cosmos-1.v4.6'
    cat_cosmos = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'../'+MAIN_OUTPUT_FILE+'.cat')
    cat_cosmos.kmag = 25.-2.5*np.log10(cat_cosmos.field('ktot'))
    zout_cosmos = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/'+MAIN_OUTPUT_FILE+'.zout')
    fout_cosmos = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'/../cosmos-1.bc03.v4.6.fout')
    
    #### GOODS-N
    MAIN_OUTPUT_FILE = 'photz'
    OUTPUT_DIRECTORY = '/research/HST/GRISM/3DHST/GOODS-N/MODS/EAZY/OUTPUT/'
    
    MODS = '/research/HST/GRISM/3DHST/GOODS-N/MODS/'
    cat_goodsn = catIO.Readfile(MODS+'/FAST/mods.cat')
    cat_goodsn.addColumn('kmag',25-2.5*np.log10(cat_goodsn.ks_ap))
    zout_goodsn = catIO.Readfile(MODS+'/FAST/mods.zout')
    fout_goodsn = catIO.ReadASCIICat(MODS+'/FAST/mods.bc03.fout')
    
    #### GOODS-S
    FIREWORKS = '/research/HST/GRISM/3DHST/GOODS-S-SN/FIREWORKS/'
    cat_fw = threedhst.analysis.catIO.ReadASCIICat(FIREWORKS+'/FIREWORKS_phot.cat')        
    cat_fw.kmag = 23.86-2.5*np.log10(cat_fw.field('Ks_totf'))
    zout_fw = threedhst.analysis.catIO.ReadASCIICat(FIREWORKS+'/fireworks.zout')
    fout_fw = threedhst.analysis.catIO.ReadASCIICat(FIREWORKS+'/fireworks.fout')
    
    #### Some K limit
    K_limit = 22.8
    q_cosmos = (cat_cosmos.kmag < K_limit) & (cat_cosmos.star_flag == 0) & (cat_cosmos.nchild <= 0)
    
    q_goodsn = cat_goodsn.kmag < K_limit
    q_fw = cat_fw.kmag < K_limit
    
    plt.plot(fout_cosmos.z[q_cosmos], fout_cosmos.lmass[q_cosmos], color='red', alpha=0.1, linestyle='None', marker='o')

    plt.plot(fout_goodsn.z[q_goodsn], fout_goodsn.lmass[q_goodsn], color='blue', alpha=0.3, linestyle='None', marker='o')

    plt.plot(fout_fw.z[q_fw], fout_fw.lmass[q_fw], color='green', alpha=0.3, linestyle='None', marker='o')
    
    plt.ylim(6,13)
    plt.xlim(0,4)

def go_field_area():
    import threedhst.analysis
    import threedhst.catIO as catIO
    
    #### COSMOS
    OUTPUT_DIRECTORY = '/Users/gbrammer/research/drg/PHOTZ/EAZY/NEWFIRM/v4.6/OUTPUT_KATE/'
    MAIN_OUTPUT_FILE = 'cosmos-1.v4.6'
    MAIN_OUTPUT_FILE = 'aegis-n2.v4.6'
    cat_cosmos = catIO.ReadASCIICat(OUTPUT_DIRECTORY+'../'+MAIN_OUTPUT_FILE+'.cat')
    
    q = cat_cosmos.wmin > 0.3
    area = threedhst.analysis.field_area(cat_cosmos.ra[q], cat_cosmos.dec[q])
    print area
    
def field_area(ra_in, dec_in):
    from shapely.geometry import MultiPoint, Point, MultiPolygon
    import numpy as np
    import matplotlib.pyplot as plt
    from descartes import PolygonPatch
    from shapely.ops import cascaded_union
    
    ra = ra_in*np.cos(np.mean(dec_in)/360.*2*np.pi)
    
    coords = []
    for i in xrange(len(ra)):
        coords.append((ra[i], dec_in[i]))
    
    env = MultiPoint(coords).envelope 
    env_area = env.area
    
    hull = MultiPoint(coords).convex_hull 
    hull_area = hull.area
    
    fig = plt.figure(1, figsize=(5,5), dpi=90)
    ax = fig.add_subplot(111)
    
    ax.plot(ra, dec_in, marker='.', linestyle='None', alpha=0.2, color='blue')
    x_e, y_e = env.boundary.xy
    ax.plot(x_e, y_e, color='red')

    x_h, y_h = hull.boundary.xy
    ax.plot(x_h, y_h, color='yellow')
        
    # plt.show()
    
    return hull_area

def LAE_morphologies():
    """
    Have a catalog of GOODS-S LAEs from Lucia via Harold.  Cut out thumbnails from my CANDELS mosaic.
    """
    import pywcs
    import matplotlib.pyplot as plt
    
    os.chdir("/research/HST/GRISM/3DHST/ANALYSIS/LAE")
    print 'Reading images...'
    f125w = pyfits.open('/research/HST/CANDELS/GOODS-S/PREP_FLT/goods-s_f125w_drz.fits.gz')
    f160w = pyfits.open('/research/HST/CANDELS/GOODS-S/PREP_FLT/goods-s_f160w_drz.fits')
    print 'Done.'
    
    lae = catIO.Readfile('GOODS-S/LuciaALLsample_34red_NB.txt')
    
    wcs = pywcs.WCS(header=f125w[1].header)
    
    xy = wcs.wcs_sky2pix(lae.ra, lae.dec, 1)
    NSUB = 16
    i = np.arange(lae.N)[lae.n20 == 135][0]
    
    plt.gray()
    
    for i in range(lae.N):
        if (xy[0][i] > NSUB) & (xy[0][i] < (wcs.naxis1-NSUB)) & (xy[1][i] > NSUB) & (xy[1][i] < (wcs.naxis2-NSUB)):
             
            print 'LAE_N20-%d.png' %lae.n20[i]
                       
            xi = xy[0][i]
            yi = xy[1][i]
            sub125 = f125w[1].data[yi-NSUB:yi+NSUB,xi-NSUB:xi+NSUB]
            sub160 = f160w[1].data[yi-NSUB:yi+NSUB,xi-NSUB:xi+NSUB]
            
            if np.std(sub125) == 0:
                continue
                
            fig = plt.figure(figsize=[6,3.4],dpi=100)
            fig.subplots_adjust(wspace=0.02,hspace=0.02,left=0.02,
                                bottom=0.02,right=0.98,top=0.86)
            
            scale = 0.128254
            size=NSUB
            asec_pix = 1./scale
            nasec = int(size/asec_pix)

            ax = fig.add_subplot(121)
            ax.imshow(sub125, vmin=-0.05,vmax=0.1, interpolation='nearest')
            
            ax.set_yticklabels([])
            xtick = ax.set_xticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size)
            ax.set_xticklabels([])
            ytick = ax.set_yticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size)
            
            plt.title('F125W')
            
            ax = fig.add_subplot(122)
            ax.imshow(sub160, vmin=-0.05,vmax=0.1, interpolation='nearest')

            ax.set_yticklabels([])
            xtick = ax.set_xticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size)
            ax.set_xticklabels([])
            ytick = ax.set_yticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size)
            
            plt.title('F160W')
            
            plt.text(-0.01, 1.12,'LAE N20-%d' %lae.n20[i],
                 horizontalalignment='center',
                 verticalalignment='center',
                 transform = ax.transAxes, fontsize=16)
            
            plt.savefig('LAE_N20-%d.png' %lae.n20[i])
            plt.close()
            