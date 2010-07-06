"""
3DHST.plotting

Utilities for plotting grism spectra.

"""
# $URL$
# $Rev$
# $Author$
# $Date$

import pyfits

import numpy as np
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.ticker as ticker
import pylab

def defaultPlotParameters():
    """
    defaultPlotParameters()
    """
    pyplot.rcParams['font.family'] = 'serif'
    pyplot.rcParams['font.serif'] = ['Times']
    pyplot.rcParams['ps.useafm'] = True
    pyplot.rcParams['patch.linewidth'] = 0.
    pyplot.rcParams['patch.edgecolor'] = 'black'
    pyplot.rcParams['text.usetex'] = True
    pyplot.rcParams['text.latex.preamble'] = ''

def plotThumb(object_number, mySexCat, in_image = None, size = 20, scale=0.128, 
              outfile='/tmp/thumb.png', close_window=False):
    """
    plotThumb(object_number, mySexCat, in_image = None, size = 20, scale=0.128,
              outfile='/tmp/thumb.png', close_window=False)
    """
    
    if in_image is not None:
        im = pyfits.open(mySexCat.filename.split('.cat')[0]+'.fits')
        in_image = im[1].data
    
    obj_str = str(object_number)
    idx = -1
    for i, number in enumerate(mySexCat.columns[mySexCat.searchcol('NUMBER')].entry):
        if number == obj_str:
            idx = i
            break
    
    if idx < 0:
        print 'Object \'%s\' not found in SExtractor catalog, %s.\n' %(obj_str, mySexCat.filename)
        return False
    
    ##### Figure out why X,Y are swapped in mySexCat.
    ##### Maybe the orientation of in_image is rotated w.r.t catalog?
    x0 = np.round(np.float(mySexCat.columns[mySexCat.searchcol('X_IMAGE')].entry[idx]))
    y0 = np.round(np.float(mySexCat.columns[mySexCat.searchcol('Y_IMAGE')].entry[idx]))
    sub = in_image[y0-size:y0+size, x0-size:x0+size]
        
    interp = 'nearest'
    asp = 'auto'
    sub_max = np.max(sub)
    if sub_max > 0:
        #### Minmax scaling
        vmin = -1.1*sub_max
        vmax = 0.1*sub_max
    else:
        vmin = -0.5*0.8
        vmax = 0.1*0.8
    
    flux = np.round(np.float(mySexCat.columns[mySexCat.searchcol('FLUX_AUTO')].entry[idx]))
    vmin = -0.03*flux
    vmax = 0.003*flux
    
    defaultPlotParameters()

    fig = pyplot.figure(figsize=[3,3],dpi=100)
    fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.02,bottom=0.02,right=0.98,top=0.98)
    ax = fig.add_subplot(111)
    ax.imshow(0-sub, interpolation=interp,aspect=asp,vmin=vmin,vmax=vmax)    
    
    asec_pix = 1./scale
    nasec = int(size/asec_pix)
    
    ax.set_yticklabels([])
    xtick = ax.set_xticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size)
    ax.set_xticklabels([])
    ytick = ax.set_yticks(np.arange(-1*nasec,nasec+1,1)*asec_pix+size)
    
    ### Save to PNG
    if outfile:
        fig.savefig(outfile,dpi=100,transparent=False)
    
    if close_window:
        status = pyplot.close()
    
    #pyplot.show()
    
def makeThumbs(SPCFile, mySexCat, path='./HTML/'):
    """
    makeThumbs(SPCFile, mySexCat,path='./HTML')
    
    Run plotThumbs for each object in SPCFile
    """
    import os
    noNewLine = '\x1b[1A\x1b[1M'
    im = pyfits.open(mySexCat.filename.split('.cat')[0]+'.fits')
    dat = im[1].data
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    for id in SPCFile._ext_map:
        idstr = '%04d' %id
        print noNewLine+'plotting.makeThumbs: %s_%s_thumb.png' %(root, idstr)
        plotThumb(id, mySexCat, in_image=dat, size= 20, scale =0.128,
                  outfile=path+root+'_'+idstr+'_thumb.png',close_window=True)

def plot2Dspec(SPCFile, object_number, outfile='/tmp/spec2D.png', close_window=False):
    """
    plot2Dspec(SPCFile, object_number, outfile='/tmp/spec2D.png', close_window=False)
    """
    import os
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    
    mef = pyfits.open('../DRIZZLE_G141/'+root+'_mef_ID'+str(object_number)+'.fits')
    
    head = mef['SCI'].header
    lmin = 10800
    lmax = 16800
    xmin = (lmin-head['CRVAL1'])/head['CDELT1']+head['CRPIX1']
    xmax = (lmax-head['CRVAL1'])/head['CDELT1']+head['CRPIX1']
    
    defaultPlotParameters()
    fig = pyplot.figure(figsize=[6,4],dpi=100)
    fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.06,bottom=0.13,right=0.98,top=0.98)
        
    interp = 'nearest'
    asp = 'auto'
    vmin = -0.6 
    vmax = 0.1
    
    pyplot.gray()
    
    mod_max = np.max(mef['MOD'].data)
    if mod_max > 0:
        vmin = -1*mod_max
        vmax = 0.1*mod_max
    else:
        vmin = -0.5*0.8
        vmax = 0.1*0.8
        
    ax = fig.add_subplot(311)
    ax.imshow(0-mef['SCI'].data, interpolation=interp,aspect=asp,vmin=vmin,vmax=vmax)    
    ax.set_xlim(xmin,xmax)
    ax.set_yticklabels([])
    ax.set_xticks((np.arange(np.ceil(lmin/1000.)*1000,np.ceil(lmax/1000.)*1000,1000)-head['CRVAL1'])/head['CDELT1']+head['CRPIX1'])
    ax.set_xticklabels([])
    pyplot.ylabel('G141')
    
    ax = fig.add_subplot(312)
    ax.imshow(0-mef['MOD'].data, interpolation=interp,aspect=asp,vmin=vmin,vmax=vmax)
    ax.set_xlim(xmin,xmax)
    ax.set_yticklabels([])
    ax.set_xticks((np.arange(np.ceil(lmin/1000.)*1000,np.ceil(lmax/1000.)*1000,1000)-head['CRVAL1'])/head['CDELT1']+head['CRPIX1'])
    ax.set_xticklabels([])
    pyplot.ylabel('Model')

    ax = fig.add_subplot(313)
    ax.imshow(0-mef['CON'].data, interpolation=interp,aspect=asp,vmin=vmin,vmax=vmax)
    ax.set_xlim(xmin,xmax)
    ax.set_yticklabels([])
    ax.set_xticks((np.arange(np.ceil(lmin/1000.)*1000,np.ceil(lmax/1000.)*1000,1000)-head['CRVAL1'])/head['CDELT1']+head['CRPIX1'])
    ax.set_xticklabels(np.arange(np.ceil(lmin/1000.)*1000,np.ceil(lmax/1000.)*1000,1000)/1.e4)
    pyplot.ylabel('Cont.')
    
    pyplot.xlabel(r'$\lambda$ [$\mu$m]')
    
    #ax.set_xticklabels([])
    #pyplot.show()
    
    ### Save to PNG
    if outfile:
        fig.savefig(outfile,dpi=100,transparent=False)
    
    if close_window:
        status = pyplot.close()

def makeSpec2dImages(SPCFile, path='./HTML/'):
    """
    makeSpec2dImages(SPCFile, path='./HTML')
    
    Run plotObject for each object in SPCFile
    """
    import os
    noNewLine = '\x1b[1A\x1b[1M'
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    for id in SPCFile._ext_map:
        idstr = '%04d' %id
        print noNewLine+'plotting.makeSpecImages: %s_%s_2D.png' %(root, idstr)
        plot2Dspec(SPCFile, id, outfile=path+root+'_'+idstr+'_2D.png',close_window=True)

def plot1Dspec(SPCFile, object_number, outfile='/tmp/spec.png', close_window=False):
    """
    plot1Dspec(SPCFile, object_number, outfile='/tmp/spec.png', close_window=False)
    """
    import os
    #import threedhst.plotting as pl
    #reload pl
    #self = pl.SPCFile('ib3721050_2_opt.SPC.fits')
    
    defaultPlotParameters()
    #object_number = 42
    spec = SPCFile.getSpec(object_number)
    flux = spec.field('FLUX')
    ferr = spec.field('FERROR')
    contam = spec.field('CONTAM')
    lam  = spec.field('LAMBDA')
    
    fig = pyplot.figure(figsize=[6,4],dpi=100)
    fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.08,bottom=0.13,right=0.98,top=0.93)
    
    ### Plot Flux and Contamination
    ax = fig.add_subplot(111)
    ax.plot(lam, flux-contam, linewidth=1.0, color='blue',alpha=1)
    ax.plot(lam, contam, linewidth=1.0, color='red',alpha=1)
    ax.plot(lam, flux,
               color='red',linewidth=1.0,alpha=0.2)
    ax.errorbar(lam, flux-contam,
               yerr= ferr,ecolor='blue',
               color='blue',fmt='o',alpha=0.5)
    
    ### Axes
    xmin = 10800
    xmax = 16800
    ax.set_xlim(xmin,xmax)
    sub = np.where((lam > xmin) & (lam < xmax))[0]
    ax.set_ylim(-0.05*np.max((flux-0*contam)[sub]),1.1*np.max((flux-0*contam)[sub]))
    
    ### Labels
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    pyplot.title('%s: \#%d' %(root,object_number))
    pyplot.xlabel(r'$\lambda$[\AA]')
    pyplot.ylabel(r'$\mathit{F}_{\lambda}$')
    
    ### Save to PNG
    if outfile:
        fig.savefig(outfile,dpi=100,transparent=False)
    
    if close_window:
        pyplot.close()

def makeSpec1dImages(SPCFile, path='./HTML/'):
    """
    makeSpec1dImages(SPCFile, path='./HTML')
    
    Run plotObject for each object in SPCFile
    """
    import os
    noNewLine = '\x1b[1A\x1b[1M'
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    for id in SPCFile._ext_map:
        idstr = '%04d' %id
        print noNewLine+'plotting.makeSpec1dImages: %s_%s_1D.png' %(root, idstr)
        plot1Dspec(SPCFile, id, outfile=path+root+'_'+idstr+'_1D.png',close_window=True)
    
def makeHTML(SPCFile, mySexCat, output='./HTML/index.html', title=None):
    """
    makeHTML(SPCFile, mySexCat, output='./HTML/index.html', title=None)
    """
    import os
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    if not title:
        title = root
    
    lines = ["""
    <html>
    <head>
    <script type="text/javascript" src="scripts/jquery-1.4.2.min.js"></script> 
    <script type="text/javascript" src="scripts/jquery.tablesorter.min.js"></script> 
    <link rel="stylesheet" href="scripts/style.css" type="text/css" id="" media="print, projection, screen" /> 
    <title>%s</title>
    <script type="text/javascript" id="js">$(document).ready(function() {
        $.tablesorter.defaults.sortList = [[0,0]]; 
    	$("table").tablesorter({
    		// pass the headers argument and assing a object
    		headers: {
    			// assign the secound column (we start counting zero)
    			4: {
    				// disable it by setting the property sorter to false
    				sorter: false
    			},
    			5: {
    				// disable it by setting the property sorter to false
    				sorter: false
    			},
    			6: {
    				// disable it by setting the property sorter to false
    				sorter: false
    			},
    		}
    	});
    });</script>    
    </head>
    <body>
    """ %title]
    
    lines.append("""
    <h2>%s</h2>
    """ %(title))
    
    lines.append("""
    <table id="myTable" cellspacing="1" class="tablesorter"> 
    <thead> 
    <tr> 
        <th width=50px>ID</th>
        <th width=100px>R.A.</th> 
        <th width=100px>Dec.</th> 
        <th width=50px>Mag</th> 
        <th width=180px>Thumb</th> 
        <th width=260px>1D Spec</th> 
        <th width=260px>2D Spec</th> 
    </tr> 
    </thead> 
    <tbody> 
    """)
        
    for id in SPCFile._ext_map:
        for idx,num in enumerate(mySexCat.columns[mySexCat.searchcol('NUMBER')].entry):
            if num == str(id):
                break
        
        ra  = mySexCat.columns[mySexCat.searchcol('X_WORLD')].entry[idx]
        dec = mySexCat.columns[mySexCat.searchcol('Y_WORLD')].entry[idx]
        mag = mySexCat.columns[mySexCat.searchcol('MAG_F1392W')].entry[idx]
        
        img = '%s_%04d' %(root,id)
        lines.append("""
        <tr> 
            <td>%d</td>
            <td>%13.6f</td> 
            <td>%13.6f</td> 
            <td>%6.2f</td> 
            <td><a href='images/%s_thumb.png'><img src='images/%s_thumb.png' width=133px></a></td> 
            <td><a href='images/%s_1D.png'><img src='images/%s_1D.png' width=200px></a></td> 
            <td><a href='images/%s_2D.png'><img src='images/%s_2D.png' width=200px></a></td> 
        </tr> 
        """ %(id,np.float(ra),np.float(dec),np.float(mag),img,img,img,img,img,img))
    
    lines.append("""
    </tbody>
    </body>
    </html>""")
    
    if not output:
        output='HTML/index.html'
    
    fp = open(output,'w')
    fp.writelines(lines)
    fp.close()
    
    

class SPCFile(object):
    """
    SPCFile
    
    Class for reading and plotting spectra from an aXe
    
    """
    def _getPath(self):
        """
        _getPath()
        
        Figure out if we're in the root directory or in DATA
        """
        import os
        if os.getcwd().split('/')[-1] == 'DATA':
            self.path = '../'+self.axe_drizzle_dir+'/'
        elif os.getcwd().split('/')[-1] == self.axe_drizzle_dir:
            self.path = './'
        else:
            self.path = './'+self.axe_drizzle_dir+'/'
    
    def _mapFitsExtensions(self):
        """
        _mapFitsExtensions()
        
        Figure out which object corresponds to which extension in the SPC.fits file
        """
        for i in range(self.N_ext):
            self._ext_map[i] = np.int(self.fits[i+1].header['EXTNAME'].split('BEAM_')[1][:-1])
            
    def __init__(self, filename='ib3721050_2_opt.SPC.fits', axe_drizzle_dir='DRIZZLE_G141'):
        """
        __init__(filename='ib3721050_2_opt.SPC.fits', axe_drizzle_dir='DRIZZLE_G141')
        """
        self.filename = filename
        self.axe_drizzle_dir = axe_drizzle_dir
        self._getPath()
        self.fits = pyfits.open(self.path+filename)
        
        self.N_ext = len(self.fits)-1
        self._ext_map = np.arange(self.N_ext)*0
        self._mapFitsExtensions()
        
    def getSpec(self, object_number):
        """
        getSpec(self, object_number)
        """
        if object_number in self._ext_map:
            idx = np.where(self._ext_map == object_number)[0][0]
            return self.fits[idx+1].data
        else:
            print "Object #%d not found in %s." (object_number, self.filename)
            return False
            