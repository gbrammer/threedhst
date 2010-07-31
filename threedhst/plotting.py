"""
3DHST.plotting

Utilities for plotting grism spectra.

"""

__version__ = "$Rev$"
# $URL$
# $Author$
# $Date$

import pyfits

import numpy as np
import matplotlib
import matplotlib.pyplot as pyplot
import matplotlib.ticker as ticker
import pylab

import threedhst

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
    for i, number in enumerate(
      mySexCat.columns[mySexCat.searchcol('NUMBER')].entry):
        if number == obj_str:
            idx = i
            break
    
    if idx < 0:
        print 'Object \'%s\' not found in SExtractor catalog, %s.\n' %(obj_str,
                             mySexCat.filename)
        return False
    
    ##### Figure out why X,Y are swapped in mySexCat.
    ##### Maybe the orientation of in_image is rotated w.r.t catalog?
    x0 = \
 np.round(np.float(mySexCat.columns[mySexCat.searchcol('X_IMAGE')].entry[idx]))
    y0 = \
 np.round(np.float(mySexCat.columns[mySexCat.searchcol('Y_IMAGE')].entry[idx]))
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
    
    flux = \
np.round(np.float(mySexCat.columns[mySexCat.searchcol('FLUX_AUTO')].entry[idx]))
    vmin = -0.03*flux
    vmax = 0.003*flux
    
    defaultPlotParameters()
    
    pyplot.gray()
    
    fig = pyplot.figure(figsize=[3,3],dpi=100)
    fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.02,
                        bottom=0.02,right=0.98,top=0.98)
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
    ids = SPCFile._ext_map+0
    ids.sort()
    for id in ids:
        idstr = '%04d' %id
        print noNewLine+'plotting.makeThumbs: %s_%s_thumb.png' %(root, idstr)
        plotThumb(id, mySexCat, in_image=dat, size= 20, scale =0.128,
                  outfile=path+'/'+root+'_'+idstr+'_thumb.png',
                  close_window=True)

def plot2Dspec(SPCFile, object_number, outfile='/tmp/spec2D.png',
               close_window=False):
    """
    plot2Dspec(SPCFile, object_number, outfile='/tmp/spec2D.png', 
               close_window=False)
    """
    import os
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    
    mef = pyfits.open('../'+threedhst.options['DRIZZLE_PATH']+
                      '/'+root+'_mef_ID'+str(object_number)+'.fits')
    
    head = mef['SCI'].header
    lmin = 10800
    lmax = 16800
    xmin = (lmin-head['CRVAL1'])/head['CDELT1']+head['CRPIX1']
    xmax = (lmax-head['CRVAL1'])/head['CDELT1']+head['CRPIX1']
    
    defaultPlotParameters()
    fig = pyplot.figure(figsize=[6,4],dpi=100)
    fig.subplots_adjust(wspace=0.2,hspace=0.02,left=0.06,
                        bottom=0.12,right=0.99,top=0.99)
    # fig = pyplot.figure(figsize=[4.2,2.8]) #,dpi=100)
    # fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.07,
    #                     bottom=0.11,right=0.99,top=0.93)
        
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
    ax.imshow(0-mef['SCI'].data, interpolation=interp,aspect=asp,
              vmin=vmin,vmax=vmax)    
    ax.set_xlim(xmin,xmax)
    ax.set_yticklabels([])
    ax.set_xticks((np.arange(np.ceil(lmin/1000.)*1000,
       np.ceil(lmax/1000.)*1000,1000)-head['CRVAL1'])/head['CDELT1']
       +head['CRPIX1'])
    ax.set_xticklabels([])
    pyplot.ylabel(threedhst.options['GRISM_NAME'])
    
    ax = fig.add_subplot(312)
    ax.imshow(0-mef['MOD'].data, interpolation=interp, aspect=asp,
              vmin=vmin,vmax=vmax)
    ax.set_xlim(xmin,xmax)
    ax.set_yticklabels([])
    ax.set_xticks((np.arange(np.ceil(lmin/1000.)*1000,
        np.ceil(lmax/1000.)*1000,1000)-head['CRVAL1'])/head['CDELT1']
        +head['CRPIX1'])
    ax.set_xticklabels([])
    pyplot.ylabel('Model')

    ax = fig.add_subplot(313)
    ax.imshow(0-mef['CON'].data, interpolation=interp, aspect=asp,
              vmin=vmin,vmax=vmax)
    ax.set_xlim(xmin,xmax)
    ax.set_yticklabels([])
    ax.set_xticks((np.arange(np.ceil(lmin/1000.)*1000,
        np.ceil(lmax/1000.)*1000,1000)-head['CRVAL1'])/head['CDELT1']
        +head['CRPIX1'])
    ax.set_xticklabels(np.arange(np.ceil(lmin/1000.)*1000,
        np.ceil(lmax/1000.)*1000,1000)/1.e4)
    pyplot.ylabel('Contam.')
    
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
    ids = SPCFile._ext_map+0
    ids.sort()
    for id in ids:
        idstr = '%04d' %id
        print noNewLine+'plotting.makeSpecImages: %s_%s_2D.png' %(root, idstr)
        plot2Dspec(SPCFile, id, outfile=path+'/'+root+'_'+idstr+'_2D.png',
                   close_window=True)

def plot1Dspec(SPCFile, object_number, outfile='/tmp/spec.png',
               close_window=False, show_test_lines=False):
    """
    plot1Dspec(SPCFile, object_number, outfile='/tmp/spec.png', 
               close_window=False, show_test_lines=False)
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
        
    ### Initialize plot
    # fig = pyplot.figure(figsize=[6,4],dpi=100)
    # fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.07,
    #                     bottom=0.11,right=0.99,top=0.93)
    fig = pyplot.figure(figsize=[4.8,3.2]) #,dpi=100)
    fig.subplots_adjust(wspace=0.2,hspace=0.2,left=0.08,
                        bottom=0.125,right=0.99,top=0.93)

    ### plot window
    ax = fig.add_subplot(111)
        
    ### Plot Flux and Contamination
    ax.plot(lam, flux-contam, linewidth=1.0, color='blue',alpha=1)
    ax.plot(lam, contam, linewidth=1.0, color='red',alpha=1)
    ax.plot(lam, flux,
               color='red',linewidth=1.0,alpha=0.2)
    ax.errorbar(lam, flux-contam,
               yerr= ferr,ecolor='blue',
               color='blue',fmt='.',alpha=0.5)
    
    xmin = 10800
    xmax = 16800
    sub = np.where((lam > xmin) & (lam < xmax))[0]
    ymax = np.max((flux-0*contam)[sub])
    
    #### Search for lines
    lines = threedhst.spec1d.findLines(SPCFile, idx=object_number)
    if lines:
        sdss_lines = threedhst.spec1d.readLinelist()
        #sdss_lines.gal_weight = sdss_lines.gal_weight+1
        sdss_lines.gal_weight /= np.max(sdss_lines.gal_weight)
        sdss_lines.qso_weight /= np.max(sdss_lines.qso_weight)
        colors = ['green','orange','blue']
        weight = sdss_lines.gal_weight
        #weight = sdss_lines.qso_weight
        sdss_use = np.where(np.array(weight) > -10)[0]
        
        iline=0
        for line in lines:
            if (line.flag == 'ok') & (line.type=='em'):
                #print iline
                ax.plot(line.wave*np.array([1,1]),np.array([-1,1]),
                        color='black',linewidth=2,alpha=0.2)
                ax.plot(line.wave*np.array([1,1]),np.array([-1,1]),'--',
                        color='orange',linewidth=2,alpha=0.7)
                ### show assuming line is OII, OIII, Ha
                compare_lines = np.array([3727., 5007, 6563.])
                if show_test_lines: 
                  for iz, z_test in enumerate(line.wave/compare_lines-1):
                    print 'z: %5.2f' %z_test
                    ygauss = lam*0.
                    gsigma = 60 ## Ang
                    for l in sdss_use:
                        ygauss += np.exp(
                            -1.*( (lam - sdss_lines.wave[l]*(1+z_test) ) 
                            /gsigma)**2)*ymax*0.1*weight[l]
                    
                    scl = 0.1
                    ax.plot(lam,ygauss+scl*ymax*iz,alpha=0.7,color='white', 
                            linewidth=4)
                    ax.plot(lam,ygauss+scl*ymax*iz,alpha=1,color=colors[iline])
                
                iline += 1
                    
            if (line.flag == 'contam') & (line.type=='em'):
                ax.plot(line.wave*np.array([1,1]),np.array([-1,1]),'--',
                        color='orange',linewidth=2,alpha=0.1)
            
            if (line.flag == 'ok') & (line.type=='abs'):
                ax.plot(line.wave*np.array([1,1]),np.array([-1,1]),
                        color='black',linewidth=2,alpha=0.2)
                ax.plot(line.wave*np.array([1,1]),np.array([-1,1]),'--',
                        color='green',linewidth=2,alpha=0.7)

            if (line.flag == 'artifact') & (line.type=='abs'):
                ax.plot(line.wave*np.array([1,1]),np.array([-1,1]),'--',
                        color='green',linewidth=2,alpha=0.1)
                
    ### Show emission line wavelengths
    # zb = np.linspace(0.5,2.5,5)
    # zi = np.linspace(0.5,2.5,20)
    # yz = np.linspace(0,ymax*1.1,5)
    # lines = [3727,4862,5007,6563]
    # for line in lines:
    #     ax.plot(line*(1+zb),yz,color='black',alpha=0.1)
    #     ax.scatter(line*(1+zb),yz,marker='o',color='black',alpha=0.1)
    
    ### Axes
    #pyplot.semilogx(subsx=[11000,12500,15000])
    ax.set_xlim(xmin,xmax)
    ax.set_ylim(-0.05*ymax,1.1*ymax)
    
    ### Labels
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    pyplot.title('%s: \#%d' %(root,object_number))
    pyplot.xlabel(r'$\lambda$[\AA]')
    pyplot.ylabel(r'$\mathit{F}_{\lambda}$')

    
    ### Save to PNG
    if outfile:
        fig.savefig(outfile,dpi=75,transparent=False)
    
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
    ids = SPCFile._ext_map+0
    ids.sort()
    for id in ids:
        idstr = '%04d' %id
        print noNewLine+'plotting.makeSpec1dImages: %s_%s_1D.png' %(root, idstr)
        plot1Dspec(SPCFile, id, outfile=path+'/'+root+'_'+idstr+'_1D.png',
                   close_window=True)
    
def makeHTML(SPCFile, mySexCat, mapParams,
             output='./HTML/index.html', title=None):
    """
    makeHTML(SPCFile, mySexCat, mapParams,
             output='./HTML/index.html', title=None)
    """
    import os
    root = os.path.basename(SPCFile.filename).split('_2')[0]
    if not title:
        title = root
    
    #### Header
    lines = ["""
    <html>
    <head>
    
    <link rel="stylesheet" href="scripts/style.css" type="text/css" id="" media="print, projection, screen" /> 

    <script type="text/javascript" src="scripts/jquery-1.4.2.min.js"></script> 

    <script type="text/javascript" src="scripts/jquery.tablesorter.min.js"></script> 

    <!--  www.astro.yale.edu/brammer/ ABQIAAAAzSrfHr_4F2D2YfSuYQD2ZBSRvJBsXGI3t4UH99Pp8ZgdgIZDpRQ9Pmiw1tbadMh-1wRrE07VYIPVOg -->

    <script src="http://maps.google.com/maps?file=api&amp;v=3&amp;key=ABQIAAAA1XbMiDxx_BTCY2_FkPh06RR20YmIEbERyaW5EQEiVNF0mpNGfBSRb_rzgcy5bqzSaTV8cyi2Bgsx3g" type="text/javascript"></script> 
    """]
    
    #### Script for sorting the table
    lines.append("""
    
    <script type="text/javascript" id="js">
    
    // Add ability to sort the table
    $(document).ready(function() {
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
    	
    	switch_layout();
    	
    });
    
    var layout = -1; // Start with horizontal
    function switch_layout() {
        if (layout == 0) {
            layout = 1;
            vertical_layout();
            $("#switchbox").text("||");
            $("#switchbox").css("cursor","s-resize");  
            map.checkResize();                    
        } else {
            layout = 0;
            horizontal_layout();
            $("#switchbox").text("=");
            $("#switchbox").css("cursor","e-resize");
            map.checkResize();          
        }
    }
    
    function vertical_layout() {
    
    	$("#title").css("width",1087);
    	
    	$("#content").css("height",$(window).height()-60);
    	$("#content").css("width","840px");
    	$("#content").css("top","55px");
    	$("#content").css("left","301px");
        
    	$("#map").css("height",$(window).height()-60);	
    	$("#map").css("width","300px");	
        
        $("#coords").css("left",
    	    parseInt($("#map").css("width"))/2.-
    	    parseInt($("#coords").css("width"))/2.);
    	    
    	c = $("#centerbox");
        c.css("left",parseFloat($("#map").css("left"))+
                     parseFloat($("#map").css("width"))/2.-
                     parseFloat(c.css("width"))/2.-0);
                     
        c.css("top",parseFloat($("#map").css("top"))+
                     parseFloat($("#map").css("height"))/2.-
                     parseFloat(c.css("height"))/2.-0);
    }
    
    function horizontal_layout() {

    	$("#title").css("width",$(window).width()-12-
    	    parseFloat($("#title").css("padding-left"))+
    	    parseFloat($("#title").css("padding-right")));

    	$("#content").css("height",170);
    	$("#content").css("width",$("#title").width()+
    	    parseFloat($("#title").css("padding-left"))+
    	    parseFloat($("#title").css("padding-right")));
    	//alert(parseFloat($("#title").css("padding-left"))+20);
    	$("#content").css("top",$(window).height()-170);
        $("#content").css("left",$("#map").css("left"));
        
    	$("#map").css("height",$(window).height()-60-170);
    	$("#map").css("width",$("#title").width()+
    	    parseFloat($("#title").css("padding-left"))+
    	    parseFloat($("#title").css("padding-right")));    
    	
    	$("#coords").css("left",
    	    parseInt($("#map").css("width"))/2.-
    	    parseInt($("#coords").css("width"))/2.);    
    	
    	c = $("#centerbox");
        c.css("left",parseFloat($("#map").css("left"))+
                     parseFloat($("#map").css("width"))/2.-
                     parseFloat(c.css("width"))/2.-0);
                     
        c.css("top",parseFloat($("#map").css("top"))+
                     parseFloat($("#map").css("height"))/2.-
                     parseFloat(c.css("height"))/2.-0);
    }
    
    function show_centerbox() {
        c = $("#centerbox");
        document.getElementById("centerbox").style.display = "block";
        c.animate({ opacity: 1.0}, 300, function() { });        
        c.animate({ opacity: 0.0}, 300, function() {
            document.getElementById("centerbox").style.display = "none";
        });
    }
    
    </script>   
    """)
    
    #### Script for the Google map
    llSW = mapParams['LLSW']
    llNE = mapParams['LLNE']
    center = mapParams['LLCENTER']
    lng_offset = mapParams['LNG_OFFSET']
    
    lines.append("""
            <script type="text/javascript"> 
    
    var map = 0; // Global
    var centerLng = %f;
    var offset = %f;
    var zoomLevel = %f;
    var root = '%s';
    
    function initialize() {        
        if (GBrowserIsCompatible()) {
            map = new GMap2(document.getElementById("map"));
            map.addControl(new GScaleControl());
            // map.addControl(new GSmallMapControl());
            // map.addControl(new GMapTypeControl());
            var copyright = new GCopyright(1,
                 new GLatLngBounds(new GLatLng(%f,%f),
                                   new GLatLng(%f,%f)),
                                   zoomLevel, "3D-HST");
            var copyrightCollection = new GCopyrightCollection('Map Data:');
            copyrightCollection.addCopyright(copyright);

            // Direct image tiles
            CustomGetDirectTileUrl=function(a,b){
                return "tiles/"+root+"_d_"+a.x+"_"+a.y+"_"+b+".jpg"
            }
            var tilelayersDirect = [new GTileLayer(copyrightCollection,
                                          zoomLevel,zoomLevel)];
            tilelayersDirect[0].getTileUrl = CustomGetDirectTileUrl;
            var custommapDirect = new GMapType(tilelayersDirect, 
                   new GMercatorProjection(zoomLevel+1), "Direct");
            map.addMapType(custommapDirect);

            // Grism image tiles
            CustomGetGrismTileUrl=function(a,b){
                return "tiles/"+root+"_g_"+a.x+"_"+a.y+"_"+b+".jpg"
            }
            var tilelayersGrism = [new GTileLayer(copyrightCollection,
                                          zoomLevel,zoomLevel)];
            tilelayersGrism[0].getTileUrl = CustomGetGrismTileUrl;
            var custommapGrism = new GMapType(tilelayersGrism, 
                   new GMercatorProjection(zoomLevel+1), "Grism");
            map.addMapType(custommapGrism);
            
            // Model image tiles
            CustomGetModelTileUrl=function(a,b){
                return "tiles/"+root+"_m_"+a.x+"_"+a.y+"_"+b+".jpg"
            }
            var tilelayersModel = [new GTileLayer(copyrightCollection,
                                          zoomLevel,zoomLevel)];
            tilelayersModel[0].getTileUrl = CustomGetModelTileUrl;
            var custommapModel = new GMapType(tilelayersModel, 
                   new GMercatorProjection(zoomLevel+1), "Model");
            map.addMapType(custommapModel);
            
            // Can't remove all three for some reason
            map.removeMapType(G_NORMAL_MAP);
            map.removeMapType(G_HYBRID_MAP);
            map.removeMapType(G_SATELLITE_MAP);
            
            map.addControl(new GMapTypeControl());
            
            // Set map center
            map.setCenter(new GLatLng(%f, offset), zoomLevel,
             custommapDirect);
            
            latLng2raDec();
            
            GEvent.addListener(map, "moveend", function() {
                latLng2raDec();
                show_centerbox();
            });
            
            plotXmlObjects();
        }
    }
    
    function latLng2raDec() {
        var mapcenter = map.getCenter();
        var dec = mapcenter.lat();                       
        var dsign = "+";
        if (dec < 0) {dsign = "-";}
        dec = Math.abs(dec);
        var ded = parseInt(dec);
        var dem = parseInt((dec-ded)*60);
        var des = parseInt(((dec-ded)*60-dem)*60);
        var dess = parseInt((((dec-ded)*60-dem)*60-des)*10);
        if (ded < 10) {ded = "0"+ded;} 
        if (dem < 10) {dem = "0"+dem;} 
        if (des < 10) {des = "0"+des;} 
        var decstr = dsign+ded+":"+dem+":"+des+"."+dess;
        document.getElementById("decInput").value = decstr;
        
        var ra = (360-mapcenter.lng()+offset-centerLng)/360.*24;
        var rah = parseInt(ra);
        var ram = parseInt((ra-rah)*60);
        var ras = parseInt(((ra-rah)*60-ram)*60);
        var rass = parseInt((((ra-rah)*60-ram)*60-ras)*100);
        if (rah < 10) {rah = "0"+rah;} 
        if (ram < 10) {ram = "0"+ram;} 
        if (ras < 10) {ras = "0"+ras;} 
        if (rass < 10) {rass = "0"+rass;} 
        var rastr = rah+":"+ram+":"+ras+"."+rass;
        document.getElementById("raInput").value = rastr;        
    }
    
    function centerOnInput() {
        var rastr = document.getElementById("raInput").value;
        var rsplit = rastr.split(":");
        
        var ra = (parseInt(rsplit[0])+
              parseInt(rsplit[1])/60.+
              parseFloat(rsplit[2])/3600.)/24.*360;
    
        var decstr = document.getElementById("decInput").value;
        var dsplit = decstr.split(":");
        var dec = Math.abs(parseInt(dsplit[0])+
              Math.abs(parseInt(dsplit[1])/60.)+
              Math.abs(parseFloat(dsplit[2])/3600.));

        /// Don't know why, but need to do this twice
        dec = Math.abs(parseInt(dsplit[0]))+
        parseInt(dsplit[1])/60.+
        parseFloat(dsplit[2])/3600.;
        
        if (parseFloat(dsplit[0]) < 0) {
            dec *= -1;
        }
        
        recenter(ra,dec);
        latLng2raDec();
    }
    
    // Globals
    var myIcon = new GIcon();
    myIcon.iconSize = new GSize(30, 25);
    myIcon.iconAnchor = new GPoint(14.5, 14.5);
    var lats = [];
    var lngs = [];
    var ids = [];
    var nObject = 0;

    // Read objects from XML file and plot regions
    function plotXmlObjects() {
        GDownloadUrl(root+".xml", function(data) {
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
    
    function recenter(ra,dec) {
        var lat = dec;
        var lng = (360-ra)-centerLng+offset
        map.setCenter(new GLatLng(lat,lng), zoomLevel);
    }
    
            </script>
        """ %(center[1],lng_offset,mapParams['ZOOMLEVEL'],
              threedhst.currentRun['root_direct'],
              llSW[0],llSW[1]-center[1]+lng_offset,
              llNE[0],llNE[1]-center[1]+lng_offset,
              center[0]))
    
    #### HTML Body   
    lines.append("""

    </head>
    <body onload="initialize()" onunload="GUnload()">
    
    <div id="map"></div>
    
    <div id="title">
        %s
    </div>
    <img src="scripts/3dHST.png" id="logo">
    <div onclick="javascript:switch_layout()" id="switchbox">
        ||
    </div>
    
    <div id="centerbox"></div>
    
    <div id="coords">
        <form>
        <input type="text" value="00:00:00.00" class="cinput" id="raInput" maxlength="11" onchange="centerOnInput()"/>
        <input type="text" value="+00:00:00.0" class="cinput" id="decInput" maxlength="11" onchange="centerOnInput()"/>
        </form>
    </div>
    
    """ %(title))
    
    lines.append("""
    
    <div id="content" style="width:840px;height:100%%;overflow:auto">
    
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
        for idx,num in enumerate(
          mySexCat.columns[mySexCat.searchcol('NUMBER')].entry):
            if num == str(id):
                break
        
        ra  = mySexCat.columns[mySexCat.searchcol('X_WORLD')].entry[idx]
        dec = mySexCat.columns[mySexCat.searchcol('Y_WORLD')].entry[idx]
        mag = mySexCat.columns[mySexCat.searchcol('MAG_F1392W')].entry[idx]
        
        img = '%s_%04d' %(root,id)
        lines.append("""
        <tr> 
            <td id="i%d" onclick="javascript:recenter(%13.6f,%13.6f)">
                  %d
            </td>
            <td>%13.6f</td> 
            <td>%13.6f</td> 
            <td>%6.2f</td> 
            <td><a href='images/%s_thumb.png'><img src='images/%s_thumb.png' width=133px></a></td> 
            <td><a href='ascii/%s.dat.gz'><img src='images/%s_1D.png' width=200px title='ascii'></a></td> 
            <td><a href='images/%s_2D.png'><img src='images/%s_2D.png' width=200px></a></td> 
        </tr> 
        """ %(id,np.float(ra),np.float(dec),id,
              np.float(ra),np.float(dec),np.float(mag),
              img,img,img,img,img,img))
    
    lines.append("""
    </tbody>
    </div> <!-- content -->
    </body>
    </html>""")
    
    if not output:
        output='HTML/index.html'
    
    fp = open(output,'w')
    fp.writelines(lines)
    fp.close()
    
def asciiSpec(SPCFile, root="spec", path="../HTML/ascii"):
    """
asciiSpec(SPCFile, root="spec", path="../HTML/ascii")
    
    Put ASCII spectra in HTML/ascii.
    """
    import os
    import gzip
    import tarfile
    
    try:
        os.mkdir(path)
    except:
        pass
    
    ids = SPCFile._ext_map+0
    ids.sort()
    noNewLine = '\x1b[1A\x1b[1M'
    
    for id in ids:
        spec = SPCFile.getSpec(id)
        lam  = spec.field('LAMBDA')
        flux = spec.field('FLUX')
        ferr = spec.field('FERROR')
        contam = spec.field('CONTAM')
    
        out = root+'_%04d.dat' %id
        print noNewLine+out
        
        fp = gzip.open(path+'/'+out+'.gz','wb')
        NL = len(lam)
        for i in range(NL):
            fp.write('%11.5e %10.3e %10.3e %10.3e\n' %(lam[i],flux[i],ferr[i],contam[i]))
        fp.close()
        
    ### make a tarfile
    oldwd = os.getcwd()
    os.chdir(path)
    fptar = tarfile.open(root+'_spec.tar.gz','w|gz')
    for id in ids:
        out = root+'_%04d.dat.gz' %id
        fptar.add(out)
    fptar.close()
    os.chdir(oldwd)
    
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
            self._ext_map[i] = np.int(
               self.fits[i+1].header['EXTNAME'].split('BEAM_')[1][:-1])
            
    def __init__(self, filename='ib3721050_2_opt.SPC.fits',
                 axe_drizzle_dir='DRIZZLE_G141'):
        """
__init__(filename='ib3721050_2_opt.SPC.fits',
         axe_drizzle_dir='DRIZZLE_G141')
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
            