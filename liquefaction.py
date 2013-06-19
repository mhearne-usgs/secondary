#!/usr/bin/env python

#stdlib imports
import ConfigParser
import os.path
import sys
import glob
import copy
import math
import warnings
from optparse import OptionParser
import glob

#turn off all warnings...
warnings.filterwarnings('ignore')

#third party imports
import numpy
from matplotlib import pyplot
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib import cm
from matplotlib.colors import LightSource
from matplotlib.table import Table
import gdal
import osr
from pylab import find

#local imports
from losspager.io.esri import EsriGrid
from losspager.io.shake import ShakeGrid
from losspager.io.grid import Grid
from losspager.io.gmt import GMTGrid
from losspager.io.shapefile import PagerShapeFile
from losspager.util.text import decToRoman,roundToNearest,ceilToNearest,floorToNearest,commify
from losspager.map import poly

CONFIGFILE = 'config.ini'
LQCOEFF = ['b0','bpga','bcti','bvs30']
LSCOEFF = ['b0','bpga','bslope','bcohesion','bpgaslope']
DEG2M = 111191
ALPHA = 0.7

def getShapes(config):
    sections = config.sections()
    shapedict = {}
    reqfields = ['filename','color']
    optfields = ['linewidth']
    for section in sections:
        if section.lower().find('shapefile_') > -1:
            t,shapename = section.lower().split('_')
            shapedict[shapename] = {}
            options = config.options(section)
            missing = []
            for field in reqfields:
                if field not in options:
                    missing.append(field)
            if len(missing):
                raise Exception,'Missing fields %s in %s' % (','.join(missing),section)
            for option in options:
                if option in reqfields or option in optfields:
                    shapedict[shapename][option.lower()] = config.get(section,option)
    return shapedict
    

def cosd(angle):
    """
    Return cosine of angle expressed in degrees.
    @param angle: Input angle in degrees.
    @return: Cosine of input angle.
    """
    return math.cos(angle*(math.pi/180))

def saveTiff(lqgrid,fname,isFloat=False):
    format = "GTiff"
    driver = gdal.GetDriverByName( format )
    metadata = driver.GetMetadata()
    nrows,ncols = lqgrid.griddata.shape
    if not isFloat:
        raster = numpy.array(numpy.round(lqgrid.griddata*100),dtype=numpy.int8)
        outdata = driver.Create(fname, ncols, nrows, 1, gdal.GDT_Byte )
    else:
        raster = lqgrid.griddata.astype(numpy.float32)
        outdata = driver.Create(fname, ncols, nrows, 1, gdal.GDT_Float32 )
    ulxmap = lqgrid.geodict['xmin']
    ulymap = lqgrid.geodict['ymax']
    xdim = lqgrid.geodict['xdim']
    ydim = lqgrid.geodict['ydim']
    outdata.SetGeoTransform([ulxmap,xdim,0,ulymap,0,-1*ydim])
    srs = osr.SpatialReference()
    srs.SetGeogCS( "My geographic coordinate system",
               "WGS_1984", 
               "My WGS84 Spheroid", 
               osr.SRS_WGS84_SEMIMAJOR,
               osr.SRS_WGS84_INVFLATTENING, 
               "Greenwich", 0.0,
        "degree", float(osr.SRS_UA_DEGREE_CONV))
    outdata.SetProjection(srs.ExportToWkt())
    outdata.GetRasterBand(1).WriteArray(raster)
    outdata = None

def latstr(parallel):
    if parallel < 0:
        parstr = str(-1*parallel) + '$\degree$ S'
    else:
        parstr = str(parallel) + '$\degree$ N'
    return parstr

def lonstr(meridian):
    if meridian < 0:
        merstr = str(-1*meridian) + '$\degree$ W'
    else:
        merstr = str(meridian) + '$\degree$ E'
    return merstr

def makeLandslideMap(topogrid,lsgrid,title=None,eventid=None,roads=None,mode='landslide'):
    #topo and liquefaction grids are already resampled to same extent and resolution
    bounds = topogrid.getRange()
    xmin,xmax,ymin,ymax = bounds
    cx = xmin + (xmax-xmin)/2.0
    cy = ymin + (ymax-ymin)/2.0
    figx = 11
    figy = 8.5
    fig = pyplot.figure(figsize=(figx,figy))
    axleft = 0.27
    axbottom = 0.21
    axwidth = 0.45
    axheight = 0.59
    ax = fig.add_axes([axleft,axbottom,axwidth,axheight])
    # setup of basemap ('lcc' = lambert conformal conic).
    # use major and minor sphere radii from WGS84 ellipsoid.
    m = Basemap(llcrnrlon=xmin,llcrnrlat=ymin,urcrnrlon=xmax,urcrnrlat=ymax,
                rsphere=(6378137.00,6356752.3142),suppress_ticks=False,
        resolution='h',projection='tmerc',lat_0=cy,lon_0=cx,ax=ax,fix_aspect=False)

    #keep the same resolution in the mapped data
    topodict = topogrid.getGeoDict()
    xdim_m = topodict['xdim'] * DEG2M * cosd(cy)
    ydim_m = topodict['ydim'] * DEG2M
    topodat = numpy.flipud(topogrid.getData().copy())

    #draw the ocean in light blue
    water_color = [.47,.60,.81]
    #water_color = [.0,.73,1.0]
    palette = cm.binary
    i = numpy.where(topodat < 0)
    topodat[i] = -1
    topodatm = numpy.ma.masked_values(topodat, -1)
    palette.set_bad(water_color,1.0)

    zdict = {0:0.00000898,
             10:0.00000912,
             20:0.00000956,
             30:0.00001036,
             40:0.00001171,
             50:0.00001395,
             60:0.00001792,
             70:0.00002619,
             80:0.00005156}
    mlat = abs(int(round(cy/10)*10))
    zfactor = zdict[mlat]
    ls = LightSource(azdeg = 90, altdeg = 20)
    #pyplot.set_cmap('bone')
    
    rgb = ls.shade(topodatm*zfactor,cmap=palette)
    im = m.imshow(rgb)
    pyplot.hold(True)

    #draw the secondary hazard probabilities (anything above 1%)
    lsdat = numpy.flipud(lsgrid.getData().copy()) * 100.0
    clear_color = [0,0,0,0.0]
    if mode == 'landslide':
        palette = cm.cool
    else:
        palette = cm.autumn_r
    i = numpy.where(lsdat < 2.0)
    lsdat[i] = 0
    lsdatm = numpy.ma.masked_equal(lsdat, 0)
    palette.set_bad(clear_color,alpha=0.0)
    lsmin = numpy.nanmin(lsdat)
    lsmax = numpy.nanmax(lsdat)
    lsprobhandle = m.imshow(lsdatm,cmap=palette,vmin=lsmin,vmax=lsmax,alpha=ALPHA)
    m.colorbar(mappable=lsprobhandle)

    #optionally draw the road networks
    if roads is not None:
        roadcolor = '#6E6E6E'
        psf = PagerShapeFile(roads)
        shapes = psf.getShapesByBoundingBox((xmin,xmax,ymin,ymax))
        for shape in shapes:
            mapxroad,mapyroad = m(shape['x'],shape['y'])
            m.plot(mapxroad,mapyroad,color=roadcolor)

    meridians = getMapLines(bounds[0],bounds[1])
    parallels = getMapLines(bounds[2],bounds[3])
    
    xmap_range = m.xmax-m.xmin
    ymap_range = m.ymax-m.ymin
    xoff = -0.09*(xmap_range)
    yoff = -0.04*(ymap_range)

    #do tick stuff
    fig.canvas.draw() #have to do this for tick stuff to show
    xlabels = [lonstr(mer) for mer in meridians]
    ylabels = [latstr(par) for par in parallels]
    xticks = []
    yticks = []
    for i in range(0,len(meridians)):
        lat = ymin
        lon = meridians[i]
        x,y = m(lon,lat)
        xticks.append(x)
    for i in range(0,len(parallels)):
        lon = xmin
        lat = parallels[i]
        x,y = m(lon,lat)
        yticks.append(y)
    pyplot.xticks(xticks,xlabels,size=6)
    pyplot.yticks(yticks,ylabels,size=6)
    pyplot.tick_params(axis='both',direction='out',right='on')
    m.drawmapboundary(color='k',linewidth=2.0)

    if title is not None:
        ax.set_title(title)
    if eventid is not None:
        pyplot.savefig(eventid+'_%s.pdf' % mode)
    else:
        pyplot.savefig('output.pdf')

  
        
def makeMatMap(topogrid,lqgrid,lsgrid,coastshapefile,riverfolder,isScenario=False,
               roads=None,borderfile=None,shapedict=None,title=None,legdict=None,
               epicenter=None,eventfolder=None):
    bounds = topogrid.getRange()
    xmin,xmax,ymin,ymax = bounds
    cx = xmin + (xmax-xmin)/2.0
    cy = ymin + (ymax-ymin)/2.0
    figx = 11
    figy = 8.5
    fig = pyplot.figure(figsize=(figx,figy))
    rect = fig.patch
    rect.set_facecolor((1,1,1,1))
    
    axleft = 0.27
    axbottom = 0.21
    axwidth = 0.45
    axheight = 0.59
    ax = fig.add_axes([axleft,axbottom,axwidth,axheight])
    #keep the same resolution in the mapped data
    topodict = topogrid.getGeoDict()
    xdim_m = topodict['xdim'] * DEG2M * cosd(cy)
    ydim_m = topodict['ydim'] * DEG2M
    topodat = topogrid.getData().copy()
    #flag the regions where topography is less than 0 (we'll color this ocean later)
    i = numpy.where(topodat < 0)

    #do shaded relief stuff
    #keys are latitude values
    #values are multiplication factor
    zdict = {0:0.00000898,
             10:0.00000912,
             20:0.00000956,
             30:0.00001036,
             40:0.00001171,
             50:0.00001395,
             60:0.00001792,
             70:0.00002619,
             80:0.00005156}
    #find the mean latitude of the map we're making, and use that to come up with a zfactor
    mlat = abs(int(round(cy/10)*10))
    zfactor = zdict[mlat]
    ls = LightSource(azdeg = 90, altdeg = 20)

    #draw the ocean in light blue
    water_color = [.47,.60,.81]
    palette1 = copy.deepcopy(cm.binary)
        
    #draw the light shaded-topography
    rgbdata = topodat*zfactor #apply the latitude specific zfactor correction
    rgb = ls.shade(rgbdata,cmap=palette1) #apply the light shading to our corrected topography

    #this is an rgb data set now, so masking the pixels won't work, but explicitly setting all of the
    #"bad" pixels to our water color will
    red = rgb[:,:,0]
    green = rgb[:,:,1]
    blue = rgb[:,:,2]
    red[i] = water_color[0]
    green[i] = water_color[1]
    blue[i] = water_color[2]
    rgb[:,:,0] =  red
    rgb[:,:,1] =  green
    rgb[:,:,2] =  blue
    
    im = pyplot.imshow(rgb,origin='upper',extent=(xmin,xmax,ymin,ymax),cmap=palette1)
    axlimits = pyplot.axis()
    pyplot.hold(True)

    print 'Drawing rivers...'
    #draw the rivers
    #first, find the shapefile that has the rivers for our map in it!
    allshapes = []
    riverfiles = glob.glob(os.path.join(riverfolder,'*river*.shp'))
    for riverfile in riverfiles:
        shpobj  = PagerShapeFile(riverfile)
        shapes = shpobj.getShapesByBoundingBox((xmin,xmax,ymin,ymax))
        allshapes = allshapes + shapes
    for shape in allshapes:
        pyplot.plot(shape['x'],shape['y'],color=water_color,lw=1)
    pyplot.axis(axlimits)

    #draw the epicenter as a black star
    if epicenter is not None:
        elat,elon = epicenter
        pyplot.plot(elon,elat,'*',markeredgecolor='k',mfc='None',mew=1.5,ms=24)

    print 'Drawing borders...'
    #draw the country borders
    if borderfile is not None:
        bordercolor = '#00FF00'
        bordershape = PagerShapeFile(borderfile)
        shapes = bordershape.getShapesByBoundingBox((xmin,xmax,ymin,ymax))
        for shape in shapes:
            pyplot.plot(shape['x'],shape['y'],color=bordercolor,lw=1.5)

    print 'Drawing coastlines...'
    #load the coastlines shapefile - just get the shapes that are in this extent
    shpobj  = PagerShapeFile(coastshapefile)
    shapes = shpobj.getShapesByBoundingBox((xmin,xmax,ymin,ymax))
    coastcolor = 'k'
    for shape in shapes:
        pyplot.plot(shape['x'],shape['y'],color=coastcolor,lw=1)
    pyplot.axis(axlimits)

    #get the polygon (if any) for a shapefile called "boundary"
    bx = None
    by = None
    for shapename,shapeinfo in shapedict.iteritems():
        if shapename.find('boundary') > -1:
            fname = shapeinfo['filename']
            psf = PagerShapeFile(fname)
            boundshape = psf.getShape(0)
            bx = boundshape['x']
            by = boundshape['y']

    print 'Drawing roads...'
    #optionally draw the road networks
    if roads is not None:
        roadcolor = '#6E6E6E'
        psf = PagerShapeFile(roads)
        shapes = psf.getShapesByBoundingBox((xmin,xmax,ymin,ymax))
        if bx is not None:
            pp = poly.PagerPolygon(bx,by)
        print 'Found %i road segments...' % len(shapes)
        for shape in shapes:
            sx = shape['x']
            sy = shape['y']
            if bx is not None:
                outside = numpy.logical_not(pp.containsPoints(sx,sy))
                sx = sx[outside]
                sy = sy[outside]
            pyplot.plot(sx,sy,color=roadcolor)

    #draw the arbitary shapefiles the user may have specified in the config file
    for shapename,shapeinfo in shapedict.iteritems():
        fname = shapeinfo['filename']
        color = shapeinfo['color']
        if shapeinfo.has_key('linewidth'):
            lw = shapeinfo['linewidth']
        else:
            lw = 1
        psf = PagerShapeFile(fname)
        shapes = psf.getShapes()
        for shape in shapes:
            pyplot.plot(shape['x'],shape['y'],color=color,linewidth=lw)

    print 'Rendering liquefaction...'
    #draw the liquefaction probabilities (anything above 1%) in an orange-red scale (OrRd)
    lqdat = lqgrid.getData().copy() * 100.0
    clear_color = [0,0,0,0.0]
    palette = cm.autumn_r
    i = numpy.where(lqdat < 2.0)
    lqdat[i] = 0
    lqdatm = numpy.ma.masked_equal(lqdat, 0)
    palette.set_bad(clear_color,alpha=0.0)
    probhandle = pyplot.imshow(lqdatm,cmap=palette,vmin=2.0,vmax=20.0,alpha=ALPHA,origin='upper',extent=(xmin,xmax,ymin,ymax))

    print 'Rendering landslide...'
    #draw the landslide probabilities (anything above 1%) in an orange-red scale (OrRd)
    lsdat = lsgrid.getData().copy() * 100.0
    clear_color = [0,0,0,0.0]
    palette = cm.cool
    i = numpy.where(lsdat < 2.0)
    lsdat[i] = 0
    lsdatm = numpy.ma.masked_equal(lsdat, 0)
    palette.set_bad(clear_color,alpha=0.0)
    lsprobhandle = pyplot.imshow(lsdatm,cmap=palette,vmin=2.0,vmax=20.0,alpha=ALPHA,origin='upper',extent=(xmin,xmax,ymin,ymax))
            
    #get the meridians we want to see
    meridians = getMapLines(bounds[0],bounds[1])
    parallels = getMapLines(bounds[2],bounds[3])

    xmap_range = xmax-xmin
    ymap_range = ymax-ymin
    xoff = -0.09*(xmap_range)
    yoff = -0.04*(ymap_range)

    #do tick stuff
    fig.canvas.draw() #have to do this for tick stuff to show
    xlabels = [lonstr(mer) for mer in meridians]
    ylabels = [latstr(par) for par in parallels]
    xticks = []
    yticks = []
    for i in range(0,len(meridians)):
        lat = ymin
        lon = meridians[i]
        xticks.append(lon)
    for i in range(0,len(parallels)):
        lon = xmin
        lat = parallels[i]
        yticks.append(lat)
    pyplot.xticks(xticks,xlabels,size=6)
    pyplot.yticks(yticks,ylabels,size=6)
    pyplot.tick_params(axis='both',direction='out',right='on')

    #draw the map boundary
    pyplot.plot([xmin,xmin],[ymin,ymax],lw=2,color='k')
    pyplot.plot([xmin,xmax],[ymin,ymin],lw=2,color='k')
    pyplot.plot([xmax,xmax],[ymin,ymax],lw=2,color='k')
    pyplot.plot([xmin,xmax],[ymax,ymax],lw=2,color='k')

    #Establish where the corners of the map are in display space (pixels)
    #this is the common coordinate system for the transformation functions of axes, figure, etc.
    ulx_pixels,uly_pixels = ax.transData.transform((xmin,ymax))
    llx_pixels,lly_pixels = ax.transData.transform((xmin,ymin))
    urx_pixels,ury_pixels = ax.transData.transform((xmax,ymax))
    lrx_pixels,lry_pixels = ax.transData.transform((xmax,ymin))

    #establish how far (in inches) we want all of the elements to be
    landslide_title_x_offset = 1.125 #to the left of the map
    landslide_title_y_offset = 0.15 #above the map
    liquefaction_title_x_offset = 0.375 #to the right of the map
    liquefaction_title_y_offset = 0.15 #to the right of the map
    liquefaction_ticks_x_offset = 0.795 #all colorbar annotation on the right side of the map
    landslide_annotation_x_offset = 0.825 #all colorbar annotation on the left side of the map
    landslide_colorbar_axis_x_offset = 0.75 #how far is the landslide colorbar axis to the left of the map
    liquefaction_colorbar_axis_x_offset = 0.5 #how far is the liquefaction colorbar axis to the right of the map
    landslide_colorbar_tick_width = 0.0625 #how wide should the tick marks on the landslide colorbar be?
    colorbar_width = 0.25 #how wide should the colorbar be?

    #how high should the colorbars be (in figure units)
    colorbar_height_pixels = uly_pixels - lly_pixels
    colorbar_width_pixels = colorbar_width * fig.dpi
    figure_width_pixels = figx * fig.dpi
    figure_height_pixels = figy * fig.dpi
    colorbar_height_figure = colorbar_height_pixels/figure_height_pixels
    colorbar_width_figure = colorbar_width_pixels/figure_width_pixels

    #draw the landslide probability colorbar title
    px = ulx_pixels - (landslide_title_x_offset * fig.dpi)
    py = uly_pixels + (landslide_title_y_offset * fig.dpi)
    xtitle,ytitle = ax.transData.inverted().transform((px,py))
    cbartitle = 'Landslide\nProbability'
    pyplot.text(xtitle,ytitle,cbartitle,multialignment='left',axes=ax)

    #draw the liquefaction probability colorbar title
    px = urx_pixels + (liquefaction_title_x_offset * fig.dpi)
    py = ury_pixels + (liquefaction_title_y_offset * fig.dpi)
    xtitle,ytitle = ax.transData.inverted().transform((px,py))
    cbartitle = 'Liquefaction\nProbability'
    pyplot.text(xtitle,ytitle,cbartitle,multialignment='left',axes=ax)

    #draw the liquefaction probability colorbar annotation
    yrange = ymax - ymin
    xrange = xmax - xmin
    cticks = [2.0,5.0,10.0,15.0,20.0]
    ctickfrac = (numpy.array(cticks) - 2)/(max(cticks) - 2)
    places = numpy.diff(ctickfrac)/2.0 + ctickfrac[0:-1]
    yplaces = places * yrange + ymin
    px = urx_pixels + (liquefaction_ticks_x_offset * fig.dpi)
    xplace,yt = ax.transData.inverted().transform((px,100))
    levels = ['Low','Medium','High','V. High']
    for i in range(0,len(yplaces)):
        pyplot.text(xplace,yplaces[i],levels[i],verticalalignment='center')

    #draw the landslide probability colorbar annotation
    px = ulx_pixels - (landslide_annotation_x_offset * fig.dpi)
    xplace,yt = ax.transData.inverted().transform((px,100))
    for i in range(0,len(yplaces)):
        pyplot.text(xplace,yplaces[i],levels[i],verticalalignment='center',horizontalalignment='right')

    #draw the tick labels on the left side of the landslide colorbar
    ctickfrac = (numpy.array(cticks) - 2)/(max(cticks) - 2)
    yplaces = ctickfrac * yrange + ymin
    px = ulx_pixels - (landslide_annotation_x_offset * fig.dpi)
    xplace,yt = ax.transData.inverted().transform((px,100))
    cticklabels = [str(int(tick))+'%' for tick in cticks]
    for i in range(0,len(yplaces)):
        pyplot.text(xplace,yplaces[i],cticklabels[i],verticalalignment='center',horizontalalignment='right')

    #draw the landslide colorbar
    px = llx_pixels - (landslide_colorbar_axis_x_offset * fig.dpi)
    axisleft,axisbottom = fig.transFigure.inverted().transform((px,lly_pixels)) #now in figure unit space
    
    cax1 = fig.add_axes([axisleft,axisbottom,colorbar_width_figure,colorbar_height_figure])
    cb1 = pyplot.colorbar(mappable=lsprobhandle,cax=cax1,ticks=[])

    #draw the tick marks on the landslide colorbar manually
    px = ulx_pixels - (landslide_colorbar_axis_x_offset * fig.dpi) #left edge of the tickmark in pixels
    px2 = px + (landslide_colorbar_tick_width * fig.dpi) #right edge of the tickmark in pixels
    xplace,yt = cax1.transData.inverted().transform((px,100)) #left edge of the tickmark in axis data space
    xplace2,yt = cax1.transData.inverted().transform((px2,100)) #right edge of the tickmark in axis data space
    for i in range(0,len(yplaces)):
        yplace = ctickfrac[i]        
        pyplot.plot([xplace,xplace2],[yplace,yplace],'k-',axes=cax1,lw=0.5)

    #draw the liquefaction colorbar
    px = lrx_pixels + (liquefaction_colorbar_axis_x_offset * fig.dpi)
    axisleft,axisbottom = fig.transFigure.inverted().transform((px,lry_pixels))
    cax2 = fig.add_axes([axisleft,axisbottom,colorbar_width_figure,colorbar_height_figure])
    cticks = [2.0,5.0,10.0,15.0,20.0]
    cb2 = pyplot.colorbar(mappable=probhandle,cax=cax2)
    cticklabels = [str(int(tick))+'%' for tick in cticks]
    cb2.set_ticks(cticks)
    cb2.set_ticklabels(cticklabels)

    #make two tables at the bottom with the secondary hazard statistics
    if legdict is not None:
        lsarea = legdict['lssum']/1e6
        if lsarea < 0.1:
            lsarea = '< 0.1'
        else:
            lsarea = commify(int(lsarea)) + ' km$^2$'
        lsncells = commify(legdict['lscells'])
        lspmax = '%.1f%%' % (legdict['lsmax']*100)
        lspmean = '%.1f%%' % (legdict['lsmean']*100)
        lspslope = '%.1f$\degree$' % legdict['slopemin']
        tcells = [['Landslide Area',lsarea],
                  ['# of cells',lsncells],
                  ['Maximum probability',lspmax],
                  ['Mean probability',lspmean],
                  ['Minimum slope',lspslope]]
        ax_lstable = fig.add_axes([0.215,0.08,0.2,0.1])
        pyplot.axis('off')
        lstable = pyplot.table(cellText=tcells,loc='lower left')

        lqarea = legdict['liqsum']/1e6
        if lqarea < 0.1:
            lqarea = '< 0.1'
        else:
            lqarea = commify(int(lqarea)) + ' km$^2$'
        lqncells = commify(legdict['liqcells'])
        lqpmax = '%.1f%%' % (legdict['liqmax']*100)
        lqpmean = '%.1f%%' % (legdict['liqmean']*100)
        lqpslope = '%.1f$\degree$' % legdict['slopemax']
        tcells = [['Liquefaction Area',lqarea],
                  ['# of cells',lqncells],
                  ['Maximum probability',lqpmax],
                  ['Mean probability',lqpmean],
                  ['Maximum slope',lqpslope]]
        ax_lqtable = fig.add_axes([0.57,0.08,0.2,0.1])
        pyplot.axis('off')
        lqtable = pyplot.table(cellText=tcells,loc='lower left')
    
    pyplot.sca(ax)
    #put the title in
    if title is not None:
        pyplot.title(title)
    
    if isScenario:
        pyplot.text(cx,cy,'SCENARIO',rotation=45,alpha=0.25,size=72,ha='center',va='center',color='red')
    
    if eventfolder is None:
        outfile = os.path.join(os.getcwd(),eventid+'.pdf')
        outfilepng = os.path.join(os.getcwd(),eventid+'.png')
    else:
        outfile = os.path.join(eventfolder,'secondary_hazards.pdf')
        outfilepng = os.path.join(eventfolder,'secondary_hazards.png')
        
    pyplot.savefig(outfile,facecolor=fig.get_facecolor())
    pyplot.savefig(outfilepng,facecolor=fig.get_facecolor())
    return outfile
    
    
def makeMap(topogrid,lqgrid,lsgrid,title=None,eventid=None,legdict=None,roads=None,isScenario=False,shapes=None):
    #topo and liquefaction grids are already resampled to same extent and resolution
    bounds = topogrid.getRange()
    xmin,xmax,ymin,ymax = bounds
    cx = xmin + (xmax-xmin)/2.0
    cy = ymin + (ymax-ymin)/2.0
    figx = 11
    figy = 8.5
    fig = pyplot.figure(figsize=(figx,figy))
    axleft = 0.27
    axbottom = 0.21
    axwidth = 0.45
    axheight = 0.59
    ax = fig.add_axes([axleft,axbottom,axwidth,axheight])
    # setup of basemap ('lcc' = lambert conformal conic).
    # use major and minor sphere radii from WGS84 ellipsoid.
    m = Basemap(llcrnrlon=xmin,llcrnrlat=ymin,urcrnrlon=xmax,urcrnrlat=ymax,
                rsphere=(6378137.00,6356752.3142),suppress_ticks=False,
        resolution='h',projection='tmerc',lat_0=cy,lon_0=cx,ax=ax,fix_aspect=False)

    #keep the same resolution in the mapped data
    topodict = topogrid.getGeoDict()
    xdim_m = topodict['xdim'] * DEG2M * cosd(cy)
    ydim_m = topodict['ydim'] * DEG2M
    topodat = numpy.flipud(topogrid.getData().copy())

    #draw the ocean in light blue
    water_color = [.47,.60,.81]
    #water_color = [.0,.73,1.0]
    palette = cm.binary
    i = numpy.where(topodat < 0)
    topodat[i] = -1
    topodatm = numpy.ma.masked_values(topodat, -1)
    #palette.set_bad(water_color,1.0)

    #do shaded relief stuff
    #z-factors
    # Latitude     Z-factor
    #  0           0.00000898
    # 10           0.00000912
    # 20           0.00000956
    # 30           0.00001036
    # 40           0.00001171
    # 50           0.00001395
    # 60           0.00001792
    # 70           0.00002619
    # 80           0.00005156
    zdict = {0:0.00000898,
             10:0.00000912,
             20:0.00000956,
             30:0.00001036,
             40:0.00001171,
             50:0.00001395,
             60:0.00001792,
             70:0.00002619,
             80:0.00005156}
    mlat = abs(int(round(cy/10)*10))
    zfactor = zdict[mlat]
    ls = LightSource(azdeg = 90, altdeg = 20)
    #pyplot.set_cmap('bone')
    
    rgb = ls.shade(topodatm*zfactor,cmap=palette)
    im = m.imshow(rgb)
    pyplot.hold(True)

    #draw the liquefaction probabilities (anything above 1%) in an orange-red scale (OrRd)
    lqdat = numpy.flipud(lqgrid.getData().copy()) * 100.0
    clear_color = [0,0,0,0.0]
    palette = cm.autumn_r
    i = numpy.where(lqdat < 2.0)
    lqdat[i] = 0
    lqdatm = numpy.ma.masked_equal(lqdat, 0)
    palette.set_bad(clear_color,alpha=0.0)
    probhandle = m.imshow(lqdatm,cmap=palette,vmin=2.0,vmax=20.0,alpha=ALPHA)

    #draw the landslide probabilities (anything above 1%) in an orange-red scale (OrRd)
    lsdat = numpy.flipud(lsgrid.getData().copy()) * 100.0
    clear_color = [0,0,0,0.0]
    palette = cm.cool
    i = numpy.where(lsdat < 2.0)
    lsdat[i] = 0
    lsdatm = numpy.ma.masked_equal(lsdat, 0)
    palette.set_bad(clear_color,alpha=0.0)
    lsprobhandle = m.imshow(lsdatm,cmap=palette,vmin=2.0,vmax=20.0,alpha=ALPHA)

    meridians = getMapLines(bounds[0],bounds[1])
    parallels = getMapLines(bounds[2],bounds[3])
    
    xmap_range = m.xmax-m.xmin
    ymap_range = m.ymax-m.ymin
    xoff = -0.09*(xmap_range)
    yoff = -0.04*(ymap_range)

    #do tick stuff
    fig.canvas.draw() #have to do this for tick stuff to show
    xlabels = [lonstr(mer) for mer in meridians]
    ylabels = [latstr(par) for par in parallels]
    xticks = []
    yticks = []
    for i in range(0,len(meridians)):
        lat = ymin
        lon = meridians[i]
        x,y = m(lon,lat)
        xticks.append(x)
    for i in range(0,len(parallels)):
        lon = xmin
        lat = parallels[i]
        x,y = m(lon,lat)
        yticks.append(y)
    pyplot.xticks(xticks,xlabels,size=6)
    pyplot.yticks(yticks,ylabels,size=6)
    pyplot.tick_params(axis='both',direction='out',right='on')
    m.drawmapboundary(color='k',linewidth=2.0)

    #make a second axes on the right side - this undoes the x tick labels, and these y labels aren't
    #in the same vertical location as those on the left side.  ??
    # yaxis2 = pyplot.twinx()
    # pyplot.yticks(yticks,ylabels,size=6)
    # pyplot.tick_params(axis='both',direction='out',right='on')

    #im = m.imshow(topodat,cm.GMT_haxby)
    m.drawcoastlines()
    border_color = [0.18,1.0,0.18]
    m.drawcountries(color=border_color,linewidth=2.0)
    m.drawstates(color=border_color,linewidth=2.0)
    m.drawrivers(color=water_color,linewidth=1.0)
    m.fillcontinents(color=[0,0,0,0],lake_color=water_color)

    #optionally draw the road networks
    if roads is not None:
        roadcolor = '#6E6E6E'
        psf = PagerShapeFile(roads)
        shapes = psf.getShapesByBoundingBox((xmin,xmax,ymin,ymax))
        for shape in shapes:
            pyplot.plot(shape['x'],shape['y'],color=roadcolor)

    #draw the arbitary shapefiles the user may have specified in the config file
    # for shapename,shapeinfo in shapes.iteritems():
    #     fname = shapeinfo['filename']
    #     color = shapeinfo['color']
    #     if shapeinfo.has_key('linewidth'):
    #         lw = shapeinfo['linewidth']
    #     else:
    #         lw = 1
    #     psf = PagerShapeFile(fname)
    #     shapes = psf.getShapes()
    #     for shape in shapes:
    #         m.plot(shape['x'],shape['y'],color=color,linewidth=lw)

    #draw SCENARIO in big red letters across the map if it's a scenario
    if isScenario:
        centermapx,centermapy = m(cx,cy)
        pyplot.text(centermapx,centermapy,'SCENARIO',rotation=45,alpha=0.25,size=72,ha='center',va='center',color='red')
            
    #Establish where the corners of the map are in display space (pixels)
    #this is the common coordinate system for the transformation functions of axes, figure, etc.
    ulx_pixels,uly_pixels = ax.transData.transform((m.xmin,m.ymax))
    llx_pixels,lly_pixels = ax.transData.transform((m.xmin,m.ymin))
    urx_pixels,ury_pixels = ax.transData.transform((m.xmax,m.ymax))
    lrx_pixels,lry_pixels = ax.transData.transform((m.xmax,m.ymin))

    #establish how far (in inches) we want all of the elements to be
    landslide_title_x_offset = 1.125 #to the left of the map
    landslide_title_y_offset = 0.15 #above the map
    liquefaction_title_x_offset = 0.375 #to the right of the map
    liquefaction_title_y_offset = 0.15 #to the right of the map
    liquefaction_ticks_x_offset = 0.795 #all colorbar annotation on the right side of the map
    landslide_annotation_x_offset = 0.825 #all colorbar annotation on the left side of the map
    landslide_colorbar_axis_x_offset = 0.75 #how far is the landslide colorbar axis to the left of the map
    liquefaction_colorbar_axis_x_offset = 0.5 #how far is the liquefaction colorbar axis to the right of the map
    landslide_colorbar_tick_width = 0.0625 #how wide should the tick marks on the landslide colorbar be?
    colorbar_width = 0.25 #how wide should the colorbar be?

    #how high should the colorbars be (in figure units)
    colorbar_height_pixels = uly_pixels - lly_pixels
    colorbar_width_pixels = colorbar_width * fig.dpi
    figure_width_pixels = figx * fig.dpi
    figure_height_pixels = figy * fig.dpi
    colorbar_height_figure = colorbar_height_pixels/figure_height_pixels
    colorbar_width_figure = colorbar_width_pixels/figure_width_pixels
    
    #draw the landslide probability colorbar title
    px = ulx_pixels - (landslide_title_x_offset * fig.dpi)
    py = uly_pixels + (landslide_title_y_offset * fig.dpi)
    xtitle,ytitle = ax.transData.inverted().transform((px,py))
    cbartitle = 'Landslide\nProbability'
    pyplot.text(xtitle,ytitle,cbartitle,multialignment='left',axes=ax)

    #draw the liquefaction probability colorbar title
    px = urx_pixels + (liquefaction_title_x_offset * fig.dpi)
    py = ury_pixels + (liquefaction_title_y_offset * fig.dpi)
    xtitle,ytitle = ax.transData.inverted().transform((px,py))
    cbartitle = 'Liquefaction\nProbability'
    pyplot.text(xtitle,ytitle,cbartitle,multialignment='left',axes=ax)

    #draw the liquefaction probability colorbar annotation
    yrange = m.ymax - m.ymin
    xrange = m.xmax - m.xmin
    cticks = [2.0,5.0,10.0,15.0,20.0]
    ctickfrac = (numpy.array(cticks) - 2)/(max(cticks) - 2)
    places = numpy.diff(ctickfrac)/2.0 + ctickfrac[0:-1]
    yplaces = places * yrange + m.ymin
    px = urx_pixels + (liquefaction_ticks_x_offset * fig.dpi)
    xplace,yt = ax.transData.inverted().transform((px,100))
    levels = ['Low','Medium','High','V. High']
    for i in range(0,len(yplaces)):
        pyplot.text(xplace,yplaces[i],levels[i],verticalalignment='center')

    #draw the landslide probability colorbar annotation
    px = ulx_pixels - (landslide_annotation_x_offset * fig.dpi)
    xplace,yt = ax.transData.inverted().transform((px,100))
    for i in range(0,len(yplaces)):
        pyplot.text(xplace,yplaces[i],levels[i],verticalalignment='center',horizontalalignment='right')
    
    #draw the tick labels on the left side of the landslide colorbar
    ctickfrac = (numpy.array(cticks) - 2)/(max(cticks) - 2)
    yplaces = ctickfrac * yrange + m.ymin
    px = ulx_pixels - (landslide_annotation_x_offset * fig.dpi)
    xplace,yt = ax.transData.inverted().transform((px,100))
    cticklabels = [str(int(tick))+'%' for tick in cticks]
    for i in range(0,len(yplaces)):
        pyplot.text(xplace,yplaces[i],cticklabels[i],verticalalignment='center',horizontalalignment='right')
    
    #draw the landslide colorbar
    px = llx_pixels - (landslide_colorbar_axis_x_offset * fig.dpi)
    axisleft,axisbottom = fig.transFigure.inverted().transform((px,lly_pixels)) #now in figure unit space
    
    cax1 = fig.add_axes([axisleft,axisbottom,colorbar_width_figure,colorbar_height_figure])
    cb1 = pyplot.colorbar(mappable=lsprobhandle,cax=cax1,ticks=[])

    #draw the tick marks on the landslide colorbar manually
    px = ulx_pixels - (landslide_colorbar_axis_x_offset * fig.dpi) #left edge of the tickmark in pixels
    px2 = px + (landslide_colorbar_tick_width * fig.dpi) #right edge of the tickmark in pixels
    xplace,yt = cax1.transData.inverted().transform((px,100)) #left edge of the tickmark in axis data space
    xplace2,yt = cax1.transData.inverted().transform((px2,100)) #right edge of the tickmark in axis data space
    for i in range(0,len(yplaces)):
        yplace = ctickfrac[i]        
        pyplot.plot([xplace,xplace2],[yplace,yplace],'k-',axes=cax1,lw=0.5)

    #draw the liquefaction colorbar
    px = lrx_pixels + (liquefaction_colorbar_axis_x_offset * fig.dpi)
    axisleft,axisbottom = fig.transFigure.inverted().transform((px,lry_pixels))
    cax2 = fig.add_axes([axisleft,axisbottom,colorbar_width_figure,colorbar_height_figure])
    cticks = [2.0,5.0,10.0,15.0,20.0]
    cb2 = pyplot.colorbar(mappable=probhandle,cax=cax2)
    cticklabels = [str(int(tick))+'%' for tick in cticks]
    cb2.set_ticks(cticks)
    cb2.set_ticklabels(cticklabels)
    
    #make a legend at the bottom with the secondary hazard statistics
    if legdict is not None:
        lsarea = legdict['lssum']/1e6
        if lsarea < 0.1:
            lsarea = '< 0.1'
        else:
            lsarea = commify(int(lsarea)) + ' km$^2$'
        lsncells = commify(legdict['lscells'])
        lspmax = '%.1f%%' % (legdict['lsmax']*100)
        lspmean = '%.1f%%' % (legdict['lsmean']*100)
        lspslope = '%.1f$\degree$' % legdict['slopemin']
        tcells = [['Landslide Area',lsarea],
                  ['# of cells',lsncells],
                  ['Maximum probability',lspmax],
                  ['Mean probability',lspmean],
                  ['Minimum slope',lspslope]]
        ax_lstable = fig.add_axes([0.215,0.08,0.2,0.1])
        pyplot.axis('off')
        lstable = pyplot.table(cellText=tcells,loc='lower left')

        lqarea = legdict['liqsum']/1e6
        if lqarea < 0.1:
            lqarea = '< 0.1'
        else:
            lqarea = commify(int(lqarea)) + ' km$^2$'
        lqncells = commify(legdict['liqcells'])
        lqpmax = '%.1f%%' % (legdict['liqmax']*100)
        lqpmean = '%.1f%%' % (legdict['liqmean']*100)
        lqpslope = '%.1f$\degree$' % legdict['slopemax']
        tcells = [['Liquefaction Area',lqarea],
                  ['# of cells',lqncells],
                  ['Maximum probability',lqpmax],
                  ['Mean probability',lqpmean],
                  ['Maximum slope',lqpslope]]
        ax_lqtable = fig.add_axes([0.57,0.08,0.2,0.1])
        pyplot.axis('off')
        lqtable = pyplot.table(cellText=tcells,loc='lower left')
    
    if title is not None:
        ax.set_title(title)
    if eventid is not None:
        pyplot.savefig(eventid+'.pdf')
    else:
        pyplot.savefig('output.pdf')

def getMapLines(dmin,dmax):
    NLINES = 4
    drange = dmax-dmin
    if drange > 4:
        near = 1
    else:
        if drange >= 0.5:
            near = 0.25
        else:
            near = 0.125
    inc = roundToNearest(drange/NLINES,near)
    if inc == 0:
        near = pow(10,round(math.log10(drange))) #make the increment the closest power of 10
        inc = ceilToNearest(drange/NLINES,near)
        newdmin = floorToNearest(dmin,near)
        newdmax = ceilToNearest(dmax,near)
    else:
        newdmin = ceilToNearest(dmin,near)
        newdmax = floorToNearest(dmax,near)
    darray = numpy.arange(newdmin,newdmax+inc,inc)
    if darray[-1] > dmax:
        darray = darray[0:-1]
    return darray

def slide(pgagrid,slopegrid,cohesiongrid,coeff,slopemin=3.0):
    pga = pgagrid.griddata
    slope = slopegrid.griddata
    slopemin = slopemin * 100
    cohesion = cohesiongrid.griddata/10.0
    pgaslope = pgagrid.griddata * slopegrid.griddata
    x = coeff['b0'] + coeff['bpga']*pga + coeff['bslope']*slope + coeff['bcohesion']*cohesion + coeff['bpgaslope']*pgaslope
    P = 1 / (1 + numpy.exp(-1*x))
    idx = (slopegrid.griddata < slopemin).nonzero()
    P[idx] = 0.0
    nrows,ncols = pgagrid.griddata.shape
    xdim = pgagrid.geodict['xdim']
    ydim = pgagrid.geodict['ydim']
    cy = pgagrid.geodict['ymax'] - (nrows/2.0)*ydim
    xdim_meters = xdim * DEG2M * cosd(cy)
    ydim_meters = ydim * DEG2M
    pravel = P.ravel()
    indices = find(pravel > 0.01)
    psum = 0
    pmax = 0
    pmeansum = numpy.nansum(pravel[indices])
    if numpy.isnan(pmeansum):
        pmean = 0
    else:
        pmean = pmeansum/len(indices)
    ncells = len(indices)
    if len(indices):
        for idx in indices:
            psum += pravel[idx] * xdim_meters * ydim_meters
            if pravel[idx] > pmax:
                pmax = pravel[idx]
        
    lsprob = Grid()
    lsprob.geodict = copy.deepcopy(pgagrid.geodict)
    lsprob.griddata = numpy.copy(P)
    return (lsprob,psum,ncells,pmax,pmean)
    

def liquefy(pgagrid,vs30grid,ctigrid,slopegrid,coeff,slopemax=3):
    pga = numpy.log(pgagrid.griddata/100.0)
    svel = numpy.log(vs30grid.griddata)
    cti = ctigrid.griddata/100.0
    x = coeff['b0'] + coeff['bpga']*pga + coeff['bcti']*cti + coeff['bvs30']*svel
    P = 1 / (1 + numpy.exp(-1*x))
    #slopemax = 1000.0 * numpy.tan(slopemax*(numpy.pi/180.0))
    slopemax = slopemax * 100
    idx = (slopegrid.griddata > slopemax).nonzero()
    P[idx] = 0.0
    nrows,ncols = pgagrid.griddata.shape
    xdim = pgagrid.geodict['xdim']
    ydim = pgagrid.geodict['ydim']
    cy = pgagrid.geodict['ymax'] - (nrows/2.0)*ydim
    xdim_meters = xdim * DEG2M * cosd(cy)
    ydim_meters = ydim * DEG2M
    pravel = P.ravel()
    indices = find(pravel > 0.01)
    psum = 0
    pmax = 0
    pmeansum = numpy.nansum(pravel[indices])
    if numpy.isnan(pmeansum):
        pmean = 0
    else:
        pmean = pmeansum/len(indices)
    ncells = len(indices)
    if len(indices):
        for idx in indices:
            psum += pravel[idx] * xdim_meters * ydim_meters
            if pravel[idx] > pmax:
                pmax = pravel[idx]
        
    lqprob = Grid()
    lqprob.geodict = copy.deepcopy(pgagrid.geodict)
    lqprob.griddata = numpy.copy(P)
    return (lqprob,psum,ncells,pmax,pmean)

def getcti(ctifolder,bounds):
    ctifiles = glob.glob(ctifolder+'*_latlon.bil')
    myfile = None
    for ctifile in ctifiles:
        ctigrid = EsriGrid(ctifile)
        hdr = ctigrid.getHeader()
        xmin = hdr['ulxmap']
        xmax = hdr['ulxmap'] + hdr['ncols']*hdr['xdim']
        ymax = hdr['ulymap']
        ymin = hdr['ulymap'] - hdr['nrows']*hdr['ydim']
        if bounds[0] <= xmax and bounds[1] >= xmin and bounds[2] <= ymax and bounds[3] > ymin:
            ctigrid.load(bounds=bounds)
            if numpy.isnan(numpy.nansum(ctigrid.griddata)) or numpy.nansum(ctigrid.griddata) == 0:
                continue
            else:
                myfile = ctifile
                break
    return myfile

if __name__ == '__main__':
    usage = """usage: %prog [options] grid.xml
    Run the landslide and liquefaction models defined by coefficients found in a config.ini file.
    Output is a pdf map with liquefaction/landslide results layered on topography."""
    parser = OptionParser(usage=usage)
    parser.add_option("-z", "--zoom", dest="zoomCoordinates",
                      help='zoom in map to coordinates', 
                      metavar='"xmin xmax ymin ymax"')
    parser.add_option("-r","--roads", dest="drawRoads",
                      action="store_true",default=False,help="Draw roads found inside map")
    parser.add_option("-t","--table", dest="drawTable",
                      action="store_true",default=False,help="Add table of summary statistics to figure")
    # parser.add_option("-q", "--quiet",
    #                   action="store_false", dest="verbose", default=True,
    #                   help="don't print status messages to stdout")
    (options, args) = parser.parse_args()
    #get input from the command line
    if not len(args):
        print 'Missing input shakemap grid.xml file.'
        parser.print_help()
        sys.exit(1)
       
    shakefile = args[0]
    print 'Processing %s' % shakefile

    #get the bounds, if present
    xmin = None
    xmax = None
    ymin = None
    ymax = None
    pgagrid = ShakeGrid(shakefile,variable='PGA')
    isScenario = pgagrid.getAttributes()['shakemap_grid']['shakemap_event_type'].lower() == 'scenario'
    if options.zoomCoordinates is not None:
        boundstr = options.zoomCoordinates
        bparts = boundstr.split()
        xmin = float(bparts[0])
        xmax = float(bparts[1])
        ymin = float(bparts[2])
        ymax = float(bparts[3])
        pdict = pgagrid.getGeoDict()
        pxmin = pdict['xmin']
        pxmax = pdict['xmax']
        pymin = pdict['ymin']
        pymax = pdict['ymax']
        if xmin < pxmin or xmax > pxmax or ymin < pymin or ymax > pymax:
            fmt = 'Your bounds are outside the bounds of the ShakeMap (%.3f %.3f %.3f %.3f)'
            print fmt % (pxmin,pxmax,pymin,pymax)
            sys.exit(1)
        
    #find the config file
    homedir = os.path.dirname(os.path.abspath(__file__)) #where is this script?
    configfile = os.path.join(homedir,CONFIGFILE)
    if not os.path.isfile(configfile):
        print 'Missing config file %s.  Exiting.' % configfile
        sys.exit(1)
        
    config = ConfigParser.ConfigParser()
    config.readfp(open(configfile))
    configopts = {}
    for opt in LQCOEFF:
        configopts[opt] = float(config.get('LIQUEFACTION',opt))
    ctifile = config.get('LIQUEFACTION','ctifile')
    topofile = config.get('LIQUEFACTION','topofile')
    slopefile = config.get('LANDSLIDE','maxslopefile')
    slopemin = float(config.get('LANDSLIDE','slopemin'))
    slopemax = float(config.get('LIQUEFACTION','slopemax'))

    #if the user wants to draw roads, get the name of the roads shapefile
    if options.drawRoads:
        roadshapefile,fext = os.path.splitext(config.get('LAYERS','roads'))
    else:
        roadshapefile = None

    coastshapefile = config.get('LAYERS','coasts')
    riverfolder = config.get('LAYERS','riverfolder')
    borderfile = config.get('LAYERS','borders')
        
    try:
        vs30grid = ShakeGrid(shakefile,variable='SVEL')
    except Exception,msg:
        print 'Could not process %s (missing Vs30 layer).'
        sys.exit(1)

    if xmin is not None:
        bounds = (xmin,xmax,ymin,ymax)
    else:
        bounds = pgagrid.getRange()

    bigbounds = (bounds[0]-1.0,bounds[1]+1.0,bounds[2]-1.0,bounds[3]+1.0)
    ctigrid = EsriGrid(ctifile)
    ctigrid.load(bounds=bounds)
    topogrid = EsriGrid(topofile)
    topogrid.load(bounds=bigbounds)
    slopegrid = EsriGrid(slopefile)
    slopegrid.load(bounds=bigbounds)
    slopegrid.interpolateToGrid(ctigrid.getGeoDict())
    pgagrid.interpolateToGrid(ctigrid.getGeoDict())
    vs30grid.interpolateToGrid(ctigrid.getGeoDict())
    topogrid.interpolateToGrid(ctigrid.getGeoDict())
    liqmap,psum,ncells,pmax,pmean = liquefy(pgagrid,vs30grid,ctigrid,slopegrid,configopts,slopemax=slopemax)

    #run the landslide model
    cohesionfile = config.get('LANDSLIDE','cohesionfile')
    configopts = {}
    for opt in LSCOEFF:
        configopts[opt] = float(config.get('LANDSLIDE',opt))

    slopegrid = EsriGrid(slopefile)
    slopegrid.load(bounds=bounds)
    cohesiongrid = EsriGrid(cohesionfile)
    cohesiongrid.load(bounds=bounds)
    pgagrid = ShakeGrid(shakefile,variable='PGA')
    pgagrid.interpolateToGrid(slopegrid.getGeoDict())
    lsmap,lspsum,lsncells,lspmax,lspmean = slide(pgagrid,slopegrid,cohesiongrid,configopts,slopemin=slopemin)
    
    shakeheader = pgagrid.getAttributes()
    mag = shakeheader['event']['magnitude']
    lat = shakeheader['event']['lat']
    lon = shakeheader['event']['lon']
    timestr = shakeheader['event']['event_timestamp'].strftime('%b %d %Y')
    location = shakeheader['event']['event_description']
    eventid = shakeheader['shakemap_grid']['shakemap_originator'] + shakeheader['shakemap_grid']['shakemap_id']
    if isScenario:
        title = location
    else:
        title = 'M%.1f %s\n %s' % (mag,timestr,location)
    xdim = shakeheader['grid_specification']['nominal_lon_spacing']
    ydim = shakeheader['grid_specification']['nominal_lat_spacing']
    if options.drawTable:
        legdict = {'liqsum':psum,'liqcells':ncells,
                   'liqmax':pmax,'shakexdim':xdim,
                   'shakeydim':ydim,'liqmean':pmean,
                   'lssum':lspsum,'lscells':lsncells,
                   'lsmax':lspmax,'lsmean':lspmean,
                   'slopemin':slopemin,'slopemax':slopemax}
    else:
        legdict = None
    #need to add arbitrary shapefile drawing into this somehow
    shapesdict = getShapes(config)
    lqtitle = 'M%.1f %s\n %s - Liquefaction Probability' % (mag,timestr,location)
    lstitle = 'M%.1f %s\n %s - Landslide Probability' % (mag,timestr,location)
    outroot = config.get('OUTPUT','folder')
    eventfolder = os.path.join(outroot,eventid)
    if not os.path.isdir(eventfolder):
        os.makedirs(eventfolder)
    makeMatMap(topogrid,liqmap,lsmap,coastshapefile,riverfolder,
               isScenario=isScenario,roads=roadshapefile,
               shapedict=shapesdict,title=title,borderfile=borderfile,
               legdict=legdict,eventfolder=eventfolder,epicenter=(lat,lon))
    liqimage = saveTiff(liqmap,os.path.join(eventfolder,'liquefaction.tif'),isFloat=True)
    lsimage = saveTiff(lsmap,os.path.join(eventfolder,'landslide.tif'),isFloat=True)
    topoimage = saveTiff(topogrid,os.path.join(eventfolder,'topography.tif'),isFloat=True)
    pgaimage = saveTiff(pgagrid,os.path.join(eventfolder,'pga.tif'),isFloat=True)
    cohesionimage = saveTiff(cohesiongrid,os.path.join(eventfolder,'cohesion.tif'),isFloat=True)
    ctiimage = saveTiff(ctigrid,os.path.join(eventfolder,'cti.tif'),isFloat=True)
    vs30image = saveTiff(vs30grid,os.path.join(eventfolder,'vs30.tif'),isFloat=True)
    slopeimage = saveTiff(slopegrid,os.path.join(eventfolder,'slope.tif'),isFloat=True)
    files = os.listdir(eventfolder)
    print '%i files saved to %s' % (len(files)-2,eventfolder)
    
