#!/usr/bin/env python

#stdlib imports
import ConfigParser
import os.path
import sys
import glob
import copy
import math
import warnings

#turn off all warnings...
warnings.filterwarnings('ignore')

#third party imports
import numpy
from matplotlib import pyplot
from mpl_toolkits.basemap import Basemap, shiftgrid
from matplotlib import cm
from matplotlib.colors import LightSource
import gdal
import osr
from pylab import find

#local imports
from losspager.io.esri import EsriGrid
from losspager.io.shake import ShakeGrid
from losspager.io.grid import Grid
from losspager.io.gmt import GMTGrid
from losspager.util.text import decToRoman,roundToNearest,ceilToNearest,floorToNearest,commify

CONFIGFILE = 'config.ini'
LQCOEFF = ['b0','bpga','bcti','bvs30']
LSCOEFF = ['b0','bpga','bslope','bcohesion','bpgaslope']
DEG2M = 111191

def cosd(angle):
    """
    Return cosine of angle expressed in degrees.
    @param angle: Input angle in degrees.
    @return: Cosine of input angle.
    """
    return math.cos(angle*(math.pi/180))

def saveTiff(lqgrid,fname):
    format = "GTiff"
    driver = gdal.GetDriverByName( format )
    metadata = driver.GetMetadata()
    nrows,ncols = lqgrid.griddata.shape
    outdata = driver.Create(fname, ncols, nrows, 1, gdal.GDT_Byte )
    ulxmap = lqgrid.geodict['xmin']
    ulymap = lqgrid.geodict['ymax']
    xdim = lqgrid.geodict['xdim']
    ydim = lqgrid.geodict['ydim']
    outdata.SetGeoTransform([ulxmap,xdim,0,ulymap,0,ydim])
    srs = osr.SpatialReference()
    srs.SetGeogCS( "My geographic coordinate system",
               "WGS_1984", 
               "My WGS84 Spheroid", 
               osr.SRS_WGS84_SEMIMAJOR,
               osr.SRS_WGS84_INVFLATTENING, 
               "Greenwich", 0.0,
        "degree", float(osr.SRS_UA_DEGREE_CONV))
    outdata.SetProjection(srs.ExportToWkt())
    raster = numpy.array(numpy.round(lqgrid.griddata*100),dtype=numpy.int8)
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

def makeMap(topogrid,lqgrid,lsgrid,title=None,eventid=None,legdict=None):
    #topo and liquefaction grids are already resampled to same extent and resolution
    bounds = topogrid.getRange()
    xmin,xmax,ymin,ymax = bounds
    cx = xmin + (xmax-xmin)/2.0
    cy = ymin + (ymax-ymin)/2.0
    fig = pyplot.figure(figsize=(11,8.5))
    ax = fig.add_axes([0.27,0.21,0.45,0.59])
    # setup of basemap ('lcc' = lambert conformal conic).
    # use major and minor sphere radii from WGS84 ellipsoid.
    m = Basemap(llcrnrlon=xmin,llcrnrlat=ymin,urcrnrlon=xmax,urcrnrlat=ymax,
                rsphere=(6378137.00,6356752.3142),suppress_ticks=False,
        resolution='h',projection='tmerc',lat_0=cy,lon_0=cx,ax=ax)

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

    #do shaded relief stuff
    ls = LightSource(azdeg = 90, altdeg = 20)
    #pyplot.set_cmap('bone')
    rgb = ls.shade(topodatm,cmap=palette)
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
    probhandle = m.imshow(lqdatm,cmap=palette,vmin=2.0,vmax=20.0)

    #draw the landslide probabilities (anything above 1%) in an orange-red scale (OrRd)
    lsdat = numpy.flipud(lsgrid.getData().copy()) * 100.0
    clear_color = [0,0,0,0.0]
    palette = cm.cool
    i = numpy.where(lsdat < 2.0)
    lsdat[i] = 0
    lsdatm = numpy.ma.masked_equal(lsdat, 0)
    palette.set_bad(clear_color,alpha=0.0)
    lsprobhandle = m.imshow(lsdatm,cmap=palette,vmin=2.0,vmax=20.0)

    meridians = getMapLines(bounds[0],bounds[1])
    parallels = getMapLines(bounds[2],bounds[3])
    
    xmap_range = m.xmax-m.xmin
    ymap_range = m.ymax-m.ymin
    xoff = -0.09*(xmap_range)
    yoff = -0.04*(ymap_range)

    #do tick stuff
    fig.canvas.draw()
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
    
    #im = m.imshow(topodat,cm.GMT_haxby)
    #m.drawcoastlines()
    border_color = [0.18,1.0,0.18]
    m.drawcountries(color=border_color,linewidth=2.0)
    m.drawstates(color=border_color,linewidth=2.0)
    m.drawrivers(color=water_color,linewidth=1.0)
    m.fillcontinents(color=[0,0,0,0],lake_color=water_color)

    #custom liquefaction colorbar
    
    # cticks = [2.0,5.0,10.0,15.0,20.0]
    # cb = m.colorbar(probhandle,"right", size="5%", pad='15%')
    # cticklabels = [str(int(tick))+'%' for tick in cticks]
    # cb.set_ticks(cticks)
    # cb.set_ticklabels(cticklabels)
    # xrange = m.xmax - m.xmin
    # yrange = m.ymax - m.ymin
    # cbartitle = 'Liquefaction\nProbability'
    # pyplot.text(m.xmax+(xrange/20),m.ymax+(yrange/25),cbartitle,multialignment='left')

    xrange = m.xmax - m.xmin
    yrange = m.ymax - m.ymin
    cbartitle = 'Landslide\nProbability'
    pyplot.text(m.xmin-(xrange/5),m.ymax+(yrange/25),cbartitle,multialignment='left',axes=ax)

    cbartitle = 'Liquefaction\nProbability'
    pyplot.text(m.xmax+(xrange/20),m.ymax+(yrange/25),cbartitle,multialignment='left',axes=ax)

    cticks = [2.0,5.0,10.0,15.0,20.0]
    #write text to the right of the liquefaction colorbar
    ctickfrac = (numpy.array(cticks) - 2)/(max(cticks) - 2)
    places = numpy.diff(ctickfrac)/2.0 + ctickfrac[0:-1]
    yplaces = places * yrange + m.ymin
    xplace = m.xmax + (m.xmax - m.xmin)/7
    levels = ['Low','Medium','High','V. High']
    for i in range(0,len(yplaces)):
        pyplot.text(xplace,yplaces[i],levels[i],verticalalignment='center')

    #write text to the left of the landslide colorbar
    xplace = m.xmin-(xrange/6.75)
    for i in range(0,len(yplaces)):
        pyplot.text(xplace,yplaces[i],levels[i],verticalalignment='center',horizontalalignment='right')
    
    #custom landslide colorbar
    ctickfrac = (numpy.array(cticks) - 2)/(max(cticks) - 2)
    yplaces = ctickfrac * yrange + m.ymin
    xplace = m.xmin-(xrange/7)
    cticklabels = [str(int(tick))+'%' for tick in cticks]
    for i in range(0,len(yplaces)):
        pyplot.text(xplace,yplaces[i],cticklabels[i],verticalalignment='center',horizontalalignment='right')
    cax1 = fig.add_axes([0.215,0.22,0.02,0.58])
    cb1 = pyplot.colorbar(mappable=lsprobhandle,cax=cax1,ticks=[])
    

    #custom liquefaction colorbar
    cax2 = fig.add_axes([0.755,0.22,0.02,0.58])
    cticks = [2.0,5.0,10.0,15.0,20.0]
    cb2 = pyplot.colorbar(mappable=probhandle,cax=cax2)
    cticklabels = [str(int(tick))+'%' for tick in cticks]
    cb2.set_ticks(cticks)
    cb2.set_ticklabels(cticklabels)

    
    
    #make a legend at the bottom with the secondary hazard statistics
    if legdict is not None:
        legfmt = 'Total area of predicted liquefaction: %s'
        liqsum = numpy.round(legdict['liqsum']/1e6)
        liqcells = legdict['liqcells']
        pmax = legdict['liqmax']
        pmean = legdict['liqmean']
        if liqsum < 0.1:
            legtext = legfmt % '< 0.1 km$^2$'
        else:
            legtext = legfmt % '%s km$^2$' % commify(int(liqsum))
        line2 = 'Total number of liquefaction cells: %s' % commify(liqcells)
        legtext = legtext + '\n' + line2
        line3 = 'Highest liquefaction probability: %i%%' % int(pmax*100)
        legtext = legtext + '\n' + line3
        line4 = 'Mean liquefaction probability: %.1f%%' % (pmean*100)
        legtext = legtext + '\n' + line4
        xdim = lqgrid.geodict['xdim']*111.191*cosd(cy)
        ydim = lqgrid.geodict['ydim']*111.191
        line5 = 'Cell dimensions: %.1f km by %.1f km' % (xdim,ydim)
        legtext = legtext + '\n' + line5
        legx = m.xmin
        legy = m.ymin - yrange/15
        pyplot.text(legx,legy,legtext,multialignment='left',verticalalignment='top')
    
    
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
    slopemax = 1000.0 * numpy.tan(slopemax*(numpy.pi/180.0))
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
    #get input from the command line
    shakefile = sys.argv[1]
    print 'Processing %s' % shakefile
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
    ctifolder = config.get('LIQUEFACTION','ctifolder')
    topofile = config.get('LIQUEFACTION','topofile')
    slopefile = config.get('LANDSLIDE','slopefile')
    slopemin = float(config.get('LANDSLIDE','slopemin'))
    slopemax = float(config.get('LIQUEFACTION','slopemax'))
    pgagrid = ShakeGrid(shakefile,variable='PGA')
    try:
        vs30grid = ShakeGrid(shakefile,variable='SVEL')
    except Exception,msg:
        print 'Could not process %s (missing Vs30 layer).'
        sys.exit(1)
    bounds = pgagrid.getRange()
    ctifile = getcti(ctifolder,bounds)
    if ctifile is None:
        print 'Could not find a file containing this shakemap.'
        sys.exit(1)

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
    timestr = shakeheader['event']['event_timestamp'].strftime('%b %d %Y')
    location = shakeheader['event']['event_description']
    eventid = shakeheader['shakemap_grid']['shakemap_originator'] + shakeheader['shakemap_grid']['shakemap_id']
    title = 'M%.1f %s, %s' % (mag,timestr,location)
    xdim = shakeheader['grid_specification']['nominal_lon_spacing']
    ydim = shakeheader['grid_specification']['nominal_lat_spacing']
    legdict = {'liqsum':psum,'liqcells':ncells,
               'liqmax':pmax,'shakexdim':xdim,
               'shakeydim':ydim,'liqmean':pmean}
    mapfile = makeMap(topogrid,liqmap,lsmap,title=title,eventid=eventid,legdict=legdict)
    tiffimage = saveTiff(liqmap,'output.tif')
            
    





