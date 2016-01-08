#!/usr/bin/env python

#stdlib imports
import copy
import os.path

#third party imports
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import LightSource
from matplotlib import cm
import matplotlib as mpl
import fiona
from shapely.geometry import shape,Polygon,LineString

#local imports
from neicio.gmt import GMTGrid
from neicutil.text import ceilToNearest,floorToNearest,roundToNearest,commify
import model
from neicmap.city import PagerCity

ALPHA = 0.7
AZDEFAULT=90
ALTDEFAULT=20

ROADCOLOR = '#6E6E6E'
COUNTRYCOLOR = '#177F10'

SEA_LEVEL = 0

def getGridExtent(grid,basemap):
    lonmin,lonmax,latmin,latmax = grid.getRange()
    xmin,ymin = basemap(lonmin,latmin)
    xmax,ymax = basemap(lonmax,latmax)
    extent = (xmin,xmax,ymin,ymax)
    return extent

def getTopoRGB(topogrid):
    topotmp = topogrid.getData().copy()
    #make a masked array
    topotmp = np.ma.array(topotmp)
    topodat = np.ma.masked_where(np.isnan(topotmp),topotmp)
    cy = topogrid.geodict['ymin'] + (topogrid.geodict['ymax'] - topogrid.geodict['ymin'])/2.0
    #flag the regions where topography is less than 0 (we'll color this ocean later)
    i = np.where(topodat == SEA_LEVEL)

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
    ls = LightSource(azdeg = AZDEFAULT, altdeg = ALTDEFAULT)

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

    rgb = np.flipud(rgb)

    return (rgb,palette1)


def hillshade(topogrid, azimuth, angle_altitude):
    """
    Most of this script borrowed from http://geoexamples.blogspot.com/2014/03/shaded-relief-images-using-gdal-python.html last accessed 9/2/2015
    """
    topotmp = topogrid.getData().copy()
    #make a masked array
    topotmp = np.ma.array(topotmp)
    topodat = np.ma.masked_where(np.isnan(topotmp), topotmp)
    topodat = np.ma.masked_where(topodat == SEA_LEVEL, topodat)
    x, y = np.gradient(topodat)
    slope = np.pi/2. - np.arctan(np.sqrt(x*x + y*y))
    aspect = np.arctan2(-x, y)
    azimuthrad = azimuth*np.pi / 180.
    altituderad = angle_altitude*np.pi / 180.
    shaded = np.sin(altituderad) * np.sin(slope) + np.cos(altituderad) * np.cos(slope) * np.cos(azimuthrad - aspect)
    return 255*(shaded + 1)/2


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
        near = np.power(10,round(math.log10(drange))) #make the increment the closest power of 10
        inc = ceilToNearest(drange/NLINES,near)
        newdmin = floorToNearest(dmin,near)
        newdmax = ceilToNearest(dmax,near)
    else:
        newdmin = ceilToNearest(dmin,near)
        newdmax = floorToNearest(dmax,near)
    darray = np.arange(newdmin,newdmax+inc,inc)
    if darray[-1] > dmax:
        darray = darray[0:-1]
    return darray

#here we are always using 2 decimal places because it makes it easier to 
#place the meridian/parallel labels in a consistent spot
#earlier versions of the code the strings had variable precision
#but were outside the map, so positioning didn't matter.
def latstr(parallel):
    if parallel < 0:
        parstr = '%.2f' % (-1*parallel) + '$\degree$ S'
    else:
        parstr = '%.2f' % (parallel) + '$\degree$ N'
    return parstr

def lonstr(meridian):
    
    if meridian < 0:
        merstr = '%.2f' % (-1*meridian) + '$\degree$ W'
    else:
        merstr = '%.2f' % (meridian) + '$\degree$ E'
    return merstr

def getMapTicks(m,xmin,xmax,ymin,ymax):
    meridians = getMapLines(xmin,xmax)
    parallels = getMapLines(ymin,ymax)
    #do tick stuff
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
    
    return (xticks,xlabels,yticks,ylabels)

def renderPanel(logmodel,colormaps,outfolder,edict):
    nparray = "<type 'numpy.ndarray'>"
    #first, figure out how many layers we have
    layerdict = logmodel.layerdict
    outfiles = []
    for smterm in model.SM_TERMS:
        for term in logmodel.terms.values():
            if term.find(smterm) > -1 and not isinstance(logmodel.shakedict[smterm],float):
                layerdict[smterm] = logmodel.shakedict[smterm]

    for layername,layergrid in layerdict.iteritems():
        fig = plt.figure(figsize=(8,8))
        ax = plt.gca()
        renderLayer(layername,layergrid,outfolder,edict,fig,ax,logmodel.model,colormaps)
        outfile = os.path.join(outfolder,'%s_%s.pdf' % (layername,logmodel.model))
        print 'Saving input layer %s to %s' % (layername,outfile)
        plt.savefig(outfile)
        outfiles.append(outfile)

    outfile = os.path.join(outfolder,'%s_model.pdf' % logmodel.model)
    fig = plt.figure(figsize=(8,8))
    ax = plt.gca()
    P = logmodel.calculate()
    pgrid = GMTGrid()
    pgrid.griddata = P.copy()
    pgrid.geodict = layergrid.geodict.copy()
    renderLayer(logmodel.model,pgrid,outfolder,edict,fig,ax,logmodel.model,colormaps)
    print 'Saving %s model to %s' % (logmodel.model,outfile)
    outfiles.append(outfile)
    return outfiles

def renderLayer(layername,layergrid,outfolder,edict,fig,ax,model,colormaps):
    try:
        xmin,xmax,ymin,ymax = layergrid.getRange()
    except:
        pass
    clat = ymin + (ymax-ymin)/2.0
    clon = xmin + (xmax-xmin)/2.0
    m = Basemap(llcrnrlon=xmin,llcrnrlat=ymin,urcrnrlon=xmax,urcrnrlat=ymax,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='l',area_thresh=1000.,projection='lcc',\
                lat_1=clat,lon_0=clon,ax=ax)
    water_color = [.47,.60,.81]
    m.drawmapboundary(fill_color=water_color)

    #render layer
    lsdat = layergrid.getData()
    #extent = getGridExtent(layergrid,m)
    cmap = 'jet'
    if colormaps.has_key(model):
        mdict = colormaps[model]
        if mdict.has_key(layername):
            cmap = mdict[layername]
    lons = np.arange(xmin,xmax,layergrid.getGeoDict()['xdim'])
    lons = lons[:layergrid.getGeoDict()['ncols']] #make sure it's the right length, sometimes it is one too long, sometimes it isn't because of GMTGrid issue
    lats = np.arange(ymax, ymin, -layergrid.getGeoDict()['ydim']) #backwards so it plots right side up
    lats = lats[:layergrid.getGeoDict()['nrows']] #make sure it's the right length
    #make meshgrid
    llons, llats = np.meshgrid(lons, lats)
    x, y = m(llons, llats)  #get projection coordinates
    lsprobhandle = m.pcolormesh(x, y, lsdat, lw=0, cmap=cm.get_cmap(cmap),rasterized=True)
    #lsprobhandle = plt.imshow(lsdat,origin='upper',extent=extent,cmap=cm.get_cmap(cmap))
    plt.colorbar(lsprobhandle)
    
    #this business apparently has to happen after something has been 
    #rendered on the map, which I guess makes sense.
    #draw the map ticks on outside of all edges
    fig.canvas.draw() #have to do this for tick stuff to show
    xticks,xlabels,yticks,ylabels = getMapTicks(m,xmin,xmax,ymin,ymax)
    plt.sca(ax)
    plt.tick_params(axis='both',direction='in',right='on',colors='white')
    plt.xticks(xticks,xlabels,size=6)
    plt.yticks(yticks,ylabels,size=6)
    for tick in ax.axes.yaxis.get_major_ticks():
        tick.set_pad(-33)
        tick.label2.set_horizontalalignment('right')
    for tick in ax.axes.xaxis.get_major_ticks():
        tick.set_pad(-10)
        tick.label2.set_verticalalignment('top')
    [i.set_color("white") for i in plt.gca().get_xticklabels()]
    [i.set_color("white") for i in plt.gca().get_yticklabels()]
    plt.title('%s' % (layername))
    
def makeDualMap(lqgrid,lsgrid,topogrid,slopegrid,eventdict,outfolder,isScenario=False,roadslist=[],colors={},cityfile=None):
    # create the figure and axes instances.
    fig = plt.figure()
    ax = fig.add_axes([0.1,0.1,0.8,0.8])
    # setup of basemap ('lcc' = lambert conformal conic).
    # use major and minor sphere radii from WGS84 ellipsoid.
    xmin,xmax,ymin,ymax = topogrid.getRange()
    clat = ymin + (ymax-ymin)/2.0
    clon = xmin + (xmax-xmin)/2.0
    m = Basemap(llcrnrlon=xmin,llcrnrlat=ymin,urcrnrlon=xmax,urcrnrlat=ymax,\
                rsphere=(6378137.00,6356752.3142),\
                resolution='h',area_thresh=1000.,projection='lcc',\
                lat_1=clat,lon_0=clon,ax=ax)

    clear_color = [0,0,0,0.0]
    
    #topoextent = getGridExtent(topogrid,m)
    #rgb,topopalette = getTopoRGB(topogrid)
    #color_tuple = rgb.transpose((1,0,2)).reshape((rgb.shape[0]*rgb.shape[1],rgb.shape[2]))/255.0

    hillsh = hillshade(topogrid, 315, 50) # divide by 2 to mute colors
    topogrid2 = GMTGrid()
    topogrid2.loadFromGrid(topogrid)
    topogrid2.interpolateToGrid(lqgrid.geodict)
    #define what cells are below sea level
    iwater = np.where(topogrid2.griddata == SEA_LEVEL)

    lons = np.arange(xmin,xmax,topogrid.getGeoDict()['xdim'])
    lons = lons[:topogrid.getGeoDict()['ncols']] #make sure right length
    lats = np.arange(ymax, ymin, -topogrid.getGeoDict()['ydim'])  # backwards so it plots right
    lats = lats[:topogrid.getGeoDict()['nrows']]
    llons, llats = np.meshgrid(lons, lats) #make meshgrid
    x, y = m(llons, llats)  #get projection coordinates
    im = m.pcolormesh(x, y, hillsh, cmap='Greys', lw=0, vmin=0.0, vmax=550.,rasterized=True)
    #im = m.imshow(rgb,cmap=topopalette,extent=topoextent)

    #figure out the aspect ratio of the axes
    print 'Map width: %s' % commify(int(m.xmax-m.xmin))
    print 'Map height: %s' % commify(int(m.ymax-m.ymin))
    aspect = (m.xmax-m.xmin)/(m.ymax-m.ymin)
    slope = -0.33
    intercept = 0.463
    leftx = slope*aspect + intercept
    print 'Left edge of left colorbar is at: %.2f' % leftx
    slope = 20.833
    intercept = -0.833
    titlex = slope*aspect + intercept
    print 'Title X of right colorbar is at: %.2f' % titlex
    
    #this business apparently has to happen after something has been 
    #rendered on the map, which I guess makes sense.
    #draw the map ticks on outside of all edges
    fig.canvas.draw() #have to do this for tick stuff to show
    xticks,xlabels,yticks,ylabels = getMapTicks(m,xmin,xmax,ymin,ymax)
    plt.sca(ax)
    plt.tick_params(axis='both',direction='in',right='on',colors='white')
    plt.xticks(xticks,xlabels,size=8)
    plt.yticks(yticks,ylabels,size=8)
    for tick in ax.axes.yaxis.get_major_ticks():
        tick.set_pad(-38)
        tick.label2.set_horizontalalignment('right')
    for tick in ax.axes.xaxis.get_major_ticks():
        tick.set_pad(-15)
        tick.label2.set_verticalalignment('top')
    [i.set_color("white") for i in plt.gca().get_xticklabels()]
    [i.set_color("white") for i in plt.gca().get_yticklabels()]
    
    #render liquefaction
    lqdat = lqgrid.getData().copy() * 100.0
    
    palettelq = cm.autumn_r
    i = np.where(lqdat < 2.0)
    lqdat[i] = 0
    lqdat[iwater] = 0
    lqdatm = np.ma.masked_equal(lqdat, 0)
    palettelq.set_bad(clear_color,alpha=0.0)

    xmin, xmax, ymin, ymax = lqgrid.getRange()
    lons = np.arange(xmin,xmax,lqgrid.getGeoDict()['xdim'])
    lons = lons[:lqgrid.getGeoDict()['ncols']] #make sure it's the right length, sometimes it is one too long, sometimes it isn't because of GMTGrid issue
    lats = np.arange(ymax, ymin, -lqgrid.getGeoDict()['ydim']) #backwards so it plots right side up
    lats = lats[:lqgrid.getGeoDict()['nrows']] #make sure it's the right length, sometimes it is one too long, sometimes it isn't because of GMTGrid issue
    #make meshgrid
    llons, llats = np.meshgrid(lons, lats)
    x, y = m(llons, llats)  #get projection coordinates
    lqprobhandle = m.pcolormesh(x, y, lqdatm, lw=0, cmap=palettelq,vmin=2.0,vmax=20.0,alpha=ALPHA,rasterized=True)
    #extent = getGridExtent(lqgrid,m)
    #lqprobhandle = m.imshow(lqdatm,cmap=palettelq,vmin=2.0,vmax=20.0,alpha=ALPHA,origin='upper',extent=extent)
    norm = mpl.colors.Normalize(vmin=2.0,vmax=20.0)
    cbarlq = m.colorbar(mappable=lqprobhandle,norm=norm,cmap=palettelq)
    cbarlq.solids.set_rasterized(True)
    # cbarlq.solids.set_edgecolor("face")
    # plt.draw()
    cbarlq.set_ticks([2.0,11.0,19.0])
    cbarlq.set_ticklabels(['Low','Medium','High'])# vertically oriented colorbar

    #render landslide
    topogrid2 = GMTGrid()
    topogrid2.loadFromGrid(topogrid)
    topogrid2.interpolateToGrid(lsgrid.geodict)
    lsdat = lsgrid.getData().copy() * 100.0
    clear_color = [0,0,0,0.0]
    palettels = cm.cool
    i = np.where(lsdat < 2.0)
    lsdat[i] = 0
    lsdat[iwater] = 0
    lsdatm = np.ma.masked_equal(lsdat, 0)
    palettels.set_bad(clear_color,alpha=0.0)

    xmin, xmax, ymin, ymax = lsgrid.getRange()
    lons = np.arange(xmin, xmax, lsgrid.getGeoDict()['xdim'])
    lons = lons[:lsgrid.getGeoDict()['ncols']] #make sure it's the right length, sometimes it is one too long, sometimes it isn't because of GMTGrid issue
    lats = np.arange(ymax, ymin, -lsgrid.getGeoDict()['ydim']) #backwards so it plots right side up
    #make meshgrid
    lats = lats[:lsgrid.getGeoDict()['nrows']]
    llons, llats = np.meshgrid(lons, lats)
    x, y = m(llons, llats)  #get projection coordinates
    lsprobhandle = m.pcolormesh(x, y, lsdatm, lw=0, cmap=palettels,vmin=2.0,vmax=20.0,alpha=ALPHA,rasterized=True)
    #extent = getGridExtent(lsgrid,m)
    #lsprobhandle = plt.imshow(lsdatm,cmap=palettels,vmin=2.0,vmax=20.0,alpha=ALPHA,origin='upper',extent=extent)

    #draw landslide colorbar on the left side
    axleft = fig.add_axes([leftx,0.1,0.033,0.8])
    norm = mpl.colors.Normalize(vmin=2.0,vmax=20.0)
    cbarls = mpl.colorbar.ColorbarBase(axleft, cmap=palettels,norm=norm,orientation='vertical')
    cbarls.solids.set_rasterized(True)
    cbarls.ax.yaxis.set_ticks_position('left')
    cbarls.set_ticks([2.0,11.0,19.0])
    cbarls.set_ticklabels(['Low','Medium','High'])# vertically oriented colorbar

    #draw roads on the map, if they were provided to us
    roadcolor = ROADCOLOR
    if colors.has_key('roadcolor'):
        roadcolor = '#'+colors['roadcolor']
    for road in roadslist:
        xy = list(road['geometry']['coordinates'])
        roadx,roady = zip(*xy)
        mapx,mapy = m(roadx,roady)
        m.plot(mapx,mapy,roadcolor,lw=1.0)

    #add city names to map with population >50,000 (add option later)
    if cityfile is not None:
        dmin = 0.04*(m.ymax-m.ymin)
        xyplotted = []
        cities = PagerCity(cityfile)
        #Find cities within bounding box
        boundcity = cities.findCitiesByRectangle(bounds=(xmin,xmax,ymin,ymax))
        #Just keep 7 biggest cities
        thresh = sorted([cit['pop'] for cit in boundcity])[-7]
        plotcity = [cit for cit in boundcity if cit['pop'] >= thresh]
        #For cities that are more than one xth of the xwidth apart, keep only the larger one
        pass #do later
        #Plot cities
        for cit in plotcity: #should sort so it plots them in order of population so larger cities are preferentially plotted - do later
            xi, yi = m(cit['lon'], cit['lat'])
            dist = [np.sqrt((xi-x0)**2+(yi-y0)**2) for x0, y0 in xyplotted]
            if not dist or np.min(dist) > dmin:
                m.scatter(cit['lon'], cit['lat'], c='k', latlon=True, marker='.')
                ax.text(xi, yi, cit['name'], ha='right', va='top', fontsize=8)
                xyplotted.append((xi, yi))

    #draw titles
    cbartitle_ls = 'Landslide\nProbability'
    plt.text(-1.0,1.03,cbartitle_ls,multialignment='left',axes=ax)

    cbartitle_ls = 'Liquefaction\nProbability'
    plt.text(titlex,1.03,cbartitle_ls,multialignment='left',axes=ax)

    axwidth = 20 #where can I get this from?

    #draw a map boundary, fill in oceans with water
    #water_color = [.47,.60,.81] # too similar to a blue shade in ls colorbar
    water_color = '#B8EEFF'
    m.drawmapboundary(fill_color=water_color)

    if isScenario:
        title = eventdict['loc']
    else:
        timestr = eventdict['time'].strftime('%b %d %Y')
        title = 'M%.1f %s v%i\n %s' % (eventdict['mag'],timestr,eventdict['version'],eventdict['loc'])

    #draw the title on the plot
    ax.set_title(title)

    #draw star at epicenter
    plt.sca(ax)
    elat,elon = eventdict['epicenter']
    ex,ey = m(elon,elat)
    plt.plot(ex,ey,'*',markeredgecolor='k',mfc='None',mew=1.5,ms=24)

    #fill in the lakes and rivers
    m.fillcontinents(color=clear_color,lake_color=water_color)
    m.drawrivers(color=water_color)

    #draw country boundaries
    countrycolor = COUNTRYCOLOR
    if colors.has_key('countrycolor'):
        countrycolor = '#'+colors['countrycolor']
    m.drawcountries(color=countrycolor,linewidth=1.0)

    #add map scale
    m.drawmapscale((xmax+xmin)/2,(ymin+(ymax-ymin)/150),clon,clat,np.round((((xmax-xmin)*110)/5)/10.)*10, barstyle='fancy')
    #draw coastlines
    m.drawcoastlines(color='#476C91',linewidth=0.5)

    #draw scenario watermark, if scenario
    if isScenario:
        plt.sca(ax)
        cx,cy = m(clon,clat)
        plt.text(cx,cy,'SCENARIO',rotation=45,alpha=0.10,size=72,ha='center',va='center',color='red')
    
    #plt.title(ptitle,axes=ax)
    outfile = os.path.join(outfolder,'%s.pdf' % eventdict['eventid'])
    pngfile = os.path.join(outfolder,'%s.png' % eventdict['eventid'])
    plt.savefig(outfile)
    plt.savefig(pngfile)
    print 'Saving map output to %s' % outfile
    print 'Saving map output to %s' % pngfile

if __name__ == '__main__':
    pass
