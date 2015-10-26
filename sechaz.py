#!/usr/bin/env python

#stdlib imports
import ConfigParser
import os.path
import sys
import warnings
import argparse
import glob
import urllib2
import tempfile

#turn off all warnings...
warnings.filterwarnings('ignore')

#local imports
from neicio.shake import ShakeGrid
from neicio.gmt import GMTGrid
from secondary.model import LogisticModel,getModelNames
from secondary.map import makeDualMap,renderPanel,ROADCOLOR
from secondary.shape import getRoadFile

#third party imports
from matplotlib import cm
import fiona
import numpy as np

CONFIGFILE = 'config.ini'
AZIMUTH_DEFAULT = 90
DEFAULT_CONFIG = os.path.join(os.path.expanduser('~'),'.secondary','config.ini')

def getColorMaps(configfile):
    config = ConfigParser.ConfigParser()
    config.readfp(open(configfile))
    sections = config.sections()
    colormaps = {}
    maps = [m for m in cm.datad]
    lmaps = [m.lower() for m in cm.datad]
    for section in sections:
        if section.lower().find('_model') > -1:
            modelname = section.split('_')[0].lower()
            layers = []
            for option in config.options(section):
                if option.lower().find('_layer') > -1:
                    layers.append(option.split('_')[0].lower())
            for option in config.options(section):
                if option.find('_colormap') > -1:
                    layer = option.split('_')[0].lower()
                    if layer not in layers:
                        print 'Colormap specified for non-existent layer "%s".  Skipping.'
                        continue
                    value = config.get(section,option).strip().lower()
                    if value not in lmaps:
                        cmapstr = ','.join(maps)
                        print 'Colormap specified for layer "%s" of "%s" is not valid.  The list of valid colormaps is: %s. Skipping' % cmapstr
                        continue
                    cmap = maps[lmaps.index(value)]
                    if colormaps.has_key(modelname):
                        mdict = colormaps[modelname]
                        mdict[layer] = cmap
                    else:
                        mdict = {layer:cmap}
                        colormaps[modelname] = mdict
        
    return colormaps                

def getGridURL(gridurl):
    gridfile = None
    try:
        fh = urllib2.urlopen(gridurl)
        data = fh.read()
        fd,gridfile = tempfile.mkstemp()
        os.close(fd)
        f = open(gridfile,'wt')
        f.write(data)
        f.close()
        fh.close()
    except:
        raise IOError,'Could not retrieve data from %s' % gridurl
    return gridfile

def isURL(gridurl):
    isURL = False
    try:
        fh = urllib2.urlopen(gridurl)
        isURL = True
    except:
        pass
    return isURL

def main(args):
    #define location for config file
    homedir = os.path.expanduser("~") #where is the user's home directory?
    configfile = args.configFile
    
    shakefile = args.shakefile

    if not os.path.isfile(shakefile):
        if isURL(shakefile):
            shakefile = getGridURL(shakefile) #returns a file object
        else:
            print 'Could not find "%s" as a file or a url.  Returning.' % (shakefile)
    
    shakemap = ShakeGrid(shakefile)
    #figure out the bounds that are greater than the biggest bounds
    #of any of the grids
    shakerange = shakemap.getRange()
    lonrange = shakerange[1] - shakerange[0]
    latrange = shakerange[3] - shakerange[2]
    xmin = shakerange[0] - lonrange*0.1
    xmax = shakerange[1] + lonrange*0.1
    ymin = shakerange[2] - latrange*0.1
    ymax = shakerange[3] + latrange*0.1
    bigbounds = (xmin,xmax,ymin,ymax)
    #
    shakeheader = shakemap.getAttributes()
    edict = {'mag':shakeheader['event']['magnitude'],
             'time':shakeheader['event']['event_timestamp'],
             'loc':shakeheader['event']['event_description'],
             'epicenter':(shakeheader['event']['lat'],shakeheader['event']['lon']),
             'version':int(shakeheader['shakemap_grid']['shakemap_version']),
             'eventid':shakeheader['shakemap_grid']['event_id']}
    config = ConfigParser.RawConfigParser()
    config.read(configfile)
    network = shakeheader['shakemap_grid']['shakemap_originator']
    eventcode =  shakeheader['shakemap_grid']['shakemap_id']
    if eventcode.startswith(network):
        eventid = eventcode
    else:
        eventid = network + eventcode
    outfolder = os.path.join(config.get('OUTPUT','folder'),eventid)
    if not os.path.isdir(outfolder):
        os.makedirs(outfolder)

    slopefile = config.get('MAPDATA','slope')
    slopegrid = GMTGrid(slopefile,bounds=shakemap.getRange())
    slopeout = os.path.join(outfolder,'slope.grd')

    cityfile = config.get('MAPDATA','cityfile')
    
    #get all of the colors that people want
    colors = {}
    for option in config.options('MAPDATA'):
        if option.endswith('color'):
            colors[option] = config.get('MAPDATA',option)

    #if they have roads configured, go find the appropriate roads segments
    hasRoads = config.has_option('MAPDATA','roadfolder')
    roadslist = []
    if hasRoads and args.roads:
        roadroot = config.get('MAPDATA','roadfolder')
        xmin,xmax,ymin,ymax = shakemap.getRange()
        for folder in os.listdir(roadroot):
            roadfolder = os.path.join(roadroot,folder)
            shpfiles = glob.glob(os.path.join(roadfolder,'*.shp'))
            if len(shpfiles):
                shpfile = shpfiles[0]
                f = fiona.open(shpfile)
                shapes = list(f.items(bbox=(xmin,ymin,xmax,ymax)))
                for shapeid,shapedict in shapes:
                    roadslist.append(shapedict)
                f.close()

    #get the thresholds for liquefaction/landslide model
    slopemin = float(config.get('MAPDATA','slopemin'))*100
    slopemax = float(config.get('MAPDATA','slopemax'))*100
    
    probdict = {}
    gridbounds = [999,-999,999,-999] #this will hold the smallest bounding box enclosing both models
    for model in getModelNames(configfile):
        lm = LogisticModel(configfile,shakefile,model)
        colormaps = getColorMaps(configfile)
        print 'Equation for %s model:' % model
        print 
        print lm.getEquation()
        print
        P = lm.calculate()
        probgrid = GMTGrid()
        probgrid.griddata = P.copy()
        probgrid.geodict = lm.layerdict[lm.layerdict.keys()[0]].geodict.copy()

        #resample the slope grid to model
        slopegrid2 = GMTGrid()
        slopegrid2.loadFromGrid(slopegrid)
        slopegrid2.interpolateToGrid(probgrid.geodict)

        if model == 'liquefaction':
            ithresh = slopegrid2.griddata > slopemax
        else:
            ithresh = slopegrid2.griddata < slopemin

        probgrid.griddata[ithresh] = 0.0
        
        xmin,xmax,ymin,ymax = probgrid.getRange()
        if xmin < gridbounds[0]:
            gridbounds[0] = xmin
        if xmax > gridbounds[1]:
            gridbounds[1] = xmax
        if ymin < gridbounds[2]:
            gridbounds[2] = ymin
        if ymax > gridbounds[3]:
            gridbounds[3] = ymax
        probdict[model] = probgrid
        probfile = os.path.join(outfolder,'%s.grd' % model)
        print 'Saving %s model output to %s' % (model,probfile)
        probgrid.save(probfile)
        #renderPanel(lm,colormaps,outfolder,edict)
        # for layername,layergrid in lm.layerdict.iteritems():
        #     layerfile = os.path.join(outfolder,layername+'.grd')
        #     print 'Saving input grid %s to %s...' % (layername,layerfile)
        #     layergrid.save(layerfile)
        #     renderLayer(layergrid,layername,outfolder,edict,model,colormaps)

    topofile = config.get('MAPDATA','topo')
    #bigbounds = shakemap.getRange()
    xdim = shakemap.geodict['xdim']
    ydim = shakemap.geodict['xdim']
    #bigbounds = (bigbounds[0]-xdim*4,bigbounds[1]+xdim*4,bigbounds[2]-ydim*4,bigbounds[3]+ydim*4)
    topogrid = GMTGrid(topofile,bounds=bigbounds)
    topogrid = adjustTopoGrid(topogrid,bigbounds) #make this grid as big as bigbounds if we hit an upper or lower bound
    topoout = os.path.join(outfolder,'topography.grd')
    print 'Saving topography to %s' % topoout
    topogrid.save(topoout)
    
    print 'Saving slope to %s' % slopeout
    slopegrid.save(slopeout)

    isScenario = shakeheader['shakemap_grid']['shakemap_event_type'].lower() == 'scenario'
    if args.noscenario:
        isScenario = False
    timestr = renderDate(shakeheader['event']['event_timestamp'])
    location = shakeheader['event']['event_description']
    #hillshfile = config.get('MAPDATA','hillshadefile')
    #hillshgrid = GMTGrid(hillshfile,bounds=bigbounds)
    makeDualMap(probdict['liquefaction'],probdict['landslide'],topogrid,slopegrid,edict,outfolder,isScenario=isScenario,roadslist=roadslist,colors=colors,cityfile=cityfile)

def adjustTopoGrid(topogrid,bigbounds):
    xdim = topogrid.geodict['xdim']
    ydim = topogrid.geodict['ydim']
    nrows,ncols = topogrid.griddata.shape
    xmin,xmax,ymin,ymax = topogrid.getRange()
    bxmin,bxmax,bymin,bymax = bigbounds
    if bymin <= ymin:
        nrows_bottom = int(np.ceil((ymin - bymin)/ydim))
        newymin = ymin - nrows_bottom*ydim
        bottom = np.ones((nrows_bottom,ncols))*np.nan
        topogrid.griddata = np.append(topogrid.griddata,bottom,axis=0)
        topogrid.geodict['ymin'] = newymin
    if ymax >= bymax:
        nrows_top = int(np.ceil((bymax - ymax)/ydim))
        newymax = ymax + nrows_top*ydim
        top = np.ones((nrows_top,ncols))*np.nan
        topogrid.griddata = np.append(top,topogrid.griddata,axis=0)
        topogrid.geodict['ymax'] = newymax
    
    newrows,newcols = topogrid.griddata.shape 
    topogrid.geodict['nrows'] = newrows
    return topogrid
        

def renderDate(dtime):
    months = ['Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec']
    month = months[dtime.month - 1]
    timestr = '%s %i %i' % (month,dtime.day,dtime.year)
    return timestr

if __name__ == '__main__':
    usage = """Run the landslide and liquefaction models defined by coefficients found in a config.ini file.
    Output is a pdf map with liquefaction/landslide results layered on topography."""
    parser = argparse.ArgumentParser(description=usage,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("shakefile", metavar='GRIDFILE',nargs='?',
                        help='ShakeMap grid.xml file (or url of same)')
    parser.add_argument('-c','--config', dest='configFile', 
                        default=DEFAULT_CONFIG,
                        help='Use a custom config file')
    parser.add_argument('-r','--roads', action='store_true', 
                        help='Draw roads')
    parser.add_argument('-n','--noscenario', action='store_true', #00FF00
                        default=False,
                        help='Turn off scenario watermark (use with caution)')
    args = parser.parse_args()
    main(args)
    
    
