#!/usr/bin/env python

#stdlib imports
import ConfigParser
import os.path
import sys
import warnings
import argparse

#turn off all warnings...
warnings.filterwarnings('ignore')

#local imports
from neicio.shake import ShakeGrid
from neicio.gmt import GMTGrid
from secondary.model import LogisticModel,getModelNames
from secondary.map import makeDualMap,renderLayer
from secondary.shape import getRoadFile

#third party imports
from matplotlib import cm

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
                

def main(args):
    #define location for config file
    homedir = os.path.expanduser("~") #where is the user's home directory?
    configfile = args.configFile
    
    shakefile = args.shakefile
    shakemap = ShakeGrid(shakefile)
    shakeheader = shakemap.getAttributes()
    edict = {'mag':shakeheader['event']['magnitude'],
             'time':shakeheader['event']['event_timestamp'],
             'loc':shakeheader['event']['event_description'],
             'epicenter':(shakeheader['event']['lat'],shakeheader['event']['lon']),
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

    #if they have roads configured, go find the appropriate roads file
    hasRoads = config.has_option('MAPDATA','roadfolder')
    hasContinents = config.has_option('MAPDATA','continents')
    if hasRoads and hasContinents and args.roads:
        roadfolder = config.get('MAPDATA','roadfolder')
        contfile = config.get('MAPDATA','continents')
        roadfile = getRoadFile(contfile,roadfolder,edict['epicenter'][0],edict['epicenter'][1])
    else:
        roadfile = None

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
        for layername,layergrid in lm.layerdict.iteritems():
            layerfile = os.path.join(outfolder,layername+'.grd')
            print 'Saving input grid %s to %s...' % (layername,layerfile)
            layergrid.save(layerfile)
            renderLayer(layergrid,layername,outfolder,edict,model,colormaps)

    topofile = config.get('MAPDATA','topo')
    bigbounds = shakemap.getRange()
    xdim = shakemap.geodict['xdim']
    ydim = shakemap.geodict['xdim']
    bigbounds = (bigbounds[0]-xdim*2,bigbounds[1]+xdim*2,bigbounds[2]-ydim*2,bigbounds[3]+ydim*2)
    topogrid = GMTGrid(topofile,bounds=bigbounds)
    topoout = os.path.join(outfolder,'topography.grd')
    print 'Saving topography to %s' % topoout
    topogrid.save(topoout)
    
    print 'Saving slope to %s' % slopeout
    slopegrid.save(slopeout)

    isScenario = shakeheader['shakemap_grid']['shakemap_event_type'].lower() == 'scenario'
    if args.noscenario:
        isScenario = False
    timestr = shakeheader['event']['event_timestamp'].strftime('%b %d %Y')
    location = shakeheader['event']['event_description']
    makeDualMap(probdict['liquefaction'],probdict['landslide'],topogrid,slopegrid,edict,outfolder,
                isScenario=isScenario,roadfile=roadfile)

if __name__ == '__main__':
    usage = """Run the landslide and liquefaction models defined by coefficients found in a config.ini file.
    Output is a pdf map with liquefaction/landslide results layered on topography."""
    parser = argparse.ArgumentParser(description=usage,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("shakefile", metavar='GRIDFILE',nargs='?',
                        help='ShakeMap grid.xml file')
    parser.add_argument('-c','--config', dest='configFile', 
                        default=DEFAULT_CONFIG,
                        help='Use a custom config file')
    parser.add_argument('-r','--roads', action='store_true', 
                        default=False,
                        help='Draw roads (can take significantly longer)')
    parser.add_argument('-n','--noscenario', action='store_true', 
                        default=False,
                        help='Turn off scenario watermark (use with caution)')
    args = parser.parse_args()
    main(args)
    
    
