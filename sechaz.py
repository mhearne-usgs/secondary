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
from secondary.map import makeDualMap

CONFIGFILE = 'config.ini'
AZIMUTH_DEFAULT = 90
DEFAULT_CONFIG = os.path.join(os.path.expanduser('~'),'.secondary','config.ini')

#TODO
# - Add roads
# - Add coastlines - drawcoastlines() is pretty slow at full res, and may not match topo.  Plot zero topo contour instead?

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

    #get the thresholds for liquefaction/landslide model
    slopemin = float(config.get('MAPDATA','slopemin'))*100
    slopemax = float(config.get('MAPDATA','slopemax'))*100
    
    probdict = {}
    gridbounds = [999,-999,999,-999] #this will hold the smallest bounding box enclosing both models
    for model in getModelNames(configfile):
        lm = LogisticModel(configfile,shakefile,model)
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

    topofile = config.get('MAPDATA','topo')
    topogrid = GMTGrid(topofile,bounds=shakemap.getRange())
    topoout = os.path.join(outfolder,'topography.grd')
    print 'Saving topography to %s' % topoout
    topogrid.save(topoout)
    
    print 'Saving slope to %s' % slopeout
    slopegrid.save(slopeout)

    isScenario = shakeheader['shakemap_grid']['shakemap_event_type'].lower() == 'scenario'
    timestr = shakeheader['event']['event_timestamp'].strftime('%b %d %Y')
    location = shakeheader['event']['event_description']
    makeDualMap(probdict['liquefaction'],probdict['landslide'],topogrid,slopegrid,edict,outfolder,isScenario=isScenario)

if __name__ == '__main__':
    usage = """Run the landslide and liquefaction models defined by coefficients found in a config.ini file.
    Output is a pdf map with liquefaction/landslide results layered on topography."""
    parser = argparse.ArgumentParser(description=usage,formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument("shakefile", metavar='GRIDFILE',nargs='?',
                        help='ShakeMap grid.xml file')
    parser.add_argument('-c','--config', dest='configFile', 
                        default=DEFAULT_CONFIG,
                        help='Use a custom config file')
    args = parser.parse_args()
    main(args)
    
    
