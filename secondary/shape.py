#!/usr/bin/env python

#stdlib imports
import os.path
import glob

#third party
import fiona
from shapely.geometry import LineString,Point,MultiLineString,MultiPolygon,shape,Polygon

CONTINENTS = {'Asia':['gROADS-v1-asia.shp'],
              'North America':['gROADS-v1-americas.shp'],
              'Europe':['gROADS-v1-europe.shp'],
              'Africa':['gROADS-v1-africa.shp'],
              'South America':['gROADS-v1-americas.shp'],
              'Oceania':['gROADS-v1-oceania-west.shp','gROADS-v1-oceania-east.shp'],
              'Australia':['gROADS-v1-oceania-east.shp'],
              'Antarctica':[]}

def expandRoadPath(shapedir,roadfile):
    folders = os.listdir(shapedir)
    shapefile = None
    for folder in folders:
        shapefolder = os.path.join(shapedir,folder)
        shapefile = os.path.join(shapefolder,roadfile)
        if os.path.isfile(shapefile):
            foundFile = True
            break
    return shapefile

def getRoadFile(contfile,shapedir,lat,lon):
    continents = fiona.open(contfile,'r')
    shapefolders = os.listdir(shapedir)
    point = Point(lon,lat)
    roadfiles = None
    for continent in continents:
        cshape = shape(continent['geometry'])
        cname = continent['properties']['CONTINENT']
        if cshape.contains(point):
            roadfiles = CONTINENTS[cname]
            break
    continents.close()
    if not len(roadfiles):
        return None
    if len(roadfiles) == 1:
        return expandRoadPath(shapedir,roadfiles[0])
    else:
        for roadfile in roadfiles:
            roadfile = expandRoadPath(shapedir,roadfile)
            roads = fiona.open(roadfile,'r')
            xmin,ymin,xmax,ymax = roads.bounds
            if lon >= xmin and lon <= xmax and lat >= ymin and lat <= ymax:
                roads.close()
                return roadfile
            roads.close()
    return None
    
if __name__ == '__main__':
    shapedir = '/Users/mhearne/data/roads/'
    contfile = '/Users/mhearne/secondary/data/continent.shp'
    points = {'Sinai Peninsula':(30.966717,33.607178),
              'New Zealand':(-43.582878,172.518311),
              'Maine':(44.777327,-68.760681),
              'French Polynesia':(-17.623491,-149.437408),
              'Papua New Guinea':(-4.592852,138.955078),
              'Hawaii':(19.651318,-155.462036)}
    for key,value in points.iteritems():
        lat,lon = value
        roadfile = getRoadFile(contfile,shapedir,lat,lon)
        p,f = os.path.split(roadfile)
        print 'Retrieved %s file for location in %s' % (f,key)
        xmin = lon - 2.0
        xmax = lon + 2.0
        ymin = lat - 2.0
        ymax = lat + 2.0
        box = Polygon([(xmin,ymin),(xmax,ymin),(xmax,ymax),(xmin,ymax)])
        roadsegments = fiona.open(roadfile,'r')
        foundRoad = False
        for roadseg in roadsegments:
            roadsegment = shape(roadseg['geometry'])
            shapeint = box.intersection(roadsegment)
            if isinstance(shapeint,LineString):
                foundRoad = True
                break
            else:
                if len(shapeint):
                    foundRoad = True
                    break
        if foundRoad:
            print '\tFound a road segment.'
        else:
            print '\tNo road segments found.'
    
            
