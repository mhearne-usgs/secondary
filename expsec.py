#!/usr/bin/env python

import sys
import os.path

from neicio.gmt import GMTGrid
from neicio.esri import EsriGrid
import numpy as np

if __name__ == '__main__':
    folder = sys.argv[1]
    popfile = sys.argv[2]
    slidefile = os.path.join(folder,'landslide.grd')
    slidegrid = GMTGrid(slidefile)
    slidegrid.load()
    bounds = slidegrid.getRange()
    popgrid = EsriGrid(popfile)
    popgrid.load(bounds=bounds)
    slidegrid.interpolateToGrid(popgrid.geodict)
    expsum = np.nansum(slidegrid.griddata * popgrid.griddata)
    parts = folder.split(os.sep)
    print 'Event %s exposure: %i' % (parts[-1],int(expsum))
    
