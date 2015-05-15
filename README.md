SecHaz
=====

SecHaz is a program that calculates liquefaction and landslide probabilities given ground motions
and a number of global data sets.

TODO:
-----
 * Add zoom capability
 * Add ability to turn off scenario watermark
 * Fix input panel plots - titles (I think) are reversed top-to-bottom.

Installation and Dependencies
-----------------------------

This package depends on:
 * numpy, the fundamental package for scientific computing with Python. <a href="http://www.numpy.org/">http://www.numpy.org/</a>  
 * matplotlib, a Python 2D plotting library which produces publication quality figures. <a href="<a href="http://matplotlib.org/index.html">http://matplotlib.org/index.html</a>
 * basemap, a mapping package built on top of matplotlib. <a href="http://matplotlib.org/basemap/">http://matplotlib.org/basemap/</a>
 * fiona, a package for reading vector data file formats (in this case, road shapefiles).
 * neicio, a Python library for reading/writing various spatial data formats (including ShakeMap grid.xml). 
 * neicmap, a Python library for doing various spatial calculations (distance, angle, etc.)
 * neicutil, a Python library which is a grab bag of interpolation routines, text manipulation functions, and time functions.

The best way to install numpy,matplotlib,scipy and fiona is to use one of the Python distributions described here:

<a href="http://www.scipy.org/install.html">http://www.scipy.org/install.html</a>

The Anaconda distributions have been successfully tested with sechaz.
 
You may need to open a new terminal window to ensure that the newly installed versions of python and pip
are in your path.

To install basemap:

If you are using anaconda (see above):

conda install basemap

Otherwise, see the installation instructions here:

http://matplotlib.org/basemap/users/installing.html

To install fiona:

If you are using anaconda (see above):

conda install fiona

Otherwise, see the installation instructions here:

http://toblerity.org/fiona/README.html#installation

To install neicio:

pip install git+git://github.com/usgs/neicio.git

To install neicmap:

pip install git+git://github.com/usgs/neicmap.git

To install neicutil:

pip install git+git://github.com/usgs/neicutil.git

To install this package:

pip install git+git://github.com/mhearne-usgs/secondary.git

To upgrade this package:

pip install -U git+git://github.com/mhearne-usgs/secondary.git

Data Sources
============
The base data used for this software comes from open sources,
and is processed using open source tools:

 * gdal - http://trac.osgeo.org/gdal/wiki/DownloadingGdalBinaries
 * GMT 5 - http://gmt.soest.hawaii.edu/projects/gmt/wiki/Installing

Topography:
USGS Global Multi-resolution Terrain Elevation Data 2010 (GMTED2010)
http://topotools.cr.usgs.gov/gmted_viewer/data/Grid_ZipFiles/md30_grd.zip

To convert this data to a GMT-style NetCDF 3 file, you can do the following on a Unix-like system:

 * unzip md30_grd.zip
 * gdal_translate -of AAIGrid md30_grd md30_grd.ascii
 * gmt grdreformat md30_grd.ascii=ei md30_grid_gmt5.grd=cs

Roads:
Global Roads Open Access Data Set (gROADS)
http://sedac.ciesin.columbia.edu/data/set/groads-global-roads-open-access-v1/data-download (requires registration)

Download all of the regional *shapefile* road data sets into a common data directory, and unzip them in place.
The *roadfolder* configuration option should point to that common data directory.

Automation
==========

An automation program is included with this distribution - autosec.py.

This program, when run (manually or at intervals by something like cron),
will search the NEIC one hour earthquake data feed (http://earthquake.usgs.gov/earthquakes/feed/v1.0/summary/4.5_hour.geojson) for new ShakeMaps for events
matching magnitude, mmi, or PAGER impact level criteria in a config file (see below).

When an event matching one or more of those thresholds is found, autosec.py will run sechaz.py and then 
email the PNG and PDF probability map to a configured list of email recipients.

Configuring SecHaz
==================
The liquefaction and landslide models are defined in a configuration file, which must be
located in the user's home folder in a sub-folder called ".secondary".  The file must be named
"config.ini".  As an example, for a user called "frogers" would need to have a file called:

Mac:
/Users/frogers/.secondary/config.ini

Linux:
/home/frogers/.secondary/config.ini

Windows (untested!):
c:\Documents and Settings\frogers\.secondary\config.ini

Below is a sample configuration file:
<pre>
[LIQUEFACTION_MODEL]
#layers must be in GMT NetCDF format
cti_layer = /Users/frogers/secondary/data/globalcti.grd
vs30_layer = /Users/frogers/secondary/data/global_vs30.grd

#coefficients for the model
b0 = 24.1
b1 = 2.067
b2 = 0.355
b3 = -4.784

#certain variables are available from the ShakeMap:
#MW = magnitude
#PGV = pgv layer from the grid
#PGA = pga layer from the grid
#MMI = mmi layer from the grid

#terms to be multiplied by corresponding coefficients
#certain operators are allowed: +,-,*,/,log,power,log10
b1_term = log((PGA/100)*(power(MW,2.56)/power(10,2.24)))
b2_term = cti/100
b3_term = log(vs30)

#what is the grid to which all other grids will be resampled?
baselayer = vs30

#custom settings for rendering of individual layers - value must be a valid matplotlib colormap name.
vs30_colormap = jet_r

[LANDSLIDE_MODEL]
#layers must be in GMT NetCDF format
cohesion_layer = /Users/frogers/secondary/data/cohesion_10i.grd
slope_layer = /Users/frogers/secondary/data/slope_50.grd

#coefficients for the model
b0 = -7.15
b1 = 0.0604
b2 = 0.000825
b3 = 0.0201
b4 = 1.45e-05

#certain variables are available from the ShakeMap:
#MW = magnitude
#PGV = pgv layer from the grid
#PGA = pga layer from the grid
#MMI = mmi layer from the grid

#terms to be multiplied by corresponding coefficients
#certain operators are allowed: +,-,*,/,log,power,log10
b1_term = PGA
b2_term = slope
b3_term = cohesion/10.0
b4_term = PGA*slope

#what is the grid wo which all other grids will be resampled?
baselayer = cohesion

[SHAKEMAP]
variables = PGV,PGA,MW

[OUTPUT]
folder = /Users/frogers/secondary/output

[MAPDATA]
#required layers
topo = /Users/frogers/secondary/data/etopo1_bed_g_f4.grd
slope = /Users/frogers/secondary/data/slope_50.grd

#required slope parameters - pixels with slope > slopemax will not be liquefied
slopemin = 5.0
slopemax = 5.0

#optional layers
roadfolder = /Users/mhearne/secondary/data/roads
roadcolor = 00FF00
countrycolor = 177F10
</pre>

Configuring AutoSec
===================

autosec.py looks for a file in the ~/.secondary folder called mailconfig.ini.  That file should 
look something like this:

<pre>
[THRESHOLDS]
eis = yellow
mag = 7.5
mmi = 7

[MAIL]
server = yourmailserver.yourorganization.org
sender = genericaddress@yourorganization.org
recipients = user1@yourorganization.org,user2@yourorganization.org,
</pre>

Running SecHaz
==============

<pre>
usage: sechaz.py [-h] [-c CONFIGFILE] [-r] [-n] [GRIDFILE]

Run the landslide and liquefaction models defined by coefficients found in a
config.ini file. Output is a pdf map with liquefaction/landslide results
layered on topography.

positional arguments:
  GRIDFILE              ShakeMap grid.xml file (default: None)

optional arguments:
  -h, --help            show this help message and exit
  -c CONFIGFILE, --config CONFIGFILE
                        Use a custom config file (default:
                        /Users/mhearne/.secondary/config.ini)
  -r, --roads           Draw roads (default: False)
  -n, --noscenario      Turn off scenario watermark (use with caution)
                        (default: False)
</pre>