SecHaz
=====

SecHaz is a program that calculates liquefaction and landslide probabilities given ground motions
and a number of global data sets.

TODO:
-----
 * Add roads
 * Add zoom capability
 * Add ability to turn off scenario watermark
 * Add ability to choose a particular config file
 * Add epicentral star
 * Add zoom capability

Installation and Dependencies
-----------------------------

This package depends on:
 * numpy, the fundamental package for scientific computing with Python. <a href="http://www.numpy.org/">http://www.numpy.org/</a>  
 * matplotlib, a Python 2D plotting library which produces publication quality figures. <a href="<a href="http://matplotlib.org/index.html">http://matplotlib.org/index.html</a>
 * basemap, a mapping package built on top of matplotlib. <a href="http://matplotlib.org/basemap/">http://matplotlib.org/basemap/</a>
 * neicio, a Python library for reading/writing various spatial data formats (including ShakeMap grid.xml). 
 * neicmap, a Python library for doing various spatial calculations (distance, angle, etc.)
 * neicutil, a Python library which is a grab bag of interpolation routines, text manipulation functions, and time functions.

The best way to install numpy,matplotlib,and scipy is to use one of the Python distributions described here:

<a href="http://www.scipy.org/install.html">http://www.scipy.org/install.html</a>

Anaconda and Enthought distributions have been successfully tested with sechaz.
 
You may need to open a new terminal window to ensure that the newly installed versions of python and pip
are in your path.

To install basemap:

If you are using anaconda (see above):

conda install basemap

Otherwise, see the installation instructions here:

http://matplotlib.org/basemap/users/installing.html

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


Configuring SecHaz
==================
The liquefaction and landslide models are defined in the configuration file, which must be
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
cti_layer = /Users/mhearne/secondary/data/globalcti.grd
vs30_layer = /Users/mhearne/secondary/data/global_vs30.grd

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

[LANDSLIDE_MODEL]
#layers must be in GMT NetCDF format
cohesion_layer = /Users/mhearne/secondary/data/cohesion_10i.grd
slope_layer = /Users/mhearne/secondary/data/slope_50.grd

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
folder = /Users/mhearne/secondary/output

[MAPDATA]
topo = /Users/mhearne/secondary/data/etopo1_bed_g_f4.grd
slope = /Users/mhearne/secondary/data/slope_50.grd
slopemin = 5.0
slopemax = 5.0
</pre>

Running SecHaz
==============

<pre>
usage: sechaz.py [-h] [GRIDFILE]

Run the landslide and liquefaction models defined by coefficients found in a
config.ini file. Output is a pdf map with liquefaction/landslide results
layered on topography.

positional arguments:
  GRIDFILE    ShakeMap grid.xml file

optional arguments:
  -h, --help  show this help message and exit
</pre>