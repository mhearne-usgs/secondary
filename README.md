SecHaz
=====

SecHaz is a program that calculates liquefaction and landslide probabilities given ground motions
and a number of global data sets.

Installation and Dependencies
-----------------------------

This package depends on:
 * numpy, the fundamental package for scientific computing with Python. <a href="http://www.numpy.org/">http://www.numpy.org/</a>  
 * matplotlib, a Python 2D plotting library which produces publication quality figures. <a href="<a href="http://matplotlib.org/index.html">http://matplotlib.org/index.html</a>
 * neicio, a Python library for reading/writing various spatial data formats (including ShakeMap grid.xml). 
 * neicmap, a Python library for doing various spatial calculations (distance, angle, etc.)
 * neicutil, a Python library which is a grab bag of interpolation routines, text manipulation functions, and time functions.

The best way to install numpy,matplotlib,and scipy is to use one of the Python distributions described here:

<a href="http://www.scipy.org/install.html">http://www.scipy.org/install.html</a>

Anaconda and Enthought distributions have been successfully tested with sechaz.
 
You may need to open a new terminal window to ensure that the newly installed versions of python and pip
are in your path.

To install neicio:

pip install git+git://github.com/usgs/neicio.git

To install neicmap:

pip install git+git://github.com/usgs/neicmap.git

To install neicutil:

pip install git+git://github.com/usgs/neicutil.git

Installing and Configuring SecHaz
----------------

<pre>
mkdir ~/src
cd ~/src
git clone https://github.com/mhearne-usgs/secondary.git
cd secondary
./sechaz.py -c
</pre>

Running the script with the -c (configuration) flag will prompt you for values for the configuration
file, which will be written to $HOME/.secondary/config.ini.

To configure the program, you will need the following data files:

Roads: Any of the shapefiles found here: http://sedac.ciesin.columbia.edu/data/set/groads-global-roads-open-access-v1/data-download*

Borders: WDBII_border_f_L1.shp as found here: http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.3.zip

Coasts: GSHHS_f_L1.shp as found here: http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.3.zip

Rivers: WDBII_river_f*.shp files found here: http://www.soest.hawaii.edu/pwessel/gshhg/gshhg-shp-2.3.3.zip

*At the time of this writing, there is no routine for selecting the appropriate roads data file from a list based on input map parameters.
Users will have to modify the config file by hand to choose a roads layer appropriate for their region.

Running SecHaz
=========

<pre>
usage: sechaz.py [-h] [-z xmin xmax ymin ymax] [-r] [-t] [-d] [-c]
                 [-s SHAPECONFIG]
                 [GRIDFILE]

Run the landslide and liquefaction models defined by coefficients found in a
config.ini file. Output is a pdf map with liquefaction/landslide results
layered on topography.

positional arguments:
  GRIDFILE              ShakeMap grid.xml file

optional arguments:
  -h, --help            show this help message and exit
  -z xmin xmax ymin ymax, --zoom xmin xmax ymin ymax
                        zoom in map to geographic coordinates (xmin xmax ymin
                        ymax)
  -r, --roads           Draw roads found inside map
  -t, --table           Add table of summary statistics to figure
  -d, --disable-scenario
                        Turn scenario text off
  -c, --configure       Create config file
  -s SHAPECONFIG, --shapeconfig SHAPECONFIG
                        Create config file
</pre>