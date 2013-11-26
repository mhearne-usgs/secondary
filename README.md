SecHaz
=====

SecHaz is a program that calculates liquefaction and landslide probabilities given ground motions
and a number of global data sets.

INSTALLATION
============

SecHaz has the following dependencies:
- Python 2.7+ (not 3.X!)
- numpy 1.5+
- matplotlib 1.3.0+
- gdal 1.10.0+
- pagerlib

Depending on your platform, there are a number of different ways to install these dependencies.  
The simplest way (by far) is to use a pre-packaged Python distribution that includes the 
"Scipy Stack" (a collection of Python tools useful for scientific analysis).  The distributions
that adhere to this standard are listed at 
<a href="http://www.scipy.org/install.html">http://www.scipy.org/install.html</a>.

<em>NB: The Pyzo distribution installs Python 3.3, which is incompatible with the 2.x Python
code in STREC.</em>

Final dependencies
------------------
Once you have a Python distribution installed (Canopy, Anaconda, etc.), you will likely still need to install 
obspy and possibly pytz.  The best way to do this is using <b>pip</b>.  pip comes bundled with Anaconda and pythonxy, 
but needs an extra step on Canopy.

From the command line, type:
<pre>
[sudo] easy_install pip
</pre>

(sudo may be necessary, depending on whether you have permissions to install software with your regular account, and how Python has been installed).

<pre>
git clone https://github.com/mhearne-usgs/pagerlib.git
pip install pagerlib/
</pre>

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