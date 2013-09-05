secondary
=========

Secondary Hazards code

configuration
=========

The config file should look like this:
<pre>
[LIQUEFACTION]
b0 = 23.21
bpga = 1.990
bcti = 0.113
bvs30 = -4.609
ctifile = /Users/mhearne/data/cti2/globalcti.grd
topofile = /Users/mhearne/data/etopo/etopo1_bed_g_f4.flt
slopemax = 5.0

[LANDSLIDE]
b0 = -7.15E+00
bpga = 6.04E-02
bslope = 8.25E-04
bcohesion = 2.01E-02
bpgaslope = 1.45E-05
bcti = 1.0
bmaxslope = 1.0
bfriction = 1.0
frictionfile = /Users/mhearne/data/godt/geology/friction.flt
cohesionfile = /Users/mhearne/data/godt/geology/cohesion_10i.flt
slopefile = /Users/mhearne/data/godt/slope/slope_50.flt
maxslopefile = /Users/mhearne/data/godt/slope/slope_50.flt
slopemin = 5.0
</pre>