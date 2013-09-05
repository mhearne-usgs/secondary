secondary
=========

Secondary Hazards code

configuration
=========

The config file should look like this, where [DATAHOME] should be replaced with path to data files:
<pre>
[LIQUEFACTION]
b0 = 23.21
bpga = 1.990
bcti = 0.113
bvs30 = -4.609
ctifile = [DATAHOME]/cti2/globalcti.grd
topofile = [DATAHOME]/etopo/etopo1_bed_g_f4.flt
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
frictionfile = [DATAHOME]/godt/geology/friction.flt
cohesionfile = [DATAHOME]/godt/geology/cohesion_10i.flt
slopefile = [DATAHOME]/godt/slope/slope_50.flt
maxslopefile = [DATAHOME]/godt/slope/slope_50.flt
slopemin = 5.0
</pre>