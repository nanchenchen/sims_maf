# to run:
# runDriver.py allSlicerCfg.py

# Example MAF config file which runs each type of available slicer.

import numpy as np
from lsst.sims.maf.driver.mafConfig import configureSlicer, configureMetric, makeDict

# Setup Database access.  Note:  Only the "root.XXX" variables are passed to the driver.
root.outputDir = './Allslicers'
root.dbAddress = {'dbAddress':'sqlite:///opsimblitz2_1060_sqlite.db'}
root.opsimName = 'opsimblitz2_1060'
root.figformat ='png'
root.dpi=200


root.verbose = True
# getConfig copies OpSim configuration parameters into output directory (if True).
root.getConfig = False
# Setup a list to hold all the slicers we want to run
slicerList=[]

# How many Healpix sides to use
nside=128

# List of SQL constraints.
# If multiple constraints are listed in a slicer object, they are looped over and each one is executed individualy.  
constraints = ["filter = \'%s\'"%'r', "filter = \'%s\' and night < 730"%'r']

# Configure a Healpix slicer:
# Configure 2 metrics to run on the Healpix slicer.  
m1 = configureMetric('CountMetric', kwargs={'col':'expMJD'},
                     plotDict={'percentileClip':80., 'units':'#'},
                     summaryStats={'MeanMetric':{},'RmsMetric':{}})
m2 = configureMetric('Coaddm5Metric',
                     plotDict={'zp':27., 'percentileClip':95, 'units':'Co-add m5 - %.1f'%27.})
# Combine metrics in a dictionary
metricDict = makeDict(m1,m2)
# Generate the slicer configuration, passing in the metric configurations and SQL constraints
slicer = configureSlicer('HealpixSlicer',
                          kwargs={"nside":nside,'spatialkey1':"fieldRA", 'spatialkey2':"fieldDec"},
                          metricDict=metricDict, constraints=constraints)
# Add the slicer to the list of slicers
slicerList.append(slicer)

# Run the same metrics, but now use the hexdither field positions:
# As before, but new spatialkeys and add a metadata keyword so the previous files don't get overwritten
slicer = configureSlicer('HealpixSlicer',
                          kwargs={"nside":nside,'spatialkey1':"ditheredRA", 'spatialkey2':"ditheredDec"},
                          metricDict=metricDict, constraints=constraints, metadata='dith')
# Add this slicer to the list of slicers
slicerList.append(slicer)


# Configure a OneDSlicer:
# Configure a new metric
m1 = configureMetric('CountMetric', kwargs={'col':'slewDist'}, plotDict={'logScale':True})
metricDict=makeDict(m1)
slicer = configureSlicer('OneDSlicer', kwargs={"sliceColName":'slewDist'},
                          metricDict=metricDict, constraints=constraints)
slicerList.append(slicer)


# Configure an OpsimFieldSlicer:
m1 = configureMetric('MinMetric', kwargs={'col':'airmass'},
                     plotDict={'cmap':'RdBu'})
m4 = configureMetric('MeanMetric', kwargs={'col':'normairmass'})
m3 = configureMetric('Coaddm5Metric')
m7 = configureMetric('CountMetric', kwargs={'col':'expMJD'},
                     plotDict={'units':"Number of Observations", 'percentileClip':80.})
metricDict = makeDict(m1,m3,m4,m7)
slicer = configureSlicer('OpsimFieldSlicer', metricDict=metricDict, constraints=constraints )
slicerList.append(slicer)


# Configure a UniSlicer.  Note new SQL constraints are passed
m1 = configureMetric('MeanMetric', kwargs={'col':'airmass'})
slicer = configureSlicer('UniSlicer', metricDict=makeDict(m1), constraints=['night < 750'] )
slicerList.append(slicer)

# Configure an Hourglass filter slicer/metric
m1=configureMetric('HourglassMetric')
slicer = configureSlicer('HourglassSlicer', metricDict=makeDict(m1), constraints=['night < 750',''])
slicerList.append(slicer)

# Try out the fO slicer
fOnside = 64
m1 = configureMetric('CountMetric',
                     kwargs={'col':'expMJD', 'metricName':'fO'},
                     plotDict={'units':'Number of Visits', 'xMin':0,
                               'xMax':1500},
                     summaryStats={'fOArea':{'nside':fOnside},
                                   'fONv':{'nside':fOnside}})
slicer = configureSlicer('fOSlicer', kwargs={'nside':fOnside},
                         metricDict=makeDict(m1),
                         constraints=[''])
slicerList.append(slicer)

# Save all the slicers to the config
root.slicers=makeDict(*slicerList)

# Optional comment string
root.comment = 'Example script that runs each of the slicers'
