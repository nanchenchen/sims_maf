import numpy as np
from lsst.sims.maf.driver.mafConfig import configureSlicer, configureMetric, makeDict

# A driver config for the plots on https://confluence.lsstcorp.org/display/SIM/MAF+documentation
# and https://confluence.lsstcorp.org/display/SIM/MAF%3A++Writing+a+new+metric

root.outputDir = './Doc'
root.dbAddress = {'dbAddress':'sqlite:///opsimblitz2_1060_sqlite.db'}
root.opsimName = 'opsimblitz2_1060'
root.figformat ='png'
root.dpi=200

root.figformat='png'

slicerList=[]

m1 = configureMetric('CountMetric', kwargs={'col':'slewTime'}, plotDict={'logScale':True})
slicer = configureSlicer('OneDSlicer', kwargs={"sliceColName":'slewTime', 'binsize':5.},
                         metricDict=makeDict(m1), constraints=[''] )
root.slicers=makeDict(slicer)
slicerList.append(slicer)

constraints = ["filter = \'%s\'"%'r']
m1 = configureMetric('MinMetric', kwargs={'col':'airmass'}, plotDict={'cmap':'RdBu'})
metricDict = makeDict(m1)
slicer = configureSlicer('OpsimFieldSlicer', metricDict=metricDict, constraints=constraints )
root.slicers=makeDict(slicer)
slicerList.append(slicer)

nside=64
constraints = ["filter = \'%s\'"%'r']
m2 = configureMetric('Coaddm5Metric', plotDict={'zp':27., 'percentileClip':95, 'units':'Co-add m5 - %.1f'%27.} )
metricDict = makeDict(m2)
slicer = configureSlicer('HealpixSlicer',
                          kwargs={"nside":nside,'spatialkey1':"fieldRA", 'spatialkey2':"fieldDec"},
                          metricDict = metricDict, constraints=constraints)
root.slicers=makeDict(slicer)
slicerList.append(slicer)

constraints = ["filter = \'%s\'"%'r']
m3 = configureMetric('Coaddm5Metric')
m4 = configureMetric('MeanMetric', kwargs={'col':'normairmass'})
metricDict = makeDict(m3,m4)
slicer = configureSlicer('OpsimFieldSlicer', metricDict=metricDict, constraints=constraints )
root.slicers=makeDict(slicer)
slicerList.append(slicer)


m1 = configureMetric('CountMetric', kwargs={'col':'slewDist'})
slicer = configureSlicer('OneDSlicer', kwargs={"sliceColName":'slewDist'},
                         metricDict=makeDict(m1), constraints=[''] )
root.slicers=makeDict(slicer)
slicerList.append(slicer)


m1 = configureMetric('CountMetric', kwargs={'col':'expMJD'}, plotDict={'ylabel':'Number of Exposures'})
slicer = configureSlicer('OneDSlicer', kwargs={"sliceColName":'night','binsize':1.},
                         metricDict=makeDict(m1), constraints=[''] )
slicerList.append(slicer)

m1=configureMetric('HourglassMetric')
slicer = configureSlicer('HourglassSlicer', metricDict=makeDict(m1), constraints=['night < 750'] )
slicerList.append(slicer)



# How many Healpix sides to use
nside=64

# List of SQL constraints.  If multiple constraints are listed in a slicer object,
# they are looped over and each one is executed individualy. 
constraints = ["filter = \'%s\'"%'r']
m1 = configureMetric('RmsMetric', kwargs={'col':'finSeeing'}, plotDict={'colorMin':0., 'colorMax':0.6})
m2 =  configureMetric('RobustRmsMetric', kwargs={'col':'finSeeing'}, plotDict={'colorMin':0., 'colorMax':0.6})
metricDict = makeDict(m1,m2)
slicer = configureSlicer('HealpixSlicer',
                          kwargs={"nside":nside,'spatialkey1':"fieldRA", 'spatialkey2':"fieldDec"},
                          metricDict = metricDict, constraints=constraints)
slicerList.append(slicer)



# Example of doing summary stats:
nside=64
constraints = ["filter = \'%s\'"%'r']
m2 = configureMetric('Coaddm5Metric', plotDict={'zp':27., 'percentileClip':95,
                                                'units':'Co-add m5 - %.1f'%27.},
                     kwargs={'metricName':'coadd_w_summary'},
                     summaryStats={'MeanMetric':{}, 'MinMetric':{}, 'MaxMetric':{}, 'RmsMetric':{}} )          
metricDict = makeDict(m2)
slicer = configureSlicer('HealpixSlicer',
                          kwargs={"nside":nside,'spatialkey1':"fieldRA", 'spatialkey2':"fieldDec"},
                          metricDict = metricDict, constraints=constraints)
root.slicers=makeDict(slicer)
slicerList.append(slicer)


# Example of merging histograms
filters = ['u','g','r','i','z','y']
colors={'u':'m','g':'b','r':'g','i':'y','z':'r','y':'k'}
for i,f in enumerate(filters):
    m1 = configureMetric('CountMetric', kwargs={'col':'airmass'},
                         histMerge={'histNum':1, 'legendloc':'upper right', 'color':colors[f],'label':'%s'%f} )
    slicer = configureSlicer('OneDSlicer', kwargs={"sliceColName":'airmass'},
                             metricDict=makeDict(m1), constraints=["filter = '%s'"%f])
    slicerList.append(slicer)



root.slicers=makeDict(*slicerList)



