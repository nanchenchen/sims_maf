# Here is an example of a very very simple MAF configuration driver script
# to run:
# runDriver.py very_simple_cfg.py

# This script uses the LSST pex_config.  This is executed as a python script, but only things that start with 'root.' are passed on to the driver script.

# Import MAF helper functions 
from lsst.sims.maf.driver.mafConfig import configureBinner, configureMetric, makeDict

# Set the output directory
root.outputDir = './Very_simple_out'
# Set the database to use (the example db included in the git repo)
root.dbAddress = {'dbAddress':'sqlite:///opsimblitz2_1060_sqlite.db'}
# Name of this run (filename base)
root.opsimName = 'VerySimpleExample'

# Make an empty list to hold all the binner configs
binList = []

# Make a set of SQL where constraints to only use each filter
constraints = ["filter = '%s'"%f for f in ['u','g','r','i','z','y'] ]

# Run 2 metrics, the mean seeing and the co-added 5-sigma limiting depth.
metric1 = configureMetric('MeanMetric', params=['finSeeing'])
metric2 = configureMetric('Coaddm5Metric')

# Configure a binner.  Use the Healpix binner to make sky maps and power spectra.
binner = configureBinner('HealpixBinner', metricDict=makeDict(metric1,metric2),
                          constraints=constraints)
binList.append(binner)

metric = configureMetric('MeanMetric', params=['finSeeing'],
                          summaryStats={'IdentityMetric':{}})
# Configure a binner.  Use the UniBinner to simply take all the data.
binner = configureBinner('UniBinner', metricDict=makeDict(metric),
                          constraints=constraints)
binList.append(binner)



root.binners = makeDict(*binList)

