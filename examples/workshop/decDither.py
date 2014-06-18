# Driver for running a declination only dither scheme

from lsst.sims.maf.driver.mafConfig import configureBinner, configureMetric, makeDict
root.outputDir = './DecDith'
root.dbAddress = {'dbAddress':'sqlite:///opsimblitz2_1060_sqlite.db'}
root.opsimName = 'Example'


binList = []

metric = configureMetric('Coaddm5Metric')
binner = configureBinner('HealpixBinner', metricDict=makeDict(metric),
                          constraints=['filter = "r"'])

binList.append(binner)

metric = configureMetric('Coaddm5Metric', kwargs={'metricName':'m5_decdith'})
binner = configureBinner('HealpixBinner', metricDict=makeDict(metric),
                          constraints=['filter = "r"'], kwargs={'spatialkey1':'fieldRA', 'spatialkey2':'decOnlyDither'})
binList.append(binner)



root.binners = makeDict(*binList)
