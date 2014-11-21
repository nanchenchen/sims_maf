import numpy
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.db as db

from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import *

from lsst.sims.photUtils import PhotometryStars, Variability
from lsst.sims.catalogs.measures.instance.fileMaps import defaultSpecMap

# Connect to opsim
dbAddress = 'sqlite:///../../../garage/OpSimData/opsimblitz2_1060_sqlite.db'
oo = db.OpsimDatabase(dbAddress)
colnames = ['expMJD', 'fieldRA', 'fieldDec']
sqlconstraint ='filter="r"'
# Get opsim simulation data
simdata = oo.fetchMetricData(colnames, sqlconstraint)

#initialize baseline (i.e. not variable) photometry
photObj = PhotometryStars()
bandPassList = ['u', 'g', 'r', 'i', 'z', 'y']
photObj.loadBandPassesFromFiles(bandPassList)
photObj.setupPhiArray_dict()

varObj = Variability()

rrLyraeDB = CatalogDBObject.from_objid('rrlystars')
obs_metadata = ObservationMetaData(unrefractedRA=0.0, unrefractedDec=0.0,
                                   boundType='circle', boundLength=2.0)

colNames = ['raJ2000', 'decJ2000','varParamStr','magNorm','sedFilename','distance']
dtype = numpy.dtype([('id',int),('raJ2000',float),('decJ2000',float),
                     ('varParamStr',str,256),('magNorm',float),('sedFilename',str,100),
                     ('distance',float)])

rrLyraeDB.dtype = dtype

results = rrLyraeDB.query_columns(colnames=colNames, obs_metadata=obs_metadata,
                                returnRecArray=True)

for metaData in results:
    sedNames = metaData['sedFilename']
    magNorms = metaData['magNorm']
    sedList = photObj.loadSeds(sedNames, magNorm=magNorms, specFileMap=defaultSpecMap)
    baselineMagnitudes = photObj.calculate_magnitudes(sedList)

    # Init the slicer, set 2 points
    slicer = slicers.UserPointsSlicer(ra=metaData['raJ2000'], dec=metaData['decJ2000'])
    # Setup slicer (builds kdTree)
    slicer.setupSlicer(simdata)
    
    #loop over objects
    for ii,vps in enumerate(metaData['varParamStr']):
        
        print vps
        
        outputName = 'rrly_'+str(metaData[ii]['id'])+'light_curve.txt'
        outputFile = open(outputName,'w')
        # Slice Point for index zero
        ind = slicer._sliceSimData(ii)
        expMJDs = simdata[ind['idxs']]['expMJD']
        expMJDs.sort()
        for expmjd in expMJDs:
            vv = varObj.calculate_stellar_variability(
                               u0=[baselineMagnitudes[ii][0]],
                               g0=[baselineMagnitudes[ii][1]],
                               r0=[baselineMagnitudes[ii][2]],
                               i0=[baselineMagnitudes[ii][3]],
                               z0=[baselineMagnitudes[ii][4]],
                               y0=[baselineMagnitudes[ii][5]],
                               magNorm=metaData[ii]['magNorm'],
                               varParams = [vps],
                               expmjd=expmjd
                               )
        
            outputFile.write("%.7f %f %f %f %f %f %f\n" %
                            (expmjd,vv[0][0],vv[1][0],vv[2][0],
                            vv[3][0],vv[4][0],vv[5][0]))
            
        outputFile.close()
