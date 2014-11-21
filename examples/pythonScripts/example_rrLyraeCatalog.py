import numpy
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.db as db

from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import *

from lsst.sims.photUtils import PhotometryStars
from lsst.sims.catalogs.measures.instance.fileMaps import defaultSpecMap

#initialize baseline (i.e. not variable) photometry
photObj = PhotometryStars()
bandPassList = ['u', 'g', 'r', 'i', 'z', 'y']
photObj.loadBandPassesFromFiles(bandPassList)
photObj.setupPhiArray_dict()

rrLyraeDB = CatalogDBObject.from_objid('rrlystars')
obs_metadata = ObservationMetaData(unrefractedRA=0.0, unrefractedDec=0.0,
                                   boundType='circle', boundLength=4.0)

colNames = ['raJ2000', 'decJ2000','varParamStr','magNorm','sedFilename']
dtype = numpy.dtype([('id',int),('raJ2000',float),('decJ2000',float),
                     ('varParamStr',str,256),('magNorm',float),('sedFilename',str,100)])

rrLyraeDB.dtype = dtype

results=rrLyraeDB.query_columns(colnames=colNames, obs_metadata=obs_metadata,
                                returnRecArray=True)
for chunk in results:
    sedNames = chunk['sedFilename']
    magNorms = chunk['magNorm']
    sedList = photObj.loadSeds(sedNames, magNorm=magNorms, specFileMap=defaultSpecMap)
    baselineMagnitudes = photObj.calculate_magnitudes(sedList)


# Connect to opsim
dbAddress = 'sqlite:///../../../garage/OpSimData/opsimblitz2_1060_sqlite.db'
oo = db.OpsimDatabase(dbAddress)
colnames = ['expMJD', 'fieldRA', 'fieldDec']
sqlconstraint ='filter="r"'
# Get opsim simulation data
simdata = oo.fetchMetricData(colnames, sqlconstraint)
# Init the slicer, set 2 points
slicer = slicers.UserPointsSlicer(ra=[0., .1], dec=[0., -.1])
# Setup slicer (builds kdTree)
slicer.setupSlicer(simdata)
# Slice Point for index zero
ind = slicer._sliceSimData(0)
expMJDs = simdata[ind['idxs']]['expMJD']
# Find the expMJDs for the 2nd point
ind = slicer._sliceSimData(1)
expMJDs = simdata[ind['idxs']]['expMJD']

