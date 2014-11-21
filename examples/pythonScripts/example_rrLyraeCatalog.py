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
    
    for i,vps in enumerate(metaData['varParamStr']):
        # Slice Point for index zero
        ind = slicer._sliceSimData(i)
        expMJDs = simdata[ind['idxs']]['expMJD']
        print expMJDs
        
    
    
    
    # Find the expMJDs for the 2nd point
    ind = slicer._sliceSimData(1)
    expMJDs = simdata[ind['idxs']]['expMJD']

