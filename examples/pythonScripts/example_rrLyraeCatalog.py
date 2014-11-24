import numpy
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.db as db

from lsst.sims.catalogs.generation.db import CatalogDBObject, ObservationMetaData
from lsst.sims.catUtils.baseCatalogModels import *

from lsst.sims.photUtils import PhotometryStars, Variability
from lsst.sims.catalogs.measures.instance.fileMaps import defaultSpecMap

class LightCurveGenerator(object):
    """
    This class reads in a database of variable stars and an opsim cadence simulation.
    It will write out light curves in [u,g,r,i,z,y] for all of the objects in the
    stellar database based on when the star was visible to the opsim cadence.
    """

    def __init__(self, address=None, filters=['u','g','r','i','z','y'], stellarDBtype='rrlystars'):
        """
        param [in] address is a string containing the connection string for the opsim database

        param [in] stellarDBtype is a string indicating what kind of CatalogDBobject containing
        variable stars should be used.  This should be one of the objid's declared for one of
        the CatalogDBObject classes defined in
        sims_catUtils/python/lsst/sims/catUtils/baseCatalogModels/StarModels.py

        param [in] filters is a list of filters for which to generate light curves
        """

        self.filters = filters
        if address is None:
            raise RuntimeError('must specify address for the OpSim database in LightCurveGenerator')
        self.dbAddress = address
        self.stellarDBtype = stellarDBtype

    def _getPointings(self, filterName):
        """
        Query the opsim database for all of the pointings made in the given
        filter name.

        param [in] filterName is a string indicating which filter is of interest
        """
        self.opDB = db.OpsimDatabase(self.dbAddress)
        colnames = ['expMJD', 'fieldRA', 'fieldDec']
        sqlconstraint = 'filter="'+filterName+'"'
        pointings = self.opDB.fetchMetricData(colnames, sqlconstraint)
        return pointings

    def _initializePhotometry(self):
        """
        Setup the objects that will do the photometric calculations
        """

        #initialize baseline (i.e. not variable) photometry
        self.photObj = PhotometryStars()
        self.photObj.loadBandPassesFromFiles(self.filters)
        self.photObj.setupPhiArray_dict()

        #initialize the object that will calculate variable magnitudes
        self.varObj = Variability()

        #query the database of variable stars
        stellarDB = CatalogDBObject.from_objid(self.stellarDBtype)
        obs_metadata = ObservationMetaData(unrefractedRA=0.0, unrefractedDec=0.0,
                                   boundType='circle', boundLength=2.0)

        colNames = ['raJ2000', 'decJ2000','varParamStr','magNorm','sedFilename','distance']

        #this dtype is necessary, otherwise varParamStr gets clipped too short and is useless
        dtype = numpy.dtype([('id',int),('raJ2000',float),('decJ2000',float),
                           ('varParamStr',str,256),('magNorm',float),('sedFilename',str,100),
                           ('distance',float)])

        stellarDB.dtype = dtype

        self.stellarResults = stellarDB.query_columns(colnames=colNames, obs_metadata=obs_metadata,
                                                      returnRecArray=True)


    def _writeFilter(self, iFilter, sedList=None, ra=None, dec=None, baselineMagnitudes=None, varParamStr=None, objectNames=None):
        """
        write the light curves for the specified filter

        param [in] iFilter is an int indicating which filter to write light curves for (it is an index for the list
        self.filters)

        param [in] sedList is a list of Sed objects representing the stars whose magnitudes we want

        param [in ] ra is a list of RAs for the stars whose magnitudes we want (in radians)

        param [in] dec is a list of Decs for the stars whose magnitudes we want (in radians)

        param [in] baselineMagnitudes is a list of zero-point (non-variable) magnitudes for the stars

        param [in] varParamStr is a list of strings containing the parameters for the variability model
        (this is read from the database of variable stars)

        param [in] objectNames is a list of unique identifiers for the stars

        This method will write a series of light curves light_curve_##_ff.txt where ## is an integer indexing
        the star and ff denotes the filter.
        """

        #first get all of the opsim pointings for the filter in question
        opsimData = self._getPointings(self.filters[iFilter])

        #Now we use the MAF infrastructure to set up a slicer that, for each object
        #will tell us which opsim pointings saw the object
        #
        # Init the slicer, set 2 points
        slicer = slicers.UserPointsSlicer(ra=ra, dec=dec)
        # Setup slicer (builds kdTree)
        slicer.setupSlicer(opsimData)

        #loop over objects
        for ii,vps in enumerate(varParamStr):
            name = objectNames[ii]
            outputName = 'light_curve_'+str(ii)+'_'+self.filters[iFilter]+'.txt'
            outputFile = open(outputName,'w')

            # Find the opsim pointings that saw the object
            ind = slicer._sliceSimData(ii)
            expMJDs = opsimData[ind['idxs']]['expMJD']
            expMJDs.sort()
            for expmjd in expMJDs:


                if self.cachedMJD is not None:
                    iclosest = numpy.fabs(self.cachedMJD - expmjd).argmin()
                else:
                    iclosest = None

                if iclosest is not None:
                    if numpy.fabs(self.cachedMJD[iclosest] - expmjd) > 1.0e-6:
                        iclosest = None

                if iclosest is None:
                    vv = self.varObj.calculate_stellar_variability(
                                       u0=[baselineMagnitudes[ii][0]],
                                       g0=[baselineMagnitudes[ii][1]],
                                       r0=[baselineMagnitudes[ii][2]],
                                       i0=[baselineMagnitudes[ii][3]],
                                       z0=[baselineMagnitudes[ii][4]],
                                       y0=[baselineMagnitudes[ii][5]],
                                       magNorm=[0.0],
                                       varParams = [vps],
                                       expmjd=expmjd
                                       )

                    mm = vv[iFilter][0]

                    if self.cachedMJD is not None:
                         self.cachedMJD = numpy.append(self.cachedMJD, expmjd)
                    else:
                         self.cachedMJD = numpy.array([expmjd])

                    subList = []
                    for jj in range(len(self.filters)):
                         subList.append(vv[jj][0])
                    self.cachedMagnitudes.append(subList)
                else:
                    mm = self.cachedMagnitudes[iclosest][iFilter]


                outputFile.write("%.7f %f\n" % (expmjd,mm))

            outputFile.close()


    def _writeChunk(self, stellarChunk):
        self.cachedMJD = None
        self.cachedMagnitudes = []
        sedNames = stellarChunk['sedFilename']
        magNorms = stellarChunk['magNorm']
        sedList = self.photObj.loadSeds(sedNames, magNorm=magNorms, specFileMap=defaultSpecMap)
        baselineMagnitudes = self.photObj.calculate_magnitudes(sedList)
        for iFilter in range(len(self.filters)):
            self._writeFilter(iFilter,sedList=sedList, baselineMagnitudes=baselineMagnitudes,
                              ra=stellarChunk['raJ2000'], dec=stellarChunk['decJ2000'],
                              objectNames=stellarChunk['id'], varParamStr=stellarChunk['varParamStr'])

    def writeLightCurves(self):
        self._initializePhotometry()
        for chunk in self.stellarResults:
            self._writeChunk(chunk)


myLC = LightCurveGenerator(address = 'sqlite:///../../../garage/OpSimData/opsimblitz2_1060_sqlite.db')
myLC.writeLightCurves()
