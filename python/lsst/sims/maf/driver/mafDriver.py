import os
import numpy as np
import matplotlib.pyplot as plt
from .mafConfig import config2dict, readMetricConfig, readSlicerConfig, readMixConfig

import lsst.sims.maf.db as db
import lsst.sims.maf.slicers as slicers
import lsst.sims.maf.metrics as metrics
import lsst.sims.maf.sliceMetrics as sliceMetrics
import lsst.sims.maf.utils as utils
import lsst.sims.maf.stackers as stackers
import lsst.sims.maf.maps as maps
import time
import collections


def dtime(time_prev):
   return (time.time() - time_prev, time.time())

class MafDriver(object):
    """Script for configuring and running metrics on Opsim output """

    def __init__(self, configvalues):
        """Load up the configuration and set the slicer and metric lists """
        # Configvalues passed from runDriver.py
        self.config = configvalues
        # Make sure only legitimate keys are in the dbAddress
        for key in self.config.dbAddress.keys():
           if key not in ['dbAddress','dbClass']:
              raise ValueError('key value "%s" not valid for dbAddress dict.  Must be "dbAddress" or "dbClass".'%key)
        # If a dbClass isn't specified, use OpsimDatabase
        if 'dbClass' not in self.config.dbAddress.keys():
           self.config.dbAddress['dbClass'] = 'OpsimDatabase'

        # Validate and freeze the config
        self.config.validate()
        self.config.freeze()

        # Check for output directory, make if needed.
        if not os.path.isdir(self.config.outDir):
            os.makedirs(self.config.outDir)

        self.verbose = self.config.verbose
        self.figformat = self.config.figformat
        self.dpi = self.config.dpi
        self.plotOnly = self.config.plotOnly

        # Import any additional modules specified by the user.
        utils.moduleLoader(self.config.modules)

        # Set up database connection.
        self.opsimdb = db.Database.getClass(self.config.dbAddress['dbClass'])(self.config.dbAddress['dbAddress'])

        time_prev = time.time()
        self.time_start = time.time()
        # Grab config info and write to disk.
        if self.config.getConfig:
            configSummary, configDetails = self.opsimdb.fetchConfig()
            f = open(os.path.join(self.config.outDir,'configSummary.txt'), 'w')
            utils.outputUtils.printDict(configSummary, 'Config Summary', filehandle=f)
            f.close()
            f = open(os.path.join(self.config.outDir, 'configDetails.txt'), 'w')
            utils.outputUtils.printDict(configDetails, 'Config Details', filehandle=f)
            f.close()
            if self.verbose:
                dt, time_prev = dtime(time_prev)
                print 'Got OpSim config info in %.3g s'%dt

        # Get proposal information (for OpSim databases).
        if self.config.dbAddress['dbClass'] == 'OpsimDatabase':
           self.propids, self.propTags = self.opsimdb.fetchPropInfo()
           if self.verbose:
               dt, time_prev = dtime(time_prev)
               print 'fetched PropID info in %.3g s'%dt
        else:
           self.propids, self.propTags = ({}, {})

        # Construct the slicers and metric objects
        self.slicerList = []
        self.metricList = []
        for i,slicer in self.config.slicers.iteritems():
            name, kwargs, metricDict, constraints, stackerDict, mapsDict, metadata, metadataVerbatim = \
                 readSlicerConfig(slicer)
            temp_slicer = slicers.BaseSlicer.getClass(name)(**kwargs )
            temp_slicer.constraints = slicer.constraints
            temp_slicer.table = slicer.table
            #check that constraints in slicer are unique
            if len(temp_slicer.constraints) > len(set(temp_slicer.constraints)):
                print 'Slicer %s has repeated constraints' %slicer.name
                print 'Constraints:  ', slicer.constraints
                raise Exception('Slicer constraints are not unique')
            temp_slicer.metadata = metadata
            temp_slicer.metadataVerbatim = metadataVerbatim
            temp_slicer.index = i
            stackersList = []
            for key in stackerDict.keys():
               stackername, kwargs = config2dict(stackerDict[key])
               stackersList.append(stackers.BaseStacker.getClass(stackername)(**kwargs))
            temp_slicer.stackers = stackersList
            mapsList = []
            mapsNames = []
            for key in mapsDict.keys():
               mapName, kwargs = config2dict(mapsDict[key])
               mapsList.append(maps.BaseMap.getClass(mapName)(**kwargs) )
               mapsNames.append(mapName)
            temp_slicer.mapsList = mapsList
            temp_slicer.mapsNames = mapsNames
            self.slicerList.append(temp_slicer)
            sub_metricList=[]
            for metric in slicer.metricDict.itervalues():
                name, kwargs, plotDict, summaryStats, histMerge, displayDict = readMetricConfig(metric)
                # Add plot parameters and display parameters to kwargs handed to metric.
                kwargs['plotDict'] = plotDict
                kwargs['displayDict'] = displayDict
                temp_metric = metrics.BaseMetric.getClass(name)(**kwargs)
                # Add an attribute to our metric which will describe the summary stats.
                temp_metric.summaryStats = []
                nameCheck=[]
                for key in summaryStats.keys():
                    summarykwargs = readMixConfig(summaryStats[key])
                    summaryMetric = metrics.BaseMetric.getClass(key.split(' ')[0])(col='metricdata', **summarykwargs)
                    temp_metric.summaryStats.append(summaryMetric)
                    nameCheck.append(summaryMetric.name)
                if len(list(set(nameCheck))) < len(nameCheck):
                   duplicates = [x for x, y in collections.Counter(nameCheck).items() if y > 1]
                   raise Exception('Summary metric names not unique. "%s" defined more than one with metric "%s"'
                                   %(duplicates[0], temp_metric.name))
                # If it is a UniSlicer, make sure at least the IdentityMetric is run
                if temp_slicer.slicerName == 'UniSlicer':
                    if len(summaryStats) == 0:
                        temp_metric.summaryStats.append(metrics.BaseMetric.registry['IdentityMetric']('metricdata'))
                temp_metric.histMerge = histMerge
                sub_metricList.append(temp_metric )
            self.metricList.append(sub_metricList)
        # Make a unique list of all SQL constraints
        self.constraints = []
        for b in self.slicerList:
            for c in b.constraints:
                self.constraints.append(c)
        self.constraints = list(set(self.constraints))
        # Make a unique list of tables
        self.tables = list(set([s.table for s in self.slicerList]))
        # Check that all filenames will be unique
        filenames=[]
        for i,slicer in enumerate(self.slicerList):
            for constraint in slicer.constraints:
                for metric in self.metricList[i]:
                    # Approximate what output filename will be
                    if slicer.metadataVerbatim:
                        comment = slicer.metadata
                    else:
                        comment = constraint.replace('=','').replace('filter','').replace("'",'')
                        comment = comment.replace('"', '').replace('  ',' ') + ' ' + slicer.metadata
                    filenames.append('_'.join([metric.name, comment, slicer.slicerName]))
        if len(filenames) != len(set(filenames)):
            duplicates = list(set([x for x in filenames if filenames.count(x) > 1]))
            counts = [filenames.count(x) for x in duplicates]
            print ['%s: %d versions' %(d, c) for d, c in zip(duplicates, counts)]
            raise Exception('Filenames for metrics will not be unique.  Add slicer metadata or change metric names.')

    def getData(self, constraint, colnames=[], stackersList=[], table=None):
        """Pull required data from database and calculate additional columns from stackers. """
        # Stacker_names describe the already-configured (via the config driver) stacker methods.
        stacker_names = [s.__class__.__name__ for s in stackersList ]
        dbcolnames = []
        sourceLookup = utils.getColInfo.ColInfo()
        # Go through all columns that the metrics need.
        for colname in colnames:
            source = sourceLookup.getDataSource(colname)
            # If data source of column is a stacker:
            if source != sourceLookup.defaultDataSource:
                stacker = source()
                for col in stacker.colsReq:
                    # Add column names that the stackers need.
                    dbcolnames.append(col)
                # If not already a configured stacker, instantiate one using defaults
                if stacker.__class__.__name__ not in stacker_names:
                    stackersList.append(stacker)
                    stacker_names.append(stacker.__class__.__name__)
            # Else if data source is just the usual database:
            else:
                dbcolnames.append(colname)
        # Remove duplicates from list of columns required from database.
        dbcolnames=list(set(dbcolnames))
        # Get the data from database.
        if (table is not None)  & (table != 'Summary'):
           self.data = self.opsimdb.fetchMetricData(sqlconstraint=constraint,colnames=colnames,
                                                    distinctExpMJD=False, groupBy=None,
                                                    tableName=table)
        else:
           self.data = self.opsimdb.fetchMetricData(sqlconstraint=constraint,
                                                    colnames=dbcolnames)
        # Calculate the data from stackers.
        for stacker in stackersList:
            self.data = stacker.run(self.data)
        # Done - self.data should now have all required columns.


    def getFieldData(self, slicer, sqlconstraint):
        """Given an opsim slicer, generate the FieldData """
        # Do a bunch of parsing to get the propids out of the sqlconstraint.
        if 'propID' not in sqlconstraint:
            propids = self.propids.keys()
        else:
            # example sqlconstraint: filter = r and (propid = 219 or propid = 155) and propid!= 90
            sqlconstraint = sqlconstraint.replace('=', ' = ').replace('(', '').replace(')', '')
            sqlconstraint = sqlconstraint.replace("'", '').replace('"', '')
            # Allow for choosing all but a particular proposal.
            sqlconstraint = sqlconstraint.replace('! =' , ' !=')
            sqlconstraint = sqlconstraint.replace('  ', ' ')
            sqllist = sqlconstraint.split(' ')
            propids = []
            nonpropids = []
            i = 0
            while i < len(sqllist):
                if sqllist[i].lower() == 'propid':
                    i += 1
                    if sqllist[i] == "=":
                        i += 1
                        propids.append(int(sqllist[i]))
                    elif sqllist[i] == '!=':
                        i += 1
                        nonpropids.append(int(sqllist[i]))
                i += 1
            if len(propids) == 0:
                propids = self.propids.keys()
            if len(nonpropids) > 0:
                for nonpropid in nonpropids:
                    if nonpropid in propids:
                        propids.remove(nonpropid)
        # And query the field Table.
        if 'Field' in self.opsimdb.tables:
            self.fieldData = self.opsimdb.fetchFieldsFromFieldTable(propids)
        else:
            fieldID, idx = np.unique(self.data[slicer.simDataFieldIDColName], return_index=True)
            ra = self.data[slicer.fieldRaColName][idx]
            dec = self.data[slicer.fieldDecColName][idx]
            self.fieldData = np.core.records.fromarrays([fieldID, ra, dec],
                                               names=['fieldID', 'fieldRA', 'fieldDec'])


    def run(self):
        """Loop over each slicer and calculate metrics for that slicer. """

        # Loop through all sqlconstraints, and run slicers + metrics that match the same sql constraints
        #   (so we only have to do one query of database per sql constraint).

        # XXX -- add a check here to make sure tables match.
        for table in self.tables:
           for sqlconstraint in self.constraints:
               # Find which slicers have an exactly matching constraint
               matchingSlicers=[]
               slicerNames=[]
               for b in self.slicerList:
                  if b.table == table:
                     if sqlconstraint in b.constraints:
                         matchingSlicers.append(b)
                         slicerNames.append(b.slicerName)
               if len(matchingSlicers) > 0:
                  # And for those slicers, find the data columns required.
                  colnames=[]
                  stackersList = []
                  for slicer in matchingSlicers:
                      for m in self.metricList[slicer.index]:
                          for cn in m.colNameArr:
                              colnames.append(cn)
                      for cn in slicer.columnsNeeded:
                          colnames.append(cn)
                      for stacker in slicer.stackers:
                          stackersList.append(stacker)
                          for col in stacker.colsReq:
                              colnames.append(col)
                  # Find the unique column names required.
                  colnames = list(set(colnames))
                  if not self.plotOnly:
                     print 'Querying with SQLconstraint:', sqlconstraint, ' from table:', table
                     # Get the data from the database + stacker calculations.
                     if self.verbose:
                         time_prev = time.time()
                     self.getData(sqlconstraint, colnames=colnames, stackersList=stackersList, table=table)
                     if self.verbose:
                         dt, time_prev = dtime(time_prev)
                     if len(self.data) == 0:
                         print '  No data matching constraint:   %s'%sqlconstraint
                  else:
                     # Set some dummy data if we are going to restore later
                     self.data = [0]

                  # Got data, now set up slicers.
                  if len(self.data) > 0:
                      if self.verbose:
                          print '  Found %i matching visits in %.3g s'%(len(self.data),dt)
                      else:
                          print '  Found %i matching visits' %(len(self.data))
                      # Special data requirements for opsim slicer.
                      self.fieldData = None
                      if 'OpsimFieldSlicer' in slicerNames and not self.plotOnly:
                          self.getFieldData(matchingSlicers[slicerNames.index('OpsimFieldSlicer')], sqlconstraint)
                      # Setup each slicer, and run through the slicepoints (with metrics) in baseSliceMetric
                      if self.verbose:
                          time_prev = time.time()
                      for slicer in matchingSlicers:

                          # Set up any additional maps
                          for m in self.metricList[slicer.index]:
                             for skyMap in m.maps:
                                if skyMap not in slicer.mapsNames:
                                   slicer.mapsList.append(maps.BaseMap.getClass(skyMap)())
                          gm = sliceMetrics.RunSliceMetric(figformat=self.figformat, dpi=self.dpi,
                                                           outDir=self.config.outDir)
                          gm._setSlicer(slicer)
                          gm._setMetrics(self.metricList[slicer.index])
                          # Make a more useful metadata comment.
                          if slicer.metadataVerbatim:
                              metadata = slicer.metadata
                          else:
                              metadata = sqlconstraint.replace('=','').replace('filter','').replace("'",'')
                              metadata = metadata.replace('"', '').replace('  ',' ') + ' '+ slicer.metadata

                          if self.plotOnly:
                             iids = gm.metricNames.keys()
                             newGm = sliceMetrics.RunSliceMetric(figformat=self.figformat, dpi=self.dpi,
                                                                 outDir=self.config.outDir,
                                                                 useResultsDb=False)
                             newGm._setSlicer(slicer)
                             restoredData = False
                             for iid in iids:
                                gm.simDataNames[iid] = self.config.opsimName
                                gm.metadatas[iid] = metadata
                                filename = gm._buildOutfileName(iid)
                                # Load all the metric data back in
                                fullFile = os.path.join(self.config.outDir, filename+'.npz')
                                if os.path.isfile(fullFile):
                                   print 'Restoring %s'%fullFile
                                   newGm.readMetricData(fullFile)
                                   # Set the filename as a property of each metric (for merged histograms)
                                   gm.metricObjs[iid].saveFile = fullFile
                                   # Set the slicer to the newly restored slicer
                                   newGm._setSlicer(newGm.slicers[iid], override=True)
                                   # Replace the restored plotting parameters
                                   newGm.plotDicts[iid] = gm.plotDicts[iid]
                                   newGm.displayDicts[iid] = gm.displayDicts[iid]
                                   restoredData = True
                             # Replot, note we are not saving the updated plotDicts to save time.
                             if restoredData:
                                newGm.plotAll(savefig=True, closefig=True, verbose=True)
                          else:
                             # Run through slicepoints in slicer, and calculate metric values.
                             print '    running slicerName =', slicer.slicerName, \
                            ' run metrics:', ', '.join([m.name for m in self.metricList[slicer.index]])
                             gm.runSlices(self.data, simDataName=self.config.opsimName,
                                          metadata=metadata, sqlconstraint=sqlconstraint,
                                          fieldData=self.fieldData, maps=slicer.mapsList)
                             if self.verbose:
                                dt,time_prev = dtime(time_prev)
                                print '    Computed metrics in %.3g s'%dt
                             # And run reduce methods for relevant metrics.
                             gm.reduceAll()
                             # And write metric data files to disk.
                             gm.writeAll()
                             # Add the metric filenames to the metric objects (for merged histograms).
                             for iid in gm.metricObjs:
                                filename = gm._buildOutfileName(iid)
                                # Load all the metric data back in
                                fullFile = os.path.join(self.config.outDir, filename+'.npz')
                                gm.metricObjs[iid].saveFile = fullFile
                             # And plot all metric values.
                             gm.plotAll(savefig=True, closefig=True, verbose=True)
                             if self.verbose:
                                dt,time_prev = dtime(time_prev)
                                print '    plotted metrics in %.3g s'%dt
                             # Loop through the metrics and calculate any summary statistics
                             for i, metric in enumerate(self.metricList[slicer.index]):
                                 if hasattr(metric, 'summaryStats'):
                                     for stat in metric.summaryStats:
                                         # If it's metric returning an OBJECT, run summary stats on
                                         # each reduced metric
                                         # (have to identify related reduced metric values first)
                                         if metric.metricDtype == 'object':
                                             iid = gm.getMetricObjIid(metric)[0]
                                             baseName = gm.metricNames[iid]
                                             all_names = gm.metricNames.values()
                                             matching_metrics = [x for x in all_names \
                                                                 if x[:len(baseName)] == baseName \
                                                                 and x != baseName]
                                             for mm in matching_metrics:
                                                 iid = gm.findIids(metricName=mm)[0]
                                                 summary = gm.computeSummaryStatistics(iid, stat)
                                         # Else it's a simple metric value.
                                         else:
                                             iid = gm.findIids(metricName=metric.name)[0]
                                             summary = gm.computeSummaryStatistics(iid, stat)
                             if self.verbose:
                                dt,time_prev = dtime(time_prev)
                                print '    Computed summarystats in %.3g s'%dt
                             if self.verbose:
                                dt,time_prev = dtime(time_prev)
                                print '    wrote output files in %.3g s'%dt

        # Create any 'merge' histograms that need merging.
        # Loop through all the metrics and find which histograms need to be merged
        histList = []
        for m1 in self.metricList:
            for m in m1:
                if 'histNum' in m.histMerge.keys():
                    histList.append(m.histMerge['histNum'])

        histList = list(set(histList))
        histList.sort()
        histDict={}
        for item in histList:
            histDict[item] = {}
            histDict[item]['files']=[]
            histDict[item]['plotkwargs']=[]

            for m1 in self.metricList:
                for m in m1:
                    if 'histNum' in m.histMerge.keys():
                        # Determine which merged histogram to put data into.
                        key = m.histMerge['histNum']
                        if hasattr(m,'saveFile') and key in histDict.keys():
                            histDict[key]['files'].append(m.saveFile)
                            temp_dict = m.histMerge
                            del temp_dict['histNum']
                            histDict[key]['plotkwargs'].append(temp_dict)

        if self.plotOnly:
           useResultsDb = False
        else:
           useResultsDb = True
        for key in histDict.keys():
            # Use a comparison slice metric per merged histogram. Only read relevant files.
            cbm = sliceMetrics.ComparisonSliceMetric(useResultsDb=useResultsDb, outDir=self.config.outDir,
                                                     figformat=self.figformat, dpi=self.dpi)
            if len(histDict[key]['files']) > 0:
                for filename in histDict[key]['files']:
                    if self.verbose:
                       print 'reading %s to make merged histogram'%fullfilename
                    cbm.readMetricData(filename)
                iids = cbm.metricValues.keys()
                fignum, title, histfile = cbm.plotHistograms(iids, savefig=True,
                                                            plotkwargs=histDict[key]['plotkwargs'])
                plt.close('all')
                if cbm.slicers[iids[0]].slicerName == 'HealpixSlicer':
                   fignum, title, psfile = cbm.plotPowerSpectra(iids, savefig=True,
                                                                plotkwargs=histDict[key]['plotkwargs'])
                   plt.close('all')

        if not self.plotOnly:
           today_date, versionInfo = utils.getDateVersion()
           # Open up a file and print the results of verison and date.
           datefile = open(self.config.outDir+'/'+'date_version_ran.dat','w')
           print >>datefile, 'date, version, fingerprint '
           print >>datefile, '%s,%s,%s'%(today_date,versionInfo['__version__'],
                                         versionInfo['__fingerprint__'])
           datefile.close()
           # Save the as-ran pexConfig file
           self.config.save(self.config.outDir+'/'+'maf_config_asRan.py')

        if self.verbose:
            dt,self.time_start = dtime(self.time_start)
            print 'Ran everything in %.3g seconds' %(dt)
