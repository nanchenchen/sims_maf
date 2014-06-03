# Example of more complex metric 
# Takes multiple columns of data (although 'night' could be calculable from 'expmjd')
# Returns variable length array of data
# Uses multiple reduce functions

import numpy as np
from .complexMetrics import ComplexMetric

class VisitPairsMetric(ComplexMetric):
    """Count the number of pairs of visits per night within deltaTmin and deltaTmax."""
    def __init__(self, timesCol='expMJD', nightsCol='night', metricName='VisitPairs',
                 deltaTmin=15.0/60.0/24.0, deltaTmax=90.0/60.0/24.0, nPairs=2, window=30,
                 **kwargs):
        """Instantiate metric.
        
        'timesCol' = column with the time of the visit (default expmjd),
        'nightsCol' = column with the night of the visit (default night),
        'deltaTmin' = minimum time of window,
        'deltaTmax' = maximum time of window."""
        self.times = timesCol   
        self.nights = nightsCol
        eps = 1e-10
        self.deltaTmin = deltaTmin - eps
        self.deltaTmax = deltaTmax
        self.nPairs = nPairs
        self.window = window
        super(VisitPairsMetric, self).__init__([self.times, self.nights],
                                               metricName=metricName, **kwargs)

    def run(self, dataSlice):
        # Identify the nights with any visits.
        uniquenights = np.unique(dataSlice[self.nights])
        nights = []
        visitPairs = []
        # Identify the nights with pairs of visits within time window.
        for i, n in enumerate(uniquenights):
            condition = (dataSlice[self.nights] == n)
            times = dataSlice[self.times][condition]
            pairnum = 0
            for t in times:
                dt = times - t
                condition2 = ((dt >= self.deltaTmin) & (dt <= self.deltaTmax))
                pairnum += len(dt[condition2])
            if pairnum > 0:
                visitPairs.append(pairnum)
                nights.append(n)
        # Convert to numpy arrays.
        visitPairs = np.array(visitPairs)
        nights = np.array(nights)
        if len(visitPairs) == 0:
            return self.badval
        return (visitPairs, nights)
        
    def reduceMedian(self, (pairs, nights)):
        """Reduce to median number of pairs per night."""
        return np.median(pairs)

    def reduceMean(self, (pairs, nights)):
        """Reduce to mean number of pairs per night."""
        return pairs.mean()
    
    def reduceRms(self, (pairs, nights)):
        """Reduce to std dev of number of pairs per night."""
        return pairs.std()

    def reduceNNightsWithPairs(self, (pairs, nights), nPairs=None):
        """Reduce to number of nights with more than 'nPairs' (default=1) visits."""
        if not nPairs:
            # Number of PAIRS (i.e. two visits = 1 pair)
            nPairs = 1
        condition = (pairs >= nPairs)
        return len(pairs[condition])

    def reduceNPairsInWindow(self, (pairs, nights), window=None):
        """Reduce to max number of pairs within 'window' (default=30 nights) of time."""
        if not window:
            window = self.window
        maxnpairs = 0
        for n in nights:
            condition = ((nights >= n) & (nights <= n+window))
            maxnpairs = max((pairs[condition].sum(), maxnpairs))
        return maxnpairs

    def reduceNNightsInWindow(self, (pairs, nights), window=None):
        """Reduce to max number of nights with a pair (or more) of visits, within 'window'."""
        if not window:
            window=self.window
        maxnights = 0
        for n in nights:
            condition = ((nights >=n) & (nights<=n+window))
            maxnights = max(len(nights[condition]), maxnights)
        return maxnights
        