import numpy as np
import matplotlib.pyplot as plt
import unittest
from lsst.sims.maf.slicers.uniSlicer import UniSlicer
from lsst.sims.maf.slicers.oneDSlicer import OneDSlicer

def makeDataValues(size=100, min=0., max=1., random=True):
    """Generate a simple array of numbers, evenly arranged between min/max, but (optional) random order."""    
    datavalues = np.arange(0, size, dtype='float')
    datavalues *= (float(max) - float(min)) / (datavalues.max() - datavalues.min()) 
    datavalues += min
    if random:
        randorder = np.random.rand(size)        
        randind = np.argsort(randorder)
        datavalues = datavalues[randind]
    datavalues = np.array(zip(datavalues), dtype=[('testdata', 'float')])
    return datavalues
    

class TestUniSlicerSetupAndSlice(unittest.TestCase):    
    def setUp(self):
        self.testslicer = UniSlicer()
        
    def tearDown(self):
        del self.testslicer
        self.testslicer = None

    def testSlicertype(self):
        """Test instantiation of slicer sets slicer type as expected."""        
        self.assertEqual(self.testslicer.slicerName, self.testslicer.__class__.__name__)
        self.assertEqual(self.testslicer.slicerName, 'UniSlicer')

    def testSlicerNbins(self):
        self.assertEqual(self.testslicer.nbins, 1)
        
    def testSetupSlicerIndices(self):
        """Test slicer returns correct indices (all) after setup. Note this also tests slicing."""
        self.assertRaises(NotImplementedError, self.testslicer.sliceSimData, 0)
        dvmin = 0
        dvmax = 1        
        nvalues = 1000
        dv = makeDataValues(nvalues, dvmin, dvmax, random=True)
        self.testslicer.setupSlicer(dv)
        # test slicing
        self.assertEqual(len(self.testslicer.indices), len(dv['testdata']))
        np.testing.assert_equal(dv[self.testslicer.indices], dv)


class TestUniSlicerIteration(unittest.TestCase):
    def setUp(self):
        self.testslicer = UniSlicer()

    def tearDown(self):
        del self.testslicer
        self.testslicer = None

    def testIteration(self):
        """Test iteration -- which is a one-step identity op for a unislicer."""
        dvmin = 0
        dvmax = 1
        nvalues = 1000
        dv = makeDataValues(nvalues, dvmin, dvmax, random=True)
        self.testslicer.setupSlicer(dv)
        for i, b in enumerate(self.testslicer):
            pass
        self.assertEqual(i, 0)

    def testGetItem(self):
        """Test that can return an individual indexed values of the slicer."""
        self.assertEqual(self.testslicer[0], 0)
        
class TestUniSlicerEqual(unittest.TestCase):
    def setUp(self):
        self.testslicer = UniSlicer()
        dvmin = 0
        dvmax = 1
        nvalues = 1000
        dv = makeDataValues(nvalues, dvmin, dvmax, random=True)
        self.testslicer.setupSlicer(dv)    

    def tearDown(self):
        del self.testslicer
        self.testslicer = None

    def testEquivalence(self):
        """Test equals method."""
        # Note that two uni slicers will be considered equal if they are both the same kind of
        # slicer (unislicer). They will not necessarily slice data equally though (the indices are
        #  not necessarily the same!).
        # These should be the same, even though data is not the same.
        testslicer2 = UniSlicer()
        dv2 = makeDataValues(100, 0, 1, random=True)
        testslicer2.setupSlicer(dv2)
        self.assertEqual(self.testslicer, testslicer2)
        # these will not be the same, as different slicer type.
        testslicer2 = OneDSlicer(sliceDataColName='testdata')
        testslicer2.setupSlicer(dv2, bins=10)
        self.assertNotEqual(self.testslicer, testslicer2)
            
if __name__ == "__main__":
    suitelist = [unittest.TestLoader().loadTestsFromTestCase(TestUniSlicerSetupAndSlice),]
    suitelist.append(unittest.TestLoader().loadTestsFromTestCase(TestUniSlicerIteration))
    suitelist.append(unittest.TestLoader().loadTestsFromTestCase(TestUniSlicerEqual))
    suite = unittest.TestSuite(suitelist)
    unittest.TextTestRunner(verbosity=2).run(suite)

    # slicing tested as part of setup here, and 'function' is identity function
    #  so equivalent to slicing. 