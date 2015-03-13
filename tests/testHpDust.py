import numpy as np
import unittest
from lsst.sims.photUtils import EBV
from lsst.sims.maf.maps import EBVhp


class TestDustMaps(unittest.TestCase):

    def testDust(self):
        np.random.seed(42)
        ra = np.random.rand(100)*2*np.pi
        dec = np.random.rand(100)*np.pi - np.pi/2.

        nside = 512

        hpVals = EBVhp(nside, ra=ra, dec=dec, interp=False)

        dustmap = EBV.EBVbase()
        dustmap.load_ebvMapNorth()
        dustmap.load_ebvMapSouth()
        schlegVals = dustmap.calculateEbv(np.array([ra,dec]), interp=False)



if __name__ == "__main__":
    unittest.main()
