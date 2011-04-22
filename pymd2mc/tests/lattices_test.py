'''
Created on 02-01-2011

@author: lisu
'''
import unittest

from math import sin, cos, pi

from lattices import Position, HexLattice
from structures.xyzfile import XYZAtom

class TestPosition(unittest.TestCase):
    def setUp(self):
        self.emptyPos = Position()
        self.pos = Position([1,2])
    def test_x0(self):
        self.assertAlmostEqual(self.emptyPos.x0[0], 0)
        self.assertAlmostEqual(self.emptyPos.x0[1], 0)
        self.assertAlmostEqual(self.pos.x0[0], 1)
        
        self.emptyPos.x0 = [1,1]
        self.assertAlmostEqual(self.emptyPos.x0[0], 1)
        self.assertAlmostEqual(self.emptyPos.x0[1], 1)
        
class TestHexLattice(unittest.TestCase):
    pass
    def setUp(self):
        self.lattice = HexLattice(1, (0,0), (1,1))
        
    def test_construct(self):
        self.assertEqual(len(self.lattice.positions), 1)
        self.assertAlmostEqual(2 * self.lattice.positions[0].x0[1], sin(pi / 3))
  
    def test_simpleProps(self):
        self.assertEqual(self.lattice.x0[0], 0)
        self.assertEqual(self.lattice.n, 1)
        self.assertEqual(self.lattice.dx[0], 1)
        self.assertEqual(len(self.lattice.pbcSites),6)
    def test_moving(self):
        self.lattice.x0 = (2,2)
        self.assertAlmostEqual(2 * (self.lattice.positions[0].x0[1] - 2), sin(pi / 3))
    def test_boxVectors(self):
        self.assertAlmostEqual(self.lattice.boxVectors[0][0], 1)
        self.assertAlmostEqual(self.lattice.boxVectors[0][1], 0)
        
        self.assertAlmostEqual(self.lattice.boxVectors[1][0], - cos(pi/3.) )
        self.assertAlmostEqual(self.lattice.boxVectors[1][1], sin(pi/3.))
        
        self.assertAlmostEqual(self.lattice.boxVectors[2][0],  0)
        self.assertAlmostEqual(self.lattice.boxVectors[2][1], 0)
        


if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()