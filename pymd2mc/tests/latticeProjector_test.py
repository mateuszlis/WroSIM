'''
Created on 02-01-2011

@author: lisu
'''
import unittest
from lattices import Position, HexLattice
from structures.xyzfile import XYZAtom
from xyz2hexnet.latticeProjector import LatticeProjectorSimple


class TestLatticeProjectorSimple(unittest.TestCase):


    def setUp(self):
        self.emptyPos = Position([0,0])
        self.pos = Position([1,2], XYZAtom('P'))
        self.lattice = HexLattice(4, [0,0,0], [100,100])
        self.atoms = []
        self.atoms.append(XYZAtom('N', 0,0,0))
        self.atoms.append(XYZAtom('N', 1,0,0))
        
        for i in range(5):
            self.atoms.append(XYZAtom('P', i + 5, i + 4, i + 3))
         
    def test_findPlace(self):
        self.assertEqual(LatticeProjectorSimple.findPlace(self.atoms[0],self.lattice, [0]), 1)
        self.assertEqual(LatticeProjectorSimple.findPlace(self.atoms[0],self.lattice), 0)
    def tearDown(self):
        pass



if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()