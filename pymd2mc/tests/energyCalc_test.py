'''
Created on 17-02-2011

@author: lisu
'''
import unittest

from lattices import Position, HexLattice
from structures.xyzfile import XYZAtom, XYZFrame, XYZFile
from xyz2hexnet.energyCalc import HexLatticeLoader, EnergyCalculator
from lattices import HexLattice

class TestEnergyCalculator(unittest.TestCase):
    def setUp(self):
        self.xyzFile = XYZFile('data/smallLattice.xyz')
        self.energyCalc = EnergyCalculator(self.xyzFile)
    def testGetNextEnergy(self):
        self.energyCalc.getNextEnergy("B")
        self.energyCalc.getNextEnergy("B")
        #self.assertAlmostEqual(self.energyCalc.getNextEnergy('B'), 575.11680054416524)
        #self.assertAlmostEqual(self.energyCalc.getNextEnergy('B'), -191.00144095067344)
        self.assertAlmostEqual(self.energyCalc.getNextEnergy('B')[0], -191.00144095067344)
    def tearDown(self):
        del(self.xyzFile)
        
class TestHexLatticeLoader(unittest.TestCase):
    def setUp(self):
        self.lattice = HexLattice(3)
        self.firstFrame, self.secondFrame = XYZFrame(), XYZFrame()
        for pos in self.lattice.positions:
            self.firstFrame.atoms.append(XYZAtom('A', *pos.x0))
            self.secondFrame.atoms.append(XYZAtom("B", *pos.x0))
        self.loader = HexLatticeLoader(self.firstFrame)

    def tearDown(self):
        pass


    def testGetNeighbors(self):
        
        properNeighbrsDict = {0 : [8, 6, 2, 1, 3, 4],
                          1 : [6, 7, 0, 2, 4, 5],
                          4 : [0,1,3,5,7,8],
                          6 : [5, 3, 8, 7, 0, 1]}
        for neighbor, properNeighbrs in properNeighbrsDict.items():
            neighList = self.loader.getNeighbors(self.firstFrame.atoms[neighbor])
            #print neighbor, ' ', properNeighbrs
            for neighbrNum in properNeighbrs:
                neigh = neighList.next()
                for i in range(3):
                    self.assertAlmostEqual(neigh.x0[i],self.firstFrame.atoms[neighbrNum].x0[i])
    def testUpdate(self):
        self.loader.updateState(self.secondFrame)
        neighList = self.loader.getNeighbors(self.firstFrame.atoms[0])
        for atom in neighList:
            self.assertEqual(atom.symbol, "B")
    def testCalcNeighborsCount(self):
        similar, different = self.loader.calcNeighborsCount('A')
        self.assertEqual(similar, 27)
        self.assertEqual(different, 0)
        
        self.loader.updateState(self.secondFrame)
        similar, different = self.loader.calcNeighborsCount('B')
        self.assertEqual(similar, 27)
        self.assertEqual(different, 0)
        
        self.secondFrame.atoms[4].symbol = 'A'
        self.loader.updateState(self.secondFrame)
        similar, different = self.loader.calcNeighborsCount('B')
        self.assertEqual(different, 6)
        self.assertEqual(similar, 21)
        
        
        

if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()