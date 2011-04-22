'''
Created on 18-01-2011

@author: lisu
'''
from copy import deepcopy
import unittest

from myExceptions import InputError
from structures.structurefile import StructureFile
from structures.xyzfile import XYZAtom
from structures.grofile import GroAtom

     
class TestStructureFile(unittest.TestCase):
    def setUp(self):
        self.xyz = StructureFile('data/sample_benzene.xyz')
        self.xyz2 = StructureFile('data/sample.xyz')
        self.gro = StructureFile('data/conf.gro')
        
    def test_loadXYZ(self):
        frame = self.xyz2.frames[0]
        
        self.assertAlmostEqual(0.0, frame.atoms[0].x)
        self.assertAlmostEqual(1.40272, frame.atoms[0].y)
        self.assertAlmostEqual(0.0, frame.atoms[0].z)
        self.assertEqual('C', frame.atoms[0].symbol)
        self.assertEqual(1, len(self.xyz.frames))
        
    def test_loadGro(self):
        frame = self.gro.frames[0]
        
        self.assertAlmostEqual(1.561, frame.atoms[0].x)
        self.assertAlmostEqual(3.310, frame.atoms[0].y)
        self.assertAlmostEqual(1.537, frame.atoms[0].z)
        self.assertEqual('C1', frame.atoms[0].symbol)
        self.assertEqual(1, len(self.xyz.frames))
        
        
    def test_loadEmpty(self):
        self.xyz3 = StructureFile('data/empty.xyz')
        self.assertRaises(InputError, self.xyz3.load)
        self.gro3 = StructureFile('data/empty.gro')
        self.assertRaises(InputError, self.gro3.load)
        
    def test_saveXYZ(self):
        frame = self.xyz2.frames[0]
        frame.atoms.append(XYZAtom('C', 1., 2., 3.))
        self.xyz2.filename = 'data/test.xyz'
        self.xyz2.save()
        
        self.xyz3 = StructureFile('data/test.xyz')
        newFrame = self.xyz3.frames[0]
        
        self.assertAlmostEqual(1., newFrame.atoms[-1].x)
        self.assertAlmostEqual(2., newFrame.atoms[-1].y)
        self.assertAlmostEqual(3., newFrame.atoms[-1].z)
    
    def test_saveGro(self):
        frame = self.gro.frames[0]
        frame.atoms.append(GroAtom('C', 1., 2., 3.))
        self.gro.filename = 'data/test.gro'
        self.gro.save()
        
        self.gro3 = StructureFile('data/test.gro')
        newFrame = self.xyz3.frames[0]
        
        self.assertAlmostEqual(1., newFrame.atoms[-1].x)
        self.assertAlmostEqual(2., newFrame.atoms[-1].y)
        self.assertAlmostEqual(3., newFrame.atoms[-1].z)
        
    def test_getFilename(self):
        self.assertEqual('data/sample_benzene.xyz', self.xyz.filename)
        self.assertEqual('data/sample.xyz', self.xyz2.filename)
    
    def test_nextFrame(self):
        for i in range(6):
            frame = self.xyz2.nextFrame()
        
            self.assertEqual( i,self.xyz2.getCurrentFrameNumber())
            self.assertAlmostEqual(0 + i * 0.1, frame.atoms[0].z)
        self.assertEqual(None, self.xyz2.nextFrame())
    
    def test__repr(self):
        self.assertEqual(24, len(self.xyz2.__repr__().split('\n')))
        self.assertEqual(14, len(self.xyz.__repr__().split('\n')))
    
    def test_addFrame(self):
        f = open('data/incrementialFile.xyz', 'w')
        f.close()
        self.xyz3 = StructureFile('data/incrementialFile.xyz')
        frame = deepcopy(self.xyz2.frames[0])
        self.xyz3.addFrame(frame)
        self.xyz3.addFrame(frame)
        
        self.assertAlmostEqual(self.xyz3.frames[0].atoms[0].x,self.xyz2.frames[0].atoms[0].x)
         
        self.xyz.load()
        self.xyz.addFrame(frame)
        self.assertAlmostEqual(self.xyz3.frames[0].atoms[0].x,self.xyz.frames[-1].atoms[0].x) 
    def test_getFrame(self):
        frame = self.xyz2.getFrame(4)
        
        self.assertAlmostEqual(0.0, frame.atoms[0].x)
        self.assertAlmostEqual(1.40272, frame.atoms[0].y)
        self.assertAlmostEqual(0.40, frame.atoms[0].z)
        self.assertEqual('C', frame.atoms[0].symbol)
         
if __name__ == "__main__":
    #import sys;sys.argv = ['', 'Test.testName']
    unittest.main()