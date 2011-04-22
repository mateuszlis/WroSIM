from copy import deepcopy
import unittest

from structures.xyzfile import XYZAtom
from structures.xyzfile import XYZFrame
from structures.xyzfile import XYZFile
from myExceptions import InputError

class TestXYZAtom(unittest.TestCase):
    def setUp(self):
        self.xyzAtom = XYZAtom('C')
        self.xyzAtom2 = XYZAtom('C',1,1,1)
    
    def test_setters(self):
        
        self.xyzAtom.x = 2
        self.xyzAtom.y = 2
        self.xyzAtom.z = 2
        self.xyzAtom.symbol = 'C'
        
        self.assertEqual(self.xyzAtom.symbol, 'C')
        self.assertEqual(self.xyzAtom.x, 2)
        self.assertEqual(self.xyzAtom.y, 2)
        self.assertEqual(self.xyzAtom.z, 2)
    
    def test_constructor(self):
        self.xyzAtom = XYZAtom('  C        0.00000        1.40272        0.00000')
       
        self.assertAlmostEqual(self.xyzAtom.x, 0.0)
        self.assertAlmostEqual(self.xyzAtom.y, 1.40272)
        self.assertAlmostEqual(self.xyzAtom.z, 0.000)
        
    def test_repr(self):
        desc = self.xyzAtom.__repr__()
        symbol, x, y, z = desc.split()
        x, y, z = float(x), float(y), float(z)
        self.assertEqual(self.xyzAtom.symbol, symbol)
        self.assertAlmostEqual(self.xyzAtom.x, x)
        self.assertAlmostEqual(self.xyzAtom.y, y)
        self.assertAlmostEqual(self.xyzAtom.z, z)
    
    def test_x0(self):
        self.assertEquals(self.xyzAtom.x0[0], 0)
        self.assertEquals(self.xyzAtom2.x0[2], 1)
        self.xyzAtom.x0 = [1,2,3]
        self.assertEquals(self.xyzAtom.x0[2], 3)
    
class TestXYZFrame(unittest.TestCase):
    def setUp(self):
        self.frameFile = open('data/sample_benzene.xyz', 'r')
    
    def compareAtoms(self, atom1, atom2):
        self.assertEqual(atom1.symbol, atom2.symbol)
        self.assertAlmostEqual(atom1.x, atom2.x)
        self.assertAlmostEqual(atom1.y, atom2.y)
        self.assertAlmostEqual(atom1.z, atom2.z)
        
    def test_construct(self):
        frameStr = ''.join(self.frameFile.readlines())
        frame = XYZFrame(frameStr)
        
        atomList = frame.atoms
        self.compareAtoms(atomList[0], XYZAtom('C        0.00000        1.40272        0.00000'))
        self.compareAtoms(atomList[1], XYZAtom('H        0.00000        2.49029        0.00000'))
        self.assertEqual(12, len(XYZFrame(frameStr).atoms))
    
    def test_addAtom(self):
        frame = XYZFrame(''.join(self.frameFile.readlines()))
        frame.atoms.append(XYZAtom('H        0.00000        2.49029        1.00000'))
        self.assertEqual(13, len(frame.atoms))
        
    def test_editAtom(self):
        frame = XYZFrame(''.join(self.frameFile.readlines()))
        frame.atoms[5] = XYZAtom('C',1,1,1)
        self.compareAtoms(frame.atoms[5], XYZAtom('C',1,1,1))
        
    def test_boxSize(self):
        def helper(frame):
            """Nature of assertRaises makes me have to do this"""
            _ = frame.boxSize
        frameList = self.frameFile.readlines()
        frameStr = ''.join(frameList)
        frame = XYZFrame(frameStr)
        
        self.assertRaises(InputError, helper, frame)

        frameList[1] = '10,10,01 \n'
        frameStr = ''.join(frameList)
        frame = XYZFrame(frameStr)
        
        self.assertAlmostEqual(frame.boxSize[0], 10.)
    
    def test_repr(self):
        frameStr = ''.join(self.frameFile.readlines())
        frame = XYZFrame(frameStr)
        reprStr = frame.__repr__()
        self.assertEqual(14, len(reprStr.split('\n')))
             
    def tearDown(self):
        self.frameFile.close()
          
class TestXYZFile(unittest.TestCase):
    def setUp(self):
        self.xyz = XYZFile('data/sample_benzene.xyz')
        self.xyz2 = XYZFile('data/sample.xyz')
        
    def test_load(self):
        frame = self.xyz2.frames[0]
        
        self.assertAlmostEqual(0.0, frame.atoms[0].x)
        self.assertAlmostEqual(1.40272, frame.atoms[0].y)
        self.assertAlmostEqual(0.0, frame.atoms[0].z)
        self.assertEqual('C', frame.atoms[0].symbol)
        self.assertEqual(1, len(self.xyz.frames))
    def test_loadEmpty(self):
        self.xyz3 = XYZFile('data/empty.xyz')
        self.assertRaises(InputError, self.xyz3.load)
        
    def test_save(self):
        frame = self.xyz2.frames[0]
        frame.atoms.append(XYZAtom('C', 1., 2., 3.))
        self.xyz2.filename = 'data/test.xyz'
        self.xyz2.save()
        
        self.xyz3 = XYZFile('data/test.xyz')
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
        self.xyz3 = XYZFile('data/incrementialFile.xyz')
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
         
    

suite1 = unittest.TestLoader().loadTestsFromTestCase(TestXYZAtom)
suite2 = unittest.TestLoader().loadTestsFromTestCase(TestXYZFrame)
suite3 = unittest.TestLoader().loadTestsFromTestCase(TestXYZFile)
XYZSuite = unittest.TestSuite([suite1])
XYZSuite2 = unittest.TestSuite([suite2])
XYZSuite3 = unittest.TestSuite([suite3])
unittest.TextTestRunner(verbosity=2).run(XYZSuite)
unittest.TextTestRunner(verbosity=2).run(XYZSuite2)
unittest.TextTestRunner(verbosity=2).run(XYZSuite3)

