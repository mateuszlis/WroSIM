'''
Created on 16-01-2011

@author: lisu
'''
from structures.baseStructs import Atom, Frame
from myExceptions import InputError

from utils import concatenate

class GroAtom(Atom):
    '''
    "Class represents atom stored in gro file"""
    '''


    def __init__(self, symbol, number = 1, x = 0., y = 0., z = 0., vx = 0., vy = 0., vz = 0.):
        if len(symbol) > 10:
            self._construct(symbol)
        else:
            Atom.__init__(self, symbol, x, y, z)
            self._vx, self._vy, self._vz, self._number = vx, vy, vz, number
            
    @property
    def vx(self):
        return self._vx
    @vx.setter
    def vx(self, vx):
        self._vx = vx
    @property
    def vy(self):
        return self._vy
    @vy.setter
    def vy(self, vy):
        self._vy = vy
    @property
    def vz(self):
        return self._vz
    @vz.setter
    def vz(self, vz):
        self._vz = vz
    @property
    def number(self):
        return self._number
    @number.setter
    def number(self, number):
        self._number = int(number)


    @property
    def v(self):
        return [self.vx, self.vy, self.vz]
    @v.setter
    def v(self, v):
        self.vx, self.vy, self.vz = v
    

    def _construct(self, groAtomStr):
        coordList = groAtomStr.split()
        if not (len(coordList) == 6 or len(coordList) == 9) :
            raise InputError(groAtomStr, 'Wrong number of arguments while \
                                    constructing the groAtom')
        self.v = [None for i in range(3)]
        if len(coordList) == 9:
            self.v = [float(i) for i in coordList[6:]]
        self.symbol, self.number = coordList[1:3]
        self.x0 = [float(i) for i in coordList[3:6]]

class Molecule(object):
    """Class hold molecule object which contains group of atoms and a residue name"""
    def __init__(self, resname, atoms = []):
        self._resname = resname
        self._atoms = atoms
    
    @property
    def resname(self):
        return self._resname
    @resname.setter
    def resname(self, resname):
        self._resname = resname
    @property 
    def atoms(self):
        return self._atoms
    @atoms.setter
    def atoms(self, atoms):
        self._atoms = atoms
        
class GroFrame(Frame):
    """Class representing frame in trajectory file"""
    def __init__(self, frameStr = None):
        Frame.__init__(self)
        self._boxSize = [0, 0, 0]
        self._molecules = []
        if self._isProper(frameStr):
            self._construct(frameStr)
    def __repr__(self):
        lines = []
        lines.append(self.comment)
        lines.append(str(len(self.atoms)))
        i, j = 1, 1
        for mol in self.molecules:
            for atom in mol.atoms:
                if atom.vx is None:
                    lines.append("%5d%5s%5s%5d%8.3f%8.3f%8.3f" % ( 
                            i, mol.resname, atom.symbol, j, atom.x, atom.y,
                            atom.z))
                else:
                    lines.append("%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f" % 
                                 (i, mol.resname, atom.symbol, j, atom.x, atom.y,
                                  atom.z, atom.vx, atom.vy, atom.vz))
                atom.number = j
                j += 1
            i += 1
        lines.append('%f %f %f' % (self.boxSize[0], self.boxSize[1], self.boxSize[2])) 
        return '\n'.join(lines)
    @property
    def molecules(self):
        return self._molecules
    @molecules.setter
    def molecules(self, molecules):
        self._molecules = molecules
    @property
    def boxSize(self):
        return self._boxSize
    @boxSize.setter
    def boxSize(self, boxSize):
        self._boxSize = boxSize
    @property
    def atoms(self):
        return concatenate([mol.atoms for mol in self.molecules])
            
        
    def _isProper(self, frameStr):
        if not frameStr is None:
            lines = frameStr.split('\n')
            if lines[-1] == '': del(lines[-1])
            if len(lines) < 3:
                raise InputError(frameStr, 'Wrong number of arguments while \
                                    constructing the XYZFrame')
            return True
        return False
    
    def _construct(self, frameStr):
        try:
            lines = frameStr.split('\n')
            length = int(lines[1])
            self._comment = lines[0]
            lastResNum = 0
            for line in lines[ 2 : (2 + length)]:
                resNum, resname= int(line[:5]), line[5:10]
                if lastResNum != resNum:
                    lastResNum = resNum
                    self.molecules.append(Molecule(resname, []))
                self.molecules[-1].atoms.append(GroAtom(line))
            self.boxSize = [float(val) for val in lines[-1].split()]   
        except:
            raise 
        