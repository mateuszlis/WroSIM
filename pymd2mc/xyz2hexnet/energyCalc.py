"""
This class calculates energy value from a state of hexagonal lattice.
Currently it handles only two types of particles in lattice file.
"""
from myExceptions import InputError
__author__ =  'Mateusz Lis'
__version__=  '0.1'

from math import log
from operator import attrgetter

class EnergyCalculator(object):
    '''
    This class calculates energy value from a state of hexagonal lattice.
    Currently it handles only two types of particles in lattice file.
    '''
    def __init__(self, xyzFile, R = 1.986, T = 310):
        self.xyzFile = xyzFile
        self.frame = self.xyzFile.nextFrame()
        self.R, self.T = R, T
        if self.frame is None:
            raise InputError(xyzFile, "no frames in file")
        self.hexLatticeLoader = HexLatticeLoader(self.frame)
    def getNextEnergy(self, symbol = 'A'):
        """Returns value of energy for next frame"""
        if self.frame is None:
            self.frame = self.xyzFile.nextFrame()
            
        if self.frame is None:
            return None
        self.hexLatticeLoader.updateState(self.frame, symbol)
        
        N_AA, N_BB, N_AB = self.hexLatticeLoader.calcNeighborsCount(symbol)
        N_A = self.hexLatticeLoader.N_A
        N_B = self.hexLatticeLoader.N_B
        PwAB = float(N_AB) / (N_AA * (N_B/float(N_A)) + N_BB * float(N_A)/N_B) 
        print "PwAB %s sim %s diff %s, %s" % (PwAB, N_AA, N_BB, self.R)
        if PwAB < 0.01:
            wAB = -10**10
            print '!!"'
        else:
            wAB =  - self.R * self.T * log(PwAB)
        self.frame = None
        return (wAB, N_AA, N_AB)

class HexLatticeLoader(object):
    def __init__(self, frame):
        self._atoms = sorted(frame.atoms, key=attrgetter('y','x'))
        self._n = int(len(self._atoms) ** (1/2.))
        self._neighSites =[(-1,-1), (-1, 0), (0,-1), (0, 1), (1, 0), (1, 1)]
        self._createPositionsDict()
        self._N_A = 0
        self._N_B = 0
        
    def updateState(self, frame, symbol):
        for atom in frame.atoms:
            self._atoms[self._indexForPosition[(atom.x, atom.y)]] = atom
            if atom.symbol == symbol:
                self._N_A += 1
            else:
                self._N_B += 1
    
    def getNeighbors(self, atom):
        atomInd = self._indexForPosition[(atom.x, atom.y)]
        for neighbor in self._neighSites:
            row = (int(atomInd / 3) + neighbor[0]) % self._n
            col = (int(atomInd % 3) + neighbor[1]) % self._n
            yield self._atoms[int(row * self._n + col)]
    def _createPositionsDict(self):
        self._indexForPosition = {}
        for atomInd in range(len(self._atoms)):
            key = (self._atoms[atomInd].x, self._atoms[atomInd].y)
            self._indexForPosition[key] = atomInd
            
    def calcNeighborsCount(self, symbol):
        N_AB, N_AA, N_BB = 0, 0, 0
        for atom in self._atoms:
            for neighAtom in self.getNeighbors(atom):
                if  atom.symbol == neighAtom.symbol:
                    if atom.symbol == symbol:
                        N_AA += 1
                    else:
                        N_BB += 1
                else:
                    #print 'numb %s symb1 %s symb2 %s' % (different, atom.symbol, neighAtom.symbol)
                    #print atom
                    #print neighAtom 
                    N_AB += 1
        return (N_AA / 2, N_BB / 2, N_AB / 2)
    @property
    def N_A(self):
        return self._N_A
    @property
    def N_B(self):
        return self._N_B