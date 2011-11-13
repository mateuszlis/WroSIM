"""
This classes are responsible of projecting continous variables trajectory,
stored in xyzFile, to lattice trajectory (with discrete fields).
"""
__author__ =  'Mateusz Lis'
__version__=  '0.1'

from math import sqrt, pi

from utils import distance, posmax, posmin, concatenate
from structures.xyzfile import XYZAtom
from structures.xyzfile import XYZFrame

class ProjectedFrame(object):
    "Class hold results of LatticeProjectorSimple.next() method."
    def __init__(self, frame, refFrame, errors):
        self._frame = frame
        self._refFrame = refFrame
        self._errors = errors
    @property
    def frame(self):
        return self._frame
    @property
    def refFrame(self):
        return self._refFrame
    @property
    def errors(self):
        return self._errors
    
class LatticeProjectorSimple(object):
    '''
    This lattice projector is capable of projecting adatoms from xyz file to
    a proper lattice. It utilizes nearest neighbour algorithm and PBC 
    multiplication to keep shape of lattice.
    '''


    def __init__(self, xyzFile, lattice):
        '''
        Constructor
        '''
        self._xyzFile = xyzFile
        self._lattice = lattice
    
    @property
    def xyzFile(self):
        """xyzFile handler"""
        return self._xyzFile
    @property
    def lattice(self):
        """lattice on which contents of xyzFile are projected"""
        return self._lattice
    def next(self):
        """
        Returns lattice containing projection of the next frame of xyzFile.
        If file has finnished it returns None value
        """
        frame = self.xyzFile.nextFrame()
        if frame is None: return None
        
        newFrame = XYZFrame()
        newFrame.boxVectors = self.lattice.boxVectors
        refFrame = XYZFrame()
        refFrame.boxVectors = self.lattice.boxVectors
        atomsLists = self.propagateAtomsThroughPbc(frame.atoms, frame.boxSize)
       
        allAtoms = concatenate(atomsLists)  
        posCount = len(atomsLists[0])
        
        match, referenceMatch, errors = self.match(atomsLists)   
        for atomIndex in range(posCount):
            newFrame.atoms.append(XYZAtom(atomsLists[0][atomIndex].symbol
                                , *self.lattice.positions[match[atomIndex]].x0))
        
        for atomIndex in range(posCount):
            refFrame.atoms.append(XYZAtom(allAtoms[referenceMatch[atomIndex]].__repr__()))    
            refFrame.atoms[-1].x += 15
            
        for atomIndex in range(len(allAtoms)):
            refFrame.atoms.append(XYZAtom(allAtoms[atomIndex].__repr__()))    
            refFrame.atoms[-1].x += 30
        
        return ProjectedFrame(newFrame, refFrame, errors)

    def propagateAtomsThroughPbc(self, atoms, boxSize):
        """Function propagates positions of atoms to other sites (as in
        periodic boundary conditions
        @return: list of lists of atoms for each propagated site"""
        atomsLists = []
        for coordinates in self.lattice.pbcSites:
            _atoms = []
            for atom in atoms:
                _atoms.append(XYZAtom(atom.symbol, atom.x + coordinates[0] * boxSize[0]
                                     , atom.y + coordinates[1] * boxSize[1]
                                     , atom.z + coordinates[2] * boxSize[2]))
                
            atomsLists.append(_atoms)
        return atomsLists
    def match(self, atomsLists):
        """Function matches atoms from atomsLists to self.lattice.positions
        @return: match, referenceMatch, errors
        match is list of indices of places for atoms match[atomIndex] = placeIndex
        referenceMatch[atomIndex] = atomIndex on allAtomsList
        errors is a list of (max error, average error)"""
        errors = []
        allAtoms = concatenate(atomsLists)      
        forbiddenPosList = []
        forbiddenAtomsList = [] #we'll need both of this to know which atom is banned
        posCount = len(atomsLists[0])
        match = [0 for _ in range(posCount)]
        referenceMatch = [0 for _ in range(posCount)]
        self._initPosSearch(len(allAtoms))
        
        for _ in range(posCount):
            
            atomIndex, placeIndex = self.matchPos(allAtoms, len(atomsLists)
                                                  , forbiddenPosList
                                                  , forbiddenAtomsList)
            if not placeIndex is None: 
                forbiddenPosList.append(placeIndex)
            forbiddenAtomsList += [(atomIndex + n * len(atomsLists[0])) % len(allAtoms) 
                                   for n in range(len(atomsLists))]
            if not placeIndex is None:
                pos = self.lattice.positions[placeIndex] 
            errors.append(distance(pos.x0, allAtoms[atomIndex].x0))            
            if not placeIndex is None:
                match[atomIndex % len(atomsLists[0])] = placeIndex
            referenceMatch[atomIndex % len(atomsLists[0])] = atomIndex
        return match, referenceMatch, errors
    def matchPos(self, allAtoms, pbcCount, forbiddenPosList, forbiddenAtomsList):
        """Matches positions for allAtoms except for forbiddenAtoms and forbidden
        Positions"""
        length = len(allAtoms) / pbcCount
        for atomIndex in range(len(allAtoms)):
            atom = allAtoms[atomIndex]
            if atomIndex in forbiddenAtomsList:
                self._places[atomIndex] = (-1)
                self._distances[atomIndex] = (-1000)
            else:
                if self._places[atomIndex] in forbiddenPosList or self._places[atomIndex] == -1:
                    self._places[atomIndex] = (self.findPlace(atom, self.lattice
                                          , forbiddenPosList))
                    if not self._places[atomIndex] is None:
                        self._distances[atomIndex] = distance(
                                            self.lattice.positions[self._places[atomIndex]].x0
                                                                , [atom.x, atom.y])
               
                
        minDists = [min([self._distances[i + n * length] for n in range(pbcCount)]) for i in range(length)]
        posMinDists = [posmin([self._distances[i + n * length] for n in range(pbcCount)]) for i in range(length)]
        indexOfAtom = posmax(minDists)
        indexOfAtom += posMinDists[indexOfAtom] * length    
        return (indexOfAtom, self._places[indexOfAtom])
   
    @staticmethod
    def findPlace(atom, lattice, forbiddenList = []):
        availablePos = range(len(lattice.positions))
        for pos in forbiddenList:
            availablePos.remove(pos)
        if len(availablePos) > 0:
            nn = availablePos[0]
        else: return None
        dist = distance([atom.x, atom.y], lattice.positions[nn].x0)
        
        for posIndex in availablePos:
            pos = lattice.positions[posIndex]
            if distance([atom.x, atom.y], pos.x0) < dist:
                nn = posIndex
                dist = distance([atom.x, atom.y], pos.x0)
        return nn

    def _initPosSearch(self, n):
        self._distances = [-1000 for _ in range(n)]
        self._places = [-1 for _ in range(n)]

class NearestNeighborLatticeProj(LatticeProjectorSimple):
    def __init__(self, xyzFile, lattice):
        '''
        Constructor
        '''
        self._xyzFile = xyzFile
        self._lattice = lattice
	apl = ( 0.5 * (lattice.dx[0] * lattice.dx[1]) * (3**(1./2))) / lattice.n**2
	self._radius = (( apl / pi )**(1./2)) 
	print self._radius, "radius"
    def next(self):
        """
        Returns lattice containing projection of the next frame of xyzFile.
        If file has finnished it returns None value
        """
        frame = self.xyzFile.nextFrame()
	if frame is None: return None
        
        newFrame = XYZFrame()
        newFrame.boxVectors = self.lattice.boxVectors
        refFrame = XYZFrame()
        refFrame.boxVectors = self.lattice.boxVectors
	atomsLists = self.propagateAtomsThroughPbc(frame.atoms, frame.boxSize)

        allAtoms = concatenate(atomsLists)  
       
        #posCount = len(atomsLists[0])
        errors = []
	differences = []
        #match, referenceMatch, errors = self.match(atomsLists)   
	atoms = []
#	print self._radius
	for pos in self.lattice.positions:
		ACounter = 0
		BCounter = 0
		for atom in allAtoms:
			dist = distance( atom.x0, pos.x0 )
			atoms.append(atom)
			if atom.symbol == "DLPC": 
				try:
					ACounter += 1. / ( dist**2 )
				except:
					ACounter = 1e308
			else:
				if atom.symbol == "DOPC": 
					try:
						BCounter += 1. / ( dist**2 )
					except:
						BCounter = 1e308
		differences.append( (ACounter - BCounter, pos.x0) )
#		if ACounter > BCounter: symbol = "DSPC"
#		else: symbol = "A"
#		newFrame.atoms.append(XYZAtom(symbol, *pos.x0))
 		errors.append(0.0)
	num = 0
	for atom in atoms:
		refFrame.atoms.append(XYZAtom(atom.symbol, *atom.x0))
		refFrame.atoms[-1].x += 17
	for atom in allAtoms:
		refFrame.atoms.append(XYZAtom(atom.symbol, *atom.x0))
		refFrame.atoms[-1].x += 48
	
	for diff, pos in sorted(differences, key = lambda tup: tup[0]):
		num += 1
		if num < 129:
			symbol = "DOPC"
		else: symbol = "DPPC"
		newFrame.atoms.append(XYZAtom(symbol, *pos))
        #for atomIndex in range(posCount):
        #    newFrame.atoms.append(XYZAtom(atomsLists[0][atomIndex].symbol
        #                        , *self.lattice.positions[match[atomIndex]].x0))
        return ProjectedFrame(newFrame, refFrame, errors)




