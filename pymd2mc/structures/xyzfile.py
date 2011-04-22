# pymd2mc.xyzfile
"""Module responsible for handling the xyz file format. xyz file format stores 
   trajectories for simulations. For further information see 
   http://en.wikipedia.org/wiki/XYZ_file_format"""
__author__ =  'Mateusz Lis'
__version__=  '0.1'

from structures.baseStructs import Atom, Frame
from myExceptions import InputError

 

class XYZAtom(Atom):
    "Class represents atom stored in XYZ file"""
    def __init__(self, symbol, x = 0.0, y = 0.0, z = 0.0):
        if len(symbol) > 10:
            self._construct(symbol)
        else:
            Atom.__init__(self, symbol, x, y, z)
        
    def __repr__(self):
        return "%s\t%.8f\t%.8f\t%.8f" % (self.symbol, self.x, self.y, self.z)
    
    def _construct(self, atomStr):
        atomList = atomStr.split()
        if len(atomList) != 4:
            raise InputError(atomStr, 'Wrong number of arguments while \
                                    constructing the XYZAtom')
        self.symbol, self.x = atomList[0].strip(), float(atomList[1])
        self.y, self.z = float(atomList[2]), float(atomList[3])
        
    
    
class XYZFrame(Frame):
    """Class representing frame in trajectory file"""
    def __init__(self, frameStr = None):
        """ 
        If second parameter is not given, constructs empty frame, otherwise
        constructs frame from the string which satisfies following pattern
        <number of atoms>
        comment line
        atom_symbol11 x-coord11 y-coord11 z-coord11
        atom_symbol12 x-coord12 y-coord11 z-coord12
        ...
        """
        Frame.__init__(self)
        if not frameStr is None:
            lines = frameStr.split('\n')
            if lines[-1] == '': del(lines[-1])
            if len(lines) < 3:
                raise InputError(frameStr, 'Wrong number of arguments while \
                                    constructing the XYZFrame')
            length = int(lines[0])
            self._comment = lines[1]
            for line in lines[2:]:
                self.atoms.append(XYZAtom(line))
            if len(self.atoms) != length:
                raise InputError(frameStr, "Wrong number of lines in file: length=%d and len(self.atoms) %d" % (length, len(self.atoms)))

    
    def __repr__(self):
        """Produces frame structure string described in constructor doc"""
        reprStr = "%d\n%s\n%s" % (len(self.atoms), self.comment,
                                  '\n'.join([repr(atom) for atom in self.atoms]))
        return reprStr
    
    @property
    def boxSize(self):
        """
        Returns size of the box from frame.
        This property is only present in specific kind of xyz files (stored in
        comment). 
        If boxSize is not present, raises input exception
        """
        try:
            return  [float(size) for size in self.comment.split(",")]
        except:
            raise(InputError(self.comment, "Cannot read size of box from comment."))
    @property
    def boxVectors(self):
        try:
            vList = [float(size) for size in self.comment.split(",")]
            for _ in range(len(vList), 9):
                vList.append(0)
            v1 = [vList[0], vList[3], vList[4]]
            v2 = [vList[5], vList[1], vList[6]]
            v3 = [vList[7], vList[8], vList[2]]
            return [v1,v2,v3]
        except:
            raise(InputError(self.comment, "Cannot read size of box from comment."))
    @boxVectors.setter
    def boxVectors(self, boxSizeVectors):
        v1,v2,v3 = boxSizeVectors
        self.comment = "%f, %f, %f, %f, %f, %f, %f, %f, %f" % (v1[0], 
                                                               v2[1], 
                                                               v3[2], 
                                                               v1[1], 
                                                               v1[2], 
                                                               v2[0], 
                                                               v2[2], 
                                                               v3[0], 
                                                               v3[1])
            

class XYZFile(object):
    """Class responsible for handling xyz file format"""
    
    def __init__(self, filename):
        """
        
        """
        self._filename = filename
        self._frames = None
        self.currentFrame = None
        self.file = None
        
    @property
    def filename(self):
        return self._filename
    @filename.setter
    def filename(self, value):
        self._filename = value
        
    
    def load(self):
        """ Loads whole file into memory."""
        if self._frames is None:
            self._frames = []
            try:
                xyzFile = open(self.filename, 'r')
            except:
                self._frames = []
                return
            lines = xyzFile.readlines()
            if len(lines) < 2: raise(InputError(lines,'Not enough lines in frame file'))
            atomsCount = int(lines[0])
            step = atomsCount + 2
            i = 0 
            while (i + step) < len(lines) + 1:
                step = int(lines[i]) + 2
                self._frames.append(XYZFrame(''.join(lines[i : i + step])[:-1]))               
                i += step
            self.currentFrame = len(self.frames) - 1
    
    def save(self):
        """Method saves current state of object into xyz file"""
        if self._frames is None:
            self.load()
        file  = open(self.filename, 'w')
        file.write(repr(self))
        file.close()
        
    def getFrame(self, number):
        """
        Returns frame of the specified number. Loads whole file into 
        memory
         """
        if self._frames is None:
            self.load()
        return self._frames[number]
    
    @property
    def frames(self):
        """
        Returns all frames from the file. Loads whole file into 
        memory
        """
        if self._frames is None:
            self.load()
        return self._frames
    def nextFrame(self):
        """ 
        Returns next frame from the file. If file has ended returns None 
        value
        """
        if self.currentFrame is None:
            self.currentFrame = 0
            self.file = open(self.filename, 'r')
        else:
            self.currentFrame += 1
        frameLst = [self.file.readline()]
        if frameLst[0] == '':
            self.file.close()
            return None
        atomsCount = int(frameLst[0])
        for _ in range(atomsCount + 1):
            lineStr = self.file.readline()
            if lineStr == '':
                print frameLst
                self.file.close()
                raise InputError('lineStr == ''',"Wrong number of lines in file")
            frameLst.append(lineStr)
        return XYZFrame(''.join(frameLst))

    
    def addFrame(self, frame):
        """
        Method adds frame after current. If file is not loaded,
        it opens the file and saves frame in proper place
        """
        if self.currentFrame is None:
            file = open(self.filename, 'a')
            file.write("%s\n" % repr(frame))
            file.close()
        else:
            self.frames.insert(self.currentFrame, frame)
        pass  
    def getCurrentFrameNumber(self):
        """
        If whole file is not loaded into memory (using just nextFrame)
        this points frame number on which nextFrame lastly pointed
        """
        return self.currentFrame
   
    def __repr__(self):
        """
        Returns xyz file string from the currently stored structure"""
        if self._frames is None:
            self.load()
        return '\n'.join([repr(frame) for frame in self._frames])
    def __del__(self):
        if not self.file is None:
            self.file.close()
