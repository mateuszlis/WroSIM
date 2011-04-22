'''
Created on 18-01-2011

@author: lisu
'''
from myExceptions import InputError
from structures.xyzfile import XYZFrame
from grofile import GroFrame

class StructureFile(object):
    """Class responsible for handling structure file format"""
    
    def __init__(self, filename):
        """
        
        """
        self._filename = filename
        if self._filename[-3:] == 'gro':
            self._frameCreator = self._createGroFrame
            self._countNum = 1
        elif self._filename[-3:] == 'xyz':
            self._frameCreator = self._createXYZFrame
            self._countNum = 0
        else:
            raise InputError(self._filename, "Wrong filename type")    
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
        #FIXME: write loaders for both files
        if self._frames is None:
            self._frames = []
            try:
                structFile = open(self.filename, 'r')
            except:
                self._frames = []
                return
            lines = structFile.readlines()
            if len(lines) < 2: raise(InputError(lines,'Not enough lines in frame file'))
            atomsCount = int(lines[self._countNum])
            step = atomsCount + 2 
            i = 0 
            while (i + step) < len(lines) + 1:
                step = int(lines[i + self._countNum]) + 2
                self._frames.append(self._frameCreator(''.join(lines[i : i + step + self._countNum])[:-1]))               
                i += step
            self.currentFrame = len(self.frames) - 1
    
    def save(self):
        """Method saves current state of object into struct file"""
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
        atomsCount = int(frameLst[self._countNum])
        for _ in range(atomsCount + 1):
            lineStr = self.file.readline()
            if lineStr == '':
                print frameLst
                self.file.close()
                raise InputError('lineStr == ''',"Wrong number of lines in file")
            frameLst.append(lineStr)
        return self._frameCreator(''.join(frameLst))

    
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
        Returns struct file string from the currently stored structure"""
        if self._frames is None:
            self.load()
        return '\n'.join([repr(frame) for frame in self._frames])
    def __del__(self):
        if not self.file is None:
            self.file.close()
            
    def _createXYZFrame(self, frameStr):
        return XYZFrame(frameStr)
    def _createGroFrame(self, frameStr):
        return GroFrame(frameStr)