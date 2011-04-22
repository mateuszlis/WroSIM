#!/usr/bin/python
"""
This program converts the xtc file to xyz file using atoms in given index file
"""
__author__ =  'Mateusz Lis'
__version__=  '0.1'

from optparse import OptionParser
import sys
from time import time

from pymacs import Model
from pymacs.ndx import IndexFile
from pymacs import openXTC, readXTCFrame


from utils import delLine, clearFile
from structures.xyzfile import *

def main(argv = None):
    
    options = parseCommandLine()

    #===========================================================================
    # opening files and so on
    #===========================================================================  
    try:
        model = Model().read(options.groFilename)
        fp = openXTC(options.xtcFilename)
        clearFile(options.xyzFilename)
        xyzFile = XYZFile(options.xyzFilename)
    except:
        raise(InputError(options.xyzFilename, "error while opening gro, xtc or xyz file"))
    
    if not options.ndxFilename is None:
        try:
            ndxFile = IndexFile(fname = options.ndxFilename)
        except:
            raise(InputError("error while opening ndx file"))
        
        atomsIndexes = []
        for (_, atomList) in ndxFile.dic.items():
            atomsIndexes += atomList
    else:
        atomsIndexes = range(1, len(model.atoms) + 1)
        
    createXYZFile(fp, model, atomsIndexes, xyzFile, options)
    
    
def createXYZFile(fp, model, atomsIndexes, xyzFile, options):
    """
    Function utilizes frames, model and atomsIndexes to create xyzFile 
    """
    startTime = time()
    while True:                     # loop over trajectory
        xyzFrame = XYZFrame()
        
        xyzFrame.comment = "%.3f,%.3f,%.3f" %  (model.box[0][0], 
                                          model.box[1][1], model.box[2][2]) 
        #comment = boxSize
        
        frame = readXTCFrame(fp)
        
        if not frame:
            break               # end of trajectory or error
        model.update(frame)      # update the model data with coords, box
        
        
        for atomIndex in atomsIndexes:
                model.atoms[atomIndex - 1].x[2] = 0.
                xyzAtom = XYZAtom(model.atoms[atomIndex - 1].name, 
                                    *model.atoms[atomIndex - 1].x)
                xyzFrame.atoms.append(xyzAtom)
        xyzFile.addFrame(xyzFrame)
        if options.verbose:
            delLine()
            print "time = %g step = %g" % (frame['time'],frame['step']),

    xyzFile.save()
    print "Done. Execution time=%f" % (time() - startTime)    
    
def parseCommandLine():
    """
    Sets up command line arguments and parses them
    """
    parser = OptionParser(usage="%prog ", version="%prog " + __version__,
                          description='''
    This program converts the xtc file to xyz file using given index file''')
    parser.add_option("-s", "--grofile", dest="groFilename",default="conf.gro",
    help="gro file (default: conf.gro)", metavar="GROFILE")
    parser.add_option("-f", "--traj", dest="xtcFilename",default="traj.xtc",
    help="xtc trajectory file (default traj.xtc)", metavar="XTCFILE")
    parser.add_option("-n", "--index", dest="ndxFilename",
    help="Index file with atoms specified to be included in xyz file", metavar="NDXFILE")
    parser.add_option("-o", "--output", dest="xyzFilename",default="conf.xyz",
    help="output xyz file. WARNING: it will be overriden", metavar="XYZFILE")
    
    parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose", default=True,
    help="don't print status messages to stdout")

    (options, _) = parser.parse_args()

    return options
    

if __name__ == "__main__":
    sys.exit(main())