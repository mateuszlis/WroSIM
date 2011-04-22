#!/usr/bin/python
# pymd2mc.xyzfile
"""
This program converts raw xyz file to xyz file with hexagonal network of particles

"""

__author__ =  'Mateusz Lis'
__version__=  '0.1'


from optparse import OptionParser
import sys
from time import time

from utils import delLine, clearFile
from latticeProjector import LatticeProjectorSimple
from lattices import HexLattice
from structures.xyzfile import XYZFile
 

  
def main():
    
    options = parseCommandLine()
    avgArea = getAvgSize(options.inXyzFilename)
    lattBoxSize = ((2 * avgArea) / (3**(1/2.)) )**(1/2.)
    inFile = XYZFile(options.inXyzFilename)
    
    clearFile(options.outXyzFilename)
    outFile = XYZFile(options.outXyzFilename)
    i = 0
    startTime = time()
    errors = []
    latt = HexLattice(11, (0,0), (lattBoxSize, lattBoxSize))
    lattProj = LatticeProjectorSimple(inFile, latt)
    while True:
        i += 1
        if options.verbose:
            delLine()
            print i, 
        projLattice = lattProj.next()
        if projLattice is None:
            break
        frame = projLattice.frame
        if options.reference:
            frame.atoms.extend(projLattice.refFrame.atoms)
        err = (i, max(projLattice.errors), sum(projLattice.errors) / len(projLattice.errors))
        errors.append("{0} {1} {2}".format(*err))
        
        
        outFile.addFrame(frame)
        
    if not options.errFile is None:
        with open(options.errFile, 'w') as errFile:
            errFile.write('\n'.join(errors))
        
    outFile.save()
    if options.verbose: 
        print "Done. Execution time=%f" % (time() - startTime)    

def getAvgSize(xyzFilename):
    xyzFile = XYZFile(xyzFilename)
    i, avgArea = 0, 0
    atomsLength = 0
    while True:
        frame = xyzFile.nextFrame()
        if frame is None:
            break
        atomsLength = len(frame.atoms)
        avgArea += frame.boxSize[0] * frame.boxSize[1]
        i += 1
    
    avgArea /= i
    print avgArea
    print "Avg APL = %f" % (avgArea / atomsLength)
    del(xyzFile)
    return avgArea
    
              
def parseCommandLine():
    """
    Sets up command line arguments and parses them
    """
    parser = OptionParser(usage="%prog ", version="%prog " + __version__,
                          description='''
    This program converts raw xyz file to xyz file with hexagonal network of particles''')
    parser.add_option("-f", "--traj", dest="inXyzFilename",default = "traj.xyz",
    help="xyz input trajectory file (default traj.xyz)", metavar="INXYZFILE")
    parser.add_option("-o", "--output", dest="outXyzFilename", default="hexTraj.xyz",
    help="output xyz file. WARNING: it will be overriden", metavar="OUTXYZFILE")
    parser.add_option("-e", "--error", dest="errFile", default="errors.dat",
    help="Error output file. WARNING: it will be overriden", metavar="ERRFILE")
    parser.add_option("-r", "--reference",
    action="store_true", dest="reference", default=False,
    help="write reference frames to xyz output file")

    parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose", default=True,
    help="don't print status messages to stdout")

    (options, _) = parser.parse_args()

    return options 
 


if __name__ == '__main__':
    sys.exit(main())