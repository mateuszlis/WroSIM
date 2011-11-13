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
from latticeProjector import *
from lattices import HexLattice
from structures.xyzfile import XYZFile
 

  
def main():
    
    options = parseCommandLine()
    avgArea = getAvgSize(options.inXyzFilename)
    print avgArea
    lattBoxSize = ((2 * avgArea) / (3**(1/2.)) )**(1/2.) #area of parallelogram
    print lattBoxSize
    inFile = XYZFile(options.inXyzFilename)
    latt = HexLattice(int(options.lattSize), (0,0), (lattBoxSize, lattBoxSize))

    lattProj = NearestNeighborLatticeProj(inFile, latt)
    clearFile(options.outXyzFilename)
    outFile = XYZFile(options.outXyzFilename)
    i = 0
    startTime = time()
    errors = []
    #lattProj = LatticeProjectorSimple(inFile, latt)
    while True:
        if options.verbose:
	    delLine()
	    if i > 0: midTime = (time() - startTime) / i
	    else: midTime = 0
	    print i, "Avg time per frame: ", midTime,
	    sys.stdout.flush()
        projLattice = lattProj.next()
        if projLattice is None:
            break
        frame = projLattice.frame
        if options.reference:
            frame.atoms.extend(projLattice.refFrame.atoms)
	    symb = "DOPC"
	    length = len(frame.atoms)
	    num = 0
	    for atom in frame.atoms:
                if atom.symbol ==  symb: num += 1
            for j in range(15000 - num):
	        frame.atoms.append(XYZAtom("DOPC", -10., 0., 0.))
	    for j in range(15000 - (length - num)):
		frame.atoms.append(XYZAtom("DPPC", -10., 0., 0.))
	    atoms = sorted(frame.atoms, key=lambda w: w.symbol)
	    frame.atoms = atoms
	    #frame.atoms = []#sorted(frame.atoms, key=lambda w: w.symbol)
        i += 1

        err = (i, max(projLattice.errors), sum(projLattice.errors) / len(projLattice.errors))
        errors.append("{0} {1} {2}".format(*err))
        
        
        outFile.addFrame(frame)
 	if i > 10: break
        
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
    #print avgArea
    #print "Avg APL = %f" % (avgArea / atomsLength)
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
    parser.add_option( "--lattice-size",
     dest="lattSize", default=16,
    help="Number denoting size of the lattice")


    parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose", default=True,
    help="don't print status messages to stdout")

    (options, _) = parser.parse_args()

    return options 
 


if __name__ == '__main__':
    sys.exit(main())
