#!/usr/bin/python
"""
This program calculates radial distribution function from xyz file 
"""

__author__ =  'Mateusz Lis'
__version__=  '0.1'

from optparse import OptionParser
import sys
from time import time

from myExceptions import InputError
from utils import bruteForce, delLine, distance, vectorLength
from structures.xyzfile import XYZFile, XYZFrame, XYZAtom

def main(argv = None):
    
    options = parseCommandLine()

    xyzFile, datFile = openFiles(options)
    frame = xyzFile.nextFrame()


    frameCounter= 1
    binSize = 0.02
    boxSize = max([vectorLength(v) for v in frame.boxVectors])
    hist = calcHist(frame, propagateThroughPBC(frame.atoms, frame.boxVectors), binSize, boxSize, options)
    
    frame = xyzFile.nextFrame()
    startTime = time()
    while not frame is None:
        if options.verbose:
            delLine()
            print frameCounter,
        newHist = calcHist(frame, propagateThroughPBC(frame.atoms, frame.boxVectors), binSize, boxSize, options)
        for i in range(len(hist)):
            hist[i] += newHist[i]
            
        frameCounter += 1
        frame = xyzFile.nextFrame()
    for i in range(len(hist)):
        datFile.write('%f %f\n' %(i * binSize, hist[i] / float(frameCounter)))
    datFile.close()
    print time() - startTime
#    atoms = [XYZAtom('N')]
#    for atomList in propagateThroughPBC(atoms, [[1,0,0],[0,1,0],[0,0,0]]):
#        for i in atomList:
#            print i.x0
        
    
def calculateRDF(xyzFile):
    pass

def calcHist(frame, extendedAtomsIter, binSize, boxSize, options):
    bins = [0 for _ in range(int(boxSize / (2 * binSize)))]
    extendedAtoms = [atom for atom in extendedAtomsIter]
    for atom in frame.atoms:
        if atom.symbol == options.name1 or options.name1 == 'all':
            for atomLst in extendedAtoms:
                for atom2 in atomLst:
                    if options.name2 == "all" or options.name2 == atom2.symbol:
                        if distance(atom.x0, atom2.x0) < boxSize / 2 and distance(atom.x0, atom2.x0) != 0:
                            for i in range(len(bins)):
                                if i * binSize > distance(atom.x0, atom2.x0):
                                    bins[i] += 1
    return bins
        

def propagateThroughPBC(atoms, boxVectors):
    
    v = boxVectors
    
    for pbcSite in bruteForce(2, [0, 1, -1]):
        newAtoms = []
        for atom in atoms:
            newAtoms.append(XYZAtom(atom.symbol, *atom.x0))
            for atomCoord in range(2):
                for vectorCoord in range(2):
                    x0 = newAtoms[-1].x0
                    x0[atomCoord] += v[vectorCoord][atomCoord] * pbcSite[vectorCoord]
                    newAtoms[-1].x0 = x0
        yield newAtoms


def openFiles(options):
    try:
        xyzFile = XYZFile(options.xyzFilename)
        datFile = open(options.datFilename, 'w')
    except:
        raise(InputError(options.xyzFilename, "error while opening xyz or output file"))
    return xyzFile, datFile


def parseCommandLine():
    """
    Sets up command line arguments and parses them
    """
    parser = OptionParser(usage="%prog ", version="%prog " + __version__,
                          description='''
    This program calculates radial distribution function from xyz file ''')
    parser.add_option("-f", "--traj", dest="xyzFilename",default="traj.xyz",
    help="xyz trajectory file (default traj.xyz)", metavar="XYZFILE")
    parser.add_option("-o", "--output", dest="datFilename",default="rdf.dat",
    help="output dat file. with rdf vals WARNING: it will be overriden", metavar="DATFILE")
    parser.add_option("--n1", "--name1", dest="name1",default="all",
    help="Name of atom to be counted distance from", metavar="Name2")
    parser.add_option("--n2", "--name2", dest="name2",default="all",
    help="Name of atom to be counted distance to", metavar="Name2")
    
    parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose", default=True,
    help="don't print status messages to stdout")

    (options, _) = parser.parse_args()

    return options
    

if __name__ == "__main__":
    sys.exit(main())