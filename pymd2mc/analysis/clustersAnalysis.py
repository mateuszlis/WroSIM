#!/usr/bin/python
"""
This program calculates clusters size histogram from xyz file 
"""

__author__ = 'Mateusz Lis'
__version__ = '0.1'

from optparse import OptionParser
import sys
from time import time

from myExceptions import InputError
from structures.xyzfile import XYZFile
from xyz2hexnet.energyCalc import HexLatticeLoader

def main(argv=None):
    options = parseCommandLine()

    xyzFile, datFile = openFiles(options)
    frame = xyzFile.nextFrame()
    startTime = time()
    clusterHists = []
    frameCounter = 0
    latticeLoader = HexLatticeLoader(frame)
    sys.setrecursionlimit(max(1000, len(frame.atoms))) #this is cheating
    while not frame is None:
        latticeLoader.updateState(frame, "")
        clusterHists.append(calcClusterHist(latticeLoader, frame))
        frame = xyzFile.nextFrame()
        if options.verbose: 
            print frameCounter, clusterHists[-1] 
        frameCounter += 1
    finalHist = {}
    for histDict in clusterHists:
        for atomName, hist in histDict.items():
            if atomName == options.name or options.name == "all":
                for size, count in hist.items():
                    if finalHist.has_key(size):
                        finalHist[size] += count
                    else: 
                        finalHist[size] = count
    for size in finalHist.keys():
        finalHist[size] /= float(frameCounter)
    for size, count in finalHist.items():
        datFile.write("%s\t%s\n" % (size, count))
    datFile.close()
    if options.verbose:
        print "Execution time", time() - startTime

def calcClusterHist(latticeLoader, frame):
    atomDict = {}
    done = {}
    for atom in frame.atoms:
        #print len(done.keys())
        if not done.has_key(atom):
            size = checkCluster(latticeLoader, frame, atom, done)
            #print "done", atom, size
            if not atomDict.has_key(atom.symbol):
                atomDict[atom.symbol] = {size : 1}
            elif atomDict[atom.symbol].has_key(size):
                atomDict[atom.symbol][size] += 1
            else: atomDict[atom.symbol][size] = 1
    return atomDict

def checkCluster(latticeLoader, frame, atom, done):
    clusterSize = 1
    done[atom] = 1
    #print atom
    neighbrs = latticeLoader.getNeighbors(atom)
    for neighb in neighbrs:
        if (not done.has_key(neighb)) and neighb.symbol == atom.symbol:
            clusterSize += checkCluster(latticeLoader, frame, neighb, done)
    return clusterSize
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
    This program calculates histogram of cluster sizes in xyz file ''')
    parser.add_option("-f", "--traj", dest="xyzFilename", default="traj.xyz",
    help="xyz trajectory file (default traj.xyz)", metavar="XYZFILE")
    parser.add_option("-o", "--output", dest="datFilename", default="clusHist.dat",
    help='output dat file. with "size\t freq " format WARNING: it will be overriden', metavar="DATFILE")
    parser.add_option("-n", "--name", dest="name", default="all",
    help="Name of reference atom (if you set this, clusters of one type will be calculated)", metavar="Name")
    parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose", default=True,
    help="don't print status messages to stdout")

    (options, _) = parser.parse_args()

    return options
    

if __name__ == "__main__":
    sys.exit(main())
