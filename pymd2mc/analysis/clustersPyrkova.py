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
from structures.xyzfile import XYZFile, XYZFrame, XYZAtom
from xyz2hexnet.energyCalc import HexLatticeLoader
from utils import distance, clearFile, delLine

def main(argv=None):
    global options 
    options = parseCommandLine()

    xyzFile, xyzOutputFile = openFiles(options)
    frame = xyzFile.nextFrame()
    startTime = time()
    frameCounter = 0
    neighborPairs = {}
    while not frame is None:
        findNeighborPairs( frame, frameCounter, neighborPairs )
        clustered = findClusteredAtoms( neighborPairs, frameCounter, int( options.framesThr ) )
        if options.verbose: 
            delLine()
            print frameCounter,
            #l = neighborPairs.keys()
            #l.sort()
            #print "Neighbor pairs", [ (str(i) + ":" + str(neighborPairs[i])) for i in l ]
            #print  "Clusters", identifyClusters( neighborPairs, frameCounter, 1, frame ) 
            #print "Clustered atoms", findClusteredAtoms( neighborPairs, frameCounter, 1 )
        frameCounter += 1
        newFrame = XYZFrame()

        clusteredAtoms = []
        interAtoms = []
        nonClustAtoms = []
        for atomInd in range( len( frame.atoms ) ):
            if atomInd in clustered['clustered']:
                clusteredAtoms.append( XYZAtom( "Clustered", *frame.atoms[atomInd].x0 ) )
            elif atomInd in clustered['interchanging']:
                interAtoms.append( XYZAtom( "Inter", *frame.atoms[atomInd].x0 ) )
            else:
                nonClustAtoms.append( XYZAtom( "NonCl", *frame.atoms[atomInd].x0 ) )
        maxNum = 64
        for i in range( maxNum - len( clusteredAtoms ) ):
            clusteredAtoms.append( XYZAtom( "Clustered", 0, 0, 0, ) )
        for i in range( maxNum - len( interAtoms ) ):
            interAtoms.append( XYZAtom( "Inter", 0, 0, 0, ) )
        for i in range( maxNum - len( nonClustAtoms ) ):
            nonClustAtoms.append( XYZAtom( "Clustered", 0, 0, 0, ) )
        newFrame.atoms = clusteredAtoms + interAtoms +  nonClustAtoms 
        xyzOutputFile.addFrame( newFrame )
        frame = xyzFile.nextFrame()
    if options.verbose:
        print "Execution time", time() - startTime
    xyzOutputFile.save()
    



def findNeighborPairs( frame, frameNum, neighborPairs):
    """ Find Neighbor Pairs:
    @param neighborPairs has the following format:
    { (atomIndex1, atomIndex2) -> [ numOfFrames, lastFrameNum ] }
    denoting how long two atoms have been together and when was the last frame 
    """
    doneInd = []
    for atomInd in range( len( frame.atoms ) ):
        doneInd.append( atomInd )
        for neigh in getNeighbors( atomInd, frame ):
            if neighborPairs.has_key( ( atomInd, neigh ) ):
                updateState( neighborPairs[ ( atomInd, neigh ) ], frameNum )
            elif not neigh in doneInd:
                neighborPairs[ ( atomInd, neigh ) ] = [ 1, frameNum ]
def getNeighbors( atomInd, frame ):
    for otherAtom in range( len( frame.atoms ) ):
        if distance( frame.atoms[ atomInd ].x0, frame.atoms[ otherAtom ].x0 ) < float( options.distance ):
           yield otherAtom
def updateState( state, frameNum ):
    if state[1] + 1 < frameNum: #if atoms were not together in last frame
        state[0] = 1 #reset time period
    else:
        state[0] += 1 #else add another frame
    state[1] = frameNum

def identifyClusters( neighborPairs, frameNum, minNumberOfFrames, frame ):
    clusters = []
    for atomsPair, value in neighborPairs.iteritems():
        atomInd1, atomInd2 = atomsPair
        numOfFrames, lastFrameNum = value
        areIn = False
        if lastFrameNum == frameNum and minNumberOfFrames < numOfFrames:
            for cluster in clusters:
                if atomInd1 in cluster or atomInd2 in clusters:
                    areIn = True
                    cluster |= set( [ atomInd1, atomInd2 ] )
            if not areIn:
                clusters.append( set( [ atomInd1, atomInd2 ] ) )
    return integrateClusters( frame, clusters )

def integrateClusters( frame, clusters ):
    for atomInd in range( len( frame.atoms ) ):
        integratedClustersList = []
        for clusterInd in range( len( clusters ) ):
            if atomInd in clusters[clusterInd]:
                integratedClustersList.append( clusterInd )
        for intClusters in integratedClustersList[1:]:
            clusters[integratedClustersList[0]] |= clusters[intClusters]
            del( clusters[intClusters] )
    return clusters

def findClusteredAtoms( neighborPairs, frameNum, minNumberOfFrames ):
    clusteredIndices = set()
    interchangingIndices = set() 
    for atomsPair, val in neighborPairs.iteritems():
        atomInd1, atomInd2 = atomsPair
        numOfFrames, lastFrameNum = val
        if lastFrameNum == frameNum and minNumberOfFrames < numOfFrames:
            clusteredIndices |= set( atomsPair )
        if lastFrameNum == frameNum and minNumberOfFrames >= numOfFrames:
            interchangingIndices |= set( atomsPair )
        if lastFrameNum < frameNum and frameNum - lastFrameNum < minNumberOfFrames:
            interchangingIndices |= set( atomsPair )
    interchangingIndices &= clusteredIndices
    clusteredIndices -= interchangingIndices
    return {"clustered" : clusteredIndices, "interchanging" : interchangingIndices }

        
def openFiles(options):
    xyzFile = XYZFile(options.xyzFilename)
    clearFile( options.xyzOutput )
    outFile = XYZFile( options.xyzOutput )
    return xyzFile, outFile


def parseCommandLine():
    """
    Sets up command line arguments and parses them
    """
    parser = OptionParser(usage="%prog ", version="%prog " + __version__,
                          description='''
    This program calculates histogram of cluster sizes in xyz file ''')
    parser.add_option("-f", "--traj", dest="xyzFilename", default="traj.xyz",
    help="xyz trajectory file (default traj.xyz)", metavar="XYZFILE")
    parser.add_option("-d", "--distance", dest="distance", default="1.0",
    help='Threshold distance of the neighboring atoms', metavar="DISTANCE")
    parser.add_option("-t", "--time", dest="framesThr", default="1",
    help='Threshold time (frames number) for clusters formation', metavar="TIME")
    parser.add_option("-o", "--output", dest="xyzOutput", default="clustersTraj.xyz",
    help='output xyz file. with "size\t freq " format WARNING: it will be overriden', metavar="OUTPUTFILE")
    parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose", default=True,
    help="don't print status messages to stdout")

    (options, _) = parser.parse_args()

    return options
    

if __name__ == "__main__":
    sys.exit(main())
