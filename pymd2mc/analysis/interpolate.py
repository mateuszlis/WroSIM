#!/usr/bin/python
# pymd2mc.
"""
This program converts xyz trajectory file into serie of interpolated images with
separated domains of given atoms.
"""

__author__ =  'Mateusz Lis'
__version__=  '0.1'

import Image
from itertools import product
from math import sqrt
from optparse import OptionParser
import sys
from time import time

from utils import delLine, clearFile, distance
from structures.xyzfile import XYZFile
 

  
def main():
    options = parseCommandLine()
    symbols, atomNum = identifyAtomTypes( options.inXyzFilename, options )
    print symbols 
    inFile = XYZFile( options.inXyzFilename )
    clearFile( options.outPngFileTemp )
    i = 0
    startTime = time()
    
    frame = inFile.nextFrame()

    colours = [ 
            ( 255, 255, 255 ) #white
            , ( 150, 150, 150 ) # grey
            , ( 0, 0, 0 ) # black
            , ( 150, 0, 0) # blue
            , ( 0, 150, 0 ) # green
            ]
    sizeX, sizeY = 300, 300
    countX, countY = 25, 25
    while not frame is None:
        stdoutStep( options, i )
        i += 1

        #prepare image
        im = Image.new( "RGBA", ( countX, countY ) )
        pix = im.load()
        
        rescaleFrameToImage( frame, countX, countY )
        for coords in product( range( countX ), range( countY ) ):
            pix[coords[0], coords[1]] =  colours[symbols.index( findClosestAtom( frame,
                coords ) )]
            #print findClosestAtom( frame, coords), coords
        for atom in frame.atoms:
            print atom.x, atom.y
            pix[ int( atom.x ), int( atom.y )] = colours[-1]

        im = im.resize( ( sizeX, sizeY ), Image.ANTIALIAS )
        im = im.transpose( Image.FLIP_TOP_BOTTOM )
        #im.show()
        im.save( options.outPngFileTemp + str(i) + ".png" )
        frame = inFile.nextFrame()
        
        
    
    if options.verbose: 
        print "Done. Execution time=%f" % (time() - startTime)    
      

def findClosestAtom( frame, coords ):
    extendedCoords = ( coords[0], coords[1], 0 )
    minDist = distance(  extendedCoords , frame.atoms[0].x0 )
    minName = frame.atoms[0].symbol

    for atom in frame.atoms:
        if distance( extendedCoords, atom.x0 ) < minDist:
            minDist = distance( extendedCoords, atom.x0 ) 
            minName = atom.symbol
    return minName

        
def rescaleFrameToImage( frame, sizeX, sizeY ):
    maxX = frame.atoms[0].x
    minX = maxX
    maxY = frame.atoms[0].y
    minY = maxY
    for atom in frame.atoms:
        if atom.x < minX: minX = atom.x
        if atom.x > maxX: maxX = atom.x
        if atom.y < minY: minY = atom.y
        if atom.y > maxY: maxY = atom.y
    for atom in frame.atoms: # move atoms to appropriate positions
        atom.x -= minX
        atom.y -= minY
    print maxX, maxY, ( ( sizeX ) / maxX )
    scaleFactorX = ( sizeX - 1 ) / ( maxX - minX )
    scaleFactorY = ( sizeY - 1 ) / ( maxY - minY )
    for atom in frame.atoms: # resize atom positions to fit image
        atom.x *= scaleFactorX
        atom.y *= scaleFactorY

def identifyAtomTypes( xyzFilename, options ):
    file = XYZFile( xyzFilename )
    frame = file.nextFrame()
    symbols = set()
    atomNum = len( frame.atoms )
    frameCntr = 0
    if options.verbose:
        print "Preparations for image creation"
    while not frame is None:
        if  frameCntr % options.freq == 0:
            for atom in frame.atoms:
                l = [ atom.symbol ] # cause set( atom.symbol ) would split into
                                    # letters
                symbols |= set( l )
            frameCntr += 1
            stdoutStep( options, frameCntr )
            if atomNum < len( frame.atoms ):
                atomNum = len( frame.atoms )
        frame = file.nextFrame()
    symbols = list( symbols )
    symbols.sort()
    return symbols, atomNum

###############################################################################


def stdoutStep( options, i ):
    """ Outputs progress information if options.verbose is set"""
    if options.verbose:
        delLine()
        print i, 
        sys.stdout.flush()
        
def parseCommandLine():
    """
    Sets up command line arguments and parses them
    """
    parser = OptionParser(usage="%prog ", version="%prog " + __version__,
                          description='''
    This program creates png images out of trajectory written in xyz file''')
    parser.add_option("-f", "--traj", dest="inXyzFilename",default = "traj.xyz",
    help="xyz input trajectory file (default traj.xyz)", metavar="INXYZFILE")
    parser.add_option("-o", "--output", dest="outPngFileTemp", default="inter",
    help="Frequency of output ", metavar="FREQ")
    parser.add_option("-p", "--freq", dest="freq", default=1,
    help="pattern for output png files. WARNING: it will be overriden",
    metavar="OUTPNGFILE")
    parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose", default=True,
    help="don't print status messages to stdout")

    (options, _) = parser.parse_args()

    return options 
 


if __name__ == '__main__':
    sys.exit(main())
