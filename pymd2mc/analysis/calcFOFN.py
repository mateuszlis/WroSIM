#!/usr/bin/python
"""
This program calculates fraction of first neighbors on trajectory images
"""

__author__ = 'Mateusz Lis'
__version__ = '0.1'

from optparse import OptionParser
import sys
from time import time

from myExceptions import InputError
from PIL import Image
import os
from random import randint
from structures.xyzfile import XYZFile
from sys import stdout
from utils import distance, clearFile, delLine
from xyz2hexnet.energyCalc import HexLatticeLoader

def main(argv=None):
    options = parseCommandLine()
    
    startTime = time()
    begin = int( options.begin )
    end = int( options.end )
    clearFile( options.datFilename )
    with open( options.datFilename, 'w' ) as output:
        try:
            interval = int( options.interval )
            for frameCounter in range( begin, end, interval ):
                if options.verbose: 
                    delLine()
                    print frameCounter,
                    stdout.flush()
                path = options.Filename + str( frameCounter ) + ".png" #TODO: png should be configurable
                im = Image.open( path )
                width, height = im.size
                correctionFactor = ( float( width - height ) / height )
                pixMap = im.load()
                correction = correctionFactor
                neighbors = [ ( 1, 0 ), ( 0, 1 ), ( -1, 0 ), ( 0, -1 ) ]
                allPixelsCount = 0
                similarPixelsCount = 0
                for h in range( 1, height - 1 ):
                    for w in range( 1, height - 1 ):
                        neighbor = neighbors[ randint( 0, len( neighbors ) - 1 ) ]
                        x = int( w + correction )
                        y = h
                        nX = x + neighbor[ 0 ]
                        nY = y + neighbor[ 1 ]
                        if pixMap[ x, y ] - pixMap[ nX, nY ] < 35:
                            similarPixelsCount += 1
                        allPixelsCount += 1
                    correction += correctionFactor
                fraction = float( similarPixelsCount ) / allPixelsCount
                output.write( str( frameCounter ) + " " + str( fraction )  + "\n" )
        except:
            pass
    
    if options.verbose:
        print "Execution time", time() - startTime



def parseCommandLine():
    """
    Sets up command line arguments and parses them
    """
    parser = OptionParser(usage="%prog ", version="%prog " + __version__,
                          description='''
    This program calculates fraction of first neighbors in image files ''')
    parser.add_option("-f", "--traj", dest="Filename", default="image_",
    help="image file name pattern", metavar="IMAGE_FILE")
    parser.add_option("-o", "--output", dest="datFilename",
            default="fraction_of_first_neighbors_image.dat",
    help='output dat file. WARNING: it will be overriden', metavar="DATFILE")
    parser.add_option("-b", "--begin", dest="begin", default="0",
    help="Begin number", metavar="start")
    parser.add_option("-e", "--end", dest="end", default="0",
    help="End number", metavar="start")
    parser.add_option("-i", "--interval", dest="interval", default="10000",
    help="End number", metavar="start")
    parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose", default=True,
    help="don't print status messages to stdout")

    (options, _) = parser.parse_args()

    return options
    

if __name__ == "__main__":
    sys.exit(main())
