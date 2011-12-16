#!/usr/bin/python
"""
This script calculates dG (difference of energy) from rdf
"""

__author__ = 'Mateusz Lis'
__version__ = '0.1'

from constants import R, T
import math
from optparse import OptionParser
import sys
from time import time


def main(argv=None):
    global options 
    options = parseCommandLine()

    startTime = time()

    fIn = options.rdfFilename
    fOut = options.datFilename
    outList = []

    with open( fIn, 'r') as rdfFile:
        for line in rdfFile.readlines():
            if not ( "#" in line or "@" in line ) :
                r, rdf = line.split()
                r, rdf = float( r ), float( rdf )
                if ( rdf != 0 ):
                    dG = - R * T * math.log( rdf )
                else:
                    dG = 0
                outList.append( ( r, dG ) )

    with open( fOut, 'w' ) as datFile:
        for r, dG in outList:
            datFile.write( "%s\t%s\n" % ( r, dG ) )

    print "Minimum energy:", min( [ x[1] for x in outList ] )
    if options.verbose:
        print "Execution time", time() - startTime
    
def parseCommandLine():
    """
    Sets up command line arguments and parses them
    """
    parser = OptionParser(usage="%prog ", version="%prog " + __version__,
                          description='''
This script calculates dG (difference of energy) from rdf
    ''')
    parser.add_option("-f", "--rdf", dest="rdfFilename", default="rdf.xvg",
    help="RDF trajectory (default rdf.xvg)", metavar="RDFFILE")
    parser.add_option("-o", "--output", dest="datFilename", default="dG.dat",
    help='Output file name', metavar="DATFILE")
    parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose", default=True,
    help="don't print status messages to stdout")

    (options, _) = parser.parse_args()

    return options
    

if __name__ == "__main__":
    sys.exit(main())
