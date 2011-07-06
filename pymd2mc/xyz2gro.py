#!/usr/bin/python
"""
This program converts the xyz file to gromacs gro file
"""
__author__ =  'Mateusz Lis'
__version__=  '0.1'

from copy import deepcopy
from optparse import OptionParser
import sys
from time import time

from structures.xyzfile import XYZFile
from utils import delLine
def main():
    options = parseCommandLine()
    xyzFile = XYZFile(options.xyzFilename)
    for frame in xyzFile.frames:
        print ''
        print len(frame.atoms)
        i = 1
        for atom in frame.atoms:
            print ("%5d%5s%5s%5d%8.3f%8.3f%8.3f" % (i,'Mol  ',atom.symbol,i,atom.x, atom.y, atom.z))
            i += 1
#            if options.verbose:
#                delLine()
#                print i, 
        #print '3.061 5.304 0.2'       
        for size in [13.755, 13.755, 6.625]: print size,
        print ''
        
    #with open("x.gro", 'w') as f:
     #   f.write("%5d%5s%5s%5d%8.3f%8.3f%8.3f%8.4f%8.4f%8.4f" % (1,'water','P',8,1.22,2.22,3.44,5.33,5.11,6.22))


def parseCommandLine():
    """
    Sets up command line arguments and parses them
    """
    parser = OptionParser(usage="%prog ", version="%prog " + __version__,
                          description='''
    This program converts the xyz file to xtc and gro files''')
    parser.add_option("-o", "--output", dest="groFilename",default="conf.gro",
    help="output gro file (default: conf.gro)", metavar="GROFILE")
    parser.add_option("-f", "--xyz", dest="xyzFilename",default="hexTraj.xyz",
    help="xyz trajectory file (default hexTraj.xyz)", metavar="XYZFILE")
    
    parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose", default=True,
    help="don't print status messages to stdout")

    (options, _) = parser.parse_args()

    return options

if __name__ == "__main__":
    sys.exit(main())