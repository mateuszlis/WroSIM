#!/usr/bin/python
"""
This program calculates moving average of a data file
"""

__author__ =  'Mateusz Lis'
__version__=  '0.1'

from optparse import OptionParser
import sys
from time import time

from myExceptions import InputError
from utils import bruteForce, delLine, distance, vectorLength, movingAverage
from structures.xyzfile import XYZFile, XYZFrame, XYZAtom

def main(argv = None):
    
    options = parseCommandLine()
    options.c, options.n = int(options.c), int(options.n)
    
    datFile, outFile = openFiles(options)
    
    valsList = []
    newFileList = []
    for line in datFile.readlines():
        lineContents = line.split()
        val = float(lineContents[options.c])
        valsList.append(val)
        newFileList.append(lineContents)
    timeDist = float(newFileList[1][0]) * options.n / 2
    newFileIter = iter(newFileList)    
    for i in movingAverage(valsList, options.n):
        line = newFileIter.next()
        line[0] = str(float(line[0]) + timeDist)
        line[options.c] = "%f" % (i)
        outFile.write("\t".join(line))
        outFile.write('\n')
        
         
    datFile.close()
    outFile.close()


    
        
    
def openFiles(options):
    try:
        outFile = open(options.outFilename, 'w')
        datFile = open(options.datFilename)
    except:
        raise(InputError(options.outFilename, "error while opening output file"))
    return  datFile, outFile


def parseCommandLine():
    """
    Sets up command line arguments and parses them
    """
    parser = OptionParser(usage="%prog ", version="%prog " + __version__,
                          description='''
    This program calculates moving average of a data file ''')
    parser.add_option("-f", "--file", dest="datFilename",default="data.dat",
    help="Data file (with numbers in columns) (default data.dat)", metavar="DATFILE")
    parser.add_option("-o", "--output", dest="outFilename",default="mvAvg.dat",
    help="output dat file. with moving average data WARNING: it will be overriden", metavar="DATFILE")
    parser.add_option("-c", "--col", dest="c",default="1",
    help="Column number to be averaged (default 1, indexed from 1)", metavar="COLNUM")
    parser.add_option("-n", "--window-size", dest="n",default="10",
    help="Size of the averaging window", metavar="AVGWINDOW")
    
    parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose", default=True,
    help="don't print status messages to stdout")

    (options, _) = parser.parse_args()

    return options
    

if __name__ == "__main__":
    sys.exit(main())