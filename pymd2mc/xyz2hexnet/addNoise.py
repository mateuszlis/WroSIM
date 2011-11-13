#!/usr/bin/python
# pymd2mc.xyzfile
"""

"""

__author__ =  'Mateusz Lis'
__version__=  '0.1'


from optparse import OptionParser
from random import random
import sys
from time import time

from structures.xyzfile import XYZFile, XYZFrame, XYZAtom
from utils import delLine, clearFile
 

  
def main():
    
    options = parseCommandLine()
    inFile = XYZFile(options.inXyzFilename)
    
    clearFile(options.outDatFilename)
    outXYZ = XYZFile(options.outDatFilename)
    i = 0
    startTime = time()
    dist = 1
    num = int(options.number)
    while True:
	frame = inFile.nextFrame()
	if i == 0:
		x, y = [], []
		for atom in frame.atoms:
			x.append(atom.x)
			y.append(atom.y)
		fieldX = ( max(x) - min(x) ) / 1.73
		fieldY = max(y) - min(y)
	newFrame = XYZFrame()
	if frame is None:
		break
	for atom in frame.atoms:
		newFrame.atoms.append(atom)
		for k in range(num):
			x0 = atom.x0
			for coord in range(len(x0)):
				x0[coord] += (random() - 0.4) * dist * 2
			newFrame.atoms.append(XYZAtom(atom.symbol, *x0))
	newFrame.comment = "%s, %s, 0.00" % (fieldX, fieldY)	
	outXYZ.addFrame(newFrame)	
        i += 1
        if options.verbose:
            delLine()
            print i, 

     
    if options.verbose: 
        print "Done. Execution time=%f" % (time() - startTime)    
def parseCommandLine():
    """
    Sets up command line arguments and parses them
    """
    parser = OptionParser(usage="%prog ", version="%prog " + __version__,
                          description='''
    This program adds noise atoms to xyz trajectory.''') 
    parser.add_option("-f", "--traj", dest="inXyzFilename",default = "hexTraj.xyz",
    help="xyz input trajectory file (default traj.xyz)", metavar="INXYZFILE")
    parser.add_option("-n", "--number", dest="number",default = "10",
    help="Number of atoms tha will be added for every atom", metavar="NUMBER")
    parser.add_option("-o", "--output", dest="outDatFilename", default="noiseTraj.xyz",
    help="WARNING: it will be overriden", metavar="OUTXYZFILE")
    
    parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose", default=True,
    help="don't print status messages to stdout")

    (options, _) = parser.parse_args()

    return options 
 


if __name__ == '__main__':
    sys.exit(main())
