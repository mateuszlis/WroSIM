#!/usr/bin/python
# pymd2mc.xyzfile
"""

"""

__author__ =  'Mateusz Lis'
__version__=  '0.1'


from optparse import OptionParser
import sys
from time import time

from constants import R, T
from energyCalc import EnergyCalculator
from latticeProjector import LatticeProjectorSimple
from lattices import HexLattice
from structures.xyzfile import XYZFile
from utils import delLine, clearFile
 

  
def main():
    
    options = parseCommandLine()
    inFile = XYZFile(options.inXyzFilename)
    
    clearFile(options.outDatFilename)
    outFile = open(options.outDatFilename, 'w')
    i = 0
    startTime = time()
    omegas = []

    sumOmegas = 0L
    calc = EnergyCalculator(inFile, R, T)
    
    while True:
        i += 1
        if options.verbose:
            delLine()
            print i, 
        omega = calc.getNextEnergy(options.symbol)
        if omega is None:
            break
        omega , sim, diff = omega
        if omega > -10**4 and omega < 10**10: 
            omegas.append(omega)
            sumOmegas += omega
            outFile.write("%d  %f %f %f \n" % (i, omega, sim, diff))

        
    outFile.close()
    if options.verbose: 
        print "Done. Execution time=%f" % (time() - startTime)    
    print "omegas" ,sumOmegas, (sum(omegas))
    lenOmegas = len(omegas)
    midOmega = (sum(omegas)/len(omegas))
    print "Result omegaAB = %f" % midOmega
    sd = 0
    for omega in omegas:
        sd += (midOmega - omega)**2
    sd /= len(omegas)
    sd **= (1./2.)
    print "Standard deviation = %f" % sd
def parseCommandLine():
    """
    Sets up command line arguments and parses them
    """
    parser = OptionParser(usage="%prog ", version="%prog " + __version__,
                          description='''
    This program calculates omegaAB value from a hexagonal lattice trajectory
    stored in xyz file (see for more details)''')
    parser.add_option("-f", "--traj", dest="inXyzFilename",default = "hexTraj.xyz",
    help="xyz input trajectory file (default traj.xyz)", metavar="INXYZFILE")
    parser.add_option("-r", "--reference", dest="symbol",default = "P11",
    help="reference particle name", metavar="ADATOM")
    parser.add_option("-o", "--output", dest="outDatFilename", default="omega.dat",
    help="output dat file with omega values for each frame. WARNING: it will be overriden", metavar="OUTXYZFILE")
    
    parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose", default=True,
    help="don't print status messages to stdout")

    (options, _) = parser.parse_args()

    return options 
 


if __name__ == '__main__':
    sys.exit(main())
