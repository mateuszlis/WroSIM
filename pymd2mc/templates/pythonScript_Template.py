#!/usr/bin/python
"""
Program Docstring
"""

__author__ = 'Mateusz Lis'
__version__ = '0.1'

from optparse import OptionParser
import sys
from time import time


def main(argv=None):
    global options 
    options = parseCommandLine()

    startTime = time()
    if options.verbose:
        print "Execution time", time() - startTime
    
def parseCommandLine():
    """
    Sets up command line arguments and parses them
    """
    parser = OptionParser(usage="%prog ", version="%prog " + __version__,
                          description='''
    Program Docstring ''')
    parser.add_option("-q", "--quiet",
    action="store_false", dest="verbose", default=True,
    help="don't print status messages to stdout")

    (options, _) = parser.parse_args()

    return options
    

if __name__ == "__main__":
    sys.exit(main())
