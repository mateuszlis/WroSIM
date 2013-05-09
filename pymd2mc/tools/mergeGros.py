#!/usr/bin/env python
"""
This script merges two gro files into one
usage:
    script.py first.gro second.gro > merged.gro
"""

from tools.changeLipids import Atom, Molecule, System
import sys

class SystemMerger(System):
    def merge(self, another_system):
        self.molecules += another_system.molecules
        self.box_x = max(self.box_x, another_system.box_x)
        self.box_y = max(self.box_y, another_system.box_y)
        self.box_z = max(self.box_z, another_system.box_z)

def main():
    
    if len(sys.argv) < 3:
        print __doc__
        return -1
    system1 = SystemMerger(sys.argv[1])
    system2 = System(sys.argv[2])
    system1.merge(system2)
    system1.printgro()

if __name__=="__main__":
    main()
