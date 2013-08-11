#!/usr/bin/env python
"""
This script merges two gro files into one
usage:
    script.py first.gro second.gro > merged.gro
"""

from tools.changeLipids import Atom, Molecule, System
import sys
import random
from copy import deepcopy
import time
from itertools import product
from numpy import arange as range

class SystemMerger(System):
    def merge(self, another_system):
        self.box_x = max(self.box_x, another_system.box_x) 
        self.box_y = max(self.box_y, another_system.box_y)  
        self.box_z = max(self.box_z, another_system.box_z) + 1.3
        

        mols = deepcopy(another_system.molecules)
        diff_z = (self.box_z - 0.4) - self.molecules[0].atoms[0].z
        iterations = int(sys.argv[3])
        if sys.argv[4] == "1":
            source = deepcopy(self.molecules)
            self.molecules = []
        else:
            source = deepcopy(mols)
        base_y = self.box_y * 0.7
        base_x = self.box_z * 0.7
        base_list_y = range(0, base_y, base_y / iterations**(1./2))
        base_list_x = range(0, base_x, base_x / iterations**(1./2))
        for diff_y, diff_x in list(product(base_list_y, base_list_x))[:iterations]:
            #diff_y = random.random() * self.box_y * 0.7
            #diff_z = random.random() * self.box_z  * 0.7

            mult = deepcopy(source)
            for mol in mult:
                for atom in mol.atoms:
                    atom.x += diff_x
                    atom.y += diff_y
                    atom.z += diff_z
            self.molecules += mult

        self.molecules += deepcopy(mols)

def main():
    
    random.seed(time.time())
    if len(sys.argv) < 4:
        print __doc__
        return -1
    system1 = SystemMerger(sys.argv[1])
    system2 = System(sys.argv[2])
    system1.merge(system2)
    system1.printgro()

if __name__=="__main__":
    main()
