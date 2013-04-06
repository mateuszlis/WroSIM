#!/usr/bin/env python

"""
grocut.py, by L. Cwiklik, 2009, cwiklik[at]gmail.com
upgraded by M. Lis, 2009 mateusz.lis[at]gmail.com

Oxidizedes DOPC phospholipids producing D1O1+DOO1
usage:
grocut.py GRO_FILE OX_FREQUENCY
Input parameters:
GRO_FILE - input gro file
OX_FREQUENCU - number between 0 and 1 with the ratio of DOPC
residues to be oxidized (random oxidation will be performed
so the number of product can vary)
Prints the resulting gro file.
"""

import sys, copy, random
class Atom:
    def __init__(self, symbol, number, x, y, z, vx, vy, vz):
        self.symbol = symbol
        self.number = number
        self.x = x
        self.y = y
        self.z = z
        self.vx = vx
        self.vy = vy
        self.vz = vz
    def __cmp__(self, other):
        return self.z - other.z

class Molecule:
    def __init__(self):
        self.resname = ""
        self.atoms = []

    def move(self, dx, dy, dz):
        for atom in self.atoms:
            atom.x += dx
            atom.y += dy
            atom.z += dz
    def __cmp__(self, other):
        return self.atoms[0] < other.atoms[0]

class System:
    def __init__(self, filename):
        self.title = ""
        self.n_atoms = 0
        self.molecules = None
        self.box_x = 0.0
        self.box_y = 0.0
        self.box_z = 0.0
        self.DOPCsNum = 0.0
        self.__readgro__(filename)

    def __readgro__(self, filename):

        f = open(sys.argv[1])
        lines = f.readlines()
        title = lines[0].strip()
        n_atoms = int(lines[1])
        box_x, box_y, box_z = lines[-1].split()
        box_x = float(box_x)
        box_y = float(box_y)
        box_z = float(box_z)

        molecules = []

        n_residue_previous = -1
        atom_n = -1
        for i in range(2, len(lines)-1):

            n_residue = int(lines[i][:5])

            #check if a new Molecule object should be created
            if n_residue != n_residue_previous:
                molecule = Molecule()
                resname = lines[i][5:10].strip()
                molecule.resname = resname
                if molecule.resname == 'DOPC': 
                    self.DOPCsNum += 1
            spl = lines[i].split()
            atom_symbol = lines[i][10:15].strip()
            atom_n = int(lines[i][15:20].strip())
            atom_x = float(lines[i][20:28].strip())
            atom_y = float(lines[i][28:36].strip())
            atom_z = float(lines[i][36:44].strip())
            if len(spl) > 6:
                atom_vx = float(lines[i][44:52].strip())
                atom_vy = float(lines[i][52:60].strip())
                atom_vz = float(lines[i][60:68].strip())
            else:
                atom_vx = atom_vy = atom_vz = 0.0

            atom = Atom(atom_symbol, atom_n, atom_x, atom_y, atom_z,
                        atom_vx, atom_vy, atom_vz)

            molecule.atoms.append(atom)

            #checking again
            if n_residue != n_residue_previous:
		#print molecule.resname
                molecules.append(molecule)
                n_residue_previous = n_residue

        f.close()

        self.title = title
        self.n_atoms = n_atoms
        self.molecules = molecules
        self.box_x = box_x
        self.box_y = box_y
        self.box_z = box_z
	
    def put_first(self, resname):
        self.molecules.sort(cmp=lambda x,y: int(x.resname!=resname and y.resname==resname)-1)

    def order(self, order):
        new_molecules = []
        for i in range(len(order)):
            new_molecules.append(0)
        for i in range(len(self.molecules)):
            mol = self.molecules[i]
            resid = i+1
            if resid in order:
                new_molecules[order.index(resid)] = mol
            else:
                new_molecules.append(mol)
        self.molecules = new_molecules

    def popc2dopc(self, popcMol):
        newmol = Molecule()
        newmol.resname = "DOPC"
        newmol.atoms = copy.deepcopy(popcMol.atoms[:8])
        newmol.atoms[6].symbol = "D3A"
        newmol.atoms.append(Atom("C5A",0,newmol.atoms[7].x, newmol.atoms[7].y, 2* newmol.atoms[7].z - newmol.atoms[6].z, 0, 0, 0))
        newmol.atoms += copy.deepcopy(popcMol.atoms[8:])
        return newmol
    def dlpc2dspc(self, dlpcMol):
        newmol = Molecule()
        newmol.resname = "DSPC"
        newmol.atoms = copy.deepcopy(dlpcMol.atoms[:7])
        last = newmol.atoms[-1]
        dz = last.z - newmol.atoms[-2].z

    def dppc2dmpc(self, dppcMol):
        newmol = Molecule()
        newmol.resname = "DMPC"
        newmol.atoms = copy.deepcopy(dppcMol.atoms[:29])
        for atom in dppcMol.atoms[33:]:
            num = int(atom.symbol[1:])
            newAtom = copy.deepcopy(atom)
            newAtom.symbol = atom.symbol[0] + str(num - 4)
            newmol.atoms.append( newAtom )
        newmol.atoms[32].symbol = "O33"
        return newmol

    def card2pope(self, cardMol):
        newmol = Molecule()
        newmol.resname = "POPE"
        newmol.atoms = copy.deepcopy(cardMol.atoms[:13])
        newmol.atoms[0].symbol = "NH3"
        newmol.atoms[1].symbol = "PO4"
        newmol.atoms[2].symbol = "GL1"
        newmol.atoms[3].symbol = "GL2"
        return newmol

	while abs(dz) < 0.2:
		dz *= 2 
        newmol.atoms.append(Atom("C4A",0, last.x, last.y, last.z + dz, 0, 0, 0))
        newmol.atoms.append(Atom("C5A",0, last.x, last.y, last.z + 2 * dz, 0, 0, 0))
        newmol.atoms += copy.deepcopy(dlpcMol.atoms[7:10])
        last = newmol.atoms[-1]
        dz = last.z - newmol.atoms[-2].z
	while abs(dz) < 0.2:
		dz *= 2 
	newmol.atoms.append(Atom("C4B",0, last.x, last.y, last.z + dz, 0, 0, 0))
        newmol.atoms.append(Atom("C5B",0, last.x, last.y, last.z + 2 * dz, 0, 0, 0))
        newmol.atoms += copy.deepcopy(dlpcMol.atoms[10:])
        return newmol
    def dspc2dppc(self, dspcMol):
        newmol = Molecule()
        newmol.resname = "DPPC"
	newmol.atoms = copy.deepcopy(dspcMol.atoms[:8])
	newmol.atoms += copy.deepcopy(dspcMol.atoms[9:13])
	return newmol
  

    def is_top(self, mol):
        atom = mol.atoms[0]
        return atom.z < 4

    def transform(self):
        molLst = []
        top_counter = 0
        bottom_counter = 0
        for m in self.molecules:
            if not m.resname in ["CARD"]:
                molLst.append(m)
                continue
            if m.resname == "CARD":
                if self.is_top(m):
                    top_counter += 1
                    counter = top_counter
                else:
                    bottom_counter += 1
                    counter = bottom_counter
                if counter < 94:
                    molLst.append(self.card2pope(m))
                else:
                    molLst.append(m)
        self.molecules = molLst

    def sort(self):
        d = {}
        for mol in self.molecules:
            if mol.resname in d:
                d[mol.resname].append(mol)
            else:
                d[mol.resname] = [mol]
        mols = []
        mols += sorted(d["CARD"])
        mols += sorted(d["POPE"])
        mols += d["W"]
        del d["CARD"]
        del d["POPE"]
        del d["W"]
        for mol_lst in d.values():
            mols += mol_lst
        self.molecules = mols






    def populate(self):
        molLst = []
        dx = self.box_x + 0.01 #proved to be perfect magic number
        dy = self.box_y + 0.01
        for m in self.molecules:
            molLst.append(m)
            m2 = copy.deepcopy(m)
            m2.move(dx, 0, 0)
            molLst.append(m2)
            m3 = copy.deepcopy(m)
            m3.move(0, dy, 0)
            molLst.append(m3)
            m4 = copy.deepcopy(m)
            m4.move(dx, dy, 0)
            molLst.append(m4)
        self.molecules = molLst
        self.box_x *= 2
        self.box_y *= 2


        
    def printgro(self):
        print self.title
        n_atoms = 0
        for m in self.molecules:
            n_atoms += len(m.atoms)
        print n_atoms
        n_atom = 0
        n_molecule = 0
        out_line = ""
        for m in self.molecules:
            n_molecule += 1
            
            for a in m.atoms:
                n_atom += 1
                out_line = str(n_molecule).rjust(5)
                out_line += m.resname.ljust(5)
                out_line += a.symbol.rjust(5)
                out_line += str(n_atom).rjust(5)
                out_line += str('%.3f'%a.x).rjust(8)
                out_line += str('%.3f'%a.y).rjust(8)
                out_line += str('%.3f'%a.z).rjust(8)
                out_line += str('%.4f'%a.vx).rjust(8)
                out_line += str('%.4f'%a.vy).rjust(8)
                out_line += str('%.4f'%a.vz).rjust(8)

                print out_line

        print self.box_x, self.box_y, self.box_z

    def trim(self, x_min, x_max, y_min, y_max, z_min, z_max):
    
        mol_to_remove = []
        for i in range(len(self.molecules)):
            m = self.molecules[i]
            remove = False
            for atom in m.atoms:
                if atom.x < x_min or atom.x > x_max or \
                       atom.y < y_min or atom.y > y_max or \
                       atom.z < z_min or atom.z > z_max:
                    remove = True
            if remove:
                mol_to_remove.append(i)
                
        new_molecules = []
        for i in range(len(self.molecules)):
            if i not in mol_to_remove:
                new_molecules.append(self.molecules[i])
        self.molecules = new_molecules

        # find new box sizes
        x_left = x_max # left<->max sic!
        x_right = x_min
        y_left = y_max
        y_right = y_min
        z_left = z_max
        z_right = z_min
        for m in self.molecules:
            for a in m.atoms:
                if a.x < x_left:
                    x_left = a.x
                if a.x > x_right:
                    x_right = a.x
                if a.y < y_left:
                    y_left = a.y
                if a.y > y_right:
                    y_right = a.y
                if a.z < z_left:
                    z_left = a.z
                if a.z > z_right:
                    z_right = a.z
            
        self.box_x = x_right - x_left
        self.box_y = y_right - y_left
        self.box_z = z_right - z_left
        

def main():

    if len(sys.argv) < 3:
        print __doc__
        return -1

    filename = sys.argv[1]
    frequency = float(sys.argv[2])

    system = System(filename)
    #system.trim(x_min, x_max, y_min, y_max, z_min, z_max)
    #system.put_first(resname)
    #system.order(order)
    #system.oxidize1(frequency)
    #system.put_first("D1O1")
    #system.put_first("DOO1")

    #system.printgro()
    system.transform()
	
    #order = []
    #for i in range(256):
    #    if (i+1)%2:
    #        order.append(i+1)
    #for i in range(256):
    #    if not (i+1)%2:
    #        order.append(i+1)
    #system.order(order)
	#system.transform()
    #system.populate()
    system.sort()
    system.printgro()

if __name__ == '__main__':
    main()
