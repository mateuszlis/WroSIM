#pragma once
// std
#include <iostream>
#include <string>
#include <iomanip>
#include <cstdlib>


// stl
#include <map>
#include <list>

// project local
#include "TriangularLattice.h"
#include "types.h"

class TriangularLattice;

using namespace std;
/**
 * @brief class provides functionality of moving lipids and other objects on the lattice
 *
 * Moved out from TriangularLattice to make it more flexible and loose coupled.
 **/
class LattExchanger
{
    public: // typedefs

    public: // functions
        /**
         * TODO: document
         **/
        LattExchanger( TriangularLattice* latt ) 
            : mpLatt( latt ) {}

        /**
         * TODO: document
         **/
        virtual void exchangeSites( lattIndex pos1, lattIndex pos2 ) const;

    protected: // fields
        TriangularLattice* mpLatt;

};

