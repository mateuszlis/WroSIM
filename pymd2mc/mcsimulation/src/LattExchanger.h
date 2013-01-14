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
struct vectorDist
{
    vectorDist() : row( 0 ), col( 0 ) {};
    vectorDist( lattIndex argRow, lattIndex argCol ) : row( argRow ), col( argCol ) {};
    vectorDist operator+( vectorDist added ) { return vectorDist( row + added.row, col + added.col ); }
    vectorDist operator-( vectorDist added ) { return vectorDist( row - added.row, col - added.col ); }
    vectorDist operator+=( vectorDist added ) 
    { 
        col += added.col;
        row += added.row;
        return *this;
    }
    lattIndex squareDisp() { return row*row + col*col; }
    lattIndex row;
    lattIndex col;
}; // struct vectorDist
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
        virtual void exchangeSites( lattIndex pos1, lattIndex pos2 );

        /**
         * TODO: document
         **/
        virtual bool hasMsd() { return false; }

        /**
         * TODO: document
         **/
        virtual double calcStat() { return 0; };

        /**
         * @brief Calculates distance between two lattice sites in the means of rows and columns
         *
         * @return number of rows and cols between sites pos1 and pos2
         *
         **/
        vectorDist calcDist( lattIndex pos1, lattIndex pos2 );


        virtual ~LattExchanger() {};

    protected: // fields
        TriangularLattice* mpLatt;

};

