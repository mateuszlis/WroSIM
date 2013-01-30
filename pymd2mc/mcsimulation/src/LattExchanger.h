#pragma once
// std
#include <iostream>
#include <string>
#include <iomanip>
#include <cstdlib>
#include <vector>


// stl
#include <map>
#include <list>

// project local
#include "TriangularLattice.h"
#include "types.h"

class TriangularLattice;

using namespace std;

/**
 * @brief Struct used to represent rows and cols distance between lattice sites
 *
 **/
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
    long long row;
    long long col;
}; // struct vectorDist

/**
 * @brief Outputs vector class contents. Used for debug purposes
 **/
ostream &operator<<( ostream &stream, vectorDist & dist );

/**
 * @brief class provides functionality of moving lipids and other objects on the lattice
 *
 * Moved out from TriangularLattice to make it more flexible and loose coupled.
 **/
class LattExchanger
{
    public: // typedefs

    public: // fields
        lattIndex mProteinSites[ 7 ];

    public: // functions
        /**
         * TODO: document
         **/
        LattExchanger( TriangularLattice* latt, lattIndex *proteins = NULL, lattIndex proteinsCnt = 0   );

        /**
         * @brief: TODO document
         **/
        virtual void setProteins( lattIndex* proteins, lattIndex proteinCnt ) 
        {
            mProteins = proteins;
            mProteinCnt = proteinCnt;
        }

        /**
         * TODO: document
         **/
        virtual void exchangeSites( lattIndex pos1, lattIndex pos2 );
        /**
         * @brief: TODO document
         **/
        virtual void moveProtein( lattIndex site, lattIndex destination );

        /**
         * TODO: document
         **/
        virtual bool hasMsd() const { return false; }

        /**
         * TODO: document
         **/
        virtual void calcStat( double & /* param1 */, double & /*param2*/) { ; }

        /**
         * @brief Calculates distance between two lattice sites in the means of rows and columns
         *
         * @return number of rows and cols between sites pos1 and pos2
         *
         **/
        vectorDist calcDist( lattIndex pos1, lattIndex pos2 );


        virtual ~LattExchanger() {};

    protected: // functions 

        /**
         * @brief Moves protein right
         **/
        virtual void moveProteinRight( lattIndex site );

        /**
         * @brief Moves protein left
         **/
        virtual void moveProteinLeft( lattIndex site );

        /**
         * @brief Moves protein to the Top 
         **/
        virtual void moveProteinUp( lattIndex site );
        /**
         * @brief Moves protein down 
         **/
        virtual void moveProteinDown( lattIndex site );

        virtual bool isSpaceToMove( lattIndex site, const std::vector< int >& sitesMoved ) const;

        /**
         * @brief Functions created to facilitate moving of larger objects
         **/
        virtual void initPushAndPop( int site );
        virtual void pushAndPop( lattIndex site, int &value );
        virtual void push( int site, int & value );
        virtual void updateProteinArray( lattIndex site, lattIndex newPos );
        
    protected: // fields
        TriangularLattice* mpLatt;
        const int RIGHT_NEIGH_OFFSET;
        const int RIGHT_TOP_NEIGH_OFFSET;
        const int RIGHT_BOTTOM_NEIGH_OFFSET;
        const int LEFT_BOTTOM_NEIGH_OFFSET;
        const int LEFT_TOP_NEIGH_OFFSET;
        const int LEFT_NEIGH_OFFSET;

        lattIndex* mProteins;
        lattIndex mProteinCnt;

        std::vector< int > mSitesMovedMoveRight;
        std::vector< int > mSitesMovedMoveLeft;
        std::vector< int > mSitesMovedMoveUp;
        std::vector< int > mSitesMovedMoveDown;
}; // class LattExchanger

