/*
 * MPKK.h
 *
 *  Created on: 01-09-2012
 *      Author: lisu
 */

#pragma once

#include <fstream>
#include <string>


// std
#include <algorithm>
#include <assert.h>
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <vector>
using namespace std;

// project related
#include "../Metropolis.h"
#include "../TriangularLattice.h"
#include "../types.h"

class MPKK: public Metropolis
{
    public: // functions

        /**
         * @brief Constructor, see Metropolis::Metropolis
         **/
        MPKK( TriangularLattice* latt, double omegaAB, int T, int equilibSteps
                , unsigned int imageRes = 300 );

        /**
         * @brief Perform a simulation for a number of iterations
         *
         * @param steps - number of iterations to be performed
         **/
        virtual void run( int steps );

        /**
         * @brief default dtor
         **/
        virtual ~MPKK();

    protected: // fields
        static const int PRECALCULATED_ENERGIES_COUNT = 36; ///< size of array which contains precalculated energy values
        float energies[ PRECALCULATED_ENERGIES_COUNT ]; 
        static const int mStepSize = 7; ///< count of neighbors on triangular latt
        const unsigned int IMAGE_RESOLUTION;

    protected: // helper member functions
        void checkLattice();
        unsigned int calcBlockSize();
        void precalculateEnergies();
        void generateImage( TriangularLattice::lattMember *d_latt, int stepNum );
        void initDevice( TriangularLattice::lattMember * &d_latt
                       , TriangularLattice::lattMember * &h_latt
                       , unsigned int &memSizeLatt
                       , float * &d_energies
                       , unsigned int &blockSize
                       , unsigned int &blocksCount );

    public: // debug member functions - not used in production code

        /**
         * @brief helper debug function - calculates sum of lattice
         *
         **/
        long calcSum( TriangularLattice::lattMember *data, int size );


        /**
         * @brief helper debug function - prints lattice to stdout
         * 
         **/
        void printLatt( TriangularLattice::lattMember *data, int rowSize, int rowsCount );
};
