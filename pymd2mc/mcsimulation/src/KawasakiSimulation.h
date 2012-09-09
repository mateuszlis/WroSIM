/*
 * Metropolis.h
 *
 *  Created on: 10-09-2012
 *      Author: lisu
 */

#pragma once

// std
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <cstdlib>

// project related
#include "Metropolis.h"
#include "TriangularLattice.h"

using namespace std;

/**
 * @brief enum for sampling type choice
 */
enum Sampler
{
    Kawasaki, Almeida, MassiveParallelKawasaki
};

/**
 * @brief Class performs Kawasak simulations on triangular lattices with two types of lipids
 **/
class KawasakiSimulation : public Metropolis
{
    public: // member functions
        /**
         * @brief Creates instance of Kawasaki simulator.
         * 
         * @param latt - lattice to perform simulation
         * @param omegaAB - energy of interaction between lipids
         * @param T - temperature
         **/
        KawasakiSimulation( TriangularLattice* latt, double omegaAB = 0.0, int T = 325, int equilibSteps = 0 );

        /**
         * @brief Starts a simulation for defined number of iterations
         *
         * @param steps - number of iterations
         **/
        virtual void run( int steps );

        /**
         * @brief Since it is possible to use various functions to chose lattice sites to exchange - this sets which one should be used
         *
         * @param s - Sampler to be used
         **/
        virtual void setSampler( Sampler s );

        /**
         * @brief Default dtor
         **/
        virtual ~KawasakiSimulation();

    protected: // fields
        // we use raw pointers instead of std::tr1::function because 
        // we are extremely concerned about the performance
        void (*mpSampler)( TriangularLattice *latt, int &, int & ); ///< function that choses lattice sites to exchange
        
    protected: // member functions
        /**
         * @brief Performs single step of metropolis algorithm for triangular lattice
         */
        void metropolisStep(); 
};
