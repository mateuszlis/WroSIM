
#pragma once
// project related
#include "TriangularLattice.h"

class ProteinTriangularLattice : public TriangularLattice
{
    public: // functions
        /**
         * @brief Creates lattice and distributes lipids and proteins in desired number
         **/
        ProteinTriangularLattice( lattIndex latticeSize
                                , lattIndex rowSize
                                , lattIndex firstTypeParticlesCnt
                                , lattIndex proteinCnt = 0
                                , bool distributeRandomly = true );

        /**
         * @brief Distributes Proteins in lattice
         **/
        void distributeProteins();

        /**
         * @brief Distributes Proteins in lattice randomly
         **/
        void distributeProteinsRandomly();

        /**
         * @brief Moves protein from given lattice site one step to given lattice site
         **/
        void moveProtein( lattIndex site, lattIndex destination );

        /**
         * @brief Moves protein right
         **/
        void moveProteinRight( lattIndex site );

        /**
         * @brief Moves protein left
         **/
        void moveProteinLeft( lattIndex site );

        /**
         * @brief Moves protein to the Top 
         **/
        void moveProteinUp( lattIndex site );
        /**
         * @brief Moves protein down 
         **/
        void moveProteinDown( lattIndex site );
        
        /**
         * @brief Get array of protein positions
         *
         * @return pointer to array of protein positions. Size of the array is number of proteins
         **/
        lattIndex* getProteins() { return mProteins; }

        /**
         * @brief Get number of proteins in lattice
         *
         * @return number of proteins.
         **/
        lattIndex getProteinCnt() { return mProteinCnt; }

        /**
         * @brief destructor
         **/
        virtual ~ProteinTriangularLattice();

    public: // fields
        static const unsigned int mProteinSize = 7;
        lattIndex mProteinCnt; ///< count of the proteins in lattice
        lattIndex* mProteins; ///< array of protein center positions

    protected: // fields 
        lattIndex mProteinSites[ mProteinSize ];

    protected: // functions
        virtual bool isFree( lattIndex pos );
        void putProtein( lattIndex pos, LATTICE_FIELD_NAMES protein = PROTEIN_A );
        void pushAndPop( lattIndex site, lattMember &value );
        void updateProteinArray( lattIndex site, lattIndex newPos );
        bool isSpaceToMove( lattMember movedSites[], lattIndex movedSitesSize );
}; // class ProteinTriangularLattice
