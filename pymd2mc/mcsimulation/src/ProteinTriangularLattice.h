
#pragma once

#include <vector> 

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

        int calcOtherLipidNeighbors( lattIndex pos )
        {
            lattIndex sum( 0 );
            for ( int i( 0 ) ; i < mNeighbCnt ; ++i )
            {
                lattIndex currentNeigh( getNeighbIndex( pos, i ) );
                sum += ( mpLattice[ pos ] != mpLattice[ currentNeigh ] 
                        && mpLattice[ currentNeigh ] != PROTEIN_A
                        && mpLattice[ currentNeigh ] != PROTEIN_B ?
                        1 : 0 );
            }
            return sum;
        }

        int calcNeighbors( lattIndex pos, lattMember neighborType )
        {
            lattIndex sum( 0 );
            for ( int i( 0 ) ; i < mNeighbCnt ; ++i )
            {
                lattIndex currentNeigh( getNeighbIndex( pos, i ) );
                sum += ( mpLattice[ currentNeigh ] == neighborType ) ? 1 : 0;
            }
            return sum;
        }
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

        const int RIGHT_NEIGH_OFFSET;
        const int RIGHT_TOP_NEIGH_OFFSET;
        const int RIGHT_BOTTOM_NEIGH_OFFSET;
        const int LEFT_BOTTOM_NEIGH_OFFSET;
        const int LEFT_TOP_NEIGH_OFFSET;
        const int LEFT_NEIGH_OFFSET;

        std::vector< int > mSitesMovedMoveRight;
        std::vector< int > mSitesMovedMoveLeft;
        std::vector< int > mSitesMovedMoveUp;
        std::vector< int > mSitesMovedMoveDown;

    protected: // functions
        virtual bool isFree( lattIndex pos );
        void putProtein( lattIndex pos, LATTICE_FIELD_NAMES protein = PROTEIN_A );
        void pushAndPop( lattIndex site, lattMember &value );
        void updateProteinArray( lattIndex site, lattIndex newPos );
        bool isSpaceToMove( lattIndex site, lattMember movedSites[], lattIndex movedSitesSize );
        bool isSpaceToMove( lattIndex site, const std::vector< int >& sitesMoved );
                        
}; // class ProteinTriangularLattice
