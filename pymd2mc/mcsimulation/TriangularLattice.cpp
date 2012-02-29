/*
 * TriangularLattice.cpp
 *
 *  Created on: 23-05-2011
 *      Author: lisu
 */

#include "TriangularLattice.h"

typedef TriangularLattice::lattIndex lattIndex;

void TriangularLattice::clearArr()
{
    for ( lattIndex i = 0; i < mLatticeSize; i++ )
    {
        mpLattice[i] = 0;
    }
}
void TriangularLattice::distributeParticles( lattIndex firstTypeParticlesCnt )
{
    for ( lattIndex i = 0; i < firstTypeParticlesCnt; i++ )
    {
        lattIndex pos = 0;
        while ( this->mpLattice[pos = rand() % mLatticeSize] )
        {
        };
        this->mpLattice[pos] = 1;
    }
}

void TriangularLattice::pushNeighborsToQueue( std::list< lattIndex > & queue, lattIndex siteInd )
{
    for( lattIndex i = 0 ; i < mNeighbCnt ; ++i )
    {
        queue.push_back( getNeighbIndex( siteInd, i ) );
    }
}
// TODO Copy constructor and assignment operator

TriangularLattice::TriangularLattice( string /*filename*/ )
    : mpLattice( NULL )
      , mLatticeSize( 0 )
      , mRowSize( 0 )
{
    for ( lattIndex i = 0; i < mNeighbCnt; ++i )
    {
        mNeighb[i] = 0;
    }
    //mNeighb = { 0, 0, 0, 0, 0, 0 };
}
TriangularLattice::TriangularLattice( lattIndex latticeSize, lattIndex rowSize, lattIndex firstTypeParticlesCnt )
{
    mpLattice = new lattMember[latticeSize];

    mLatticeSize = latticeSize;
    mRowSize = rowSize;
    clearArr();
    lattIndex neighb[] = { 1, -1, rowSize, -rowSize, rowSize - 1, -rowSize + 1 }; // bloody helical boundary conditions
    for ( lattIndex i = 0; i < mNeighbCnt; i++ )
    {
        mNeighb[i] = neighb[i];
    }

    this->distributeParticles( firstTypeParticlesCnt);

}

TriangularLattice::lattMember TriangularLattice::operator[]( lattIndex index ) const
{
    return mpLattice[index];
}

lattIndex TriangularLattice::getLatticeSize() const
{
    return mLatticeSize;
}

lattIndex TriangularLattice::getRowSize() const
{
    return mRowSize;
}
void TriangularLattice::exchangeSites( lattIndex pos1, lattIndex pos2 )
{
    // fast variable values exchange trick (pos1 != pos2)
	mpLattice[pos1] ^= mpLattice[pos2];
    mpLattice[pos2] ^= mpLattice[pos1];
    mpLattice[pos1] ^= mpLattice[pos2];

}
lattIndex TriangularLattice::simNeighbCount( lattIndex pos )
{
    lattIndex sum = 0;
    for ( int i = 0; i < 6; i++ )
    {
        lattIndex currentNeigh = getNeighbIndex(pos, i);
        sum += ( mpLattice[pos] == mpLattice[currentNeigh] ? 1 : 0 );
    }
    return sum;
}

lattIndex TriangularLattice::getNeighbIndex( lattIndex pos, int neighborNum ) const
{
    lattIndex translationIndex( mNeighb[ neighborNum ] );
	lattIndex neigh = ( pos + translationIndex );
	if ( neigh >= mLatticeSize )
    {
		neigh -= mLatticeSize;
    }
	if ( neigh < 0 )
    {
	    neigh += mLatticeSize;
    }
    return  neigh;
}
unsigned int TriangularLattice::getNeighborsCnt() const
{
    return mNeighbCnt;
}

void TriangularLattice::calculateClusters( TriangularLattice::clustersMap& map )
{
    const lattMember kind = 1;

    doneMap doneSites;
    for( lattIndex startPos = 0 ; startPos < getLatticeSize() ; ++startPos )
    {
        lattIndex clusterSize = 1;
        if( ! doneSites[ startPos ] > 0 && mpLattice[ startPos ] == kind )
        {
            lattIndex currentSite = startPos;
            doneSites[ currentSite ] = 1;
            std::list< lattIndex > queue;
            pushNeighborsToQueue( queue, currentSite  );
            while( queue.size() )
            {
                currentSite = queue.back();
                queue.pop_back();
                if( ! doneSites[ currentSite ] && mpLattice[ currentSite ] == kind )
                {
                    clusterSize++;
                    pushNeighborsToQueue( queue, currentSite );
                    doneSites[ currentSite ] = 1;
                }
            }
        map[ clusterSize ]++;
        }
    }
}

        
TriangularLattice::~TriangularLattice()
{

    delete[] this->mpLattice;
}

ostream &operator<<( ostream &stream, TriangularLattice &latt )
{
    stream << latt.getLatticeSize() << endl;
    stream << "Simulation" << endl;
    for ( lattIndex i = 0; i < latt.getLatticeSize(); i++ )
        if ( latt[i] )
        {
	    lattIndex line( i / latt.getRowSize() );
            double y = line * 0.866025;
            double x = ( i % latt.getRowSize() ) - line * 0.5;
            stream << "A\t" << setprecision( 8) << x << "\t" << y << "\t0.00000000" << endl;
        }
    for ( lattIndex i = 0; i < latt.getLatticeSize(); i++ )
    {
        if ( ! ( latt.mpLattice[i] ) )
        {
	    lattIndex line( i / latt.getRowSize() );
            double y = line * 0.866025;
            double x = ( i % latt.getRowSize() ) - ( line * 0.5 ); 
            stream << "B\t" << setprecision( 8) << x << "\t" << y << "\t0.00000000" << endl;
        }
    }
    return stream;
}

