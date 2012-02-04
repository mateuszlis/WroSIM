/*
 * TriangularLattice.cpp
 *
 *  Created on: 23-05-2011
 *      Author: lisu
 */

#include "TriangularLattice.h"

void TriangularLattice::clearArr()
{
    for ( int i = 0; i < mLatticeSize; i++ )
    {
        mpLattice[i] = 0;
    }
}
void TriangularLattice::distributeParticles( int firstTypeParticlesCnt )
{
    for ( int i = 0; i < firstTypeParticlesCnt; i++ )
    {
        int pos = 0;
        while ( this->mpLattice[pos = rand() % mLatticeSize] )
        {
        };
        this->mpLattice[pos] = 1;
    }
}
// TODO Copy constructor and assignment operator

TriangularLattice::TriangularLattice( string /*filename*/ )
    : mpLattice( NULL )
      , mLatticeSize( 0 )
      , mRowSize( 0 )
{
    for ( int i = 0; i < mNeighbCnt; ++i )
    {
        mNeighb[i] = 0;
    }
    //mNeighb = { 0, 0, 0, 0, 0, 0 };
}
TriangularLattice::TriangularLattice( int latticeSize, int rowSize, int firstTypeParticlesCnt )
{
    mpLattice = new int[latticeSize];

    mLatticeSize = latticeSize;
    mRowSize = rowSize;
    clearArr();
    int neighb[] = { 1, -1, rowSize, -rowSize, rowSize - 1, -rowSize + 1 }; // bloody helical boundary conditions
    int neighbL[] = { 1, rowSize - 1, rowSize, -rowSize, 2 * rowSize - 1, -rowSize + 1 }; // bloody periodic boundary conditions
    int neighbR[] = { -rowSize + 1, -1, rowSize, -rowSize, rowSize - 1, - 2 * rowSize + 1 }; // bloody helical boundary conditions
    for ( int i = 0; i < mNeighbCnt; i++ )
    {
        mNeighb[i] = neighb[i];
        mNeighbLeft[i] = neighbL[i];
        mNeighbRight[i] = neighbR[i];
    }

    this->distributeParticles( firstTypeParticlesCnt);

}

int TriangularLattice::operator[]( int index ) const
{
    return mpLattice[index];
}

int TriangularLattice::getLatticeSize() const
{
    return mLatticeSize;
}

int TriangularLattice::getRowSize() const
{
    return mRowSize;
}
void TriangularLattice::exchangeSites( int pos1, int pos2 )
{
    // fast variable values exchange trick (pos1 != pos2)
	mpLattice[pos1] ^= mpLattice[pos2];
    mpLattice[pos2] ^= mpLattice[pos1];
    mpLattice[pos1] ^= mpLattice[pos2];

}
int TriangularLattice::simNeighbCount( int pos )
{
    int sum = 0;
    for ( int i = 0; i < 6; i++ )
    {
        int currentNeigh = getNeighbIndex(pos, i);
        sum += ( mpLattice[pos] == mpLattice[currentNeigh] ? 1 : 0 );
    }
    return sum;
}

int TriangularLattice::getNeighbIndex( int pos, int neighborNum ) const
{
    int translationIndex( 0 );
    if ( pos % mRowSize == 0 )
    {
       translationIndex = mNeighbLeft[neighborNum];
    }
    else if ( pos % mRowSize == ( mRowSize - 1 ) )
    {
        translationIndex = mNeighbRight[neighborNum];
    }
    else
    {
        translationIndex = mNeighb[neighborNum];
    }
	int neigh = ( pos + translationIndex );
	if ( neigh >= mLatticeSize )
		neigh -= mLatticeSize;
	if ( neigh < 0 )
	    neigh += mLatticeSize;
    return neigh;
}
int TriangularLattice::getNeighborsCnt() const
{
    return mNeighbCnt;
}

TriangularLattice::~TriangularLattice()
{

    delete[] this->mpLattice;
}

ostream &operator<<( ostream &stream, TriangularLattice &latt )
{
    stream << latt.getLatticeSize() << endl;
    stream << "Simulation" << endl;
    for ( int i = 0; i < latt.getLatticeSize(); i++ )
        if ( latt[i] )
        {
	    int line( i / latt.getRowSize() );
            double y = line * 0.866025;
            double x = ( i % latt.getRowSize() ) - line * 0.5;
            stream << "A\t" << setprecision( 8) << x << "\t" << y << "\t0.00000000" << endl;
        }
    for ( int i = 0; i < latt.getLatticeSize(); i++ )
    {
        if ( ! ( latt.mpLattice[i] ) )
        {
	    int line( i / latt.getRowSize() );
            double y = line * 0.866025;
            double x = ( i % latt.getRowSize() ) - ( line * 0.5 ); 
            stream << "B\t" << setprecision( 8) << x << "\t" << y << "\t0.00000000" << endl;
        }
    }
    return stream;
}

