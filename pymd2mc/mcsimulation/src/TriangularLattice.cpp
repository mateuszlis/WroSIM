/*
 * TriangularLattice.cpp
 *
 *  Created on: 23-05-2011
 *      Author: lisu
 */

#include "TriangularLattice.h"
#include "LattExchanger.h"
#include <string.h>
#include <list>
#include <stdexcept>

typedef TriangularLattice::lattIndex lattIndex;

void TriangularLattice::clearArr()
{
    for ( lattIndex i = 0; i < mLatticeSize; i++ )
    {
        mpLattice[i] = LIPID_A;
    }
}
void TriangularLattice::distributeParticlesRandomly( lattIndex firstTypeParticlesCnt )
{
    for ( lattIndex i = 0; i < firstTypeParticlesCnt; ++i )
    {
        lattIndex pos = rand() % getLatticeSize(); //FIXME: if RAND_MAX is too small this causes bad distribution
        while ( !isFree( pos ) )
        {
				pos = rand() % getLatticeSize();
        };
        this->mpLattice[pos] = LIPID_B;
    }
}
bool TriangularLattice::isFree( lattIndex pos ) const
{
    return mpLattice[ pos ] == LIPID_A;
}

TriangularLattice::lattMember TriangularLattice::get( int index ) const
{ 
    return mpLattice[ cutToPos( index ) ];
}

void TriangularLattice::distributeParticles( lattIndex firstTypeParticlesCnt )
{
    for ( lattIndex pos = 0; pos < firstTypeParticlesCnt; ++pos  )
    {
        this->mpLattice[pos] = LIPID_B;
    }
}

// TODO Copy constructor and assignment operator

TriangularLattice::TriangularLattice( string /*filename*/ )
    : mpLattice( NULL )
      , mLatticeSize( 0 )
      , mRowSize( 0 )
      , selfLattExchanger( true )
{
    for ( lattIndex i = 0; i < mNeighbCnt; ++i )
    {
        mNeighb[i] = 0;
    }
    //mNeighb = { 0, 0, 0, 0, 0, 0 };
    mpExchanger = new LattExchanger( this );
}
TriangularLattice::TriangularLattice( lattIndex latticeSize
                                    , lattIndex rowSize
                                    , lattIndex firstTypeParticlesCnt
                                    , bool distributeRandomly )
      : selfLattExchanger( true )
{
    if ( firstTypeParticlesCnt > latticeSize )
    {
        throw InputParametersException( "Number of lipids of this type cannot be larger than lattice size" );
    }

    mpLattice = new lattMember[ latticeSize ];


    mLatticeSize = latticeSize;
    this->mRowSize = rowSize;
    clearArr();
    lattIndex neighb[] = { 1, rowSize, rowSize - 1, -1, -rowSize, -rowSize + 1 }; // bloody helical boundary conditions
    for ( lattIndex i = 0; i < mNeighbCnt; i++ )
    {
        mNeighb[i] = neighb[i];
    }
    if ( distributeRandomly )
    {
        this->distributeParticlesRandomly( firstTypeParticlesCnt );
    }
    else 
    {
        this->distributeParticles( firstTypeParticlesCnt );
    }
    
    mpExchanger = new LattExchanger( this );

}

TriangularLattice::TriangularLattice( const TriangularLattice &  )
      : mpLattice( nullptr )
        , mLatticeSize( 0 )
        , mRowSize( 0 )
        , mpExchanger( nullptr )
        , selfLattExchanger( true )
        , mNeighb{0, 0, 0, 0, 0, 0}

{
    std::cout << "Undesired copy instruction" << std::endl;
    throw std::runtime_error( "Undesired copy in Triangular lattice" );
}

TriangularLattice::lattMember &TriangularLattice::operator[]( lattIndex index ) const
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
    mpExchanger->exchangeSites( pos1, pos2 );
}

lattIndex TriangularLattice::simNeighbCount( lattIndex pos )
{
    lattIndex sum = 0;
    for ( int i = 0; i < mNeighbCnt; ++i )
    {
        lattIndex currentNeigh = getNeighbIndex(pos, i);
        sum += ( mpLattice[pos] == mpLattice[currentNeigh] ? 1 : 0 );
    }
    return sum;
}

lattIndex TriangularLattice::getNeighbIndex( lattIndex pos, int neighborNum ) const
{
    lattIndex translationIndex( mNeighb[ neighborNum ] );
	return cutToPos( pos + translationIndex );
}

lattIndex TriangularLattice::cutToPos( int pos ) const
{
	if ( pos >= mLatticeSize )
    {
		pos -= mLatticeSize;
    }
	if ( pos < 0 )
    {
	    pos += mLatticeSize;
    }
    return  pos;
}
int TriangularLattice::getNeighborsCnt() const
{
    return mNeighbCnt;
}

bool TriangularLattice::gotDifferentNeighbors(list<int> neighLabels, int currentLabel)
{
	for(list<int>::iterator iterator = neighLabels.begin(); iterator != neighLabels.end(); ++iterator)
		if(*iterator != currentLabel) return true;

	return false;
}

int TriangularLattice::findAncestor(int currentLabel, TriangularLattice::clustersMap& map)
{
	int label = currentLabel;

	while(map[label] < 0 && label >= 0)
	{
		label += map[label];
	}
	return label;

}
void TriangularLattice::calculateClusters( TriangularLattice::clustersMap& map )
{
    const lattMember kind = LIPID_B;
	
	clustersMap clMap;

	int * labels = new int[getLatticeSize()];
	memset(labels, 0, getLatticeSize() * sizeof(int));

	list<int> neighLabels;

	int currentLabel = 0;

	for(int i = 0; i < getLatticeSize(); i++)
	{
		int index = i % getLatticeSize();
		if(mpLattice[index] == kind)
		{
			neighLabels.clear();
			for(int n = 0; n < 6; ++n)
				if(labels[getNeighbIndex(index,n)] > 0 ) neighLabels.push_back(findAncestor(labels[getNeighbIndex(index,n)], clMap));

			neighLabels.sort();
			neighLabels.unique();

			if(neighLabels.empty())
			{
				labels[index] = ++currentLabel;
				clMap[currentLabel] = 1;
			}
			else if(!gotDifferentNeighbors(neighLabels, currentLabel))
			{
				labels[index] = currentLabel;
				clMap[findAncestor(currentLabel, clMap)]++;
			}
			else
			{
				int ancestorLabel = neighLabels.front();
				neighLabels.pop_front();

				labels[index] = ancestorLabel;
				clMap[ancestorLabel]++;

				for(list<int>::iterator it = neighLabels.begin(); it != neighLabels.end(); ++it)
				{
					clMap[ancestorLabel] += clMap[neighLabels.front()];
					clMap[neighLabels.front()] = ancestorLabel - neighLabels.front();
				}
			}

		}
	}

	for(std::map< lattIndex, lattIndex >::iterator mapIt = clMap.begin(); mapIt != clMap.end(); ++mapIt)
	{
		if((*mapIt).second > 0)
			map[(*mapIt).second]++;
	}

	delete[] labels;
}

void TriangularLattice::clearExchanger()
{
    if ( selfLattExchanger )
    {
        delete mpExchanger;
        selfLattExchanger = false;
    }
}

void TriangularLattice::setExchanger( LattExchanger* exchanger )
{
    clearExchanger();
    mpExchanger = exchanger;
}

LattExchanger* TriangularLattice::getExchanger()
{
    return mpExchanger;
}

        
TriangularLattice::~TriangularLattice()
{
    clearExchanger();
    delete[] this->mpLattice;
}

ostream &operator<<( ostream &stream, TriangularLattice &latt )
{
    stream << latt.getLatticeSize() << endl;
    stream << "Simulation" << endl;
    for ( lattIndex i = 0; i < latt.getLatticeSize(); i++ )
        if ( latt[i] == LIPID_A )
        {
	    lattIndex line( ( latt.getLatticeSize() - i ) / latt.getRowSize() );
            double y = line * 0.866025;
            double x = ( i % latt.getRowSize() ) - line * 0.5;
            stream << "A\t" << setprecision( 8) << x << "\t" << y << "\t0.00000000" << endl;
        }
    for ( lattIndex i = 0; i < latt.getLatticeSize(); i++ )
    {
        if ( latt[i] == LIPID_B )
        {
	    lattIndex line( ( latt.getLatticeSize() - i ) / latt.getRowSize() );
            double y = line * 0.866025;
            double x = ( i % latt.getRowSize() ) - ( line * 0.5 ); 
            stream << "B\t" << setprecision( 8) << x << "\t" << y << "\t0.00000000" << endl;
        }
    }
    for ( lattIndex i = 0; i < latt.getLatticeSize(); i++ )
    {
        if ( latt[i] == PROTEIN_A )
        {
	    lattIndex line( ( latt.getLatticeSize() - i ) / latt.getRowSize() );
            double y = line * 0.866025;
            double x = ( i % latt.getRowSize() ) - ( line * 0.5 ); 
            stream << "PROTEIN\t" << setprecision( 8) << x << "\t" << y << "\t0.00000000" << endl;
        }
    }
    for ( lattIndex i = 0; i < latt.getLatticeSize(); i++ )
    {
        if ( latt[i] == PROTEIN_B )
        {
	    lattIndex line( ( latt.getLatticeSize() - i ) / latt.getRowSize() );
            double y = line * 0.866025;
            double x = ( i % latt.getRowSize() ) - ( line * 0.5 ); 
            stream << "PROTEIN_B\t" << setprecision( 8) << x << "\t" << y << "\t0.00000000" << endl;
        }
    }
    return stream;
}

