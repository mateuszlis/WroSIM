#include "testUtils.h"

namespace testUtils
{
    void printPermutation( lattIndex *latt, int rowSize, int rowsCount )
    {
        for ( int i = 0 ; i < rowsCount ; ++i )
        {
            cout << setw( 2 * ( i + 1 ) ) << "  " ;
            for( int j = 0 ; j < rowSize; ++j )
            {   
                cout << setw(2) << static_cast< int >( latt[ i * rowSize + j ] ) << "  ";
            }
            cout << endl;
        }
    }
    void printLatt( lattMember *latt, int rowSize, int rowsCount )
    {
        for ( int i = 0 ; i < rowsCount ; ++i )
        {
            cout << setw( 2 * ( i + 1 ) ) << "  " ;
            for( int j = 0 ; j < rowSize; ++j )
            {   
                cout << setw(2) << static_cast< int >( latt[ i * rowSize + j ] ) << "  ";
            }
            cout << endl;
        }
    }
    int calcSum( TriangularLattice *latt )
    {
        int sum = 0;
        for ( int i = 0; i < latt->getLatticeSize(); i++ )
        {
            sum += ( *latt )[i];
        }
        return sum;
    }

    std::set< int > getNeighborsOf( int site, TriangularLattice* latt )
    {
        std::set< int > neighbors;
        for( int neighNum = 0 ; neighNum < 6 ; ++neighNum )
        {
            neighbors.insert( latt->getNeighbIndex( site, neighNum ) );
        }
        return neighbors;
    }

    bool checkIfIsNeighborOf( const int& neighbor, const std::set< int >& elements, TriangularLattice* latt )
    {
        for ( std::set< int >::const_iterator it = elements.begin() 
                ; it != elements.end() 
                ; ++it )
        {
            std::set< int > elementsNeighbors = getNeighborsOf( *it, latt );
            if ( elementsNeighbors.find( neighbor ) == elementsNeighbors.end() )
            {
                return false;
            }
        }
        return true;
    }

} // namespace testUtils
