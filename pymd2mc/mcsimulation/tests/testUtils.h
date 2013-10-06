
// stl
#include <set>

// project related
#include "types.h"
#include "TriangularLattice.h"

namespace testUtils
{
    void printPermutation( lattIndex *latt, int rowSize, int rowsCount );

    template < typename T >
    void printLatt( T *latt, int rowSize, int rowsCount )
    {
        for ( int i = 0 ; i < rowsCount ; ++i )
        {
            cout << setw( 2 * ( i + 1 ) ) << "  " ;
            for( int j = 0 ; j < rowSize; ++j )
            {
                cout << setw(2) << static_cast< double >( latt[ i * rowSize + j ] ) << "  ";
            }
            cout << endl;
        }
    }
    int calcSum( TriangularLattice* );
    std::set< int > getNeighborsOf( int site, TriangularLattice* latt );

    /**
     * @brief for each element from elements, checks if neighbor is its neighbor on
     * the lattice latt
     */
    bool checkIfIsNeighborOf( const int& neighbor, const std::set< int >& elements, TriangularLattice* latt );
}; // namespace testUtils
