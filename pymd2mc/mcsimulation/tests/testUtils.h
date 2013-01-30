
// stl
#include <set>

// project related
#include "types.h"
#include "TriangularLattice.h"

namespace testUtils
{
    void printPermutation( lattIndex *latt, int rowSize, int rowsCount );
    void printLatt( lattMember *latt, int rowSize, int rowsCount );
    int calcSum( TriangularLattice* );
    std::set< int > getNeighborsOf( int site, TriangularLattice* latt );

    /**
     * @brief for each element from elements, checks if neighbor is its neighbor on
     * the lattice latt
     */
    bool checkIfIsNeighborOf( const int& neighbor, const std::set< int >& elements, TriangularLattice* latt );
}; // namespace testUtils
