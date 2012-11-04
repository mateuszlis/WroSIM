#include "MsdEnabledLattEx.h"

bool const MsdEnabledLattEx::isNeighbor[3][3] = {
    { 0, 1, 1 },
    { 1, 1, 1 },
    { 1, 1, 0 }
};

ostream &operator<<( ostream &stream, vectorDist & dist )
{
    stream << "row=" << dist.row << " col=" << dist.col << ". ";
    return stream;
}
