#include "Atom.h"


bool compareDistances(const Distance& d1, const Distance& d2)
{
    if (d1.d < d2.d)
        return true;
    else return false;
}


