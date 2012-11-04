// Google test tools
#include "gtest/gtest.h"
#include "gmock/gmock.h"

// project-local
#include "MsdEnabledLattEx.h"
#include "TriangularLattice.h"

TEST( MsdEnabledLattEx, Construct )
{
    TriangularLattice *latt = new TriangularLattice( 100, 10, 20);
    MsdEnabledLattEx ex( latt );
    delete latt;
}

TEST( TriangularLattice, ConstructException )
{
    TriangularLattice *latt;
    try
    {
        latt = new TriangularLattice( 100, 10, 101 );
    }
    catch ( InputParametersException& e )
    {
        SUCCEED();
        return;
    }
    FAIL();
}

