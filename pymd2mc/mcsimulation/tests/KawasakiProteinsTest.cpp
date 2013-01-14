// std
#include <iostream>
#include <assert.h>
#include <tr1/memory>

// stl
#include <set>

// Google test tools
#include "gtest/gtest.h"
#include "gmock/gmock.h"

// project-local
#include "ProteinTriangularLattice.h"

#define protected public // to test hidden features
#include "KawasakiProteins.h"

// tests
#include "testUtils.h"
using namespace testUtils;


TEST( KawasakiProteins, calcEnergy )
{
    static const int LIPID_B_COUNT( 1 );
    static const int PROTEIN_COUNT( 1 );
    std::tr1::shared_ptr< ProteinTriangularLattice > latt( new ProteinTriangularLattice( 9, 3, LIPID_B_COUNT, PROTEIN_COUNT, false ) );
    printLatt( latt->getLattice(), 3, 3 );
    ProteinTriangularLattice* localPtr( latt.get() ); // constructor in line below requires raw ptr
    std::tr1::shared_ptr< KawasakiProteins > simulation( new KawasakiProteins( localPtr, 0, 310, 0 ) );
    EXPECT_EQ( simulation->calcEnergy(), 5 * -400 + 5 * 800 );
}

#undef protected
