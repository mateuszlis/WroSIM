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
#include "InteractionsTracker.h"

// tests
#include "testUtils.h"
using namespace testUtils;


TEST( KawasakiProteins, calcEnergy )
{
    static const int LIPID_B_COUNT( 12 );
    static const int PROTEIN_COUNT( 1 );
    std::tr1::shared_ptr< ProteinTriangularLattice > latt( new ProteinTriangularLattice( 64, 8, LIPID_B_COUNT, PROTEIN_COUNT, false ) );
    printLatt( latt->getLattice(), 8, 8 );
    ProteinTriangularLattice* localPtr( latt.get() ); // constructor in line below requires raw ptr
    std::tr1::shared_ptr< KawasakiProteins > simulation( new KawasakiProteins( localPtr, 0, -400, 800 ) );
    EXPECT_DOUBLE_EQ( simulation->calcEnergy(), 18.499995 * -400 + 6.499999 * 800 );
    InteractionsTracker tracker( latt.get(), latt->mProteins, PROTEIN_COUNT );
    testUtils::printLatt( tracker.mpProtInt.get(), 8, 8 );
}

#undef protected
