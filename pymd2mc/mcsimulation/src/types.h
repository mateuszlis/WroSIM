#pragma once

#define likely(x)      __builtin_expect(!!(x), 1)
#define unlikely(x)    __builtin_expect(!!(x), 0)

typedef unsigned char lattMember;
typedef int lattIndex;
enum LATTICE_FIELD_NAMES 
{
    LIPID_A = 0,
    PROTEIN_A = 120,
    PROTEIN_B = 200,
    LIPID_B = 255
};
bool isProtein( lattMember memb );


