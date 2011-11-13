#include <iostream>
#include <string>
#include <cstdlib>

using namespace std;
#pragma once
struct Atom 
{
    public:
        string resname;
        string name;
        double x;
        double y;
        double z;
        Atom() {}
        Atom(string resname, string name, double x, double y, double z):
            resname( resname),
            name( name ),
            x( x ), y( y ), z( z )
        {}
};

struct Distance {
public:
    double d;
    int at2Ind;
    Distance(double d, int at2Ind): 
        at2Ind(at2Ind), d(d)
    {}
    ~Distance() {}
};

bool compareDistances(const Distance& d1, const Distance& d2);

