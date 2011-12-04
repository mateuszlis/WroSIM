#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>
#include <iomanip>
#include "Atom.h"
#include "ClustersAnalyzer.h"
using namespace std;

#include "analyzePyrkovaCommandLine.hxx" //command line library
const double TEMPERATURE = 333.0;
const double R = 1.985877;
const string AT1 = "DOPCPO4";
const string AT2 = "DPPCPO4";

void usage()
{
    cerr << "usage: analyzePyrkova [options]" << endl << "options:" << endl;
    options::print_usage( cerr);
}

inline std::string trim(std::string str)
{
    str.erase(0, str.find_first_not_of(' '));       //prefixing spaces
    str.erase(str.find_last_not_of(' ')+1);         //surfixing spaces
return str;
}

int main(int argc,char *argv[]) {
    float minDist;
    int minFrames;
    string file;
    try // handling command line
    {
        int end;
        options opt( argc, argv, end );

        if ( opt.help() )
        {
            usage();
            return 0;
        }
        minDist = opt.d();
        minDist *= minDist; // we don't calculate square root when calculating distances
        // so we're interested in square of threshold distance
        minFrames = opt.t();
        file = opt.f();
    }
    catch ( const cli::exception& e )
    {
        cerr << e << endl;
        usage();
        return 1;
    }
        
    int global_n_atoms;
    Atom *atoms;
    long double box_x( 0 ), box_y( 0 );

    //processing GRO file 
    ifstream ifile;
    ifile.open(file.c_str());
    string line;
    ClustersAnalyzer analyzer( minFrames, minDist ); // ClustersAnalyzer is heart of algorithm
    double omega_sumAt1( 0 ), omega_sumAt2( 0 ); // use this to calculate average
    int frameCounter( 0 );
    cout << setw(6) << "#Frame" << setw(10) << "omegaDOPC" << setw(10) << "omegaDPPC" <<  "  fraction of lipids clustered" << endl;
    if (ifile.is_open())
    {
        double box_x, box_y;//, box_z;
        int n_atoms;
        while (ifile.good() )
        {
            int i = 0;
            n_atoms = 0;
            while ( (i-3) <= n_atoms && ifile.good())
            {
                getline(ifile, line);
                ++i;
                if (i==2)
                {
                    n_atoms = atoi(line.c_str());
                    atoms = new Atom[n_atoms];
                }
                if (i < 3)
                    continue;
                if (i == n_atoms+3)
                {
                    box_x += atof(trim(line.substr(0,10)).c_str());
                    box_y += atof(trim(line.substr(10,10)).c_str());
                    //box_z += atof(trim(line.substr(20,10)).c_str());
                    break;
                }
                
                string resname = trim(line.substr(5, 5));
                string name = trim(line.substr(10, 5));
                double x = atof(trim(line.substr(20,8)).c_str());
                double y = atof(trim(line.substr(28,8)).c_str());
                double z = atof(trim(line.substr(36,8)).c_str());
                if ( i - 3 >= 0 ) 
                    atoms[i-3] = Atom(resname, name, x, y, z);


            }
            if ( ! ifile.good() )
                break;
            //end of processing GRO file
            //analyzing
            global_n_atoms = n_atoms;
            for (int i=0; i<n_atoms; ++i)
            {
                vector<Distance> distances2; //will collect all distances to j-th atoms
                for (int j=0; j<n_atoms; ++j)
                {
                    if (i==j) continue;
                    double dx = abs(atoms[i].x - atoms[j].x);
                    double dy = abs(atoms[i].y - atoms[j].y);
                    
                    if (dx > 0.5*box_x) // fix pbc
                        dx = box_x - dx;
                    if (dy > 0.5*box_y)
                        dy = box_y - dy;

                    double dist2 = dx*dx + dy*dy;
                    
                    distances2.push_back(Distance(dist2, j));
                }
                
                sort(distances2.begin(), distances2.end(), compareDistances); //the shortest distances go first
                analyzer.registerAtom( i, distances2, frameCounter, atoms );
            }
            int nonMixedAt1( 0 ), nonMixedAt2( 0 ); // store counts of atoms separately to calculate omega_AB
            // for each lipid
            int mixedAt1( 0 ), mixedAt2( 0 );
            for (int i=0; i<n_atoms; ++i)
            {
                if ( analyzer.isClustered( i, frameCounter ) )
                {
                    if ( analyzer.isInMixedCluster( i, frameCounter ) )
                    {
                        if ( atoms[i].resname + atoms[i].name  == AT1 )
                            mixedAt1++;
                        else 
                            mixedAt2++;
                    }
                    else
                    {
                        if ( atoms[i].resname + atoms[i].name == AT2 )
                            nonMixedAt1++;
                        else
                            nonMixedAt2++;
                    }
                }
            }

            double omega_abAt1( 0 ), omega_abAt2;
            if ( mixedAt1 && frameCounter > minFrames )
            {
                double PAt1;
                PAt1 = static_cast< double >( mixedAt1 ) / nonMixedAt1;
                omega_abAt1 = - R * TEMPERATURE * log(PAt1);
                omega_sumAt1 += omega_abAt1;
            }
            if ( mixedAt2 && frameCounter > minFrames )
            {
                double PAt2;
                PAt2 = static_cast< double >( mixedAt2 ) / nonMixedAt2;
                omega_abAt2 = - R * TEMPERATURE * log(PAt2);
                omega_sumAt2 += omega_abAt2;
            }

            float clusteredPercent = static_cast< double > ( nonMixedAt1 + nonMixedAt2 + mixedAt1 + mixedAt2 ) / n_atoms;
            cout << setw(6) << frameCounter <<  setw(10) << omega_abAt1 << setw(10)  << omega_abAt2 << setw(10)  << clusteredPercent << endl;
            delete[] atoms; 
            frameCounter++;
        }
            cout << "atoms " << n_atoms << endl;
    }
    else cout << "Unable to open file!" << endl;
    ifile.close();
    cout << "#final " << omega_sumAt1 / frameCounter << "\t" << omega_sumAt2 / frameCounter << endl;
    double avg_box_x( box_x / frameCounter ), avg_box_y( box_y / frameCounter );
    cout << "#APL " << ( avg_box_x * avg_box_y ) / global_n_atoms << "\t" 
        << "box size\t" << avg_box_x << " " << avg_box_y <<  endl;

    return 0;
}
