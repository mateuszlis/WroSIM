#include <iostream>
#include <fstream>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>

#include "Atom.h"
#include "ClustersAnalyzer.h"
using namespace std;

#include "analyzePyrkovaCommandLine.hxx" //command line library
const double TEMPERATURE = 333.0;
const int N_NEIGHBORS = 6;
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
    //cout << "Starting..." << endl;
    float minDist;
    int minFrames;
    string file;
    try
    {
        int end;
        options opt( argc, argv, end );

        if ( opt.help() )
        {
            usage();
            return 0;
        }
        minDist = opt.d();
        minFrames = opt.t();
        file = opt.f();
    }
    catch ( const cli::exception& e )
    {
        cerr << e << endl;
        usage();
        return 1;
    }
        
    int n_atoms;
    Atom *atoms;
    double box_x, box_y, box_z;

    //processing GRO file (first frame only)
    ifstream ifile;
    ifile.open(file.c_str());
    string line;
    ClustersAnalyzer analyzer( minFrames, minDist );
    double omega_sum = 0;
    int frameCounter( 0 );
    if (ifile.is_open())
    {
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
                    box_x = atof(trim(line.substr(0,10)).c_str());
                    box_y = atof(trim(line.substr(10,10)).c_str());
                    box_z = atof(trim(line.substr(20,10)).c_str());
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
            double n_aa = 0;
            double n_bb = 0;
            double n_ab = 0;
             
            for (int i=0; i<n_atoms; ++i)
            {
                vector<Distance> distances2; //will collect all distances to j-th atoms
                for (int j=0; j<n_atoms; ++j)
                {
                    if (i==j) continue;
                    double dx = abs(atoms[i].x - atoms[j].x);
                    double dy = abs(atoms[i].y - atoms[j].y);
                    
                    if (dx > 0.5*box_x)
                        dx = box_x - dx;
                    if (dy > 0.5*box_y)
                        dy = box_y - dy;

                    double dist2 = dx*dx + dy*dy;
                    
                    distances2.push_back(Distance(dist2, j));

                    
                }
                
                sort(distances2.begin(), distances2.end(), compareDistances); //the shortest distances go first
                analyzer.registerAtom( i, distances2, frameCounter, atoms );

                //cout << "AT=" << atoms[i].name << endl;
                
                
                for (int n=0; n<N_NEIGHBORS; ++n)
                {
                    if (atoms[i].resname+atoms[i].name == atoms[distances2[n].at2Ind].resname+atoms[distances2[n].at2Ind].name) {
                        if (atoms[i].resname+atoms[i].name == AT1)
                            ++n_aa;
                        else ++n_bb;
                    }
                    else ++n_ab;
                        
                }
                
            }
            int nonMixed( 0 );
            int mixed( 0 );
            for (int i=0; i<n_atoms; ++i)
            {
                if ( analyzer.isClustered( i, frameCounter ) )
                {
                    if ( analyzer.isInMixedCluster( i, frameCounter ) )
                    {
                        mixed++;   
                    }
                    else
                    {
                        nonMixed++;
                    }
                }
            }
            n_aa = 0.5*n_aa; n_ab = 0.5*n_ab; n_bb = 0.5*n_bb; //0.5 since we calculated each neighboring pair twice

            //cout << "N: " << n_aa << " " << n_ab << " " << n_bb << endl;
            double omega_ab( 0 );
            double P;
            if ( mixed && frameCounter > minFrames )
            {
                P = static_cast< double >( mixed ) / nonMixed;
                omega_ab = - 1.985877 * TEMPERATURE * log(P);
                omega_sum += omega_ab;
            }
            //double K = n_ab * n_ab / (n_aa * n_bb);
            //double omega_ab = -0.5 * 1.9859 * TEMPERATURE * log(K/4.0); //in cal*mol-1*K-1
            float clusteredPercent = static_cast< double > ( nonMixed + mixed ) / n_atoms;
            cout << frameCounter << "\t" << omega_ab << "\t" << P << "\t" << clusteredPercent << endl;
            delete[] atoms; 
            frameCounter++;
        }
    }
    else cout << "Unable to open file!" << endl;
    ifile.close();
    cout << "#final " << omega_sum / frameCounter << endl;
    return 0;
}
