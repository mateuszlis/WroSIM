#include <iostream>
#include <fstream>
#include <iomanip>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <vector>
#include <map>
#include <algorithm>

using namespace std;

const double TEMPERATURE = 333.0;
const int N_NEIGHBORS = 6;
const string AT1 = "DOPCP8";
const string AT2 = "POPCP8";

class Atom {
public:
    string resname;
    string name;
    double x;
    double y;
    double z;
    Atom()
        : x( 0 )
          , y( 0 )
          , z( 0 )
    {}
    Atom(string resname, string name, double x, double y, double z)
        : resname( resname ) 
        , name( name )
        , x( x ) 
        , y( y )
        , z( z )
    {}
};

        
        
class Distance {
public:
    int at2_ind;
    double d;
    Distance(double d, int at2_ind): 
        at2_ind(at2_ind), d(d)
    {}
    ~Distance() {}
};

bool compareDistances(const Distance& d1, const Distance& d2)
{
    if (d1.d < d2.d)
        return true;
    else return false;
}

inline std::string trim(std::string str)
{
    str.erase(0, str.find_first_not_of(' '));       //prefixing spaces
    str.erase(str.find_last_not_of(' ')+1);         //surfixing spaces
return str;
}

int main(int /*argc*/,char *argv[]) {

    Atom *atoms;
    int global_n_atomsAT1 = 0;
    int global_n_atomsAT2 = 0;

    //processing GRO file (first frame only)
    ifstream ifile;
    ifile.open(argv[1]);
    string line;
    double omega_sum( 0 );
    int frameCounter( 0 );
    // this will be utilized to calc neighbors histogram
    double neighHistAT1[] = { 0, 0, 0, 0, 0 ,0, 0 };
    int neighHistCalcAT1[] = { 0, 0, 0, 0, 0, 0, 0 };
    double neighHistAT2[] = { 0, 0, 0, 0, 0 ,0, 0 };
    int neighHistCalcAT2[] = { 0, 0, 0, 0, 0, 0, 0 };
    // fraction of first neighbors
    vector< double > firstNeighFr;

    cout << setw(10) << "#Frame" << setw(10) << "omega_AB" << endl;
    
    if (ifile.is_open())
    {
        while (ifile.good())
        {
            int i = 0;
            int n_atoms = 0;
            double box_x( 0 ), box_y( 0 );
            while ( ( i - 3 ) <= n_atoms && ifile.good() )
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
                    //box_z = atof(trim(line.substr(20,10)).c_str());
                    break;
                }
                
                string resname = trim(line.substr(5, 5));
                string name = trim(line.substr(10, 5));
                double x = atof(trim(line.substr(20,8)).c_str());
                double y = atof(trim(line.substr(28,8)).c_str());
                double z = atof(trim(line.substr(36,8)).c_str());

                atoms[i-3] = Atom(resname, name, x, y, z);
            }
            if ( ! ifile.good() ) 
                break;
            //end of processing GRO file

            //analyzing
            double n_aa = 0;
            double n_bb = 0;
            double n_ab = 0;
            int firstSimNeigh = 0;
            for (int i=0; i<n_atoms; ++i)
            {
                if ( atoms[i].resname + atoms[i].name == AT1 )
                {
                    ++global_n_atomsAT1;
                }
                else
                {
                    ++global_n_atomsAT2;
                }
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
                
                
                
                int simNeighbors( 0 );
                int* neighHistCalc;
                if ( atoms[i].resname + atoms[i].name == AT1 )
                {
                    neighHistCalc = neighHistCalcAT1;
                }
                else 
                {
                    neighHistCalc = neighHistCalcAT2;
                }
                if ( atoms[i].resname + atoms[i].name == atoms[distances2[0].at2_ind].resname + atoms[distances2[0].at2_ind].name )
                {
                    ++firstSimNeigh;
                }


                for (int n=0; n<N_NEIGHBORS; ++n)
                {
                    if (atoms[i].resname+atoms[i].name == atoms[distances2[n].at2_ind].resname+atoms[distances2[n].at2_ind].name) 
                    {
                        ++simNeighbors;
                        if (atoms[i].resname+atoms[i].name == AT1)
                            ++n_aa;
                        else ++n_bb;
                    }
                    else ++n_ab;
                        
                }
                ++neighHistCalc[ simNeighbors ];
                
            }
            firstNeighFr.push_back( static_cast< double >( firstSimNeigh ) / ( n_atoms ) );
            

            
            n_aa = 0.5*n_aa; n_ab = 0.5*n_ab; n_bb = 0.5*n_bb; //0.5 since we calculated each neighboring pair twice

            //cout << "N: " << n_aa << " " << n_ab << " " << n_bb << endl;

            double K = n_ab * n_ab / (n_aa * n_bb);
            double omega_ab = -0.5 * 1.9859 * TEMPERATURE * log(K/4.0); //in cal*mol-1*K-1

            cout << setw(15) << frameCounter << setw(15) <<  omega_ab << endl;
            omega_sum += omega_ab;
            delete[] atoms; 
            frameCounter++; 
        }
    }
    double finalOmega( omega_sum / frameCounter );
    cout << setw(10) << "Final" << setw(10) <<  finalOmega << endl;
    
    ofstream neighFileAT1;
    neighFileAT1.open( string( "neighbHist" + AT1 + ".dat").c_str() );
    ofstream neighFileAT2;
    neighFileAT2.open( string( "neighbHist" + AT2 + ".dat" ).c_str() );
    for ( int i = 0 ; i < 7 ; ++i )
    {
        neighHistAT1[ i ] = neighHistCalcAT1[ i ] / ( static_cast< double >( global_n_atomsAT1 ) );
        neighFileAT1 <<  i << "\t" << neighHistAT1[ i ] << endl;
        neighHistAT2[ i ] = neighHistCalcAT2[ i ] / ( static_cast< double >( global_n_atomsAT2 ) );
        neighFileAT2 <<  i << "\t" << neighHistAT2[ i ] << endl;
    }
    // write firction of first neighbors
    ofstream fraction;
    fraction.open( "firstNeighbFr.dat" );
    int i = 0;
    for ( vector< double >::const_iterator it = firstNeighFr.begin() ; it != firstNeighFr.end() ; ++it, ++i )
    {
        fraction << i << "\t" << ( *it ) << endl;
    }
    return 0;
}
