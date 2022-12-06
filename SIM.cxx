#include <sstream>
using namespace std;

#include "m.h"
class RandomFile{
public:
    ifstream r;
    RandomFile(string filename){
        r = ifstream(filename, std::ifstream::in);
        if(!r){
            cerr << "Error opening random file" << endl;
            ::exit(1);
        }
    }
    double getU(){
        double u = 0.0;
        if(!r.eof()) {
            r >> u;
        }
        else{
            cerr << "Ran out of random numbers" << endl;
            ::exit(1);
        }
        return u;
    }
    ~RandomFile(){
        r.close();
    }
};

class Welford{
public:
    int i = 0;
    double xibar = 0.0;
    double vi = 0.0;

    void addDataPoint(double xi){
        double diff = xi - xibar;
        i++;
        vi = vi + ((i - 1) / (double)i) * pow(diff, 2);
        xibar = xibar + (1 / (double)i) * diff;
    }

    double getStandardDeviation(){
        return sqrt(vi / i);
    }
};

double Uniform(double alpha, double beta, double u){
    return alpha + (beta - alpha) * u;
}

long Equilikely(double alpha, double beta, double u){
    return (alpha + (long)((beta - alpha + 1) * u));
}

int main( int argc, char* argv[] ){
	if( argc == 4 ) {
		experiment = argv[1][0];
		istringstream( argv[2] ) >> seed;
		istringstream( argv[3] ) >> threshold;
	}
	return 0;
}

