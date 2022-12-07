#include <sstream>
#include <queue>
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

class Order{
public:
    int numberOfCars;
    double deadline;
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

class Event{
public:
    string type;
    double at;

    bool operator<( const Event& rhs ) const {
        // .at is activation time of the event
        if( this->at == rhs.at) {
            return this->type < rhs.type;
        }
        // inverted! we want the least // activation time to have
        // higher priority
        return !( this->at < rhs.at );
    }
};

double Uniform(double alpha, double beta, double u){
    return alpha + (beta - alpha) * u;
}

double Triangular(double a, double c, double b, double u){
    if(u <= (c - a) / (b - a)){
        return a + sqrt((b - a) * (c - a) * u);
    }
    else{
        return b - sqrt((b - a) * (b - c) * (1.0 - u));
    }
}

long Equilikely(double alpha, double beta, double u){
    return (alpha + (long)((beta - alpha + 1) * u));
}

void runSim(Welford &w, RandomFile &r, int S, int s){
    double t = 0.0;

    priority_queue<Event, vector<Event> > eventList = priority_queue<Event, vector<Event> >();
    int inventoryLevel = S;

    while(!eventList.empty()){
        Event currentEvent = eventList.top();
        eventList.pop();
        t = currentEvent.at;

        //TODO process events here
        if(currentEvent.type == "inventoryReview"){
            if(inventoryLevel < s){
                //TODO schedule inventory restock
            }
        }
        else if(currentEvent.type == "inventoryRestock"){

        }
    }
}

int main( int argc, char* argv[] ){
	//TODO take in arguments
    int S, s;


	return 0;
}

