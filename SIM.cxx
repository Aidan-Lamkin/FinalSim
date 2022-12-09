#include <sstream>
#include <queue>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
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
    double ti = 0.0;

    int l;
    int c = 0;;

    int minInventory = 99999;
    int penalties = 0;
    int maxRumorMill = 0;
    int orders = 0;

};

class Event{
public:
    string type;
    double at;
    int numberOfCars;

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

void runSim(Welford &w, RandomFile &r, int S, int s, double start, double end){
    double t = 0.0;

    priority_queue<Event, vector<Event> > eventList = priority_queue<Event, vector<Event> >();
    int inventoryLevel = S;

    while(!eventList.empty()){
        Event currentEvent = eventList.top();
        eventList.pop();
        t = currentEvent.at;

        //TODO process events here
        if(currentEvent.type == "inventoryReview"){
            if(w.l + w.c <= s){
                Event nextRestock;
                nextRestock.type = "inventoryRestock";
                nextRestock.numberOfCars = S - (w.l + w.c);
            }

            //Scheduling next inventory review in 60 working hours
            Event nextReview;
            nextReview.type = "inventoryReview";
            nextReview.at = t + 60;
            eventList.push(nextReview);
        }
        else if(currentEvent.type == "inventoryRestock"){
            
        }
        else if(currentEvent.type == "carDemand"){
            if(w.l <= 0){
                while(Equilikely(0, 1, r.getU()) == 1){
                    double u = r.getU();
                    double h = Uniform(2, 9, u);

                    Event demand;
                    demand.at = t + h;
                    demand.type = "rumorMillCarDemand";
                    demand.numberOfCars = 1;
                    eventList.push(demand);
                }
            }
        }
    }
}

void runTriangle(RandomFile &r, double a, double b, double c){
    while(true){
        double u = r.getU();
        cout << "OUTPUT:" << Triangular(a, b, c, u);
    }
}

int main( int argc, char* argv[] ){
	//Take in the arguments from command line
    int S, s;
    double a, b, c;
    string runMode;
    double start, end;
    ifstream definitionFile;

    Welford w = Welford();

    runMode = argv[1];

    string randomFileName(argc > 2 ? argv[2] : "/dev/null");
    RandomFile r(randomFileName);

    if(runMode != "TRIANGLE"){
        definitionFile = ifstream(argc > 3 ? argv[3] : "/dev/null", std::ifstream::in);
        if(!definitionFile){
            cerr << "Error opening definition file" << endl;
                ::exit(1);
        }
    }

    a = atof(argv[4]);
    b = atof(argv[5]);
    c = atof(argv[6]);

    

    S = stoi(argv[7]);
    s = stoi(argv[8]);

    start = atof(argv[9]);
    end = atof(argv[10]);

    w.l = S;
    w.c = 0;

    if(runMode == "SIM"){
        runSim(w, r, S, s, start, end);
    }
    else if(runMode == "TRIANGLE"){
        runTriangle(r, a, b, c);
    }
	return 0;
}

