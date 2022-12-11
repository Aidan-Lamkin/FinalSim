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
    double orderPlacementTime;  // Save the order placement time to calculate final discount to calculate penalty
     bool operator<( const Order& rhs ) const {
        return !( this->deadline < rhs.deadline );
    }
};

class Welford{
public:
    int i = 0;
    double ti = 0.0;

    int l;
    int c = 0;

    double previousDeliveryTime = 0.0;

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

class Demand{
public:
    double tLeft;
    double tRight;
    double arrivalLeft;
    double arrivalRight;
    double slope;
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

double Exponential(double mew, double u){
    return (-mew * log(1.0 - u));
}
double calculateDeliveryLag(double previousTime, double q, double a, double b, double c, double u){
    // TODO: Fix this when values: A, M, Pi are provided 
    // return previousTime + ((A + q) / M) + Pi + Triangular(a, b, c, u);
    return previousTime + ((387 + q) / 13.0) + 23.7 + Triangular(a, c, b, u); 
}

void runSim(Welford &w, RandomFile &r, double a, double b, double c, int S, int s, double start, double end){
    double t = 0.0;

    priority_queue<Event, vector<Event> > eventList = priority_queue<Event, vector<Event> >();
    priority_queue<Order, vector<Order> > backOrders = priority_queue<Order, vector<Order> >();

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
                nextRestock.at = t + calculateDeliveryLag(w.previousDeliveryTime, nextRestock.numberOfCars, a, b, c, r.getU());
                eventList.push(nextRestock);
            }

            //Scheduling next inventory review in 60 working hours
            Event nextReview;
            nextReview.type = "inventoryReview";
            nextReview.at = t + 60;
            eventList.push(nextReview);
        }
        else if(currentEvent.type == "inventoryRestock"){
            w.previousDeliveryTime = t; // update the previous delivery time to current time
            w.l += currentEvent.numberOfCars;

            // Everytime we get new orders, cycle through the N new cars arriving, and calculate the penalty
            for(int i = 0; i < currentEvent.numberOfCars; i++){
                Order currOrder = backOrders.top();
                backOrders.pop();
                
                w.penalties += floor((t - currOrder.orderPlacementTime) / 10.0) * 100;
            }
        }
        else if(currentEvent.type == "carDemand"){
            if(w.l <= 0){ // no cars in stock, increase order count
                //TODO need to figure out deadline
                Order order;
                order.numberOfCars = currentEvent.numberOfCars;
                order.orderPlacementTime = t;
                backOrders.push(order);
                w.c += currentEvent.numberOfCars; // increase the number of orders placed

                while(Equilikely(0, 1, r.getU()) == 1){
                    double u = r.getU();
                    double h = Uniform(2, 9, u);

                    Event demand;
                    demand.at = t + h;
                    demand.type = "rumorMillCarDemand";
                    demand.numberOfCars = 1;
                    eventList.push(demand);
                }
            } else { // cars in stock, decrease level - don't need to schedule any new orders
                w.l -= currentEvent.numberOfCars;
            }
        }
        else if(currentEvent.type == "rumorMillCarDemand"){
            if(w.l <= 0){ // only order if there are no cars in stock
                //TODO need to figure out deadline
                Order order;
                order.numberOfCars = currentEvent.numberOfCars;
                order.orderPlacementTime = t;
                backOrders.push(order);
                w.c += currentEvent.numberOfCars; // increase the number of orders placed

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
        if(t >= start){
            if(w.l < w.minInventory){
                w.minInventory = w.l;
            }
        }
        if(t > end){
            cout << "OUTPUT MININVENTORY" << w.minInventory << endl;
            cout << "OUTPUT PENALTIES" << w.penalties << endl;
            cout << "OUTPUT MAXRUMORMILL" << w.maxRumorMill << endl;
            cout << "OUTPUT ORDERS" << w.orders << endl;
            ::exit(0);
        }
    }
}

void runTriangle(RandomFile &r, double a, double b, double c){
    while(true){
        double u = r.getU();
        cout << "OUTPUT:" << Triangular(a, b, c, u);
    }
}

double nextArrival(double previousArrival, vector<Demand> demands, RandomFile &r){
    int last = demands.size() - 1;
    double tK = demands[last].tRight;
    double cumulativeArrivals = demands[last].arrivalRight;

    double nextArrival = previousArrival - floor(previousArrival / tK) * tK;
    int ji = 0;
    while(nextArrival > demands[ji].tRight){
        ji++;
    }

    double u = r.getU();
    double A = previousArrival + Exponential(1.0 , u);
}

int main( int argc, char* argv[] ){
	//Take in the arguments from command line
    int S, s;
    double a, b, c;
    string runMode;
    double start, end;
    ifstream definitionFile;

    Welford w = Welford();
    vector<Demand> demands = vector<Demand>();

    // check that simulation has enough parameters - prevents segfault
    if(argc < 11){
        cerr << "Missing Required Parameters." << endl;
        exit(-1);
    }

    runMode = argv[1];

    string randomFileName(argc > 2 ? argv[2] : "/dev/null");
    RandomFile r(randomFileName);

    if(runMode != "TRIANGLE"){
        definitionFile = ifstream(argc > 3 ? argv[3] : "/dev/null", std::ifstream::in);
        if(!definitionFile){
            cerr << "Error opening definition file" << endl;
                ::exit(1);
        }
        else{
            //fills in demand objects
            double startT = 0.0;
            double L = 0.0;
            while(!definitionFile.eof()){
                Demand d;
                d.tLeft = startT;
                definitionFile >> d.tRight;
                definitionFile >> d.slope;
                d.arrivalLeft = L;
                double l = d.slope * (d.tRight - d.tLeft);
                L += l;
                d.arrivalRight = L;
                startT = d.tRight;
                demands.push_back(d);
            }
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
        runSim(w, r, a, b, c, S, s, start, end);
    }
    else if(runMode == "TRIANGLE"){
        runTriangle(r, a, b, c);
    }
	return 0;
}

