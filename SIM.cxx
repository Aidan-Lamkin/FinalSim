#include <sstream>
#include <queue>
#include <cmath>
#include <iostream>
#include <fstream>
#include <stdlib.h>
using namespace std;

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
            ::exit(0);
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
    double orderPlacementTime;  // Save the order placement time to calculate final discount to calculate penalty
    bool operator<( const Order& rhs ) const {
        return !( this->orderPlacementTime < rhs.orderPlacementTime );
    }
};

class Welford{
public:
    int i = 0;
    double ti = 0.0;

    int inventory;
    int onOrder = 0;

    int backorderCount = 0;

    double previousDeliveryTime = 0.0;

    int minInventory = 99999;
    long penalties = 0;
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
    if(u <= (double(c - a) / double(b - a))){
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

double nextArrival(double previousArrival, vector<Demand> demands, RandomFile &r){
    int last = demands.size() - 1;
    double tK = demands[last].tRight;
    double cumulativeArrivals = demands[last].arrivalRight;

    double ai = previousArrival - floor(double(previousArrival) / tK) * tK;

    int ji = 0;
    while(ai > demands[ji].tRight){
        ji++;
    }

    double Ai = demands[ji].arrivalLeft + demands[ji].slope * (ai - demands[ji].tLeft);

    double u = r.getU();
    double Ai1 = Ai + Exponential(1.0 , u);

    double tCycle = 0;
    if(Ai1 >= cumulativeArrivals){
        int w = floor(double(Ai1 - Ai) / cumulativeArrivals);

        Ai1 -= w * cumulativeArrivals;
        tCycle += w * cumulativeArrivals;

        Ai1 -= cumulativeArrivals;
        tCycle += tK - ai;

        ji = 0;
    }

    int ji1 = ji;
    while(Ai1 > demands[ji1].arrivalRight){
        ji1++;
    }

    double ai1 = (double(Ai1 - demands[ji1].arrivalLeft) / demands[ji1].slope) + demands[ji1].tLeft;

    if(tCycle > 0){
        return previousArrival + tCycle + ai1;
    }
    else{
        return previousArrival + (ai1 - ai);
    }
}

double calculateDeliveryLag(double t, double previousTime, double q, double a, double b, double c, double u){
    // return previousTime + ((A + q) / M) + Pi + Triangular(a, b, c, u);
    return max(previousTime, t) + (double(387 + q) / 13.0) + 23.7 + Triangular(a, c, b, u); 
}

void runSim(Welford &w, RandomFile &r, double a, double b, double c, int S, int s, double start, double end, vector<Demand> demands){
    double t = 0.0;

    //initialize event list and back order queue
    priority_queue<Event, vector<Event> > eventList = priority_queue<Event, vector<Event> >();
    priority_queue<Order, vector<Order> > backOrders = priority_queue<Order, vector<Order> >();

    //schedule first arrival
    Event firstArrival;
    double ar = nextArrival(0.0, demands,r);
    firstArrival.type = "carDemand";
    firstArrival.at = ar;
    firstArrival.numberOfCars = 1;
    eventList.push(firstArrival);

    //schedule first inventory review
    Event firstInventoryReview;
    firstInventoryReview.type = "inventoryReview";
    firstInventoryReview.at = 60;
    firstInventoryReview.numberOfCars = 0;
    eventList.push(firstInventoryReview);


    while(!eventList.empty()){

        //get and pop imminent event
        Event currentEvent = eventList.top();
        eventList.pop();
        t = currentEvent.at;

        //initialize number of rumor mill for demand generated to be 0 
        int rumorMill = 0;

        if(currentEvent.type == "inventoryReview"){

            //Scheduling next inventory review in 60 working hours
            Event nextReview;
            nextReview.type = "inventoryReview";
            nextReview.at = t + 60;
            if(nextReview.at < end){
                eventList.push(nextReview);
            }

            //If inventory is below threshold schedule restock
            if(w.inventory + w.onOrder <= s){
                Event nextRestock;
                nextRestock.type = "inventoryRestock";
                nextRestock.numberOfCars = S - (w.inventory + w.onOrder);
                nextRestock.at = t + calculateDeliveryLag(t,w.previousDeliveryTime, nextRestock.numberOfCars, a, b, c, r.getU());

                w.onOrder += nextRestock.numberOfCars;
                eventList.push(nextRestock);
                if(t >= start){
                    w.orders++;
                }
            }
        }
        else if(currentEvent.type == "inventoryRestock"){
            // cout << "DEBUG: restock of " << currentEvent.numberOfCars << " at time of " << t << endl; 

            w.onOrder -= currentEvent.numberOfCars;
            w.inventory += currentEvent.numberOfCars;

            w.previousDeliveryTime = t; // update the previous delivery time to current time
            
            // Everytime we get new orders, cycle through the N new cars arriving, and calculate the penalty
            int numberOfCarsToGive = currentEvent.numberOfCars;
            while(!backOrders.empty() && numberOfCarsToGive > 0){
                Order currOrder = backOrders.top();
                backOrders.pop();
                
                if(t > start){
                    w.penalties += floor(double(t - currOrder.orderPlacementTime) / 10.0) * 100;
                }
                numberOfCarsToGive--;
            }
        }
        else if(currentEvent.type == "carDemand"){
            // cout << "DEBUG: demand at time " << t << endl;

            //schedule next demand
            Event nextDemand;
            nextDemand.type = "carDemand";
            nextDemand.at = nextArrival(t, demands, r);
            nextDemand.numberOfCars = 1;
            eventList.push(nextDemand);


            if(w.inventory <= 0){
                // no cars in stock, increase order count
                Order order;
                order.numberOfCars = currentEvent.numberOfCars;
                order.orderPlacementTime = t;
                backOrders.push(order);

                w.inventory -= currentEvent.numberOfCars;
                w.backorderCount += currentEvent.numberOfCars;
                // increase the number of orders placed
                
                //rumormill
                while(Equilikely(0, 1, r.getU()) == 1){
                    double u = r.getU();
                    double h = Uniform(2, 9, u);

                    Event demand;
                    demand.at = t + h;
                    demand.type = "rumorMillCarDemand";
                    demand.numberOfCars = 1;
                    eventList.push(demand);
                    rumorMill++;
                }
            }
            else { 
                // cars in stock, decrease level - don't need to schedule any new orders
                w.inventory -= currentEvent.numberOfCars;
            }

           
        }
        else if(currentEvent.type == "rumorMillCarDemand"){

            if(w.inventory <= 0){ 
                // only order if there are no cars in stock
                Order order;
                order.numberOfCars = currentEvent.numberOfCars;
                order.orderPlacementTime = t;
                backOrders.push(order);

                w.inventory -= currentEvent.numberOfCars;
                w.backorderCount += currentEvent.numberOfCars;

                //rumormill
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
        //update statistics if in results window
        if(t >= start){
            if(w.inventory < w.minInventory){
                w.minInventory = w.inventory;
            }
            if(rumorMill > w.maxRumorMill){
                w.maxRumorMill = rumorMill;
            }
        }
        //if past results window print OUTPUT and exit
        if(t > end){
            // cout << "OUTPUT MININVENTORY " << w.minInventory << endl;
            // cout << "OUTPUT PENALTIES " << w.penalties << endl;
            // cout << "OUTPUT MAXRUMORMILL " << w.maxRumorMill << endl;
            // cout << "OUTPUT ORDERS " << w.orders << endl;
            return;
        }
    }
}

void runTriangle(RandomFile &r, double a, double b, double c){
    while(true){
        double u = r.getU();
        cout << "OUTPUT " << Triangular(a, c, b, u) << endl;
    }
}

void runArrival(RandomFile &r, vector<Demand> demands){
    double next = 0.0;
    while(true){
        next = nextArrival(next, demands, r);
        cout << "OUTPUT " << next << endl;
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
    vector<Demand> demands = vector<Demand>();
    // check that simulation has enough parameters - prevents segfault
    if(argc < 11){
        cerr << "Missing Required Parameters." << endl;
        ::exit(-1);
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

    w.inventory = S;
    w.onOrder = 0;

    if(runMode == "RESULTS"){
        for(int newS = 0; newS <= 99; newS++){
            Welford newW = Welford();
            newW.inventory = S;
            newW.onOrder = 0;
            for(int i = 0; i < 5; i++){
                runSim(newW, r, a, b, c, S, newS, start, end, demands);
            }
            cout << (newW.penalties / 5.0) << endl;
        }
    }
    else if(runMode == "TRIANGLE"){
        runTriangle(r, a, b, c);
    } else if (runMode == "ARRIVALS"){
        runArrival(r, demands);
    }

	return -1;
}