#include <vector>
#include <limits>

#include "ConVRPSol.h"

using namespace std;

/**
 * Default constructor
 * Updated: 07/18/2016
 */
ConVRPSol::ConVRPSol() {

    initMemVars();
    createDailyRoutes();
    computeTravDist();
}

/**
 * Destructor
 * Updated: 07/18/2016
 */
ConVRPSol::~ConVRPSol() {

}

/**
 * Initialize member variables
 * Updated: 07/16/2016
 */
void ConVRPSol::initMemVars() {

    //allocate space for vehicle routes on each day
    for (size_t d = 0; d < myNumDays; d++) {
        vector<Customer> myvec;
        myvec.reserve(myNumNode);
        for (size_t i = 0; i < myNumNode; i++) {
            myvec.push_back(Customer(i, d));
        }
        myPlan.push_back(myvec);
    }

    //set indicator for each customer on each day
    for (size_t d = 0; d < myNumDays; d++) {
        for (size_t i = 0; i < myNumNode; i++) {
            if (myCusDemand[i][d] > 0) {
                myPlan[d][i].needVisit = true;
            } else {
                myPlan[d][i].needVisit = false;
            }
        }
    }

    //allocate space for vehicle summary information
    for (size_t d = 0; d < myNumDays; d++) {
        Route myRt;
        vector<Route> myvec(myNumVeh, myRt);
        myRoutes.push_back(myvec);
    }

    //initialize objective value
    myTravDist = 0.0;
}

/**
 * This function creates vehicle routes on each day
 * Updated: 07/16/2016
 */
void ConVRPSol::createDailyRoutes() {

    //start with a random customer and do the insertion
    size_t index = rand() % myNumCus;
    size_t iter = index;
    size_t cus; //the customer to be inserted
    vector<UIntDbl> myVehList;
    while (true) {
        cus = myCusVecAngleSorted[iter];

        //for each vehicle, check whether it is feasible to do the insertion and its associated cost
        myVehList.clear();
        bool hasEmptyVeh = false;
        for (size_t v = 0; v < myNumVeh; v++) {
            //check whether this vehicle is empty or not
            size_t sumStops = 0;
            for (size_t d = 0; d < myNumDays; d++) {
                sumStops += myRoutes[d][v].myNumStops;
            }
            if (sumStops < 1) {
                if (hasEmptyVeh == false) {
                    //insert this vehicle to the list
                    double distInc = myCusSrvFreq[cus] * 2 * myNodeDistMat[0][cus];
                    myVehList.push_back(UIntDbl(v, distInc));
                    hasEmptyVeh = true;
                }
            } else {
                bool isfeasible = true;
                //check vehicle capacity
                for (size_t d = 0; d < myNumDays; d++) {
                    if (myCusDemand[cus][d] > 0) {
                        if (myRoutes[d][v].myLoad + myCusDemand[cus][d] > myVehCap) {
                            isfeasible = false;
                            break;
                        }
                    }
                }

                //check route duration
                if (isfeasible == true) {
                    for (size_t d = 0; d < myNumDays; d++) {
                        if (myCusDemand[cus][d] > 0) {
                            //compute route duration increase after insertion
                            double duraInc = 0.0;
                            if (myRoutes[d][v].myNumStops < 1) {
                                duraInc = 2 * myNodeDistMat[0][cus] + myCusSrvTime[cus][d];
                            } else {
                                duraInc = myNodeDistMat[myRoutes[d][v].myTail][cus]
                                        + myCusSrvTime[cus][d]
                                        + myNodeDistMat[cus][0]
                                        - myNodeDistMat[myRoutes[d][v].myTail][0];
                            }
                            //check duration constraint
                            if (myRoutes[d][v].myDuration + duraInc > myMaxDuration) {
                                isfeasible = false;
                                break;
                            }
                        }
                    }
                }

                //check time consistency
                if (isfeasible == true && myCusSrvFreq[cus] > 1) {
                    double minarrvt = numeric_limits<double>::max();
                    double maxarrvt = -numeric_limits<double>::max();
                    for (size_t d = 0; d < myNumDays; d++) {
                        if (myCusDemand[cus][d] > 0) {
                            //compute the arrival time at the customer
                            double arrvt = 0.0;
                            if (myRoutes[d][v].myNumStops < 1) {
                                arrvt = myNodeDistMat[0][cus];
                            } else {
                                myPlan[d][myRoutes[d][v].myTail].getLeaveTime() + myNodeDistMat[myRoutes[d][v].myTail][cus];
                            }
                            if (arrvt < minarrvt) {
                                minarrvt = arrvt;
                            }
                            if (arrvt > maxarrvt) {
                                maxarrvt = arrvt;
                            }
                        }
                    }
                    if (maxarrvt > 0) {
                        if (maxarrvt - minarrvt > myMaxTimeDiff) {
                            isfeasible = false;
                        }
                    }
                }

                //compute insertion cost
                if (isfeasible == true) {
                    double distInc = 0.0;
                    for (size_t d = 0; d < myNumDays; d++) {
                        if (myCusDemand[cus][d] > 0) {
                            if (myRoutes[d][v].myNumStops < 1) {
                                distInc += 2 * myNodeDistMat[0][cus];
                            } else {
                                distInc += myNodeDistMat[myRoutes[d][v].myTail][cus]
                                        + myNodeDistMat[cus][0]
                                        - myNodeDistMat[myRoutes[d][v].myTail][0];
                            }
                        }
                    }
                    myVehList.push_back(UIntDbl(v, distInc));
                }
            }
        }//end for each vehicle computing the insertion cost

        //make sure the feasible insertion vehicle list is not empty
        if (myVehList.empty() == true) {
            cout << "ConVRPSol::createDailyRoutes(): the vehicle insertion list is empty!" << endl;
            exit(-1);
        }

        //sort vehicle list
        sort(myVehList.begin(), myVehList.end(), UIntDbl::sortByDblInc);

        //do insertion
        size_t veh = myVehList.front().myUInt;
        for (size_t d = 0; d < myNumDays; d++) {
            if (myCusDemand[cus][d] > 0) {
                addCusToRouteTail(d, veh, cus);
            }
        }

        //continue to the next customer
        iter++;
        if (iter == myNumCus) {
            iter = 0;
        }
        if (iter == index) {
            break;
        }
    }//end insertion of each customer
}

/**
 * Check whether a solution is feasible or not regarding capacity, duration and
 * time consistency constraints
 * @return true if feasible, false otherwise
 * Updated: 07/18/2016
 */
bool ConVRPSol::isFeasible() {

    bool isfeasible = true;
    //check capacity constraint
    for (size_t d = 0; d < myNumDays; d++) {
        for (size_t v = 0; v < myNumVeh; v++) {
            if (myRoutes[d][v].myLoad > myVehCap) {
                isfeasible = false;
                goto CapCheckStop;
            }
        }
    }

CapCheckStop:
    for (size_t d = 0; d < myNumDays; d++) {
        for (size_t v = 0; v < myNumVeh; v++) {
            if (myRoutes[d][v].myDuration > myMaxDuration) {
                isfeasible = false;
                goto DuraCheckStop;
            }
        }
    }

DuraCheckStop:
    for (size_t i = 1; i < myNumNode; i++) {
        if (myCusSrvFreq[i] > 1) {
            double minarrvt = numeric_limits<double>::max();
            double maxarrvt = -numeric_limits<double>::max();
            for (size_t d = 0; d < myNumDays; d++) {
                if (myCusDemand[i][d] > 0) {
                    if (myPlan[d][i].myArrvTime < minarrvt) {
                        minarrvt = myPlan[d][i].myArrvTime;
                    }
                    if (myPlan[d][i].myArrvTime > maxarrvt) {
                        maxarrvt = myPlan[d][i].myArrvTime;
                    }
                }
            }
            if (maxarrvt > 0) {
                double diff = maxarrvt - minarrvt;
                if (diff > myMaxTimeDiff) {
                    isfeasible = false;
                    break;
                }
            }
        }
    }

    return isfeasible;
}

/**
 * Attach a customer to the tail of a specified route
 * @param d: the day on which the insertion happens
 * @param v: the vehicle in which the customer is inserted
 * @param cus: the customer to be inserted
 * Updated: 07/18/2016
 */
void ConVRPSol::addCusToRouteTail(size_t d, size_t v, size_t cus) {

    if (myPlan[d][cus].needVisit == false) {
        cout << "ConVRPSol::addCusToRouteTail(): this customer does not need visit!" << endl;
        return;
    }

    if (myRoutes[d][v].myNumStops < 1) {
        //in this case, there is no customer existing on the route
        //set the starting customer
        myRoutes[d][v].myNumStops = 1;
        myRoutes[d][v].myHead = cus;
        myRoutes[d][v].myTail = cus;

        //update route
        myPlan[d][cus].isVisited = true;
        myPlan[d][cus].myDriver = v;
        myPlan[d][cus].myPrevCus = 0;
        myPlan[d][cus].myNextCus = 0;
        myPlan[d][cus].myArrvTime = myNodeDistMat[0][cus];

        myRoutes[d][v].myLoad = myCusDemand[cus][d];
        myRoutes[d][v].myTravDist = 2 * myNodeDistMat[0][cus];
        myRoutes[d][v].myDuration = myRoutes[d][v].myTravDist + myCusSrvTime[cus][d];
    } else {
        //in this case, there are customers existing on the route
        //update route
        myPlan[d][cus].isVisited = true;
        myPlan[d][cus].myDriver = v;
        myPlan[d][cus].myPrevCus = myRoutes[d][v].myTail;
        myPlan[d][cus].myNextCus = 0;
        myPlan[d][cus].myArrvTime = myPlan[d][myRoutes[d][v].myTail].getLeaveTime()
                + myNodeDistMat[myRoutes[d][v].myTail][cus];
        myPlan[d][myRoutes[d][v].myTail].myNextCus = cus;

        myRoutes[d][v].myNumStops++;
        myRoutes[d][v].myTail = cus;
        myRoutes[d][v].myLoad += myCusDemand[cus][d];
        double distinc = myNodeDistMat[myPlan[d][cus].myPrevCus][cus]
                + myNodeDistMat[cus][0]
                - myNodeDistMat[myPlan[d][cus].myPrevCus][0];
        myRoutes[d][v].myTravDist += distinc;
        myRoutes[d][v].myDuration += distinc + myCusSrvTime[cus][d];
    }
}

/**
 * This function summarizes the travel distance of all vehicles
 * Updated: 07/18/2016
 */
void ConVRPSol::computeTravDist() {

    myTravDist = 0.0;
    for (size_t d = 0; d < myNumDays; d++) {
        for (size_t v = 0; v < myNumVeh; v++) {
            myTravDist += myRoutes[d][v].myTravDist;
        }
    }
}

/**
 * This function prints the solution representation
 * Updated: 07/20/2016
 */
void ConVRPSol::printRep() {
}

void ConVRPSol::printSol() {

}