#include "Route.h"

Route::Route() {

    myNumStops = 0;
    myHead = 0;
    myTail = 0;

    myLoad = 0;
    myTravDist = 0.0;
    myDuration = 0.0;
}

Route::~Route() {

}

double Route::getPenalizedTD() {

    double myPenalizedTD = myTravDist;

    if (myLoad > myVehCap) {
        myPenalizedTD += (myLoad - myVehCap) * 1.0;
    }

    if (myDuration > myMaxDuration) {
        myPenalizedTD += (myDuration - myMaxDuration) * 1.0;
    }

    return myPenalizedTD;
}

void Route::resetState() {

    myNumStops = 0;
    myHead = 0;
    myTail = 0;
    myLoad = 0;
    myTravDist = 0.0;
    myDuration = 0.0;
}

