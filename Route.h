/* 
 * File:   Route.h
 * Author: klian
 *
 * Created on July 11, 2016, 11:20 AM
 */

#ifndef ROUTE_H
#define	ROUTE_H
#include "Data.h"

//define a vehicle route
class Route {    
public:
    size_t myNumStops;
    size_t myHead;
    size_t myTail;

    size_t myLoad;
    double myTravDist;
    double myDuration;

public:
    Route();
    ~Route();

    void resetState();
    double getPenalizedTD();
};


#endif	/* ROUTE_H */

