/* 
 * File:   ConVRPSol.h
 * Author: klian
 *
 * Created on July 16, 2016, 4:23 PM
 */

#ifndef CONVRPSOL_H
#define	CONVRPSOL_H
#include "Data.h"
#include "AssistFunc.h"
#include "Customer.h"
#include "Route.h"

class ConVRPSol {
public:
    std::vector<std::vector<Customer> > myPlan;
    std::vector<std::vector<Route> > myRoutes;
    
    double myTravDist;
    
public:
    ConVRPSol();
    ~ConVRPSol();
    
    bool isFeasible();
    
    void printRep();
    void printSol();
    
    void computeTravDist();
    
    void addCusToRouteTail(size_t d, size_t v, size_t cus);
    
private:
    void initMemVars();
    void createDailyRoutes();
};


#endif	/* CONVRPSOL_H */

