/* 
 * File:   Data.h
 * Author: klian
 * Comment: This file reads instance data from files. It is responsible for filling out
 * all instance information associated with an instance. This data file is essential for
 * all other functions to work. It should be included in other files.
 * 
 * Created on January 12, 2016, 10:07 AM
 */

#ifndef DATA_H
#define	DATA_H
#include "UserStruct.h"

const double myEPS = 1.0e-4;
const double myPrecision = 1.0e4;

extern size_t myInstID; //instance id
extern size_t myNumNode; //number of nodes, including depot and customers
extern size_t myNumCus; //number of customers
extern size_t myNumDays; //number of days in the planning horizon
extern size_t myNumVeh; //number of available vehicles
extern size_t myVehCap; //vehicle capacity
extern double myMaxDuration; //maximum vehicle duration
extern double myMaxTimeDiff; //maximum arrival time differential

extern double **myNodeCoord; //node coordinates of size myNumNode * 2
extern double **myCusSrvTime; //customer service time of size 1 * myNumNode
extern size_t **myCusDemand; //customer demand of size myNumNode * myNumNode
extern size_t **mySrvIndict; //service indicator: w_{id}
extern double **myNodeDistMat; //node distance matrix of size myNumNode * myNumNode
extern double *myCusDepotAngle;
extern size_t *myCusSrvFreq;
extern std::vector<size_t> myCusVec; //
extern std::vector<size_t> myCusVecAngleSorted; //customer list sorted by angle with the depot

void readSmallInst(size_t instsize, size_t instid);
void readNewInst(size_t instsize, size_t numdays, std::string dpttype, size_t instid);
void readNewTestInst(size_t instsize, size_t locid, size_t numdays, std::string dpttype, size_t instid);
void createOtherVars();
double roundVal(double x);
double getCusDepotAngle(double dx, double dy, double cx, double cy);

#endif	/* DATA_H */

