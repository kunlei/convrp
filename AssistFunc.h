/* 
 * File:   AssistFunc.h
 * Author: klian
 * Comment: this file contains facilitating functions that will be used by other files
 *
 * Created on July 11, 2016, 10:26 AM
 */

#ifndef ASSISTFUNC_H
#define	ASSISTFUNC_H
#include "Data.h"

//time recording functions
double get_wall_time();
double get_cpu_time();

//general printing functions
void print(std::vector<size_t> &myvec);
void print(std::vector<double> &myvec);
void print(std::vector<std::vector<int> > &mymat);
void print(std::vector<std::vector<size_t> > &mymat);
void print(std::vector<std::vector<double> > &mymat);
void print(std::vector<std::vector<std::vector<double> > > &mymat);
void print(std::vector<std::vector<std::vector<size_t> > > &mymat);
void print(std::vector<std::vector<std::vector<int> > > &mymat);

//instance printing functions
void printNodeCoord();
void printCusDemand();
void printSrvIndict();
void printCusSrvTime();
void printNodeDistMat();
void printCusDepotAngle();
void printCusVecAngleSorted();
void printCusSrvFreq();

//TSP-related functions
size_t getOrient(size_t orig, size_t p1, size_t p2);
double getNodeDist(size_t node1, size_t node2);
void getConvexHull(std::vector<size_t> &myvec, std::vector<size_t> &myhull);
void createTSPTourGreedy(std::vector<size_t> &myvec, std::string alg);
size_t getTotalServiceTime();

//create new test instances
double triangular(double a, double b, double c);
void createNewInstance(size_t numInst,
        size_t instsize, size_t numDays,
        std::string dpttype, std::string custype,
        double xlim, double ylim,
        double minDmd, double maxDmd, double dmdProb,
        double minSrvT, double maxSrvT);
void createNewInstance();
void createCusLocationFiles(size_t numrep, size_t numcus, std::string dpttype, std::string custype, double xlim, double ylim);
void createNewInstUseFixedLoc(size_t locid, size_t numinst, size_t numcus, size_t numdays,
        double mindmd, double maxdmd, double dmdprob,
        double minsrvt, double maxsrvt);
void verifyParaAgainstGoldenInst(size_t instsize, size_t instid);

//create run script
void createShellScriptForOne(size_t instsize, size_t numdays, std::string dpttype, size_t instid);
void createShellScriptForAll();

void getScriptForConVRPModelOneSmallInst(size_t instsize, size_t instid, size_t modelid);
void getScriptForConVRPModelAllSmallInst(size_t modelid);
void getScriptForConVRPModelOneLargeInst(size_t instsize, size_t numdays, std::string dpttype, size_t instid, size_t modelid);
void getScriptForConVRPModelAllLargeInst(size_t modelid);

#endif	/* ASSISTFUNC_H */

