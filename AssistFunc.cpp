#include "AssistFunc.h"

using namespace std;

/**
 * get the current wall time
 * @return
 */
double get_wall_time() {

    struct timeval time;
    if (gettimeofday(&time, NULL)) {
        //  Handle error
        return 0;
    }
    return (double) time.tv_sec + (double) time.tv_usec * .000001;
}

/**
 * get current CPU time
 * @return
 */
double get_cpu_time() {

    return (double) clock() / CLOCKS_PER_SEC;
}

/**
 * print a vector of size_t numbers
 * @param myvec
 */
void print(vector<size_t> &myvec) {

    if (myvec.empty()) return;
    for (size_t i = 0; i < myvec.size(); i++) {
        cout << setw(4) << myvec[i];
    }
    cout << endl;
}

/**
 * print a vector of double numbers
 * @param myvec
 */
void print(vector<double> &myvec) {

    if (myvec.empty()) return;
    for (size_t i = 0; i < myvec.size(); i++) {
        cout << setw(12) << myvec[i];
    }
    cout << endl;
}

/**
 * print a matrix of integer values
 * @param mymat
 */
void print(std::vector<std::vector<int> >& mymat) {

    for (size_t i = 0; i < mymat.size(); i++) {
        for (size_t j = 0; j < mymat[i].size(); j++) {
            cout << setw(8) << mymat[i][j];
        }
        cout << endl;
    }
}

/**
 * print a matrix of double values
 * @param mymat
 */
void print(std::vector<std::vector<double> >& mymat) {

    for (size_t i = 0; i < mymat.size(); i++) {
        for (size_t j = 0; j < mymat[i].size(); j++) {
            cout << setw(12) << mymat[i][j];
        }
        cout << endl;
    }
}

/**
 * print a matrix of size_t values
 * @param mymat
 */
void print(vector<vector<size_t> > &mymat) {

    for (size_t i = 0; i < mymat.size(); i++) {
        for (size_t j = 0; j < mymat[i].size(); j++) {
            cout << setw(8) << mymat[i][j];
        }
        cout << endl;
    }
}

/**
 * print a matrix of double values
 * @param mymat
 */
void print(vector<vector<vector<double> > > &mymat) {

    for (size_t d = 0; d < mymat.size(); d++) {
        cout << "d = " << d << endl;
        for (size_t i = 0; i < mymat[d].size(); i++) {
            for (size_t j = 0; j < mymat[d][i].size(); j++) {
                cout << setw(8) << mymat[d][i][j];
            }
            cout << endl;
        }
        cout << endl;
    }
}

/**
 * print a matrix of size_t values
 * @param mymat
 */
void print(vector<vector<vector<size_t> > > &mymat) {

    for (size_t d = 0; d < mymat.size(); d++) {
        cout << "Day = " << d << endl;
        cout << setw(8) << "0";
        for (size_t i = 1; i < mymat[d].size(); i++) {
            cout << setw(4) << i;
        }
        cout << endl;
        for (size_t i = 0; i < mymat[d].size(); i++) {
            cout << setw(4) << i;
            for (size_t j = 0; j < mymat[d][i].size(); j++) {
                cout << setw(4) << mymat[d][i][j];
            }
            cout << endl;
        }
    }
}

/**
 * print a matrix of integer values
 * @param mymat
 */
void print(vector<vector<vector<int> > > &mymat) {

    for (size_t d = 0; d < mymat.size(); d++) {
        cout << "Day = " << d << endl;
        cout << setw(8) << "0";
        for (size_t i = 1; i < mymat[d].size(); i++) {
            cout << setw(4) << i;
        }
        cout << endl;
        for (size_t i = 0; i < mymat[d].size(); i++) {
            cout << setw(4) << i;
            for (size_t j = 0; j < mymat[d][i].size(); j++) {
                cout << setw(4) << mymat[d][i][j];
            }
            cout << endl;
        }
    }
}

/**
 * Print node distance matrix
 * Updated: 04/20/2016
 */
void printNodeDistMat() {

    cout << setw(12) << "0";
    for (size_t i = 1; i < myNumNode; i++) {
        cout << setw(8) << i;
    }
    cout << endl;
    for (size_t i = 0; i < myNumNode; i++) {
        cout << setw(4) << i << setw(8);
        for (size_t j = 0; j < myNumNode; j++) {
            cout << myNodeDistMat[i][j] << setw(8);
        }
        cout << endl << endl;
    }
    cout << endl;
}

/**
 * Print node coordinates
 * Updated: 04/20/2016
 */
void printNodeCoord() {

    for (size_t i = 0; i < myNumNode; i++) {
        cout << setw(2) << i << ":" << setw(8) << myNodeCoord[i][0]
                << setw(8) << myNodeCoord[i][1]
                << endl;
    }
}

/**
 * Print customer demand
 * Updated: 04/20/2016
 */
void printCusDemand() {

    for (size_t i = 0; i < myNumNode; i++) {
        cout << setw(2) << i << ":";
        for (size_t j = 0; j < myNumDays; j++) {
            cout << setw(12) << myCusDemand[i][j];
        }
        cout << endl;
    }
}

/**
 * This function prints service indicator
 * Updated: 05/08/2016
 */
void printSrvIndict() {

    for (size_t i = 0; i < myNumNode; i++) {
        cout << setw(2) << i << ":";
        for (size_t j = 0; j < myNumDays; j++) {
            cout << setw(12) << mySrvIndict[i][j];
        }
        cout << endl;
    }
}

/**
 * Print customer service time
 * Updated: 04/20/2016
 */
void printCusSrvTime() {

    for (size_t i = 0; i < myNumNode; i++) {
        cout << setw(2) << i << ":";
        for (size_t j = 0; j < myNumDays; j++) {
            cout << setw(8) << myCusSrvTime[i][j];
        }
        cout << endl;
    }
}

/**
 * Print customer depot angle
 * Updated: 04/20/2016
 */
void printCusDepotAngle() {

    cout << "Print customer depot angle: " << endl;
    for (size_t i = 0; i < myNumNode; i++) {
        cout << setw(12) << myCusDepotAngle[i];
    }
    cout << endl;
}

/**
 * print the list of all customers sorted by the depot they make with the depot
 */
void printCusVecAngleSorted() {

    cout << "Print customer vector angle sorted: " << endl;
    for (size_t i = 0; i < myCusVecAngleSorted.size(); i++) {
        cout << setw(5) << myCusVecAngleSorted[i];
    }
    cout << endl;
}

/**
 * print customer service frequency
 */
void printCusSrvFreq() {

    cout << "Print customer service frequency: " << endl;
    cout << setw(8) << "CusID:";
    for (size_t i = 0; i < myNumNode; i++) {
        cout << setw(4) << i;
    }
    cout << endl;
    cout << setw(8) << "SrvFreq:";
    for (size_t i = 0; i < myNumNode; i++) {
        cout << setw(4) << myCusSrvFreq[i];
    }
    cout << endl;
}

/**
 * get distance between two nodes
 * @param node1
 * @param node2
 * @return
 */
double getNodeDist(size_t node1, size_t node2) {

    return myNodeDistMat[node1][node2];
}

/**
 * facilitating function for getConvexHull()
 * @param orig
 * @param p1
 * @param p2
 * @return
 */
size_t getOrient(size_t orig, size_t p1, size_t p2) {

    double prod = (myNodeCoord[p1][0] - myNodeCoord[orig][0]) * (myNodeCoord[p2][1] - myNodeCoord[orig][1])
            - (myNodeCoord[p2][0] - myNodeCoord[orig][0]) * (myNodeCoord[p1][1] - myNodeCoord[orig][1]);
    if (prod < myEPS) {
        return 0;
    } else {
        return (prod > 0) ? 1 : 2;
    }
}

/**
 * find convex hull of a given set of customers
 * @param myvec: the input set of customers, note: non-convex hull customers are also saved in this vector
 * @param myhull: saves convex hull points
 */
void getConvexHull(vector<size_t>& myvec, vector<size_t>& myhull) {

    if (myvec.size() > 3) {
        size_t node;
        //determine the bottom-most point
        size_t p0 = myvec.front();
        double x0 = myNodeCoord[p0][0];
        double y0 = myNodeCoord[p0][1];
        vector<size_t>::iterator it = myvec.begin();
        ++it;
        for (; it != myvec.end(); it++) {
            node = *it;
            if (myNodeCoord[node][1] < y0) {
                p0 = node;
                x0 = myNodeCoord[node][0];
                y0 = myNodeCoord[node][1];
            } else if (myNodeCoord[node][1] == y0) {
                if (myNodeCoord[node][0] < x0) {
                    p0 = node;
                    x0 = myNodeCoord[node][0];
                    y0 = myNodeCoord[node][1];
                }
            }
        }

        //sort the remaining points according to the angle they made with p0
        for (it = myvec.begin(); it != myvec.end(); it++) {
            if (*it == p0) {
                break;
            }
        }
        myvec.erase(it); //delete the bottom-most point
        vector<NPoint> myPointVec;
        for (it = myvec.begin(); it != myvec.end(); it++) {
            node = *it;
            myPointVec.push_back(NPoint(node, myNodeCoord[node][0] - x0, myNodeCoord[node][1] - y0));
        }
        myvec.clear();
        sort(myPointVec);

        myhull.push_back(p0);
        for (size_t i = 0; i < 2; i++) {
            myhull.push_back(myPointVec.at(i).myCusID);
        }
        size_t top, nexttotop;
        size_t orient;
        size_t vsize = myPointVec.size();
        for (size_t i = 2; i < vsize; i++) {
            if (myPointVec.at(i).flag == true) {
                top = myhull.back();
                nexttotop = myhull.at(myhull.size() - 2);
                orient = getOrient(nexttotop, top, myPointVec.at(i).myCusID);
                while (orient != 1) {
                    myvec.push_back(myhull.back());
                    myhull.pop_back();
                    if (myhull.size() <= 1) {
                        break;
                    } else {
                        top = myhull.back();
                        nexttotop = myhull.at(myhull.size() - 2);
                        orient = getOrient(nexttotop, top, myPointVec.at(i).myCusID);
                    }
                }
                myhull.push_back(myPointVec.at(i).myCusID);
            } else {
                myvec.push_back(myPointVec.at(i).myCusID);
            }
        }
    }
}

/**
 * create a TSP tour based on a given set of customers and a specified algorithm
 * @param myvec
 * @param alg
 */
void createTSPTourGreedy(std::vector<size_t>& myvec, std::string alg) {

    if (myvec.size() <= 3) return;

    if (alg == "FarthestInsertionMaxMin" ||
            alg == "ConvexHullFarthestInsertionMaxMin") {
        //farthest insertion heuristic for tsp
        list<size_t> mylist; //stores inserted nodes
        list<UIntDbl> uninscuslist; //nodes to be inserted
        if (alg == "FarthestInsertionMaxMin") {
            //FI heuristic starting with random three nodes
            for (size_t i = 0; i < 3; i++) {
                mylist.push_back(myvec.back());
                myvec.pop_back();
            }
            for (vector<size_t>::iterator it = myvec.begin(); it != myvec.end(); it++) {
                uninscuslist.push_back(UIntDbl(*it, 0.0));
            }
            myvec.clear();
        } else if (alg == "ConvexHullFarthestInsertionMaxMin") {
            //FI heuristic starting with convex hull
            vector<size_t> myhull;
            getConvexHull(myvec, myhull);
            for (vector<size_t>::iterator it = myhull.begin(); it != myhull.end(); it++) {
                mylist.push_back(*it);
            }
            for (vector<size_t>::iterator it = myvec.begin(); it != myvec.end(); it++) {
                uninscuslist.push_back(UIntDbl(*it, 0.0));
            }
            myvec.clear();
        }

        size_t cus;
        double dist, distinc;
        list<size_t>::iterator it, bestit, prev, last;
        list<UIntDbl>::iterator uit, bestuit;
        //find the minimal distance of inserted nodes to un-inserted nodes
        for (uit = uninscuslist.begin(); uit != uninscuslist.end(); uit++) {
            it = mylist.begin();
            (*uit).myDbl = getNodeDist((*uit).myUInt, *it);
            ++it;
            for (; it != mylist.end(); it++) {
                dist = getNodeDist((*uit).myUInt, *it);
                if (dist < (*uit).myDbl) {
                    (*uit).myDbl = dist;
                }
            }
        }
        //insert remaining nodes
        while (uninscuslist.empty() == false) {
            //find the node with maximal distance value
            uit = uninscuslist.begin();
            bestuit = uit;
            ++uit;
            for (; uit != uninscuslist.end(); uit++) {
                if ((*uit).myDbl > (*bestuit).myDbl) {
                    bestuit = uit;
                }
            }
            //insert the chosen node into mylist
            cus = (*bestuit).myUInt;
            uninscuslist.erase(bestuit);
            last = mylist.end();
            distinc = getNodeDist(mylist.back(), cus)
                    + getNodeDist(cus, mylist.front())
                    - getNodeDist(mylist.front(), mylist.back());
            bestit = last;
            it = mylist.begin();
            ++it;
            for (; it != last; it++) {
                prev = it;
                --prev;
                dist = getNodeDist(*prev, cus)
                        + getNodeDist(cus, *it)
                        - getNodeDist(*prev, *it);
                if (dist < distinc) {
                    bestit = it;
                    distinc = dist;
                }
            }
            mylist.insert(bestit, cus);
            //update distance after insertion
            for (uit = uninscuslist.begin(); uit != uninscuslist.end(); uit++) {
                dist = getNodeDist((*uit).myUInt, cus);
                if (dist < (*uit).myDbl) {
                    (*uit).myDbl = dist;
                }
            }
        }
        for (it = mylist.begin(); it != mylist.end(); it++) {
            myvec.push_back(*it);
        }
    }
}

/**
 * get total service time required by all customers in an instance
 * @return
 * Updated: 07/11/2016
 */
size_t getTotalServiceTime() {

    size_t totalTime = 0;
    for (size_t i = 1; i < myNumNode; i++) {
        for (size_t d = 0; d < myNumDays; d++) {
            if (myCusDemand[i][d] > 0) {
                totalTime += 1;
            }
        }
    }

    return totalTime;
}

/**
 * Generate a triangular distribution number
 * Updated: 06/28/2016
 **/
double triangular(double a, double b, double c) {
    double U = (double) rand() / (double) RAND_MAX;
    double F = (c - a) / (b - a);
    if (U <= F)
        return a + sqrt(U * (b - a) * (c - a));
    else
        return b - sqrt((1 - U) * (b - a) * (b - c));
}

/**
 * This function creates new instances
 * @numInst: number of instances to be created for this parameter combination
 * @instsize: instance size
 * @num: number of instances with this size
 * @xlim: the x-coordinate range
 * @ylim: the y-coordinate range
 * @dpttype: depot location type: center, origin or random
 * @custypeL customer location type: random
 */
void createNewInstance(size_t numInst,
        size_t instsize, size_t numDays,
        string dpttype, string custype,
        double xlim, double ylim,
        double minDmd, double maxDmd, double dmdProb,
        double minSrvT, double maxSrvT) {

    for (size_t n = 0; n < numInst; n++) {
        size_t numNode = instsize + 1;
        size_t numCus = instsize;

        //create customer and depot locations
        vector<vector<int> > nodeCoord;
        for (size_t j = 0; j < numNode; j++) {
            vector<int> vec(2, 0);
            nodeCoord.push_back(vec);
        }
        if (dpttype == "center") {
            nodeCoord[0][0] = (int) (xlim / 2);
            nodeCoord[0][1] = (int) (ylim / 2);
        } else if (dpttype == "eccentric") {
            nodeCoord[0][0] = 0;
            nodeCoord[0][1] = 0;
        } else if (dpttype == "random") {
            nodeCoord[0][0] = (int) (xlim * (double) rand() / RAND_MAX);
            nodeCoord[0][1] = (int) (ylim * (double) rand() / RAND_MAX);
        } else {
            cout << "createNewInstance(): wrong depot location type!" << endl;
        }
        for (size_t j = 1; j < numNode; j++) {
            nodeCoord[j][0] = (int) (xlim * (double) rand() / RAND_MAX);
            nodeCoord[j][1] = (int) (ylim * (double) rand() / RAND_MAX);
            bool flag = false;
            for (size_t k = 0; k < j; k++) {
                if (nodeCoord[j][0] == nodeCoord[k][0] &&
                        nodeCoord[j][1] == nodeCoord[k][1]) {
                    flag = true;
                    break;
                }
            }
            while (flag == true) {
                nodeCoord[j][0] = (int) (xlim * (double) rand() / RAND_MAX);
                nodeCoord[j][1] = (int) (ylim * (double) rand() / RAND_MAX);
                flag = false;
                for (size_t k = 0; k < j; k++) {
                    if (nodeCoord[j][0] == nodeCoord[k][0] &&
                            nodeCoord[j][1] == nodeCoord[k][1]) {
                        flag = true;
                        break;
                    }
                }
            }
        }

        //compute node distance matrix
        vector<vector<int> > nodeDistMat;
        for (size_t i = 0; i < numNode; i++) {
            vector<int> vec(numNode, 0);
            nodeDistMat.push_back(vec);
        }
        for (size_t i = 0; i < numNode; i++) {
            for (size_t j = 0; j < numNode; j++) {
                nodeDistMat[i][j] = (int) (sqrt(pow(nodeCoord[j][0] - nodeCoord[i][0], 2.0) +
                        pow(nodeCoord[j][1] - nodeCoord[i][1], 2.0)));
            }
        }

        //create demand on each day
        vector<vector<int> > cusDemand;
        for (size_t i = 0; i < numCus; i++) {
            vector<int> vec(numDays, 0);
            size_t ct = 0;
            for (size_t d = 0; d < numDays; d++) {
                if ((double) rand() / RAND_MAX < dmdProb) {
                    ct++;
                    vec[d] = minDmd + (int) ((maxDmd - minDmd) * (double) rand() / RAND_MAX);
                } else {
                    vec[d] = -1;
                }
            }
            //make sure a customer requires service on at least one day
            if (ct < 1) {
                vec[rand() % numDays] = minDmd + (int) ((maxDmd - minDmd) * (double) rand() / RAND_MAX);
            }
            cusDemand.push_back(vec);
        }

        //create service time on each day
        vector<vector<int> > cusSrvTime;
        for (size_t i = 0; i < numCus; i++) {
            vector<int> vec(numDays, 0);
            cusSrvTime.push_back(vec);
        }
        for (size_t i = 0; i < numCus; i++) {
            for (size_t d = 0; d < numDays; d++) {
                if (cusDemand[i][d] > 0) {
                    cusSrvTime[i][d] = minSrvT + (int) ((maxSrvT - minSrvT) * (double) rand() / RAND_MAX);
                } else {
                    cusSrvTime[i][d] = -1;
                }
            }
        }

        //determine vehicle capacity
        size_t vehCap = 0;
        double aval = 5.0;
        double bval = 0.70 * numCus;
        double cval = 0.50 * numCus;
        double r = triangular(aval, bval, cval);
        int numVisit = 0, totalDemand = 0;
        for (size_t i = 0; i < numCus; i++) {
            for (size_t d = 0; d < numDays; d++) {
                if (cusDemand[i][d] > 0) {
                    totalDemand += cusDemand[i][d];
                    numVisit++;
                }
            }
        }
        vehCap = (size_t) (r * ((double) totalDemand / numVisit));

        //determine maximum route duration
        size_t vehDura = 0;
        int numArcs = 0, totalDist = 0;
        for (size_t i = 0; i < numNode; i++) {
            for (size_t j = 0; j < numNode; j++) {
                if (nodeDistMat[i][j] > 0) {
                    numArcs++;
                    totalDist += nodeDistMat[i][j];
                }
            }
        }
        double avgArcDist = (double) totalDist / (double) numArcs;
        int totalTime = 0;
        for (size_t i = 0; i < numCus; i++) {
            for (size_t d = 0; d < numDays; d++) {
                if (cusDemand[i][d] > 0) {
                    totalTime += cusSrvTime[i][d];
                }
            }
        }
        double avgSrvTime = (double) totalTime / (double) numVisit;
        aval = 0.40 * numCus;
        bval = 0.70 * numCus;
        cval = 0.60 * numCus;
        r = triangular(aval, bval, cval);
        vehDura = (size_t) (r * (avgArcDist + avgSrvTime));

        //create instance filename
        string filename = "convrp-n";
        stringstream convert;
        convert << numCus;
        filename.append(convert.str());
        convert.str("");
        filename.append("-d");
        convert << numDays;
        filename.append(convert.str());
        convert.str("");
        filename.append("-");
        if (dpttype == "center") {
            filename.append("c");
        } else if (dpttype == "eccentric") {
            filename.append("e");
        } else if (dpttype == "random") {
            filename.append("r");
        } else {
            cout << "createNewInstance(): wrong tyep!" << endl;
        }
        convert << n + 1;
        filename.append(convert.str());
        filename.append(".vrp");
        cout << "filename = " << filename << endl;

        //write instance data to file
        ofstream ofile(filename.c_str(), ios::out);
        ofile << "NAME: " << filename << endl;
        ofile << "TYPE: " << "ConVRP" << endl;
        ofile << "DIMENSION: " << numNode << endl;
        ofile << "NUM_DAYS: " << numDays << endl;
        ofile << "CAPACITY: " << vehCap << endl;
        ofile << "DISTANCE: " << vehDura << endl;
        ofile << "EDGE_WEIGHT_TYPE: " << "EXPLICIT" << endl;
        ofile << "EDGE_WEIGHT_FORMAT: " << "EUC_2D" << endl;
        ofile << "NODE_COORD_TYPE: " << "TWOD_COORDS" << endl;
        ofile << "NODE_COORD_SECTION" << endl;
        for (size_t i = 1; i < numNode; i++) {
            ofile << i << setw(8) << nodeCoord[i][0]
                    << setw(8) << nodeCoord[i][1] << endl;
        }
        ofile << "EDGE_WEIGHT_SECTION" << endl;
        for (size_t i = 0; i < numNode; i++) {
            ofile << setw(5) << nodeDistMat[i][0];
            for (size_t j = 1; j < numNode; j++) {
                ofile << setw(12) << nodeDistMat[i][j];
            }
            ofile << endl;
        }
        ofile << "DEMAND_SECTION" << endl;
        for (size_t i = 0; i < numCus; i++) {
            ofile << i + 1;
            for (size_t d = 0; d < numDays; d++) {
                ofile << setw(8) << cusDemand[i][d];
            }
            ofile << endl;
        }
        ofile << "SRV_TIME_SECTION" << endl;
        for (size_t i = 0; i < numCus; i++) {
            ofile << i + 1;
            for (size_t d = 0; d < numDays; d++) {
                ofile << setw(8) << cusSrvTime[i][d];
            }
            ofile << endl;
        }
        ofile << "DEPOT_SECTION" << endl;
        ofile << nodeCoord[0][0]
                << setw(8) << nodeCoord[0][1] << endl;
        ofile << "-1" << endl;
        ofile << "EOF";

        ofile.close();
    }
}

void createNewInstance() {

    //size = 10
    createNewInstance(5, 10, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 10, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 10, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 12
    createNewInstance(5, 12, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 12, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 12, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 15
    createNewInstance(5, 15, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 15, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 15, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 18
    createNewInstance(5, 18, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 18, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 18, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 20
    createNewInstance(5, 20, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 20, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 20, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 22
    createNewInstance(5, 22, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 22, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 22, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 25
    createNewInstance(5, 25, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 25, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 25, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 28
    createNewInstance(5, 28, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 28, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 28, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 30
    createNewInstance(5, 30, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 30, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 30, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 32
    createNewInstance(5, 32, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 32, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 32, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 35
    createNewInstance(5, 35, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 35, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 35, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 38
    createNewInstance(5, 38, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 38, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 38, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 40
    createNewInstance(5, 40, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 40, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 40, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 42
    createNewInstance(5, 42, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 42, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 42, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 45
    createNewInstance(5, 45, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 45, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 45, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 48
    createNewInstance(5, 48, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 48, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 48, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);

    //size = 50
    createNewInstance(5, 50, 3, "center", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 50, 3, "eccentric", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
    createNewInstance(5, 50, 3, "random", "random", 1000.0, 1000.0, 1.0, 10.0, 0.70, 5.0, 10.0);
}

/**
 * generate customer location files
 * @param numrep: number of files to be generated
 * @param numcus: number of customers in the instance
 * @param xlim: x-coordinate limit
 * @param ylim: y-coordinate limit
 * Updated: 07/11/2016
 */
void createCusLocationFiles(size_t numrep, size_t numcus,
        string dpttype, string custype,
        double xlim, double ylim) {

    size_t numNode = numcus + 1;
    for (size_t f = 0; f < numrep; f++) {
        vector<vector<int> > nodeCoord;
        for (size_t i = 0; i < numNode; i++) {
            vector<int> vec(2, 0);
            nodeCoord.push_back(vec);
        }

        if (dpttype == "c") {
            nodeCoord[0][0] = (int) (xlim / 2);
            nodeCoord[0][1] = (int) (ylim / 2);
        } else if (dpttype == "e") {
            nodeCoord[0][0] = 0;
            nodeCoord[0][1] = 0;
        } else if (dpttype == "r") {
            nodeCoord[0][0] = (int) (xlim * (double) rand() / RAND_MAX);
            nodeCoord[0][1] = (int) (ylim * (double) rand() / RAND_MAX);
        } else {
            cout << "createCusLocationFiles(): wrong depot type parameter!" << endl;
        }

        for (size_t i = 1; i < numNode; i++) {
            nodeCoord[i][0] = (int) (xlim * (double) rand() / RAND_MAX);
            nodeCoord[i][1] = (int) (ylim * (double) rand() / RAND_MAX);
            //make sure no two points share the same location
            bool flag = false;
            for (size_t j = 0; j < i; j++) {
                if (nodeCoord[j][0] == nodeCoord[i][0] &&
                        nodeCoord[j][1] == nodeCoord[i][1]) {
                    flag = true;
                    break;
                }
            }
            while (flag == true) {
                nodeCoord[i][0] = (int) (xlim * (double) rand() / RAND_MAX);
                nodeCoord[i][1] = (int) (ylim * (double) rand() / RAND_MAX);
                flag = false;
                for (size_t j = 0; j < i; j++) {
                    if (nodeCoord[j][0] == nodeCoord[i][0] &&
                            nodeCoord[j][1] == nodeCoord[i][1]) {
                        flag = true;
                        break;
                    }
                }
            }
        }//end generating location for each node

        //save node location to file
        string filename = "CusLoc-n";
        stringstream convert;
        convert << numcus;
        filename.append(convert.str());
        filename.append("-");
        convert.str("");
        convert << f + 1;
        filename.append(convert.str());
        filename.append(".txt");

        fstream pfile(filename.c_str(), ios::out);
        if (!pfile) {
            cout << "createCusLocationFiles(): cannot open file for writing customer location!" << endl;
            exit(-1);
        } else {
            for (size_t i = 0; i < numNode; i++) {
                pfile << setw(4) << nodeCoord[i][0]
                        << setw(4) << nodeCoord[i][1] << endl;
            }
        }
        pfile.close();
    }
}

/**
 * create new instances using input customer location
 * @param locid: input file id
 * @param numinst: number of instances to be generated
 * @param numdays: number of days in the planning horizon
 * @param mindmd: minimum demand
 * @param maxdmd: maximum demand
 * @param dmdprob: demand probability
 * @param minsrvt: minimum service time
 * @param maxsrvt: maximum service time
 * Updated: 07/11/2016
 */
void createNewInstUseFixedLoc(size_t locid, size_t numinst, size_t numcus, size_t numdays,
        double mindmd, double maxdmd, double dmdprob,
        double minsrvt, double maxsrvt) {

    //create filename of customer location file
    string cuslocfilename = "Instance/CusLocFile/CusLoc-n";
    stringstream convert;
    convert << numcus;
    cuslocfilename.append(convert.str());
    cuslocfilename.append("-");
    convert.str("");
    convert << locid;
    cuslocfilename.append(convert.str());
    cuslocfilename.append(".txt");

    size_t numNode = numcus + 1;
    vector<vector<int> > nodeCoord;
    for (size_t i = 0; i < numNode; i++) {
        vector<int> vec(2, 0);
        nodeCoord.push_back(vec);
    }
    fstream pfile(cuslocfilename.c_str(), ios::in);
    if (!pfile) {
        cout << "createNewInstUseFixedLoc(): cannot open file for reading customer locations" << endl;
        exit(-1);
    } else {
        for (size_t i = 0; i < numNode; i++) {
            pfile >> nodeCoord[i][0] >> nodeCoord[i][1];
        }
    }
    //    ::print(nodeCoord);

    //compute node distance matrix
    vector<vector<int> > nodeDistMat;
    for (size_t i = 0; i < numNode; i++) {
        vector<int> vec(numNode, 0);
        nodeDistMat.push_back(vec);
    }
    for (size_t i = 0; i < numNode; i++) {
        for (size_t j = 0; j < numNode; j++) {
            nodeDistMat[i][j] = (int) (sqrt(pow(nodeCoord[j][0] - nodeCoord[i][0], 2.0) +
                    pow(nodeCoord[j][1] - nodeCoord[i][1], 2.0)));
        }
    }

    //create instance files
    for (size_t n = 0; n < numinst; n++) {
        //create demand on each day
        vector<vector<int> > cusDemand;
        for (size_t i = 0; i < numcus; i++) {
            vector<int> vec(numdays, 0);
            size_t ct = 0;
            for (size_t d = 0; d < numdays; d++) {
                if ((double) rand() / RAND_MAX < dmdprob) {
                    ct++;
                    vec[d] = mindmd + (int) ((maxdmd - mindmd) * (double) rand() / RAND_MAX);
                } else {
                    vec[d] = -1;
                }
            }
            //make sure a customer needs service on at least one day
            if (ct < 1) {
                vec[rand() % numdays] = mindmd + (int) ((maxdmd - mindmd) * (double) rand() / RAND_MAX);
            }
            cusDemand.push_back(vec);
        }

        //create service time on each day
        vector<vector<int> > cusSrvTime;
        for (size_t i = 0; i < numcus; i++) {
            vector<int> vec(numdays, 0);
            cusSrvTime.push_back(vec);
        }
        for (size_t i = 0; i < numcus; i++) {
            for (size_t d = 0; d < numdays; d++) {
                if (cusDemand[i][d] > 0) {
                    cusSrvTime[i][d] = minsrvt + (int) ((maxsrvt - minsrvt) * (double) rand() / RAND_MAX);
                } else {
                    cusSrvTime[i][d] = -1;
                }
            }
        }

        //determine vehicle capacity
        size_t vehCap = 0;
        double val = 0.5;
        double aval = val * numcus;
        double bval = (val + 0.2) * numcus;
        double cval = (val + 0.1) * numcus;
        double r = triangular(aval, bval, cval);
        int numVisit = 0, totalDemand = 0;
        for (size_t i = 0; i < numcus; i++) {
            for (size_t d = 0; d < numdays; d++) {
                if (cusDemand[i][d] > 0) {
                    totalDemand += cusDemand[i][d];
                    numVisit++;
                }
            }
        }
        vehCap = (size_t) (r * ((double) totalDemand / numVisit));
        cout << "vehCap = " << vehCap << endl;

        //determine maximum route duration
        size_t vehDura = 0;
        int numArcs = 0, totalDist = 0;
        for (size_t i = 0; i < numNode; i++) {
            for (size_t j = 0; j < numNode; j++) {
                if (nodeDistMat[i][j] > 0) {
                    numArcs++;
                    totalDist += nodeDistMat[i][j];
                }
            }
        }
        double avgArcDist = (double) totalDist / (double) numArcs;
        //        cout << "avgArcDist = " << avgArcDist << endl;
        int totalTime = 0;
        for (size_t i = 0; i < numcus; i++) {
            for (size_t d = 0; d < numdays; d++) {
                if (cusDemand[i][d] > 0) {
                    totalTime += cusSrvTime[i][d];
                }
            }
        }
        double avgSrvTime = (double) totalTime / (double) numVisit;
        //        cout << "avgSrvTime = " << avgSrvTime << endl;
        val = 0.70;
        aval = val * numcus;
        bval = (val + 0.20) * numcus;
        cval = (val + 0.10) * numcus;
        r = triangular(aval, bval, cval);
        vehDura = (size_t) (r * (avgArcDist + avgSrvTime));
        cout << "vehDura = " << vehDura << endl;

        //determine maximum arrival time differential
        size_t maxTimeDiff = (size_t) (avgArcDist * 0.20);
        cout << "maxTimeDiff = " << maxTimeDiff << endl;

        //create instance filename
        string filename = "convrp-n";
        stringstream convert;
        convert << numcus;
        filename.append(convert.str());
        convert.str("");
        filename.append("-l");
        convert << locid;
        filename.append(convert.str());
        convert.str("");
        filename.append("-d");
        convert << numdays;
        filename.append(convert.str());
        convert.str("");
        filename.append("-");
        filename.append("c");
        convert << n + 1;
        filename.append(convert.str());
        filename.append(".vrp");
        //        cout << "filename = " << filename << endl;

        //write instance data to file
        ofstream ofile(filename.c_str(), ios::out);
        ofile << "NAME: " << filename << endl;
        ofile << "TYPE: " << "ConVRP" << endl;
        ofile << "DIMENSION: " << numNode << endl;
        ofile << "NUM_DAYS: " << numdays << endl;
        ofile << "CAPACITY: " << vehCap << endl;
        ofile << "DISTANCE: " << vehDura << endl;
        ofile << "TIMEDIFFLIMIT: " << maxTimeDiff << endl;
        ofile << "EDGE_WEIGHT_TYPE: " << "EXPLICIT" << endl;
        ofile << "EDGE_WEIGHT_FORMAT: " << "EUC_2D" << endl;
        ofile << "NODE_COORD_TYPE: " << "TWOD_COORDS" << endl;
        ofile << "NODE_COORD_SECTION" << endl;
        for (size_t i = 1; i < numNode; i++) {
            ofile << i << setw(8) << nodeCoord[i][0]
                    << setw(8) << nodeCoord[i][1] << endl;
        }
        ofile << "EDGE_WEIGHT_SECTION" << endl;
        for (size_t i = 0; i < numNode; i++) {
            ofile << setw(5) << nodeDistMat[i][0];
            for (size_t j = 1; j < numNode; j++) {
                ofile << setw(12) << nodeDistMat[i][j];
            }
            ofile << endl;
        }
        ofile << "DEMAND_SECTION" << endl;
        for (size_t i = 0; i < numcus; i++) {
            ofile << i + 1;
            for (size_t d = 0; d < numdays; d++) {
                ofile << setw(8) << cusDemand[i][d];
            }
            ofile << endl;
        }
        ofile << "SRV_TIME_SECTION" << endl;
        for (size_t i = 0; i < numcus; i++) {
            ofile << i + 1;
            for (size_t d = 0; d < numdays; d++) {
                ofile << setw(8) << cusSrvTime[i][d];
            }
            ofile << endl;
        }
        ofile << "DEPOT_SECTION" << endl;
        ofile << nodeCoord[0][0]
                << setw(8) << nodeCoord[0][1] << endl;
        ofile << "-1" << endl;
        ofile << "EOF";

        ofile.close();
    }//end for each instance
}

/**
 * This function reads in the Golden instances and computes the capacity, duration and
 * arrival time differential parameters using the procedure described in the new instance generation
 * Updated: 07/15/2016
 */
void verifyParaAgainstGoldenInst(size_t instsize, size_t instid) {

    //read in small instances
    readSmallInst(instsize, instid);
    //    ::printNodeCoord();
    //    ::printNodeDistMat();
    //    ::printCusDemand();

    //determine vehicle capacity
    double vehCap = 0.0;
    double aval = 0.6 * myNumCus;
    double bval = 0.8 * myNumCus;
    double cval = 0.7 * myNumCus;
    double r = triangular(aval, bval, cval);
    size_t numVisit = 0, totalDemand = 0;
    for (size_t i = 1; i < myNumNode; i++) {
        for (size_t d = 0; d < myNumDays; d++) {
            if (myCusDemand[i][d] > 0) {
                numVisit++;
                totalDemand += myCusDemand[i][d];
            }
        }
    }
    vehCap = r * (double) totalDemand / (double) numVisit;
    cout << "vehCap = " << vehCap << endl;

    //determine maximum route duration
    double vehDura = 0.0;
    size_t numArcs = 0;
    double totalDist = 0.0;
    for (size_t i = 0; i < myNumNode; i++) {
        for (size_t j = 0; j < myNumNode; j++) {
            if (myNodeDistMat[i][j] > 0) {
                numArcs++;
                totalDist += myNodeDistMat[i][j];
            }
        }
    }
    double avgDist = (double) totalDist / (double) numArcs;
    double totalTime = 0.0;
    for (size_t i = 1; i < myNumNode; i++) {
        for (size_t d = 0; d < myNumDays; d++) {
            if (myCusDemand[i][d] > 0) {
                totalTime += myCusSrvTime[i][d];
            }
        }
    }
    double avgSrvTime = (double) totalTime / (double) numVisit;
    aval = 0.4 * myNumCus;
    bval = 0.7 * myNumCus;
    cval = 0.6 * myNumCus;
    r = triangular(aval, bval, cval);
    vehDura = r * (double) (avgDist + avgSrvTime);
    cout << "vehDura = " << vehDura << endl;

    //determine maximum arrival time differential
    double maxTimeDiff = avgDist * 0.9;
    //    cout << "avgArcDist = " << avgDist << endl;
    cout << "maxTimeDiff = " << maxTimeDiff << endl;
}

/**
 * This function is used to generate bash scripts
 * @param instsize
 * @param numdays
 * @param dpttype
 * @param instid
 * Updated: 06/29/2016
 */
void createShellScriptForOne(size_t instsize, size_t numdays, std::string dpttype, size_t instid) {

    string pbsFileName = "pbs-"; //name for pbs file
    string logFileName = "cplex-"; //name for log file
    string proFileName = "convrp-"; //name for running process

    string common = "n";
    stringstream convert;
    convert << instsize;
    common.append(convert.str());
    convert.str("");
    common.append("-d");
    convert << numdays;
    common.append(convert.str());
    convert.str("");
    common.append("-");
    common.append(dpttype);
    convert << instid;
    common.append(convert.str());

    //create filename for pbs file
    pbsFileName.append(common);
    pbsFileName.append(".sh");
    ofstream ofile(pbsFileName.c_str(), ios::out);
    if (!ofile) {
        cout << "solveByBranchAndPriceXF(): cannot open output file!" << endl;
    }

    //create filename for log file
    logFileName.append(common);
    logFileName.append(".log");

    //create filename for process
    proFileName.append(common);

    ofile << "#!/bin/bash" << endl;
    ofile << "#PBS -N " << proFileName << endl;
    ofile << "#PBS -j oe" << endl;
    ofile << "#PBS -m ae" << endl;
    ofile << "#PBS -o ~/run/bp1/logfile/zzz.$PBS_JOBID.script" << endl;
    ofile << "#PBS -l nodes=1:ppn=12,walltime=72:00:00" << endl;
    ofile << "#PBS -q med12core" << endl;
    ofile << "cd $PBS_O_WORKDIR" << endl;
    ofile << "~/run/bp1/mainfile/main-" << common << " >" << "~/run/bp1/logfile/" << logFileName << endl;
    ofile.close();
}

/**
 * Create shell scripts for all processes
 * Updated: 06/29/2016
 */
void createShellScriptForAll() {

    for (size_t i = 1; i <= 5; i++) {
        //instsize = 10
        size_t instsize = 10;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 12;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 15;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 18;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 20;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 22;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 25;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 28;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 30;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 32;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 35;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 38;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 40;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 42;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 45;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 48;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);

        instsize = 50;
        createShellScriptForOne(instsize, 3, "c", i);
        createShellScriptForOne(instsize, 3, "e", i);
        createShellScriptForOne(instsize, 3, "r", i);
    }
}

/**
 * generate bash script for ConVRP model to solve 10 small instances
 * Updated: 07/03/2016
 */
void getScriptForConVRPModelOneSmallInst(size_t instsize, size_t instid, size_t modelid) {

    string pbsFileName = "pbs-"; //name for pbs file
    string logFileName = "cplex-"; //name for log file
    string proFileName = "convrp-"; //name for running process

    string common = "";
    stringstream convert;
    convert << instsize;
    common.append(convert.str());
    convert.str("");
    common.append("-");
    convert << instid;
    common.append(convert.str());

    //create filename for pbs file
    pbsFileName.append(common);
    pbsFileName.append(".sh");
    ofstream ofile(pbsFileName.c_str(), ios::out);
    if (!ofile) {
        cout << "solveByBranchAndPriceXF(): cannot open output file!" << endl;
    }

    //create filename for log file
    logFileName.append(common);
    logFileName.append(".log");

    //create filename for process
    proFileName.append(common);

    ofile << "#!/bin/bash" << endl;
    ofile << "#PBS -N " << proFileName << endl;
    ofile << "#PBS -j oe" << endl;
    ofile << "#PBS -m ae" << endl;
    ofile << "#PBS -o ~/run/convrp/convrp" << modelid << "/logfile/zzz.$PBS_JOBID.script" << endl;
    ofile << "#PBS -l nodes=1:ppn=12,walltime=06:00:00" << endl;
    ofile << "#PBS -q tiny12core" << endl;
    ofile << "cd $PBS_O_WORKDIR" << endl;
    ofile << "~/run/convrp/convrp" << modelid << "/mainfile/main-" << common << " >"
            << "~/run/convrp/convrp" << modelid << "/logfile/" << logFileName << endl;
    ofile.close();
}

/**
 * generate bash scripts for ConVRP model to solve all ten small instances
 * Updated: 07/03/2016
 */
void getScriptForConVRPModelAllSmallInst(size_t modelid) {

    for (size_t i = 1; i <= 5; i++) {
        //instsize = 10
        size_t instsize = 10;
        getScriptForConVRPModelOneSmallInst(instsize, i, modelid);

        //instsize = 12
        instsize = 12;
        getScriptForConVRPModelOneSmallInst(instsize, i, modelid);
    }
}

void getScriptForConVRPModelOneLargeInst(size_t instsize, size_t numdays, std::string dpttype, size_t instid, size_t modelid) {

    string pbsFileName = "pbs-"; //name for pbs file
    string logFileName = "cplex-"; //name for log file
    string proFileName = "convrp-"; //name for running process

    string common = "n";
    stringstream convert;
    convert << instsize;
    common.append(convert.str());
    convert.str("");
    common.append("-d");
    convert << numdays;
    common.append(convert.str());
    convert.str("");
    common.append("-");
    common.append(dpttype);
    convert << instid;
    common.append(convert.str());

    //create filename for pbs file
    pbsFileName.append(common);
    pbsFileName.append(".sh");
    ofstream ofile(pbsFileName.c_str(), ios::out);
    if (!ofile) {
        cout << "solveByBranchAndPriceXF(): cannot open output file!" << endl;
    }

    //create filename for log file
    logFileName.append(common);
    logFileName.append(".log");

    //create filename for process
    proFileName.append(common);

    ofile << "#!/bin/bash" << endl;
    ofile << "#PBS -N " << proFileName << endl;
    ofile << "#PBS -j oe" << endl;
    ofile << "#PBS -m ae" << endl;
    ofile << "#PBS -o ~/run/convrp/convrp" << modelid << "/logfile/zzz.$PBS_JOBID.script" << endl;
    ofile << "#PBS -l nodes=1:ppn=12,walltime=72:00:00" << endl;
    ofile << "#PBS -q med12core" << endl;
    ofile << "cd $PBS_O_WORKDIR" << endl;
    ofile << "~/run/convrp/convrp" << modelid << "/mainfile/main-" << common << " >"
            << "~/run/convrp/convrp" << modelid << "/logfile/" << logFileName << endl;
    ofile.close();
}

void getScriptForConVRPModelAllLargeInst(size_t modelid) {

    for (size_t i = 1; i <= 5; i++) {
        //instsize = 10
        size_t instsize = 10;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 12;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 15;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 18;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 20;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 22;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 25;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 28;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 30;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 32;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 35;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 38;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 40;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 42;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 45;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 48;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);

        instsize = 50;
        getScriptForConVRPModelOneLargeInst(instsize, 3, "c", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "e", i, modelid);
        getScriptForConVRPModelOneLargeInst(instsize, 3, "r", i, modelid);
    }
}
