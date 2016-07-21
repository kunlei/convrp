#include "Data.h"

using namespace std;

size_t myInstID = 0;
size_t myNumNode;
size_t myNumCus;
size_t myNumDays;
size_t myNumVeh;
size_t myVehCap;
double myMaxDuration;
double myMaxTimeDiff;

double **myNodeCoord; //myNumNode * 2
size_t **myCusDemand; //myNumNode * myNumDays
double **myCusSrvTime; //myNumNode * myNumDays
size_t **mySrvIndict; //myNumNode * myNumDays
double **myNodeDistMat; //myNumNode * myNumNode
double *myCusDepotAngle; //1 * myNumNode
size_t *myCusSrvFreq; //1 * myNumNode
vector<size_t> myCusVec; //a list consisting of all the customers
vector<size_t> myCusVecAngleSorted; //all customers sorted in angle made with depot

/**
 * read small-sized instances, there are a total of 10 such instances
 * @param instsize: instance size, either 10 or 12
 * @param instid: instance id, from 1 to 5
 * Updated: 02/27/2016, 04/20/2016
 */
void readSmallInst(size_t instsize, size_t instid) {

    if (instid < 1 || instid > 5) {
        cout << "readSmallInst(): instance id out of range!" << endl;
        exit(1);
    }

    if (myInstID > 0) {
        //In this case a previous instance was read, free the instance memory first
        //before reading the next instance file
        for (size_t i = 0; i < myNumNode; i++) {
            delete[] myNodeCoord[i];
            delete[] myCusDemand[i];
            delete[] mySrvIndict[i];
            delete[] myCusSrvTime[i];
            delete[] myNodeDistMat[i];
        }
        delete[] myNodeCoord;
        delete[] myNodeDistMat;
        delete[] myCusSrvTime;
        delete[] myCusDemand;
        delete[] mySrvIndict;
        delete[] myCusDepotAngle;
        delete[] myCusSrvFreq;
        myCusVec.clear();
        myCusVecAngleSorted.clear();
    }

    myInstID = instid;

    //determine filename
    string filename = "Instance/convrp_";
    ostringstream convert;
    convert << instsize;
    filename.append(convert.str());
    filename.append("_test_");
    convert.str("");
    convert << myInstID;
    filename.append(convert.str());
    filename.append(".vrp");
    //    cout << filename << endl;

    ifstream pfile(filename.c_str(), ios::in);
    if (!pfile) {
        cerr << "readSmallInst(): Cannot open file!" << endl;
        exit(1);
    } else {
        size_t i, j;
        //read instance description information
        string comment;
        for (i = 0; i < 10; i++) {
            pfile >> comment;
        }
        //read total number of nodes in the instance
        pfile >> comment >> myNumNode;
        myNumCus = myNumNode - 1;
        myNumVeh = myNumCus; //the number of available vehicles
        //read total number of days in the planning horizon
        pfile >> comment >> myNumDays;
        //read vehicle capacity
        pfile >> comment >> myVehCap;
        //read route duration
        pfile >> comment >> myMaxDuration;
        myMaxTimeDiff = 5.0;
        //read other descriptive information
        for (i = 0; i < 7; i++)
            pfile >> comment;

        //read node coordinates
        myNodeCoord = new double *[myNumNode];
        for (i = 0; i < myNumNode; i++)
            myNodeCoord[i] = new double[2]();
        for (i = 1; i < myNumNode; i++)
            pfile >> j >> myNodeCoord[i][0] >> myNodeCoord[i][1];

        //read customer demand
        pfile >> comment;
        myCusDemand = new size_t *[myNumNode];
        myCusDemand[0] = new size_t[myNumDays]();
        mySrvIndict = new size_t *[myNumNode];
        mySrvIndict[0] = new size_t[myNumDays]();
        int tDemand = 0;
        for (i = 1; i < myNumNode; i++) {
            myCusDemand[i] = new size_t[myNumDays]();
            mySrvIndict[i] = new size_t[myNumDays]();
            pfile >> j;
            for (j = 0; j < myNumDays; j++) {
                pfile >> tDemand;
                if (tDemand > 0) {
                    myCusDemand[i][j] = (size_t) tDemand;
                    mySrvIndict[i][j] = 1;
                } else {
                    myCusDemand[i][j] = 0;
                    mySrvIndict[i][j] = 0;
                }
            }
        }

        //read customer service time
        pfile >> comment;
        myCusSrvTime = new double *[myNumNode];
        myCusSrvTime[0] = new double[myNumDays]();
        for (i = 1; i < myNumNode; i++) {
            myCusSrvTime[i] = new double[myNumDays]();
            pfile >> j;
            for (j = 0; j < myNumDays; j++) {
                pfile >> myCusSrvTime[i][j];
            }
        }

        //read depot coordinates
        pfile >> comment;
        pfile >> myNodeCoord[0][0] >> myNodeCoord[0][1];

        //end of file
        pfile >> j >> comment;
        if (comment == "EOF") {
            pfile.close();
        }

        createOtherVars();
    }
}

/**
 * read newly-generated instances
 * @param instsize
 * @param numdays
 * @param dpttype
 * @param instid
 * Updated: 06/29/2016
 */
void readNewInst(size_t instsize, size_t numdays, string dpttype, size_t instid) {

    if (myInstID > 0) {
        //in this case, a previous instance was read, free the memory before reading another instance
        for (size_t i = 0; i < myNumNode; i++) {
            delete[] myNodeCoord[i];
            delete[] myCusDemand[i];
            delete[] mySrvIndict[i];
            delete[] myCusSrvTime[i];
            delete[] myNodeDistMat[i];
        }
        delete[] myNodeCoord;
        delete[] myNodeDistMat;
        delete[] myCusSrvTime;
        delete[] myCusDemand;
        delete[] mySrvIndict;
        delete[] myCusDepotAngle;
        delete[] myCusSrvFreq;
        myCusVec.clear();
        myCusVecAngleSorted.clear();
    }

    myInstID = instid;

    //determine filename
    string filename = "Instance/convrp-n";
    stringstream convert;
    convert << instsize;
    filename.append(convert.str());
    convert.str("");
    filename.append("-d");
    convert << numdays;
    filename.append(convert.str());
    convert.str("");
    filename.append("-");
    filename.append(dpttype);
    convert << instid;
    filename.append(convert.str());
    filename.append(".vrp");
    cout << filename << endl;

    ifstream pfile(filename.c_str(), ios::in);
    if (!pfile) {
        cerr << "readNewInst(): cannot open file!" << endl;
        exit(1);
    } else {
        size_t i, j;
        string comment;
        for (size_t i = 0; i < 4; i++) {
            pfile >> comment;
        }
        //read total number of customers
        pfile >> comment >> myNumNode;
        myNumCus = myNumNode - 1;
        myNumVeh = myNumCus;
        //read total number of days
        pfile >> comment >> myNumDays;
        //read vehicle capacity
        pfile >> comment >> myVehCap;
        //read route duration
        pfile >> comment >> myMaxDuration;
        //read arrival time differential limit
        pfile >> comment >> myMaxTimeDiff;
        //read other information
        for (size_t i = 0; i < 7; i++) {
            pfile >> comment;
        }

        //read node coordinates
        myNodeCoord = new double *[myNumNode];
        for (i = 0; i < myNumNode; i++)
            myNodeCoord[i] = new double[2]();
        for (i = 1; i < myNumNode; i++)
            pfile >> j >> myNodeCoord[i][0] >> myNodeCoord[i][1];

        //read node distance matrix
        pfile >> comment;
        myNodeDistMat = new double *[myNumNode];
        for (i = 0; i < myNumNode; i++) {
            myNodeDistMat[i] = new double[myNumNode]();
            for (j = 0; j < myNumNode; j++) {
                pfile >> myNodeDistMat[i][j];
            }
        }


        //read customer demand
        pfile >> comment;
        myCusDemand = new size_t *[myNumNode];
        myCusDemand[0] = new size_t[myNumDays]();
        mySrvIndict = new size_t *[myNumNode];
        mySrvIndict[0] = new size_t[myNumDays]();
        int tDemand = 0;
        for (i = 1; i < myNumNode; i++) {
            myCusDemand[i] = new size_t[myNumDays]();
            mySrvIndict[i] = new size_t[myNumDays]();
            pfile >> j;
            for (j = 0; j < myNumDays; j++) {
                pfile >> tDemand;
                if (tDemand > 0) {
                    myCusDemand[i][j] = (size_t) tDemand;
                    mySrvIndict[i][j] = 1;
                } else {
                    myCusDemand[i][j] = 0;
                    mySrvIndict[i][j] = 0;
                }
            }
        }

        //read customer service time
        pfile >> comment;
        myCusSrvTime = new double *[myNumNode];
        myCusSrvTime[0] = new double[myNumDays]();
        for (i = 1; i < myNumNode; i++) {
            myCusSrvTime[i] = new double[myNumDays]();
            pfile >> j;
            for (j = 0; j < myNumDays; j++) {
                pfile >> myCusSrvTime[i][j];
            }
        }

        //read depot coordinates
        pfile >> comment;
        pfile >> myNodeCoord[0][0] >> myNodeCoord[0][1];

        //end of file
        pfile >> j >> comment;
        if (comment == "EOF") {
            pfile.close();
        }

        //create other variables
        //compute customer depot angles
        myCusDepotAngle = new double[myNumNode]();
        for (size_t i = 0; i < myNumNode; i++) {
            myCusDepotAngle[i] = getCusDepotAngle(myNodeCoord[0][0],
                    myNodeCoord[0][1],
                    myNodeCoord[i][0],
                    myNodeCoord[i][1]);
        }

        //compute customer service frequency
        myCusSrvFreq = new size_t[myNumNode]();
        for (size_t i = 1; i < myNumNode; i++) {
            for (size_t d = 0; d < myNumDays; d++) {
                if (myCusDemand[i][d] > 0) {
                    myCusSrvFreq[i]++;
                }
            }
        }

        //save all customers
        myCusVec.clear();
        myCusVecAngleSorted.clear();
        vector<UIntDbl> myvec;
        for (size_t i = 1; i < myNumNode; i++) {
            myCusVec.push_back(i);
            myvec.push_back(UIntDbl(i, myCusDepotAngle[i]));
        }
        sort(myvec.begin(), myvec.end(), UIntDbl::sortByDblInc);
        for (vector<UIntDbl>::iterator it = myvec.begin(); it != myvec.end(); it++) {
            myCusVecAngleSorted.push_back((*it).myUInt);
        }
    }
}

/**
 * This function is used to read instances that are generated to study the parameters when generating the instances
 * @param instsize: number of customers in the instance
 * @param numdays: number of days in the instance
 * @param dpttype: depot location type in the instance
 * @param instid: instance id
 * Updated: 07/12/2016
 */
void readNewTestInst(size_t instsize, size_t locid, size_t numdays, std::string dpttype, size_t instid) {

    if (myInstID > 0) {
        //in this case, a previous instance was read, free the memory before reading another instance
        for (size_t i = 0; i < myNumNode; i++) {
            delete[] myNodeCoord[i];
            delete[] myCusDemand[i];
            delete[] mySrvIndict[i];
            delete[] myCusSrvTime[i];
            delete[] myNodeDistMat[i];
        }
        delete[] myNodeCoord;
        delete[] myNodeDistMat;
        delete[] myCusSrvTime;
        delete[] myCusDemand;
        delete[] mySrvIndict;
        delete[] myCusDepotAngle;
        delete[] myCusSrvFreq;
        myCusVec.clear();
        myCusVecAngleSorted.clear();
    }

    myInstID = instid;

    //determine filename
    string filename = "Instance/TestNewInst/convrp-n";
    stringstream convert;
    convert << instsize;
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
    filename.append(dpttype);
    convert << instid;
    filename.append(convert.str());
    filename.append(".vrp");
    cout << filename << endl;

    ifstream pfile(filename.c_str(), ios::in);
    if (!pfile) {
        cerr << "readNewInst(): cannot open file!" << endl;
        exit(1);
    } else {
        size_t i, j;
        string comment;
        for (size_t i = 0; i < 4; i++) {
            pfile >> comment;
        }
        //read total number of customers
        pfile >> comment >> myNumNode;
        myNumCus = myNumNode - 1;
        myNumVeh = myNumCus;
        //read total number of days
        pfile >> comment >> myNumDays;
        //read vehicle capacity
        pfile >> comment >> myVehCap;
        //read route duration
        pfile >> comment >> myMaxDuration;
        //read arrival time differential limit
        pfile >> comment >> myMaxTimeDiff;
        //read other information
        for (size_t i = 0; i < 7; i++) {
            pfile >> comment;
        }

        //read node coordinates
        myNodeCoord = new double *[myNumNode];
        for (i = 0; i < myNumNode; i++)
            myNodeCoord[i] = new double[2]();
        for (i = 1; i < myNumNode; i++)
            pfile >> j >> myNodeCoord[i][0] >> myNodeCoord[i][1];

        //read node distance matrix
        pfile >> comment;
        myNodeDistMat = new double *[myNumNode];
        for (i = 0; i < myNumNode; i++) {
            myNodeDistMat[i] = new double[myNumNode]();
            for (j = 0; j < myNumNode; j++) {
                pfile >> myNodeDistMat[i][j];
            }
        }


        //read customer demand
        pfile >> comment;
        myCusDemand = new size_t *[myNumNode];
        myCusDemand[0] = new size_t[myNumDays]();
        mySrvIndict = new size_t *[myNumNode];
        mySrvIndict[0] = new size_t[myNumDays]();
        int tDemand = 0;
        for (i = 1; i < myNumNode; i++) {
            myCusDemand[i] = new size_t[myNumDays]();
            mySrvIndict[i] = new size_t[myNumDays]();
            pfile >> j;
            for (j = 0; j < myNumDays; j++) {
                pfile >> tDemand;
                if (tDemand > 0) {
                    myCusDemand[i][j] = (size_t) tDemand;
                    mySrvIndict[i][j] = 1;
                } else {
                    myCusDemand[i][j] = 0;
                    mySrvIndict[i][j] = 0;
                }
            }
        }

        //read customer service time
        pfile >> comment;
        myCusSrvTime = new double *[myNumNode];
        myCusSrvTime[0] = new double[myNumDays]();
        for (i = 1; i < myNumNode; i++) {
            myCusSrvTime[i] = new double[myNumDays]();
            pfile >> j;
            for (j = 0; j < myNumDays; j++) {
                pfile >> myCusSrvTime[i][j];
            }
        }

        //read depot coordinates
        pfile >> comment;
        pfile >> myNodeCoord[0][0] >> myNodeCoord[0][1];

        //end of file
        pfile >> j >> comment;
        if (comment == "EOF") {
            pfile.close();
        }

        //create other variables
        //compute customer depot angles
        myCusDepotAngle = new double[myNumNode]();
        for (size_t i = 0; i < myNumNode; i++) {
            myCusDepotAngle[i] = getCusDepotAngle(myNodeCoord[0][0],
                    myNodeCoord[0][1],
                    myNodeCoord[i][0],
                    myNodeCoord[i][1]);
        }

        //compute customer service frequency
        myCusSrvFreq = new size_t[myNumNode]();
        for (size_t i = 1; i < myNumNode; i++) {
            for (size_t d = 0; d < myNumDays; d++) {
                if (myCusDemand[i][d] > 0) {
                    myCusSrvFreq[i]++;
                }
            }
        }

        //save all customers
        myCusVec.clear();
        myCusVecAngleSorted.clear();
        vector<UIntDbl> myvec;
        for (size_t i = 1; i < myNumNode; i++) {
            myCusVec.push_back(i);
            myvec.push_back(UIntDbl(i, myCusDepotAngle[i]));
        }
        sort(myvec.begin(), myvec.end(), UIntDbl::sortByDblInc);
        for (vector<UIntDbl>::iterator it = myvec.begin(); it != myvec.end(); it++) {
            myCusVecAngleSorted.push_back((*it).myUInt);
        }
    }
}

/**
 * create other variables for this instance
 * Updated: 02/27/2016, 04/20/2016
 */
void createOtherVars() {

    size_t i, j;
    //compute node distance matrix
    myNodeDistMat = new double *[myNumNode];
    for (i = 0; i < myNumNode; i++) {
        myNodeDistMat[i] = new double[myNumNode]();
        for (j = 0; j < myNumNode; j++) {
            myNodeDistMat[i][j] = roundVal(sqrt(
                    pow(myNodeCoord[j][0] - myNodeCoord[i][0], 2.0) +
                    pow(myNodeCoord[j][1] - myNodeCoord[i][1], 2.0)));
        }
    }

    //compute customer depot angles
    myCusDepotAngle = new double[myNumNode]();
    for (size_t i = 0; i < myNumNode; i++) {
        myCusDepotAngle[i] = getCusDepotAngle(myNodeCoord[0][0],
                myNodeCoord[0][1],
                myNodeCoord[i][0],
                myNodeCoord[i][1]);
    }

    //compute customer service frequency
    myCusSrvFreq = new size_t[myNumNode]();
    for (size_t i = 1; i < myNumNode; i++) {
        for (size_t d = 0; d < myNumDays; d++) {
            if (myCusDemand[i][d] > 0) {
                myCusSrvFreq[i]++;
            }
        }
    }

    //save all customers
    myCusVec.clear();
    myCusVecAngleSorted.clear();
    vector<UIntDbl> myvec;
    for (size_t i = 1; i < myNumNode; i++) {
        myCusVec.push_back(i);
        myvec.push_back(UIntDbl(i, myCusDepotAngle[i]));
    }
    sort(myvec.begin(), myvec.end(), UIntDbl::sortByDblInc);
    for (vector<UIntDbl>::iterator it = myvec.begin(); it != myvec.end(); it++) {
        myCusVecAngleSorted.push_back((*it).myUInt);
    }
}

/**
 * round double variables
 * @param x
 * @return
 */
double roundVal(double x) {

    return (int(x * myPrecision + (x < 0 ? -0.5 : 0.5))) / (double) myPrecision;
}

/**
 * compute the angle between two points
 * @param dx
 * @param dy
 * @param cx
 * @param cy
 * @return
 */
double getCusDepotAngle(double dx, double dy, double cx, double cy) {

    double angle = 0;
    double x = cx - dx;
    double y = cy - dy;
    if (x != 0 || y != 0) {
        double len = sqrt(x * x + y * y);
        double dot = x;
        if (y >= 0) {
            angle = acos(dot / len) * 180 / M_PI;
        } else {
            angle = 360 - acos(dot / len) * 180 / M_PI;
        }
    } else {
        angle = 0.0;
    }

    return angle;
}


