#include "Customer.h"

using namespace std;

/**
 * Default constructor
 */
Customer::Customer() {

    myID = 0;
    myDay = 0;
    needVisit = false;
    isVisited = false;

    myDriver = 0;
    myPrevCus = 0;
    myNextCus = 0;
    myArrvTime = 0.0;
}

/**
 * User-defined constructor
 * @param cus
 * @param day
 */
Customer::Customer(size_t cus, size_t day) {

    myID = cus;
    myDay = day;
    needVisit = false;
    isVisited = false;

    myDriver = 0;
    myPrevCus = 0;
    myNextCus = 0;
    myArrvTime = 0.0;
}

Customer::~Customer() {

}

/**
 * Reset customer state
 */
void Customer::resetState() {

    needVisit = false;
    isVisited = false;
    myDriver = 0;
    myPrevCus = 0;
    myNextCus = 0;
    myArrvTime = 0.0;
}

/**
 * Return leave time from the customer
 * @return
 */
double Customer::getLeaveTime() const {

    return myArrvTime + myCusSrvTime[myID][myDay];
}

bool Customer::isEqualTo(Customer other) {

    bool isequal = true;
    if (this->myID != other.myID ||
            this->myDay != other.myDay ||
            this->needVisit != other.needVisit ||
            this->isVisited != other.isVisited ||
            this->myDriver != other.myDriver ||
            this->myPrevCus != other.myPrevCus ||
            this->myNextCus != other.myNextCus ||
            this->myArrvTime != other.myArrvTime) {
        isequal = false;
    }
    return isequal;
}
