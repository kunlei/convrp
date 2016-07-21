/* 
 * File:   Customer.h
 * Author: klian
 *
 * Created on July 11, 2016, 11:19 AM
 */

#ifndef CUSTOMER_H
#define	CUSTOMER_H
#include "Data.h"

//define a customer
class Customer {    
public:
    size_t myID;
    size_t myDay;
    bool needVisit;
    bool isVisited;

    size_t myDriver;
    size_t myPrevCus;
    size_t myNextCus;
    double myArrvTime;
    
public:
    Customer();
    Customer(size_t cus, size_t day);
    ~Customer();

    void resetState();
    double getLeaveTime() const;
    bool isEqualTo(Customer other);
};

#endif	/* CUSTOMER_H */

