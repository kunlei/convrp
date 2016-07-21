/* 
 * File:   UserStruct.h
 * Author: klian
 * Comment: This file defines structures that will be used by other functions. This is the fundamental
 * file that does not depend on any other files.
 * 
 * Created on July 11, 2016, 11:08 AM
 */

#ifndef USERSTRUCT_H
#define	USERSTRUCT_H

#define _USE_MATH_DEFINES
#include <algorithm>
#include <numeric>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <string>
#include <cstring>
#include <sstream>
#include <vector>
#include <list>
#include <deque>
#include <cmath>
#include <math.h>
#include <limits>
#include <time.h>
#include <cstdlib>
#include <stdio.h>
#include <stdlib.h>
#include <stdexcept>
#include <sys/stat.h>
#include <sys/types.h>
#include <sys/time.h>

template <typename Iter>
bool next(Iter begin, Iter end) {
    if (begin == end) // changed all digits
    { // so we are back to zero
        return false; // that was the last number
    }
    --end;
    if ((*end & 1) == 0) // even number is treated as zero
    {
        ++*end; // increase to one
        return true; // still more numbers to come
    } else // odd number is treated as one
    {
        --*end; // decrease to zero
        return next(begin, end); // RECURSE!
    }
}

//unsigned int and double
struct UIntDbl {
    size_t myUInt; //size_t, unsigned int
    double myDbl; //double

    UIntDbl() : myUInt(0), myDbl(0.0) {}
    UIntDbl(size_t i) : myUInt(i), myDbl(0.0) {}
    UIntDbl(size_t i, double d): myUInt(i), myDbl(d) {}

    static bool sortByDblInc(const UIntDbl &first, const UIntDbl &second) {
        return first.myDbl < second.myDbl;
    }

    static bool sortByDblDec(const UIntDbl &first, const UIntDbl &second) {
        return first.myDbl > second.myDbl;
    }
};
void print(std::vector<UIntDbl> &myvec);

//used in convex hull construction
struct NPoint {
    size_t myCusID;
    double x;
    double y;
    double length;
    bool flag;

    NPoint() : myCusID(0), x(0.0), y(0.0), length(0.0), flag(true) {}
    NPoint(size_t cus, double xc, double yc) : myCusID(cus), x(xc), y(yc) {
        length = x * x + y * y;
        flag = true;
    }
};
void sort(std::vector<NPoint> &myvec);

//represent an arc
struct Arc {
public:
    size_t myDay;
    size_t myHead;
    size_t myTail;
    double myValue;
    
public:
    Arc() : myDay(0), myHead(0), myTail(0), myValue(0.0) {};
    
    Arc(size_t day, size_t head, size_t tail): myDay(day), myHead(head), myTail(tail) {
        myValue = 0.0;
    }
    
    static bool sortByDblDec(const Arc &first, const Arc &second) {
        return first.myValue > second.myValue;
    }
};

#endif	/* USERSTRUCT_H */

