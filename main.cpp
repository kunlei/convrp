/*
 * File:   main.cpp
 * Author: Liya
 *
 * Created on July 18, 2016, 12:08 PM
 */

#include "Data.h"
#include "AssistFunc.h"
#include "ConVRPSol.h"

using namespace std;

/*
 *
 */
int main(int argc, char** argv) {

    size_t instsize = 10;
    size_t instid = 1;
    readSmallInst(instsize, instid);

    ConVRPSol mysol;
    if (mysol.isFeasible()) {
        cout << "Feasible!" << endl;
    } else {
        cout << "Infeasible!" << endl;
    }

    return 0;
}

