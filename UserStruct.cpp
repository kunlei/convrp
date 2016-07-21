#include "UserStruct.h"

using namespace std;

/**
 * Print a vector of user-defined UIntDbl objects
 * @param myvec: input data
 * Updated: 07/11/2016
 */
void print(vector<UIntDbl> &myvec) {

    if (myvec.empty()) return;
    for (size_t i = 0; i < myvec.size(); i++) {
        cout << setw(12) << myvec[i].myUInt;
    }
    cout << endl;
    for (size_t i = 0; i < myvec.size(); i++) {
        cout << setw(12) << myvec[i].myDbl;
    }
    cout << endl;
}

/**
 * Sort a vector of user-defined NPoint objects
 * @param myvec: input data
 * Updated: 07/11/2016
 */
void sort(vector<NPoint> &myvec) {

    size_t size = myvec.size();
    for (size_t i = 0; i < size - 1; i++) {
        for (size_t j = i + 1; j < size; j++) {
            double prod = myvec[i].x * myvec[j].y - myvec[j].x * myvec[i].y;
            if (prod < 0) {
                swap(myvec[i], myvec[j]);
            } else if (prod == 0.0) {
                if (myvec[i].length <= myvec[j].length) {
                    myvec[i].flag = false;
                } else {
                    myvec[j].flag = false;
                    swap(myvec[i], myvec[j]);
                }
            }
        }
    }
}
