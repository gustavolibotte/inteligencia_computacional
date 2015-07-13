#include <cstdlib>
#include <cmath>
#include "EVAL_HAMMER.h"
#include "ED.h"


int main(){
//    double solPSO[Dimension] = {0.72763, 0.19756, 0.76200, 2.75342, 2.90602, 4.99963, 100, 201};
//    double fitPSO = funcaoObjetivo(solPSO);
////
//    double solED_old[Dimension] = {0.7281, 0.1932, 0.7610, 2.7318, 2.8657, 4.9686, 100, 201};
//    double fitED_old = funcaoObjetivo(solED_old);
////
//    double solED[Dimension] = {0.8595113174, 0.1108324252, 0.6520854355, 2.2093628594, 2.3331228566, 4.4609036845, 100, 201};
//    double fitED = funcaoObjetivo(solED);

    srand (time(NULL));

    double* xotm = EvolucaoDiferencialDPS();

    cout << "O ponto otimo e':" << endl << endl;
    for (int i = 0; i < 8; i++){
        cout<< xotm[i] << endl;
    }
}
