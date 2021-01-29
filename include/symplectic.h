//
// Created by mattw on 28/01/2021.
//

#ifndef INC_1D_GPESOLVER_SYMPLECTIC_H
#define INC_1D_GPESOLVER_SYMPLECTIC_H

#include <complex>

namespace symplectic {
    double normalise(const std::complex<double> *t_wfn, int t_arrayLength, double t_dx);
}
#endif // INC_1D_GPESOLVER_SYMPLECTIC_H
