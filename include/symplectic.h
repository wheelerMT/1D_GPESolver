//
// Created by mattw on 28/01/2021.
//

#ifndef INC_1D_GPESOLVER_SYMPLECTIC_H
#define INC_1D_GPESOLVER_SYMPLECTIC_H

#include <complex>
#include <vector>

namespace symplectic {
    double getAtomNum(const std::complex<double> *wfn, int arrayLength, double dx);

    void potentialEvolution(std::complex<double> *wfn, const std::vector<double>& pot,
                            double g, double dt, int arrayLength);

    void kineticEvolution(std::complex<double> *wfn_k, const std::vector<double>& wvn, double dt, int arrayLength);
}
#endif // INC_1D_GPESOLVER_SYMPLECTIC_H
