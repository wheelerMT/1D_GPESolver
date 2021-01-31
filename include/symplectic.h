//
// Created by mattw on 28/01/2021.
//

#ifndef INC_1D_GPESOLVER_SYMPLECTIC_H
#define INC_1D_GPESOLVER_SYMPLECTIC_H

#include <complex>
#include <vector>

namespace symplectic {
    double normalise(std::complex<double> *wfn, int arrayLength, double dx, double N);

    void potentialEvolution(std::complex<double> *wfn, const std::vector<double> &pot,
                            double g, double dt);

    void kineticEvolution(std::complex<double> *wfn_k, const std::vector<double> &wvn, double dt);
}
#endif // INC_1D_GPESOLVER_SYMPLECTIC_H
