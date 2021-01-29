//
// Created by mattw on 28/01/2021.
//

#include "symplectic.h"

namespace symplectic {
    double normalise(const std::complex<double> *t_wfn, const int t_arrayLength, const double t_dx) {
        double sum = 0.;

        for (int i = 0; i < t_arrayLength; ++i) {
            sum += (std::abs(t_wfn[i]) * std::abs(t_wfn[i])) * t_dx;
        }

        return sum;
    }
} // namespace symplectic