//
// Created by mattw on 28/01/2021.
//

#include "symplectic.h"

namespace symplectic {
    double getAtomNum(const std::complex<double> *wfn, int arrayLength, double dx) {
        double atom_num = 0.;

        for (int i = 0; i < arrayLength; ++i) {
            atom_num += (std::abs(wfn[i]) * std::abs(wfn[i])) * dx;
        }

        return atom_num;
    }

    void potentialEvolution(std::complex<double> *wfn, const std::vector<double>& pot, double g, double dt, int arrayLength) {
        // Computes the potential evolution part of GPE
        for (int i = 0; i < arrayLength; ++i) {

            wfn[i] = wfn[i] * exp(-0.5 * dt * (pot[i] + g * abs(wfn[i]) * abs(wfn[i])));
        }
    }

    void kineticEvolution(std::complex<double> *wfn_k, const std::vector<double>& wvn, double dt, int arrayLength) {
        // Computes the kinetic evolution part of GPE
        for (int i = 0; i < arrayLength; ++i) {
            wfn_k[i] = (wfn_k[i] * exp(-0.5 * dt * wvn[i] * wvn[i])) / static_cast<double>(arrayLength);
        }
    }
} // namespace symplectic