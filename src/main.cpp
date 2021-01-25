#include <iostream>
#include <cmath>
#include "constants.h"

int main() {

    /* -------------- Controlled variables -------------- */
    /* Grid parameters */
    const int Nx = 128;
    const int Mx = Nx / 2;
    const int dx = 1.;
    const double dkx = PI / (Mx * dx);
    double x[Nx];
    double kx[Nx];
    double V[Nx];

    /* Generate 1D grid and harmonic trap */
    for (int i = 0; i < Nx; ++i) {
        x[i] = (-Mx + double(i)) * dx;
        V[i] = 0.5 * pow(x[i], 2);
    }

    /* Generate 1D shifted k-space grid */
    for (int i = 0; i < Nx; ++i) {
        if (i < Nx / 2) {
            kx[i] = double(i) * dkx;
        }
        if (i >= Nx / 2) {
            kx[i] = double(-Nx + i) * dkx;
        }
    }

    /* Condensate parameters */
    double c0 = 1.;  // Interaction strength

    /* Time parameters */
    int Nt = 10000;
    double dt = 1e-2;
    double t = 0;

    /* -------------- Initial State -------------- */

    return 0;
}
