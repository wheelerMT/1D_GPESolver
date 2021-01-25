#include <iostream>
#include "constants.h"

int main() {

    /* -------------- Controlled variables -------------- */
    /* Grid parameters */
    const int Nx =1024, Ny =1024;
    const int Mx = Nx / 2, My = Ny / 2;
    const int dx = 1, dy = 1;
    const float dkx = PI / (Mx * dx), dky = PI / (My * dy);
    double x[Nx], y[Ny];
    double kx[Nx], ky[Ny];

    /* Generate 1D grids */
    for (int i = 0; i < Nx; ++i) {
        x[i] = (-Mx + i) * dx;
        y[i] = (-My + i) * dy;
    }

    /* Generate 1D shifted k-space grids */
    for (int i = 0; i < Nx; ++i) {
        if (i < Nx / 2) {
            kx[i] = float(i) * dkx;
            ky[i] = float(i) * dky;
        }
        if (i >= Nx / 2) {
            kx[i] = float(-Nx + i) * dkx;
            ky[i] = float(-Ny + i) * dky;
        }
    }

    return 0;
}
