#include <iostream>
#include <cmath>
#include "constants.h"
#include <complex>
#include "fftw3.h"


int main() {

    /* -------------- Controlled variables -------------- */
    /* Grid parameters */
    const int Nx = 128;  // Grid pts
    const int Mx = Nx / 2;
    const int dx = 1.;   // Grid spacing
    const double dkx = PI / (Mx * dx);  // k-space spacing
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

    /* -------------- Set up FFT plans -------------- */
    fftw_plan p_forward, p_back;

    std::complex<double> psi[Nx], psi_k[Nx];

    p_forward = fftw_plan_dft_1d(Nx, reinterpret_cast<fftw_complex *>(psi), reinterpret_cast<fftw_complex *>(psi),
                                 FFTW_FORWARD, FFTW_ESTIMATE);
    p_back = fftw_plan_dft_1d(Nx, reinterpret_cast<fftw_complex *>(psi_k), reinterpret_cast<fftw_complex *>(psi),
                              FFTW_BACKWARD, FFTW_ESTIMATE);

    /* -------------- Generate initial Gaussian -------------- */
    for (int i = 0; i < Nx; ++i) {
        psi[i] = exp((-x[i] * x[i]) / 4);
    };


    return 0;
}
