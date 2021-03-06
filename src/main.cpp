#include "constants.h"
#include "fftw3.h"
#include "symplectic.h"
#include <complex>
#include <iostream>
#include <vector>
#include <fstream>

int main() {

    // -------------- Controlled variables -------------- //
    // Grid parameters //
    const int Nx = 4096; // Grid pts
    const int Mx = Nx / 2;
    const double dx = 1.;              // Grid spacing
    const double dkx = PI / (Mx * dx); // k-space spacing
    std::vector<double> x(Nx);
    std::vector<double> kx(Nx);
    std::vector<double> V(Nx);

    // Generate 1D grid and harmonic trap
    for (int i = 0; i < Nx; ++i) {
        x[i] = (-Mx + double(i)) * dx;
        V[i] = 0.5 * pow(x[i], 2);
    }

    // Generate 1D shifted k-space grid //
    for (int i = 0; i < Nx; ++i) {
        if (i < Nx / 2) {
            kx[i] = double(i) * dkx;
        }
        if (i >= Nx / 2) {
            kx[i] = double(-Nx + i) * dkx;
        }
    }

    // Condensate parameters //
    const double c0 = 10.;   // Interaction strength
    const double N = 100000.; // Atom number

    // Time parameters //
    const int Nt = 10000;   // Number of timesteps
    const double dt = 1e-2; // Timestep
    double t = 0;     // Time

    // Wavefunction //
    std::complex<double> psi[Nx], psi_k[Nx];

    // -------------- Set up FFT plans -------------- //
    fftw_plan p_forward, p_back;

    p_forward = fftw_plan_dft_1d(Nx, reinterpret_cast<fftw_complex *>(psi), reinterpret_cast<fftw_complex *>(psi_k),
                                 FFTW_FORWARD, FFTW_ESTIMATE);
    p_back = fftw_plan_dft_1d(Nx, reinterpret_cast<fftw_complex *>(psi_k), reinterpret_cast<fftw_complex *>(psi),
                              FFTW_BACKWARD, FFTW_ESTIMATE);

    // -------------- Generate initial Gaussian -------------- //
    for (int i = 0; i < Nx; ++i) {
        psi[i] = sqrt(N) / Nx * exp((-x[i] * x[i]) / 2000);
    }

    // -------------- Imaginary time evolution -------------- //
    for (int i = 0; i < Nt; ++i) {

        // Potential half-step evolution:
        symplectic::potentialEvolution(psi, V, c0, dt);

        fftw_execute(p_forward); // Forward FFT

        // Kinetic half-step evolution:
        symplectic::kineticEvolution(psi_k, kx, dt);

        fftw_execute(p_back); // Backward FFT

        // Potential half-step evolution:
        symplectic::potentialEvolution(psi, V, c0, dt);

        // Renormalise wavefunction:
        symplectic::normalise(psi, Nx, dx, N);

        // Print time
        if (i % 100 == 0) { std::cout << "t = " << t << '\n'; }

        t += dt;
    }

    std::ofstream outfile("data.txt", std::ofstream::binary);

    for (auto &i : psi) {
        outfile.write((char *) &i, sizeof(std::complex<double>));
    }

    return 0;
}
