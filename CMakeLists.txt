cmake_minimum_required(VERSION 3.17)
project(1D_GPESolver)

set(CMAKE_CXX_STANDARD 17)

INCLUDE_DIRECTORIES(C:/dev/cpp/libs/fftw-3.3.5-dll64 include/)
LINK_DIRECTORIES(C:/dev/cpp/libs/fftw-3.3.5-dll64)

find_package(OpenMP REQUIRED)

add_executable(1D_GPESolver
        src/main.cpp

        include/constants.h

        src/symplectic.cpp include/symplectic.h)

TARGET_LINK_LIBRARIES(1D_GPESolver libfftw3-3.lib OpenMP::OpenMP_CXX)