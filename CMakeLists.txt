cmake_minimum_required(VERSION 3.17)
project(1D_GPESolver)

set(CMAKE_CXX_STANDARD 14)

INCLUDE_DIRECTORIES(C:/dev/cpp/libs/fftw-3.3.5-dll64)
LINK_DIRECTORIES(C:/dev/cpp/libs/fftw-3.3.5-dll64)

add_executable(1D_GPESolver src/main.cpp)
TARGET_LINK_LIBRARIES(1D_GPESolver libfftw3-3.lib)