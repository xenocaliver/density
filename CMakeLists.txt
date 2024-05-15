cmake_minimum_required(VERSION 3.8)
project(density CXX)
# option for C++ standard
set(CMAKE_CXX_STANDARD 20)

# add current directory for include search path
set(CMAKE_INCLUDE_CURRENT_DIR ON)
# auto execution for moc
set(CMAKE_AUTOMOC ON)
set(CMAKE_AUTOUIC ON)
set(CMAKE_AUTORCC ON)

# search package Qt6
find_package(Qt6Core REQUIRED)
find_package(Qt6Widgets REQUIRED)

set(Qt6_DIR /ope/homebrew/Cellar/qt)
list(APPEND CMAKE_PREFIX_PATH ${Qt6_DIR})

set(Qwt_INCLUDE_DIR /opt/homebrew/Cellar/qwt/6.3.0/lib/qwt.framework/Versions/6/Headers/)
file(GLOB_RECURSE CPP_FILES *.cpp)

include_directories(.)

set(fftw3_PREFIX ${CMAKE_CURRENT_BINARY_DIR}/fftw3)
set(fftw3_INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/lib)
set(fftw3_BINARY_DIR ${CMAKE_CURRENT_BINARY_DIR}/fftw3/build)
set(fftw3_CMAKE_ARGS -DCMAKE_INSTALL_PREFIX=${fftw3_INSTALL_DIR})
set(fftw3_VERSION 3.3.10)

include(ExternalProject)
ExternalProject_Add(fftw3
### Set prefix
    PREFIX ${fftw3_PREFIX}
### Set URL
    URL https://www.fftw.org/fftw-${fftw3_VERSION}.tar.gz
### Set MD5
    URL_MD5 8ccbf6a5ea78a16dbc3e1306e234cc5c
### Set source directory
    SOURCE_DIR ${fftw3_PREFIX}/src
### prohibit extraction
    DOWNLOAD_NO_EXTRACT true
### Set download directory
    TMP_DIR ${fftw3_PREFIX}/tmp
    DOWNLOAD_DIR ${fftw3_PREFIX}/tarballs
    BINARY_DIR ${fftw3_BINARY_DIR}
    INSTALL_DIR ${fftw3_INSTALL_DIR}
    CMAKE_ARGS ${galois_CMAKE_ARGS}
### CONFIGURE_COMMAND
    CONFIGURE_COMMAND cd ${fftw3_PREFIX}/src && tar xvfz ${fftw3_PREFIX}/tarballs/fftw-${fftw3_VERSION}.tar.gz && cd ${fftw3_PREFIX}/src/fftw-${fftw3_VERSION} && ./configure --enable-threads && make
### BUILD_COMMAND
    BUILD_COMMAND ""
### UPDATE COMMAND
    UPDATE_COMMAND ""
### install command
    INSTALL_COMMAND cp ${fftw3_PREFIX}/src/fftw-${fftw3_VERSION}/.libs/libfftw3.a ${fftw3_INSTALL_DIR}/ && 
                    cp ${fftw3_PREFIX}/src/fftw-${fftw3_VERSION}/threads/.libs/libfftw3_threads.a ${fftw3_INSTALL_DIR}/ &&
                    mkdir ${PROJECT_BINARY_DIR}/include &&
                    cp -r ${fftw3_PREFIX}/src/fftw-${fftw3_VERSION}/api/fftw3.h ${PROJECT_BINARY_DIR}/include/
)

set(Qwt_LIBRARY_DIR /opt/homebrew/lib/qwt.framework)

add_executable(density density.cpp common.hpp plot.hpp plot.cpp load_json.cpp)
add_executable(biterror biterror.cpp common.hpp plot.hpp plot.cpp load_json.cpp)

set(CMAKE_CXX_FLAGS "-O3 -Wall -g -fPIC -fmessage-length=0 -I/usr/local/include -I/opt/homebrew/include -I${Qwt_INCLUDE_DIR} -I${PROJECT_BINARY_DIR}/include -L/opt/homebrew/lib")
target_link_libraries(density LINK_PUBLIC Qt6::Widgets Qt6::Core ${Qwt_LIBRARY_DIR}/qwt ${PROJECT_BINARY_DIR}/lib/libfftw3.a ${PROJECT_BINARY_DIR}/lib/libfftw3_threads.a m)
target_link_libraries(biterror LINK_PUBLIC Qt6::Widgets Qt6::Core ${Qwt_LIBRARY_DIR}/qwt ${PROJECT_BINARY_DIR}/lib/libfftw3.a ${PROJECT_BINARY_DIR}/lib/libfftw3_threads.a m)