# density evolution demo program

## Preparation
Install following tools.
- gcc, g++
- make
- cmake
- boost
- [Qt6](https://www.qt.io/product/qt6)
- [qwt](https://qwt.sourceforge.io/index.html)
- [fftwi3](http://fftw.org)

gcc, g++, make, cmake, boost, Qt6, and qwt are installed by package management system.

## Downdloading the demo program
Clone demo program's repository.

```sh
git clone https://github.com/xenocaliver/density.git
```
## Build the demo program
Before building the demo program, you modify CMakeLists.txt file for your system. This CMakeLists.txt is writen for mac OS and homebrew.

```sh
cd density
mkdir build
cd build
cmake ..
make
```

## Run the demo program
Run demo program as follows.

```sh
./density <degree distribution polynomial json file> <sigma>
```