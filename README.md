# density evolution demo program

## Preparation
Install following tools.
- gcc, g++
- make
- cmake
- boost
- Qt5
- qwt
- fftw

## Downdloading the demo program
Clone demo program's repository.

```sh
git clone https://github.com/xenocaliver/density.git
```
## Build the demo program
Build demo program as follows.

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