# Distributed under the OSI-approved BSD 3-Clause License.  See accompanying
# file Copyright.txt or https://cmake.org/licensing for details.

cmake_minimum_required(VERSION 3.5)

file(MAKE_DIRECTORY
  "/Users/jason/workspace/density/build/fftw3/src"
  "/Users/jason/workspace/density/build/fftw3/build"
  "/Users/jason/workspace/density/build/lib"
  "/Users/jason/workspace/density/build/fftw3/tmp"
  "/Users/jason/workspace/density/build/fftw3/src/fftw3-stamp"
  "/Users/jason/workspace/density/build/fftw3/tarballs"
  "/Users/jason/workspace/density/build/fftw3/src/fftw3-stamp"
)

set(configSubDirs )
foreach(subDir IN LISTS configSubDirs)
    file(MAKE_DIRECTORY "/Users/jason/workspace/density/build/fftw3/src/fftw3-stamp/${subDir}")
endforeach()
if(cfgdir)
  file(MAKE_DIRECTORY "/Users/jason/workspace/density/build/fftw3/src/fftw3-stamp${cfgdir}") # cfgdir has leading slash
endif()
