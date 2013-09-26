Build for linux:
    Required compilers: 
      - CMake
      - gcc-4.6 (or newer) with C++11 support

    Required libraries
     - OpenCV2.3.1 or newer (build from source)
     - Intel TBB 3.0 for Open Source 
     - gzib (zlib)
     - (optional protobuffer if compiling with java JNI support) 

    On Ubuntu 11.10 or newer run:
      sudo apt-get install gcc g++ cmake swig libtbb-dev libprotoc-dev protobuf-compiler libprotobuf-java zlib1g-dev
      <build OpenCV>
      <go to libhop root folder>
      mkdir build
      cd build
      ccmake .. 
      <IMPORTANT !! during cmake configuration press t and add '-std=c++0x' to CMAKE_CXX_FLAGS and CMAKE_C_FLAGS !!>
      make
      sudo make install

Please check the User Manual for the details of installation steps.
