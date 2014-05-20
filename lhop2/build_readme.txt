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

Build for windows:
    Required compilers: 
     - CMake
     - MS Visual Studio 2010 (SP1) with C++11 support (VS2012 may not work properly)

    Required libraries:
     - OpenCV2.3.1 or newer (build from source)
       - after building and compiling add path to binaries (DLL) to the system variable PATH
       - also add system variable OpenCV_DIR that points to build folder of OpenCV (e.g. c:\work\C\OpenCV2.3.1\build)
     - Intel TBB 3.0 for Open Source 
       - after installing make sure a system variable PATH contains path to correct version for binaries (DLL)
     - gzib (zlib)
       - zlib can be build from source http://www.zlib.net (best way is to compile static library and directly embedd it to libhop library)
     - (optional protobuffer if compiling with java JNI support) 
     
    Example of CMake variables:        
    OpenCV:
     - OpenCV_DIR: C:/work/C/OpenCV2.3.1/build
    TBB:
     - TBB_ARCHITECTURE: intel64 
     - TBB_INCLUDE_DIR: C:/Program Files/Intel/TBB/include 
     - TBB_INCLUDE_DIRS: C:/Program Files/Intel/TBB/include 
     - TBB_INSTALL_DIR: C:/Program Files/Intel/TBB;C:/Program Files (x86)/Intel/TBB 
     - TBB_LIBRARY: C:/Program Files/Intel/TBB/lib/intel64/vc10/tbb.lib 
     - TBB_LIBRARY_DEBUG: C:/Program Files/Intel/TBB/lib/intel64/vc10/tbb_debug.lib 
     - TBB_LIBRARY_DIRS: C:/Program Files/Intel/TBB/lib/intel64/vc10 
     - TBB_MALLOC_LIBRARY: C:/Program Files/Intel/TBB/lib/intel64/vc10/tbbmalloc.lib 
     - TBB_MALLOC_LIBRARY_DEBUG: C:/Program Files/Intel/TBB/lib/intel64/vc10/tbbmalloc_debug.lib 
    ZLIB: 
     - ZLIB_INCLUDE_DIR: C:/work/libhop_vs10/include 
     - ZLIB_LIBRARY: C:/work/libhop_vs10/lib/gzip_x64_vs10.lib
    OpenCL (if enabled):
     - OPENCL_INCLUDE_DIR: C:/Program Files (x86)/ATI Stream/include 
     - OPENCL_LIBRARY: C:/Program Files (x86)/ATI Stream/lib/x86/OpenCL.lib 