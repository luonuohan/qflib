
set(qflib_SOURCES
    math/interpol/piecewisepolynomial.cpp 
    math/stats/errorfunction.cpp
    math/linalg/choldcmp.cpp 
    math/linalg/eigensym.cpp
    math/linalg/spectrunc.cpp
    methods/montecarlo/pathgenerator.cpp
    methods/pde/pdebase.cpp
    methods/pde/pde1dsolver.cpp
    pricers/simplepricers.cpp
    pricers/bsmcpricer.cpp
    pricers/multiassetbsmcpricer.cpp 
    pricers/bsmcsimpledeltahedger.cpp    
    market/market.cpp
    market/yieldcurve.cpp
    market/volatilitytermstructure.cpp
)

add_library(qflib STATIC ${qflib_SOURCES})
# add_library(name STATIC paths_of_source_files)

if(APPLE)
    # for #include<>
    target_include_directories(qflib PRIVATE
        .. # parent directory of qflib/
        ${ARMADILLO_PREFIX}/include # where armadillo headers live
        
    )
else()
    target_include_directories(qflib PRIVATE
        ..
        ${THIRDPARTY_DIRECTORY}/armadillo-${ARMA_VERSION}/include
    )
endif()

# compiler and linker option adjustments for targets
if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    target_compile_options(qflib PRIVATE "/permissive-")
elseif(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" OR CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang" OR CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
    target_compile_options(qflib PRIVATE -fPIC)
endif()
