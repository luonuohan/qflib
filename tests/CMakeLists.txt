find_package(Catch2 3 REQUIRED)

add_executable(qflib_tests
    test_utils.cpp
    test_sptrmap.cpp
    test_errorfunction.cpp
    test_normaldistribution.cpp
    test_interpolation1d.cpp
    test_polynomial.cpp
    test_meanvarcalculator.cpp
    test_products.cpp
)

target_include_directories(qflib_tests PRIVATE
    ${CMAKE_SOURCE_DIR}        # for #include <qflib/...>
)

if(APPLE)
    target_include_directories(qflib_tests PRIVATE
        ${ARMADILLO_PREFIX}/include
    )
else()
    target_include_directories(qflib_tests PRIVATE
        ${THIRDPARTY_DIRECTORY}/armadillo-${ARMA_VERSION}/include
    )
endif()

target_link_libraries(qflib_tests PRIVATE qflib Catch2::Catch2WithMain)

include(CTest)
include(Catch)
catch_discover_tests(qflib_tests)
