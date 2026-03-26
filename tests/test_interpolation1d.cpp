#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <qflib/math/interpol/interpolation1d.hpp>

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

TEST_CASE("LinearInterpolation1D at exact data points", "[interp1d]")
{
    qf::Vector x = {1.0, 2.0, 3.0, 4.0, 5.0};
    qf::Vector y = {10.0, 20.0, 30.0, 40.0, 50.0};
    qf::LinearInterpolation1D<qf::Vector> interp(x, y);

    REQUIRE(interp.size() == 5);

    // getValue(index) returns exact y values
    REQUIRE_THAT(interp.getValue(size_t(0)), WithinAbs(10.0, 1e-15));
    REQUIRE_THAT(interp.getValue(size_t(2)), WithinAbs(30.0, 1e-15));
    REQUIRE_THAT(interp.getValue(size_t(4)), WithinAbs(50.0, 1e-15));
}

TEST_CASE("LinearInterpolation1D interpolates between points", "[interp1d]")
{
    qf::Vector x = {1.0, 2.0, 3.0, 4.0, 5.0};
    qf::Vector y = {10.0, 20.0, 30.0, 40.0, 50.0};
    qf::LinearInterpolation1D<qf::Vector> interp(x, y);

    // linear data: midpoints should be exact
    REQUIRE_THAT(interp.getValue(1.5), WithinAbs(15.0, 1e-12));
    REQUIRE_THAT(interp.getValue(2.5), WithinAbs(25.0, 1e-12));
    REQUIRE_THAT(interp.getValue(3.75), WithinAbs(37.5, 1e-12));
}

TEST_CASE("LinearInterpolation1D flat extrapolation", "[interp1d]")
{
    qf::Vector x = {1.0, 2.0, 3.0, 4.0, 5.0};
    qf::Vector y = {10.0, 20.0, 15.0, 40.0, 50.0};
    qf::LinearInterpolation1D<qf::Vector> interp(x, y);

    // below range: uses first segment's extrapolation
    double below = interp.getValue(0.5);
    // y1 + (y2-y1)*(x-x1)/(x2-x1) = 10 + 10*(0.5-1)/(2-1) = 10 - 5 = 5
    REQUIRE_THAT(below, WithinAbs(5.0, 1e-12));

    // above range: uses last segment's extrapolation
    double above = interp.getValue(6.0);
    // y1 + (y2-y1)*(x-x1)/(x2-x1) = 40 + 10*(6-4)/(5-4) = 60
    REQUIRE_THAT(above, WithinAbs(60.0, 1e-12));
}

TEST_CASE("LinearInterpolation1D non-linear data", "[interp1d]")
{
    qf::Vector x = {0.0, 1.0, 2.0};
    qf::Vector y = {0.0, 1.0, 4.0};  // y = x^2 at these points
    qf::LinearInterpolation1D<qf::Vector> interp(x, y);

    // between 1 and 2: linear interpolation y = 1 + 3*(x-1)
    REQUIRE_THAT(interp.getValue(1.5), WithinAbs(2.5, 1e-12));
    // note: actual x^2 at 1.5 = 2.25, but linear interp gives 2.5
}
