#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <qflib/math/stats/errorfunction.hpp>

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;
using qf::ErrorFunction;

TEST_CASE("erf known values", "[errorfunction]")
{
    // erf(0) = 0
    REQUIRE_THAT(ErrorFunction::erf(0.0), WithinAbs(0.0, 1e-15));

    // erf is odd: erf(-x) = -erf(x)
    REQUIRE_THAT(ErrorFunction::erf(-1.0), WithinRel(-ErrorFunction::erf(1.0), 1e-14));
    REQUIRE_THAT(ErrorFunction::erf(-0.5), WithinRel(-ErrorFunction::erf(0.5), 1e-14));

    // known values (from standard tables)
    REQUIRE_THAT(ErrorFunction::erf(1.0), WithinRel(0.8427007929497149, 1e-10));
    REQUIRE_THAT(ErrorFunction::erf(2.0), WithinRel(0.9953222650189527, 1e-10));

    // erf(large) -> 1
    REQUIRE_THAT(ErrorFunction::erf(5.0), WithinRel(1.0, 1e-8));
}

TEST_CASE("erfc known values", "[errorfunction]")
{
    // erfc(x) = 1 - erf(x)
    REQUIRE_THAT(ErrorFunction::erfc(0.0), WithinRel(1.0, 1e-14));
    REQUIRE_THAT(ErrorFunction::erfc(1.0), WithinRel(1.0 - ErrorFunction::erf(1.0), 1e-14));

    // erfc(large) -> 0
    REQUIRE_THAT(ErrorFunction::erfc(5.0), WithinAbs(0.0, 1e-8));

    // erfc(-large) -> 2
    REQUIRE_THAT(ErrorFunction::erfc(-5.0), WithinRel(2.0, 1e-8));
}

TEST_CASE("inverf roundtrip", "[errorfunction]")
{
    // inverf(erf(x)) == x
    // The Numerical Recipes implementation has ~1e-3 precision on the inverse
    for (double x : {-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0}) {
        double y = ErrorFunction::erf(x);
        REQUIRE_THAT(ErrorFunction::inverf(y), WithinAbs(x, 2e-3));
    }
}

TEST_CASE("inverfc roundtrip", "[errorfunction]")
{
    // inverfc(erfc(x)) == x
    for (double x : {-2.0, -1.0, 0.0, 0.5, 1.0, 2.0}) {
        double y = ErrorFunction::erfc(x);
        REQUIRE_THAT(ErrorFunction::inverfc(y), WithinAbs(x, 2e-3));
    }
}
