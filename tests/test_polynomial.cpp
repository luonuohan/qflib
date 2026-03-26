#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <qflib/math/optim/polyfunc.hpp>

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

TEST_CASE("Polynomial constant", "[polynomial]")
{
    // p(x) = 5
    qf::Vector c = {5.0};
    qf::Polynomial p(c);
    REQUIRE_THAT(p(0.0), WithinAbs(5.0, 1e-15));
    REQUIRE_THAT(p(100.0), WithinAbs(5.0, 1e-15));
    REQUIRE_THAT(p(-7.3), WithinAbs(5.0, 1e-15));
}

TEST_CASE("Polynomial linear", "[polynomial]")
{
    // p(x) = 2 + 3x
    qf::Vector c = {2.0, 3.0};
    qf::Polynomial p(c);
    REQUIRE_THAT(p(0.0), WithinAbs(2.0, 1e-15));
    REQUIRE_THAT(p(1.0), WithinAbs(5.0, 1e-15));
    REQUIRE_THAT(p(-1.0), WithinAbs(-1.0, 1e-15));
    REQUIRE_THAT(p(2.5), WithinAbs(2.0 + 3.0 * 2.5, 1e-12));
}

TEST_CASE("Polynomial quadratic", "[polynomial]")
{
    // p(x) = 1 - 2x + 3x^2
    qf::Vector c = {1.0, -2.0, 3.0};
    qf::Polynomial p(c);
    REQUIRE_THAT(p(0.0), WithinAbs(1.0, 1e-15));
    REQUIRE_THAT(p(1.0), WithinAbs(2.0, 1e-12));     // 1 - 2 + 3 = 2
    REQUIRE_THAT(p(-1.0), WithinAbs(6.0, 1e-12));    // 1 + 2 + 3 = 6
    REQUIRE_THAT(p(2.0), WithinAbs(9.0, 1e-12));     // 1 - 4 + 12 = 9
}

TEST_CASE("Polynomial cubic", "[polynomial]")
{
    // p(x) = x^3 => coeffs = {0, 0, 0, 1}
    qf::Vector c = {0.0, 0.0, 0.0, 1.0};
    qf::Polynomial p(c);
    REQUIRE_THAT(p(0.0), WithinAbs(0.0, 1e-15));
    REQUIRE_THAT(p(2.0), WithinAbs(8.0, 1e-12));
    REQUIRE_THAT(p(-3.0), WithinAbs(-27.0, 1e-12));
}

TEST_CASE("Polynomial empty coefficients throw", "[polynomial]")
{
    qf::Vector c;  // empty
    REQUIRE_THROWS_AS(qf::Polynomial(c), qf::Exception);
}
