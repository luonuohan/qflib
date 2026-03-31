#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <qflib/math/stats/normaldistribution.hpp>

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

TEST_CASE("Standard normal pdf", "[normaldist]")
{
    qf::NormalDistribution N;  // standard normal (mu=0, sigma=1)

    // pdf(0) = 1/sqrt(2*pi) ~ 0.3989422804
    REQUIRE_THAT(N.pdf(0.0), WithinRel(M_1_SQRT2PI, 1e-12));

    // pdf is symmetric: pdf(-x) == pdf(x)
    REQUIRE_THAT(N.pdf(-1.5), WithinRel(N.pdf(1.5), 1e-14));

    // pdf(1) = (1/sqrt(2*pi)) * exp(-0.5)
    double expected = M_1_SQRT2PI * std::exp(-0.5);
    REQUIRE_THAT(N.pdf(1.0), WithinRel(expected, 1e-12));
}

TEST_CASE("Standard normal cdf", "[normaldist]")
{
    qf::NormalDistribution N;

    // cdf(0) = 0.5
    REQUIRE_THAT(N.cdf(0.0), WithinRel(0.5, 1e-12));

    // cdf(-x) + cdf(x) = 1 (symmetry)
    REQUIRE_THAT(N.cdf(-1.0) + N.cdf(1.0), WithinRel(1.0, 1e-12));
    REQUIRE_THAT(N.cdf(-2.0) + N.cdf(2.0), WithinRel(1.0, 1e-12));

    // known values
    REQUIRE_THAT(N.cdf(1.0), WithinRel(0.8413447460685429, 1e-8));
    REQUIRE_THAT(N.cdf(-1.0), WithinRel(0.1586552539314571, 1e-8));
    REQUIRE_THAT(N.cdf(1.96), WithinRel(0.9750021048517796, 1e-6));
}

TEST_CASE("Standard normal invcdf", "[normaldist]")
{
    qf::NormalDistribution N;

    // invcdf(0.5) = 0
    REQUIRE_THAT(N.invcdf(0.5), WithinAbs(0.0, 1e-10));

    // roundtrip: invcdf(cdf(x)) == x
    // After inverfc bug fix, precision improved from ~1e-3 to near machine epsilon
    for (double x : {-3.0, -2.0, -1.0, 0.0, 0.5, 1.5, 2.0, 3.0}) {
        double p = N.cdf(x);
        REQUIRE_THAT(N.invcdf(p), WithinAbs(x, 1e-10));
    }

    // known quantiles (true z_0.975 = 1.959963984540054...)
    REQUIRE_THAT(N.invcdf(0.975), WithinRel(1.959963984540054, 1e-10));
    REQUIRE_THAT(N.invcdf(0.025), WithinRel(-1.959963984540054, 1e-10));
}

TEST_CASE("Non-standard normal distribution", "[normaldist]")
{
    double mu = 5.0, sigma = 2.0;
    qf::NormalDistribution N(mu, sigma);

    // cdf at the mean is 0.5
    REQUIRE_THAT(N.cdf(mu), WithinRel(0.5, 1e-12));

    // pdf at the mean = 1/(sigma * sqrt(2*pi))
    double expected_pdf = M_1_SQRT2PI / sigma;
    REQUIRE_THAT(N.pdf(mu), WithinRel(expected_pdf, 1e-12));

    // invcdf(0.5) = mu
    REQUIRE_THAT(N.invcdf(0.5), WithinAbs(mu, 1e-10));

    // cdf roundtrip
    double x = 6.0;
    REQUIRE_THAT(N.invcdf(N.cdf(x)), WithinAbs(x, 1e-10));
}

TEST_CASE("NormalDistribution sigma must be positive", "[normaldist]")
{
    REQUIRE_THROWS_AS(qf::NormalDistribution(0.0, 0.0), qf::Exception);
    REQUIRE_THROWS_AS(qf::NormalDistribution(0.0, -1.0), qf::Exception);
}
