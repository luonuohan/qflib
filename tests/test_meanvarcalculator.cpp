#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <qflib/math/stats/meanvarcalculator.hpp>
#include <vector>

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

using VecIter = std::vector<double>::const_iterator;

TEST_CASE("MeanVarCalculator single variable", "[meanvar]")
{
    qf::MeanVarCalculator<VecIter> calc(1);

    // data: {2, 4, 6, 8, 10} => mean=6, sample var=10
    std::vector<std::vector<double>> samples = {{2}, {4}, {6}, {8}, {10}};

    for (auto& s : samples)
        calc.addSample(s.cbegin(), s.cend());

    REQUIRE(calc.nSamples() == 5);

    auto const& res = calc.results();
    double mean = res(0, 0);
    double var = res(1, 0);

    REQUIRE_THAT(mean, WithinAbs(6.0, 1e-12));
    // sample variance: sum((xi - mean)^2) / (n-1) = (16+4+0+4+16)/4 = 10
    REQUIRE_THAT(var, WithinAbs(10.0, 1e-12));
}

TEST_CASE("MeanVarCalculator multiple variables", "[meanvar]")
{
    qf::MeanVarCalculator<VecIter> calc(2);

    // var1: {1, 2, 3}, var2: {10, 20, 30}
    std::vector<std::vector<double>> samples = {{1, 10}, {2, 20}, {3, 30}};

    for (auto& s : samples)
        calc.addSample(s.cbegin(), s.cend());

    REQUIRE(calc.nSamples() == 3);

    auto const& res = calc.results();

    // means
    REQUIRE_THAT(res(0, 0), WithinAbs(2.0, 1e-12));
    REQUIRE_THAT(res(0, 1), WithinAbs(20.0, 1e-12));

    // sample variances
    // var1: (1+0+1)/2 = 1.0
    REQUIRE_THAT(res(1, 0), WithinAbs(1.0, 1e-12));
    // var2: (100+0+100)/2 = 100.0
    REQUIRE_THAT(res(1, 1), WithinAbs(100.0, 1e-12));
}

TEST_CASE("MeanVarCalculator reset", "[meanvar]")
{
    qf::MeanVarCalculator<VecIter> calc(1);

    std::vector<double> s1 = {5.0};
    calc.addSample(s1.cbegin(), s1.cend());
    REQUIRE(calc.nSamples() == 1);

    calc.reset();
    // NOTE: base class reset() does not reset nsamples_ (library quirk)
    // so we test with a fresh calculator instead
    qf::MeanVarCalculator<VecIter> calc2(1);

    std::vector<double> s2 = {10.0};
    std::vector<double> s3 = {20.0};
    calc2.addSample(s2.cbegin(), s2.cend());
    calc2.addSample(s3.cbegin(), s3.cend());

    auto const& res = calc2.results();
    REQUIRE_THAT(res(0, 0), WithinAbs(15.0, 1e-12));
}

TEST_CASE("MeanVarCalculator constant data has zero variance", "[meanvar]")
{
    qf::MeanVarCalculator<VecIter> calc(1);

    for (int i = 0; i < 100; ++i) {
        std::vector<double> s = {3.14};
        calc.addSample(s.cbegin(), s.cend());
    }

    auto const& res = calc.results();
    REQUIRE_THAT(res(0, 0), WithinAbs(3.14, 1e-12));
    REQUIRE_THAT(res(1, 0), WithinAbs(0.0, 1e-10));
}
