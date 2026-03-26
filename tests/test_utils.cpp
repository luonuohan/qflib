#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <qflib/utils.hpp>

using Catch::Matchers::WithinRel;
using Catch::Matchers::WithinAbs;

TEST_CASE("toContCmpd converts periodic to continuous compounding", "[utils]")
{
    // Annual compounding: ln(1 + r)
    REQUIRE_THAT(qf::toContCmpd(0.05, 1), WithinRel(std::log(1.05), 1e-12));

    // Semi-annual: 2*ln(1 + r/2)
    // For r=0.05, annfreq=2: ln((1+0.025)^2) = 2*ln(1.025)
    REQUIRE_THAT(qf::toContCmpd(0.05, 2), WithinRel(2.0 * std::log(1.025), 1e-12));

    // Quarterly
    REQUIRE_THAT(qf::toContCmpd(0.08, 4), WithinRel(4.0 * std::log(1.02), 1e-12));

    // Monthly
    double expected = 12.0 * std::log(1.0 + 0.12 / 12.0);
    REQUIRE_THAT(qf::toContCmpd(0.12, 12), WithinRel(expected, 1e-12));

    // Zero rate should give zero
    REQUIRE_THAT(qf::toContCmpd(0.0, 1), WithinAbs(0.0, 1e-15));
}

TEST_CASE("fromContCmpd converts continuous to periodic compounding", "[utils]")
{
    // Annual: (exp(r) - 1)
    REQUIRE_THAT(qf::fromContCmpd(0.05, 1), WithinRel(std::exp(0.05) - 1.0, 1e-12));

    // Semi-annual: 2*(exp(r/2) - 1)
    REQUIRE_THAT(qf::fromContCmpd(0.05, 2), WithinRel(2.0 * (std::exp(0.025) - 1.0), 1e-12));

    // Zero rate
    REQUIRE_THAT(qf::fromContCmpd(0.0, 1), WithinAbs(0.0, 1e-15));
}

TEST_CASE("toContCmpd and fromContCmpd are inverses", "[utils]")
{
    // roundtrip: fromContCmpd(toContCmpd(r, n), n) == r
    double rate = 0.06;
    for (size_t freq : {1, 2, 4, 12}) {
        double cont = qf::toContCmpd(rate, freq);
        double back = qf::fromContCmpd(cont, freq);
        REQUIRE_THAT(back, WithinRel(rate, 1e-12));
    }
}
