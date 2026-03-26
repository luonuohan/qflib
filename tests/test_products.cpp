#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include <qflib/products/europeancallput.hpp>
#include <qflib/products/digitalcallput.hpp>

using Catch::Matchers::WithinAbs;

// --- European Call/Put ---

TEST_CASE("European call payoff", "[products]")
{
    double K = 100.0, T = 1.0;
    qf::EuropeanCallPut call(1, K, T);

    REQUIRE(call.nAssets() == 1);
    REQUIRE(call.fixTimes().size() == 1);
    REQUIRE(call.payTimes().size() == 1);

    SECTION("ITM: S > K") {
        qf::Matrix path(1, 1);
        path(0, 0) = 120.0;
        call.eval(path);
        REQUIRE_THAT(call.payAmounts()[0], WithinAbs(20.0, 1e-12));
    }

    SECTION("ATM: S == K") {
        qf::Matrix path(1, 1);
        path(0, 0) = 100.0;
        call.eval(path);
        REQUIRE_THAT(call.payAmounts()[0], WithinAbs(0.0, 1e-12));
    }

    SECTION("OTM: S < K") {
        qf::Matrix path(1, 1);
        path(0, 0) = 80.0;
        call.eval(path);
        REQUIRE_THAT(call.payAmounts()[0], WithinAbs(0.0, 1e-12));
    }
}

TEST_CASE("European put payoff", "[products]")
{
    double K = 100.0, T = 1.0;
    qf::EuropeanCallPut put(-1, K, T);

    SECTION("ITM: S < K") {
        qf::Matrix path(1, 1);
        path(0, 0) = 80.0;
        put.eval(path);
        REQUIRE_THAT(put.payAmounts()[0], WithinAbs(20.0, 1e-12));
    }

    SECTION("OTM: S > K") {
        qf::Matrix path(1, 1);
        path(0, 0) = 120.0;
        put.eval(path);
        REQUIRE_THAT(put.payAmounts()[0], WithinAbs(0.0, 1e-12));
    }

    SECTION("ATM: S == K") {
        qf::Matrix path(1, 1);
        path(0, 0) = 100.0;
        put.eval(path);
        REQUIRE_THAT(put.payAmounts()[0], WithinAbs(0.0, 1e-12));
    }
}

TEST_CASE("European call-put parity on payoff", "[products]")
{
    double K = 100.0, T = 0.5;
    qf::EuropeanCallPut call(1, K, T);
    qf::EuropeanCallPut put(-1, K, T);

    for (double S : {80.0, 90.0, 100.0, 110.0, 120.0}) {
        qf::Matrix path(1, 1);
        path(0, 0) = S;
        call.eval(path);
        put.eval(path);
        // call - put = S - K (at expiry)
        double diff = call.payAmounts()[0] - put.payAmounts()[0];
        REQUIRE_THAT(diff, WithinAbs(S - K, 1e-12));
    }
}

TEST_CASE("European option PDE eval interface", "[products]")
{
    double K = 100.0, T = 1.0;
    qf::EuropeanCallPut call(1, K, T);

    qf::Vector spots = {110.0};
    call.eval(size_t(0), spots, 0.0);
    REQUIRE_THAT(call.payAmounts()[0], WithinAbs(10.0, 1e-12));
}

TEST_CASE("European option input validation", "[products]")
{
    REQUIRE_THROWS_AS(qf::EuropeanCallPut(0, 100.0, 1.0), qf::Exception);   // bad payoff type
    REQUIRE_THROWS_AS(qf::EuropeanCallPut(1, -10.0, 1.0), qf::Exception);   // bad strike
    REQUIRE_THROWS_AS(qf::EuropeanCallPut(1, 100.0, -1.0), qf::Exception);  // bad expiry
}

// --- Digital Call/Put ---

TEST_CASE("Digital call payoff", "[products]")
{
    double K = 100.0, T = 1.0;
    qf::DigitalCallPut dcall(1, K, T);

    SECTION("S >= K pays 1") {
        qf::Matrix path(1, 1);
        path(0, 0) = 100.0;
        dcall.eval(path);
        REQUIRE_THAT(dcall.payAmounts()[0], WithinAbs(1.0, 1e-15));

        path(0, 0) = 150.0;
        dcall.eval(path);
        REQUIRE_THAT(dcall.payAmounts()[0], WithinAbs(1.0, 1e-15));
    }

    SECTION("S < K pays 0") {
        qf::Matrix path(1, 1);
        path(0, 0) = 99.99;
        dcall.eval(path);
        REQUIRE_THAT(dcall.payAmounts()[0], WithinAbs(0.0, 1e-15));
    }
}

TEST_CASE("Digital put payoff", "[products]")
{
    double K = 100.0, T = 1.0;
    qf::DigitalCallPut dput(-1, K, T);

    SECTION("S < K pays 1") {
        qf::Matrix path(1, 1);
        path(0, 0) = 99.0;
        dput.eval(path);
        REQUIRE_THAT(dput.payAmounts()[0], WithinAbs(1.0, 1e-15));
    }

    SECTION("S >= K pays 0") {
        qf::Matrix path(1, 1);
        path(0, 0) = 100.0;
        dput.eval(path);
        REQUIRE_THAT(dput.payAmounts()[0], WithinAbs(0.0, 1e-15));
    }
}

TEST_CASE("Digital call + put = 1", "[products]")
{
    double K = 100.0, T = 1.0;
    qf::DigitalCallPut dcall(1, K, T);
    qf::DigitalCallPut dput(-1, K, T);

    for (double S : {80.0, 99.99, 100.0, 100.01, 120.0}) {
        qf::Matrix path(1, 1);
        path(0, 0) = S;
        dcall.eval(path);
        dput.eval(path);
        REQUIRE_THAT(dcall.payAmounts()[0] + dput.payAmounts()[0], WithinAbs(1.0, 1e-15));
    }
}

TEST_CASE("Digital option input validation", "[products]")
{
    REQUIRE_THROWS_AS(qf::DigitalCallPut(2, 100.0, 1.0), qf::Exception);
    REQUIRE_THROWS_AS(qf::DigitalCallPut(1, 0.0, 1.0), qf::Exception);
    REQUIRE_THROWS_AS(qf::DigitalCallPut(1, 100.0, 0.0), qf::Exception);
}
