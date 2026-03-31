/**
@file  bench_fwdrate_precompute.cpp
@brief Benchmark: live fwdRate() calls vs precomputed table lookup in the MC loop.

Simulates the work pattern of BsMcSimpleDeltaHedger before and after the
precomputation refactor, without running the full hedger.

Build: linked into qflib_bench via CMakeLists.
Run:   ./qflib_bench
*/

#include <qflib/market/yieldcurve.hpp>
#include <chrono>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

using Clock = std::chrono::high_resolution_clock;
using Ms    = std::chrono::duration<double, std::milli>;

// ---------------------------------------------------------------------------
// Build a realistic piecewise yield curve
// ---------------------------------------------------------------------------
static qf::SPtrYieldCurve makeCurve()
{
  std::vector<double> mats  = {0.25, 0.5, 1.0, 2.0, 3.0, 5.0, 7.0, 10.0};
  std::vector<double> rates = {0.04, 0.042, 0.045, 0.047, 0.049, 0.05, 0.051, 0.052};
  return std::make_shared<qf::YieldCurve>(
      mats.begin(), mats.end(),
      rates.begin(), rates.end(),
      qf::YieldCurve::InputType::SPOTRATE);
}

// ---------------------------------------------------------------------------
// BEFORE: live fwdRate() call inside the MC loop (old pattern)
// ---------------------------------------------------------------------------
static double benchLive(qf::YieldCurve const& yc,
                        std::vector<double> const& times,
                        unsigned long nPaths)
{
  unsigned N = (unsigned)times.size() - 1;
  double T   = times[N];
  double acc = 0.0;

  for (unsigned long p = 0; p < nPaths; ++p) {
    // cash accrual per step + priceAndDelta rate lookup (2 calls per step)
    for (unsigned j = 0; j < N; ++j) {
      double r1 = yc.fwdRate(times[j], times[j+1]);   // cash accrual
      double r2 = yc.fwdRate(times[j], T);             // priceAndDelta
      acc += r1 + r2;
    }
    // final accrual
    acc += yc.fwdRate(times[N-1], T);
  }
  return acc; // prevent dead-code elimination
}

// ---------------------------------------------------------------------------
// AFTER: precomputed tables, MC loop does array lookups only
// ---------------------------------------------------------------------------
static double benchPrecomputed(std::vector<double> const& fwdrates,
                               std::vector<double> const& fwdratesToT,
                               unsigned long nPaths)
{
  unsigned N = (unsigned)fwdrates.size();
  double acc = 0.0;

  for (unsigned long p = 0; p < nPaths; ++p) {
    for (unsigned j = 0; j < N; ++j) {
      double r1 = fwdrates[j];     // cash accrual
      double r2 = fwdratesToT[j];  // priceAndDelta
      acc += r1 + r2;
    }
    acc += fwdrates[N-1];          // final accrual
  }
  return acc;
}

// ---------------------------------------------------------------------------
int main()
{
  auto spyc = makeCurve();
  qf::YieldCurve const& yc = *spyc;

  // Typical hedger parameters: hedgeFreq=52 (weekly), T=1yr -> N=52 steps
  const unsigned    hedgeFreq = 52;
  const double      T         = 1.0;
  const unsigned    N         = (unsigned)std::ceil(hedgeFreq * T);
  const unsigned long nPaths  = 100000UL;

  // Build time grid
  std::vector<double> times(N+1);
  for (unsigned i = 0; i <= N; ++i) times[i] = double(i) / hedgeFreq;

  // Precompute tables (the one-time cost before the MC loop)
  std::vector<double> fwdrates(N), fwdratesToT(N);
  for (unsigned i = 0; i < N; ++i) {
    fwdrates[i]    = yc.fwdRate(times[i], times[i+1]);
    fwdratesToT[i] = yc.fwdRate(times[i], T);
  }

  std::cout << "=== BsMcSimpleDeltaHedger fwdRate precomputation benchmark ===\n";
  std::cout << "  N (steps)    = " << N        << "\n";
  std::cout << "  nPaths       = " << nPaths   << "\n";
  std::cout << "  fwdRate calls before: ~" << (2*N + 1) * nPaths
            << " (" << 2*N+1 << " per path)\n";
  std::cout << "  fwdRate calls after : " << 2*N
            << " (precomputed once, then " << 2*N << " table reads per path)\n\n";

  // --- BEFORE ---
  double dummy1 = 0;
  auto t0 = Clock::now();
  dummy1 = benchLive(yc, times, nPaths);
  auto t1 = Clock::now();
  double ms_live = Ms(t1 - t0).count();

  // --- AFTER ---
  double dummy2 = 0;
  auto t2 = Clock::now();
  dummy2 = benchPrecomputed(fwdrates, fwdratesToT, nPaths);
  auto t3 = Clock::now();
  double ms_precomp = Ms(t3 - t2).count();

  // dummy use to prevent dead-code elimination
  if (dummy1 == dummy2 + 1e300) std::cout << "(should not print)\n";

  std::cout << "  BEFORE (live calls)  : " << ms_live    << " ms\n";
  std::cout << "  AFTER  (precomputed) : " << ms_precomp << " ms\n";
  std::cout << "  Speedup              : " << ms_live / ms_precomp << "x\n";

  return 0;
}
