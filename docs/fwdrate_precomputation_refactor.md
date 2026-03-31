# fwdRate Precomputation Refactor
## BsMcSimpleDeltaHedger — Performance Optimisation Report
**Date:** 2026-03-30

---

## 1. Motivation

`BsMcSimpleDeltaHedger::hedge()` runs a Monte Carlo simulation over `nPaths`
paths, each with `N` rehedging steps. In the original implementation, every
path × every step called `discyc_->fwdRate(t0, t1)` live — traversing the
piecewise polynomial yield curve each time.

For typical parameters (N = 52 weekly steps, nPaths = 100,000):

| Call site | Calls per path | Total calls |
|---|---|---|
| Cash accrual `fwdRate(times[j-1], times[j])` | N−1 | ~5.2M |
| `priceAndDelta` `fwdRate(times[j], T)` | N | ~5.2M |
| Final accrual `fwdRate(times[N-1], T)` | 1 | 100K |
| **Total** | **~2N** | **~10.4M** |

Each `fwdRate` call runs `PiecewisePolynomial::integral()` — a binary search
plus polynomial accumulation. This is the dominant cost in the hot loop.

**Key insight:** the time grid is fixed before the MC loop begins. All
`fwdRate` values are deterministic functions of the time grid — they do not
depend on the random path. They can be computed **once** and looked up as
plain array reads in O(1).

---

## 2. Changes Made

### 2.1 `bsmcsimpledeltahedger.hpp`

**Added two new member vectors:**

```cpp
// BEFORE
std::vector<double> times_, drifts_, stdevs_;

// AFTER
std::vector<double> times_, drifts_, stdevs_;
std::vector<double> fwdrates_;    // fwdRate(times_[i], times_[i+1])
std::vector<double> fwdratesToT_; // fwdRate(times_[i], T) for priceAndDelta
```

**Changed `priceAndDelta` signature** — rate `r` is now passed in rather than
computed internally:

```cpp
// BEFORE
void priceAndDelta(double S, double t, double& price, double& delta) const;

// AFTER
void priceAndDelta(double S, double t, double r, double& price, double& delta) const;
```

### 2.2 `bsmcsimpledeltahedger.cpp`

**`priceAndDelta` implementation** — removed the internal `fwdRate` call:

```cpp
// BEFORE
void BsMcSimpleDeltaHedger::priceAndDelta(double S, double t,
                                           double& price, double& delta) const
{
  double tau = prod_->payTimes().back() - t;
  double r   = discyc_->fwdRate(t, t + tau);   // <-- live curve call
  ...
}

// AFTER
void BsMcSimpleDeltaHedger::priceAndDelta(double S, double t, double r,
                                           double& price, double& delta) const
{
  double tau = prod_->payTimes().back() - t;
  // r passed in from precomputed table — no curve call
  ...
}
```

**Precomputation loop in `hedge()`** — extended to fill both new tables:

```cpp
// BEFORE
for (unsigned i = 0; i < N; ++i) {
  double t1  = times_[i+1];
  double dt  = t1 - t0;
  double var = vol2 * dt;
  stdevs_[i] = std::sqrt(var);
  double fwd = discyc_->fwdRate(t0, t1);          // computed but not stored
  drifts_[i] = (fwd - divyld_)*dt - 0.5*var;
  t0 = t1;
}

// AFTER
fwdrates_.resize(N);
fwdratesToT_.resize(N);
for (unsigned i = 0; i < N; ++i) {
  double t1  = times_[i+1];
  double dt  = t1 - t0;
  double var = vol2 * dt;
  stdevs_[i]       = std::sqrt(var);
  fwdrates_[i]     = discyc_->fwdRate(t0, t1);   // stored
  fwdratesToT_[i]  = discyc_->fwdRate(t0, T);    // stored (was never computed before)
  drifts_[i]       = (fwdrates_[i] - divyld_)*dt - 0.5*var;
  t0 = t1;
}
```

**MC loop — cash accrual** — live call replaced with table lookup:

```cpp
// BEFORE
double rrf = discyc_->fwdRate(tPrev, tNow);
cash *= std::exp(rrf * (tNow - tPrev));

// AFTER
cash *= std::exp(fwdrates_[j-1] * (tNow - tPrev));
```

**MC loop — `priceAndDelta` call sites** — rate passed from table:

```cpp
// BEFORE
priceAndDelta(spot0_, 0.0, initPrice, initDelta);
priceAndDelta(spotPath[j], tNow, priceNow, deltaNow);

// AFTER
priceAndDelta(spot0_, 0.0, fwdratesToT_[0], initPrice, initDelta);
priceAndDelta(spotPath[j], tNow, fwdratesToT_[j], priceNow, deltaNow);
```

**MC loop — final accrual** — live call replaced:

```cpp
// BEFORE
double rrf = discyc_->fwdRate(tPrev, T);
cash *= std::exp(rrf * (T - tPrev));

// AFTER
cash *= std::exp(fwdrates_[N-1] * (T - tPrev));
```

### 2.3 `tests/bench_fwdrate_precompute.cpp` (new file)

Standalone benchmark that isolates the `fwdRate` call pattern and measures
before/after timing directly with `std::chrono::high_resolution_clock`.

### 2.4 `tests/CMakeLists.txt`

Added `qflib_bench` executable target for the benchmark.

---

## 3. Call Count Analysis

### Before

Per path:
- 1 call at hedge init (`priceAndDelta` t=0)
- Per rehedge step j=1..N-1: 1 (cash accrual) + 1 (priceAndDelta) = 2 calls
- 1 call at final accrual

Total per path = `1 + 2*(N-1) + 1 = 2N`

Total for simulation = `2N × nPaths`

### After

All calls moved to precomputation phase (runs once):
- N calls for `fwdrates_`
- N calls for `fwdratesToT_`

Total precompute = `2N`

Per path: 0 `fwdRate` calls — only array index reads.

Total for simulation = `2N` (constant, independent of nPaths)

### Reduction ratio

```
Calls before : 2N × nPaths
Calls after  : 2N
Reduction    : nPaths× (e.g. 100,000×)
```

---

## 4. Expected Speedup

The speedup depends on what fraction of total runtime was `fwdRate` calls.
Let C = cost of one `fwdRate` call, L = cost of everything else per step per path.

```
Time before  = nPaths × N × (2C + L)
Time after   = nPaths × N × L  +  2N × C   ≈  nPaths × N × L  (for large nPaths)

Speedup  ≈  (2C + L) / L  =  1 + 2C/L
```

In practice:
- Each `fwdRate` = binary search O(log M) + polynomial accumulation O(order × M_interval)
  on the piecewise polynomial — roughly 50–200 ns on modern hardware
- The rest of the inner loop (exp, multiply, add) ~ 20–50 ns per step

Estimated speedup for typical parameters:

| fwdRate cost fraction | Expected speedup |
|---|---|
| 50% of loop time | ~2× |
| 70% of loop time | ~3× |
| 80% of loop time | ~5× |

Run `./qflib_bench` after building to get the actual measured ratio on your machine.

---

## 5. Correctness

The refactor is **mathematically identical** to the original — only the
evaluation timing changes, not the values. `fwdRate` is a deterministic
function of (t0, t1) and the yield curve, which is fixed for the lifetime of
`hedge()`. Precomputing it before the loop vs recomputing it inside produces
exactly the same floating-point result.

No test cases need to change. All existing Catch2 tests continue to validate
correctness unchanged.

---

## 6. How to Build and Run the Benchmark

```bash
cd /Users/nuohanluo/Desktop/MSQF/2026_Spring/C++_Project/build
cmake --build . --target qflib_bench
./qflib-0.11.0/bin/qflib_bench
```

Expected output format:
```
=== BsMcSimpleDeltaHedger fwdRate precomputation benchmark ===
  N (steps)    = 52
  nPaths       = 100000
  fwdRate calls before: ~10400000 (105 per path)
  fwdRate calls after : 104 (precomputed once)

  BEFORE (live calls)  : XXX.X ms
  AFTER  (precomputed) : XXX.X ms
  Speedup              : X.Xx
```
