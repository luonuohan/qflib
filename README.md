# qflib — Quantitative Finance Library in C++

A modular C++ library for quantitative finance featuring Monte Carlo pricing, delta hedging simulation, yield curve construction, PDE solvers, and numerical methods — with a full Python front-end via the C API.

---

## Table of contents

- [Overview](#overview)
- [Project structure](#project-structure)
- [Key components](#key-components)
- [Build setup](#build-setup)
- [Python installation](#python-installation)
- [Run tests](#run-tests)
- [Python usage example](#python-usage-example)
- [Design decisions and optimizations](#design-decisions-and-optimizations)

---

## Overview

qflib is a self-contained quantitative finance engine written in C++20. The core library (`qflib/`) implements pricing models, numerical methods, and Monte Carlo infrastructure. The Python interface (`pyqflib/`) wraps the C++ core using the Python C API — no Pybind11 or SWIG — giving full control over type conversions and object lifetime.

**Project 1: Delta Hedging Simulation**

The primary application is a Monte Carlo delta hedging study under the Black-Scholes model. The engine generates stock price paths and dynamically delta hedges European and digital options on each path. The P&L from hedging and final settlement is collected into pluggable statistics calculators (`MeanVarCalculator`, `HistogramCalculator`) and returned as the mean, standard deviation, and histogram of the hedge P&L distribution.

Two volatilities are distinguished:
- `simVol` — the true volatility driving the stock (unknown to the trader)
- `markVol` — the volatility used for pricing and delta computation (what the trader controls)

When `markVol == simVol`, the mean P&L is zero and the standard deviation scales as `O(hedgeFreq^{-0.5})`. When they diverge, the distribution shifts and widens in a payoff-type-dependent way.

---

## Project structure

```
qflib-0.11.0/
│
├── CMakeLists.txt               # Top-level CMake: configures all three subdirectories
│
├── qflib/                       # Core C++ static library
│   ├── CMakeLists.txt
│   ├── defines.hpp              # Version macros, math constants, namespace macros
│   ├── exception.hpp            # qf::Exception and QF_ASSERT macro
│   ├── sptrmap.hpp              # Named object registry (SPtrMap<T>)
│   ├── sptr.hpp                 # Shared pointer type aliases
│   ├── utils.hpp                # Compounding conversions (toContCmpd, fromContCmpd)
│   │
│   ├── math/
│   │   ├── matrix.hpp           # Armadillo type aliases (Vector, Matrix)
│   │   ├── interpol/
│   │   │   ├── interpolation1d.hpp       # LinearInterpolation1D template
│   │   │   ├── piecewisepolynomial.hpp   # PiecewisePolynomial class
│   │   │   └── piecewisepolynomial.cpp
│   │   ├── linalg/
│   │   │   ├── linalg.hpp       # Cholesky, eigenvalue, spectral truncation
│   │   │   ├── choldcmp.cpp
│   │   │   ├── eigensym.cpp
│   │   │   └── spectrunc.cpp
│   │   ├── optim/
│   │   │   ├── polyfunc.hpp     # Polynomial functor
│   │   │   └── roots.hpp        # Bracket and secant root finders
│   │   ├── random/
│   │   │   ├── normalrng.hpp    # NormalRng<URNG> template
│   │   │   └── rng.hpp          # Type aliases: NormalRngMt19937, NormalRngRanLux*, etc.
│   │   └── stats/
│   │       ├── statisticscalculator.hpp  # Abstract base class
│   │       ├── meanvarcalculator.hpp     # Running mean and variance
│   │       ├── histogramcalculator.hpp   # Frequency histogram (strategy pattern)
│   │       ├── errorfunction.hpp         # std::erf/erfc + inverfc (Halley refinement)
│   │       ├── errorfunction.cpp
│   │       ├── normaldistribution.hpp    # NormalDistribution with pdf/cdf/invcdf
│   │       └── univariatedistribution.hpp
│   │
│   ├── market/
│   │   ├── market.hpp           # Global SPtrMap registry for market objects
│   │   ├── market.cpp
│   │   ├── yieldcurve.hpp       # YieldCurve: discount, spotRate, fwdRate
│   │   ├── yieldcurve.cpp
│   │   ├── volatilitytermstructure.hpp
│   │   └── volatilitytermstructure.cpp
│   │
│   ├── products/
│   │   ├── product.hpp          # Abstract Product base class
│   │   ├── europeancallput.hpp  # European call/put payoff
│   │   ├── digitalcallput.hpp   # Digital call/put payoff
│   │   ├── americancallput.hpp  # American call/put (PDE only)
│   │   └── asianbasketcallput.hpp
│   │
│   ├── methods/
│   │   ├── montecarlo/
│   │   │   ├── pathgenerator.hpp         # Abstract PathGenerator base
│   │   │   ├── pathgenerator.cpp
│   │   │   ├── eulerpathgenerator.hpp    # EulerPathGenerator<NRNG> template
│   │   │   └── mcparams.hpp              # McParams: URNG type, path gen type
│   │   └── pde/
│   │       ├── pde1dsolver.hpp/cpp       # 1D Crank-Nicolson PDE solver
│   │       ├── pdebase.hpp/cpp
│   │       ├── pdegrid.hpp
│   │       ├── pdeparams.hpp
│   │       ├── pderesults.hpp
│   │       └── tridiagonalops1d.hpp
│   │
│   └── pricers/
│       ├── simplepricers.hpp/cpp         # Closed-form: BS European, digital, caplet, CDS, KO fwd
│       ├── bsmcpricer.hpp/cpp            # Single-asset Black-Scholes MC pricer
│       ├── multiassetbsmcpricer.hpp/cpp  # Multi-asset correlated MC pricer
│       └── bsmcsimpledeltahedger.hpp/cpp # Delta hedging simulation engine
│
├── pyqflib/                     # Python C-extension (compiled as pyqflib.so/.pyd)
│   ├── CMakeLists.txt
│   ├── pymodule.cpp             # PyMODINIT_FUNC, method table registration
│   ├── pycpp.hpp                # PY_BEGIN/PY_END macros, type converters
│   ├── pyutils.hpp              # Helper functions for Python<->C++ conversion
│   ├── pytestfunc.hpp           # Echo functions for type round-trip testing
│   ├── pyfunctions0.hpp         # Math functions: erf, normal, interpolation, etc.
│   ├── pyfunctions1.hpp         # Pricing: fwdPrice, euroBS, digiBS, koFwd, etc.
│   ├── pyfunctions2.hpp         # Market objects: yield curves, vol surfaces
│   ├── pyfunctions3.hpp         # MC pricers: euroBSMC, asianBasketBSMC
│   ├── pyfunctions4.hpp         # PDE pricers + BSHedge
│   └── qflib/                   # Python package (installed by pip)
│       ├── __init__.py          # Re-exports pyqflib symbols with docstrings
│       └── pyproject.toml       # Package metadata (name="qflib", version="0.11.0")
│
└── tests/                       # Catch2 unit tests and benchmarks
    ├── CMakeLists.txt
    ├── test_utils.cpp
    ├── test_sptrmap.cpp
    ├── test_errorfunction.cpp
    ├── test_normaldistribution.cpp
    ├── test_interpolation1d.cpp
    ├── test_polynomial.cpp
    ├── test_meanvarcalculator.cpp
    ├── test_products.cpp
    └── bench_fwdrate_precompute.cpp
```

---

## Key components

### Delta hedging simulation engine

`BsMcSimpleDeltaHedger` (`qflib/pricers/bsmcsimpledeltahedger.hpp/cpp`) is the main simulation engine. It:

1. Builds a uniform time grid with `hedgeFreq` rehedging events per year
2. Precomputes all drifts, standard deviations, and forward rates before the MC loop
3. On each path: generates a full log-normal spot path via the Euler scheme, then rebalances the delta hedge at each grid node in a self-financing manner
4. At maturity: settles the option payoff and liquidates the hedge, recording the total P&L
5. Feeds the P&L into each `StatisticsCalculator` in `statsCalcs`

### Statistics calculators (strategy pattern)

The `hedge()` method accepts any combination of `StatisticsCalculator` objects. Two are provided:

- **`MeanVarCalculator`** — accumulates running sum and sum-of-squares; computes unbiased sample mean and variance on demand at `results()`
- **`HistogramCalculator`** — bins samples into user-defined intervals using `std::upper_bound` for O(log B) lookup; returns a counts matrix (bins × variables)

Both can be passed to the same simulation run simultaneously.

### Yield curves and market objects

`YieldCurve` wraps a `PiecewisePolynomial` log-discount factor and exposes `discount(t)`, `spotRate(t)`, and `fwdRate(t0, t1)`. Market objects are stored in a global `SPtrMap` registry and looked up by name from Python.

### Random number generation

The default RNG is Mersenne Twister (`std::mt19937`) with period 2^19937 − 1, equidistributed in up to 623 dimensions. Normal deviates are produced via `std::normal_distribution<double>`. Other supported URNGs: `minstd_rand`, `ranlux24`, `ranlux48`.

### Error function

`erf`/`erfc` delegate to `std::erf`/`std::erfc` (C++11, full double precision). `inverfc` uses a rational approximation for the initial guess, then refines with 4 iterations of Halley's method, achieving ~1e-12 precision.

### PDE solver

A 1D Crank-Nicolson finite difference solver for European and American options, accessible via `qf.euroBSPDE` and `qf.amerBSPDE`.

---

## Build setup

### Prerequisites

| Dependency | Version | Notes |
|---|---|---|
| C++ compiler | C++20 | GCC 14+, AppleClang (Xcode 15+), MSVC 2022 |
| CMake | 3.25+ | |
| Armadillo | 14.2.2 | Linear algebra (BLAS/LAPACK backend) |
| Python | 3.12 | With NumPy installed |
| Catch2 | 3.x | For unit tests only |

**macOS:**
```bash
brew install armadillo catch2
```

**Linux:** Download Armadillo 14.2.2 from [arma.sourceforge.net](https://arma.sourceforge.net) and extract to `~/thirdparty/armadillo-14.2.2/`.

**Windows:** Extract Armadillo to a directory of your choice.

---

### CMakeLists.txt changes required before building

**1. Top-level `CMakeLists.txt`** — set your Armadillo path:

```cmake
# macOS (Homebrew default — change if installed elsewhere):
set(ARMADILLO_PREFIX "/opt/homebrew/opt/armadillo")

# Linux — set the directory containing armadillo-14.2.2/:
set(THIRDPARTY_DIRECTORY "$ENV{HOME}/thirdparty")
set(ARMA_VERSION 14.2.2)

# Windows — set the directory containing armadillo-14.2.2/:
set(THIRDPARTY_DIRECTORY "C:/path/to/thirdparty")
set(ARMA_VERSION 14.2.2)
```

**2. `pyqflib/CMakeLists.txt`** — set your Python/conda environment path:

```cmake
# macOS
set(PYTHON_HOME "$ENV{HOME}/anaconda3/envs/YOUR_ENV_NAME")

# Linux
set(PYTHON_HOME $ENV{HOME}/miniconda3/envs/YOUR_ENV_NAME)

# Windows
set(PYTHON_HOME "C:/path/to/conda/envs/YOUR_ENV_NAME")
```

The Python include directory, NumPy include directory, and Python library path are all derived automatically from `PYTHON_HOME`.

---

### Configure and build

The build directory should be created **one level above** `qflib-0.11.0/`:

```
C++_Project/
├── qflib-0.11.0/    ← source
└── build/           ← build output (create this)
```

**macOS / Linux:**
```bash
cd C++_Project
mkdir build && cd build
cmake ../qflib-0.11.0 -G "Unix Makefiles"
cmake --build .
```

**macOS with Ninja (faster):**
```bash
mkdir build && cd build
cmake ../qflib-0.11.0 -G Ninja
cmake --build .
```

**Windows (Visual Studio 2022):**
```bash
mkdir build && cd build
cmake ..\qflib-0.11.0 -G "Visual Studio 17 2022" -A x64
cmake --build . --config Release
```

Build outputs after a successful build:
```
qflib-0.11.0/lib/libqflib.a           # static C++ library (macOS/Linux)
qflib-0.11.0/lib/qflib.lib            # static C++ library (Windows)
qflib-0.11.0/pyqflib/qflib/pyqflib.so # Python extension (macOS/Linux)
qflib-0.11.0/pyqflib/qflib/pyqflib.pyd # Python extension (Windows)
qflib-0.11.0/bin/qflib_tests          # Catch2 test executable
qflib-0.11.0/bin/qflib_bench          # Benchmark executable
```

---

## Python installation

After building, install the Python package in editable mode from the `pyqflib/` directory:

```bash
cd qflib-0.11.0/pyqflib
pip install -e .
```

This registers the `qflib` package on Python's path without copying any files. The `pyqflib.so` compiled by CMake is already in place inside `pyqflib/qflib/`. After this, `import qflib` works in any terminal, Jupyter notebook, or script in the same environment.

To verify:
```python
import qflib
print(qflib.version())  # should print "0.11.0"
```

---

## Run tests

```bash
# Build the test executable (from the build directory)
cmake --build . --target qflib_tests

# Run all 163 assertions across 36 test cases
./qflib-0.11.0/bin/qflib_tests

# Run a specific test group
./qflib-0.11.0/bin/qflib_tests "[errorfunction]"
./qflib-0.11.0/bin/qflib_tests "[normaldist]"

# Verbose output
./qflib-0.11.0/bin/qflib_tests -v high
```

Run the forward rate precomputation benchmark:
```bash
cmake --build . --target qflib_bench
./qflib-0.11.0/bin/qflib_bench
```

---

## Python usage example

```python
import qflib as qf

# Create a yield curve (name, times, rates, compounding convention)
qf.ycCreate("USD", [0.25, 0.5, 1.0, 2.0, 5.0], [0.04, 0.042, 0.045, 0.047, 0.05], 0)

# Black-Scholes delta hedge simulation
# Returns: [mean_pnl, std_pnl, hist_counts_array, hist_edges_array]
result = qf.BSHedge(
    opttype=0,        # 0 = European, 1 = Digital
    payofftype=1,     # 1 = call, -1 = put
    strike=100.0,
    expiry=1.0,
    ycname="USD",
    divyield=0.02,
    spot=100.0,
    simvol=0.20,      # volatility driving the simulated paths
    markvol=0.20,     # volatility used for BS pricing and delta
    urngtype=1,       # 1 = MT19937
    pathgentype=0,    # 0 = Euler
    hedgefreq=52,     # rehedgings per year
    npaths=100000,
    histmin=-10.0,
    histmax=10.0,
    nbins=100
)

mean_pnl = result[0]
std_pnl  = result[1]
hist_counts = result[2]   # numpy array of length nbins
hist_edges  = result[3]   # numpy array of length nbins+1

# European option price and Greeks
# returns [price, delta, gamma, vega, theta, rho]
greeks = qf.euroBS(1, 100.0, 100.0, 1.0, 0.045, 0.02, 0.20)

# Vectorize over an array of strikes (qflib functions are scalar)
import numpy as np
euroBS_vec = np.vectorize(qf.euroBS)
strikes = np.linspace(80, 120, 41)
prices = euroBS_vec(1, 100.0, strikes, 1.0, 0.045, 0.02, 0.20)[:, 0]
```

---

## Design decisions and optimizations

### Cache-friendly path generation (`EulerPathGenerator`)

The original `next()` used a scratch array `normalDevs_` — generating random values into a temporary buffer, then copying them into `pricePath` column by column. Since Armadillo matrices are column-major, `begin_col(j)` / `end_col(j)` expose contiguous memory iterators. `NormalRng::next()` accepts any iterator, so the scratch array is eliminated entirely — random values are written directly into the destination column in one pass:

```cpp
// Before: generate into scratch → copy into matrix (two passes, extra allocation)
nrng_.next(normalDevs_.begin(), normalDevs_.end());
for (size_t i = 0; i < ntimesteps_; ++i)
    pricePath(i, j) = normalDevs_(i);

// After: write directly into column (one pass, no scratch array)
nrng_.next(pricePath.begin_col(j), pricePath.end_col(j));
```

This removes the intermediate buffer and halves the number of memory writes in the hot path.

### Performance optimization: forward rate precomputation

The original `hedge()` called `fwdRate()` on every path at every time step — ~10.4M redundant curve traversals for 52 steps × 100K paths. Since the time grid is fixed before the MC loop, all forward rates are deterministic and can be precomputed once.

**Fix:** all forward rates precomputed into `fwdrates_[]` and `fwdratesToT_[]` before the MC loop. The hot loop performs O(1) array lookups.

| Metric | Before | After |
|---|---|---|
| `fwdRate()` calls per simulation | ~10.4M | 104 (2 per grid point) |
| Reduction factor | — | 100,000× |

### `HistogramCalculator` implementation

Implemented `HistogramCalculator<ITER>` that plugs into the `StatisticsCalculator` strategy pattern. Takes user-defined bin edges, uses `std::upper_bound` for O(log B) bin lookup per sample, and returns a counts matrix (rows = bins, cols = variables). Samples outside the range are silently ignored. Both `MeanVarCalculator` and `HistogramCalculator` can be passed to the same simulation run simultaneously.

### Analytical digital option Greeks

Added closed-form delta for digital call and put options to `simplepricers.hpp/cpp`, enabling `BsMcSimpleDeltaHedger` to delta hedge digital options alongside Europeans.

### Numerical precision: `erf`/`erfc`/`inverfc`

Replaced the custom 28-coefficient Chebyshev approximation for `erf`/`erfc` with `std::erf`/`std::erfc` (C++11 standard library). The `inverfc` inverse uses a rational approximation for the initial guess, then refines with Halley's method (cubic convergence) to achieve ~1e-12 precision. `NormalDistribution::invcdf` inherits this precision directly.

### Named object registry (`SPtrMap`)

`SPtrMap<T>` is a type-safe string-to-shared-pointer registry used to store and retrieve market objects (yield curves, vol surfaces) by name from both C++ and Python. Key design details:

- **Name normalization** — `processName()` trims whitespace and uppercases every key before storage, so `"usd"`, `" USD "`, and `"USD"` all resolve to the same object. Internal blanks are rejected with `QF_ASSERT`.
- **Versioning** — every `set()` call increments a per-entry version counter, making it straightforward to detect stale references after an object is overwritten.
- **Ownership** — values are `std::shared_ptr<T>`, so the registry never outlives its objects and callers can hold their own references without dangling pointers.

### Unit test suite (Catch2)

163 assertions across 36 test cases covering: compounding conversions, `SPtrMap`, error function, normal distribution, 1D interpolation, polynomials, mean/variance calculator, and option product payoffs.

