# BsMcSimpleDeltaHedger — Code Walkthrough

## Overview

`BsMcSimpleDeltaHedger` simulates a **discrete delta-hedging strategy** for a
European or Digital option under the Black-Scholes model. For each Monte Carlo
path, it:

1. Simulates a log-normal stock price path
2. Initializes a delta-hedge portfolio at t = 0
3. Rebalances the portfolio at each hedging date
4. Computes the terminal P&L as: `payoff + cash account`

The distribution of P&L across all paths measures how well discrete hedging
replicates the option.

---

## Constructor

```cpp
BsMcSimpleDeltaHedger(prod, discountCurve, divYield, simVol, markVol,
                      spot, mcparams, optType, payoffType, strike, maturity)
```

Stores all inputs as member variables. No computation happens here — the time
grid and path generator are built lazily inside `hedge()`.

| Parameter | Meaning |
|---|---|
| `prod` | Product object (EuropeanCallPut or DigitalCallPut) |
| `discountCurve` | Yield curve for discounting / forward rates |
| `divYield` | Continuous dividend yield q |
| `simVol` | Volatility used to simulate the stock path (true vol) |
| `markVol` | Volatility used for BS pricing and delta (hedge vol) |
| `spot` | Initial stock price S₀ |
| `optType` | European or Digital |
| `payoffType` | +1 = Call, -1 = Put |
| `strike` | Strike price K |
| `maturity` | Time to expiry T |

---

## Helper: `priceAndDelta`

```cpp
void priceAndDelta(double S, double t, double r, double& price, double& delta)
```

Computes the BS **theoretical price and delta** at a given hedge time.

- `tau = T - t` — remaining time to maturity
- `r` = precomputed forward rate from `t` to `T` (passed in, not looked up live)
- Calls `europeanOptionBS` or `digitalOptionBS` from `simplepricers`
- Returns `price = V(S, t)` and `delta = ∂V/∂S`

This is called once at t = 0 and once per rehedge step per path.

---

## Main Method: `hedge(statsCalcs, hedgeFreq, nPaths)`

### Step 1 — Build the Time Grid

```
times_ = [0, 1/N_freq, 2/N_freq, ..., T]   (N+1 points, N intervals)
```

- `N = ceil(hedgeFreq × T)` total time steps
- `hedgeFreq` = rehedges per year (e.g. 52 = weekly)
- Grid is uniform with spacing `dt = 1/hedgeFreq`

---

### Step 2 — Precompute Per-Step Constants

For each interval `[t_i, t_{i+1}]`, computed **once before the MC loop**:

| Array | Formula | Purpose |
|---|---|---|
| `stdevs_[i]` | `simVol × sqrt(dt)` | Stock path diffusion per step |
| `fwdrates_[i]` | `fwdRate(t_i, t_{i+1})` | Cash accrual rate over one step |
| `fwdratesToT_[i]` | `fwdRate(t_i, T)` | Discount rate for BS pricing at `t_i` |
| `drifts_[i]` | `(r_i - q) × dt - 0.5 × simVol² × dt` | Log-normal drift per step |

**Why precompute?** `fwdRate` traverses the piecewise-polynomial yield curve
(binary search + integration) on every call. For N=52 steps and 100,000 paths,
this would be ~10M redundant curve evaluations. Precomputing reduces this to
2N lookups total.

---

### Step 3 — Create the Path Generator

```cpp
pathGen_ = std::make_unique<EulerPathGenerator<NormalRngMt19937>>(
               times_.begin(), times_.end(), 1);
```

- Uses Euler discretization with Mersenne Twister RNG
- `nfactors = 1` (single asset, no correlation matrix needed)
- `pathGen_->next(path)` fills an `(N+1) × 1` matrix with standard normal deviates Z ~ N(0,1)

---

### Step 4 — Monte Carlo Loop (per path)

#### 4a — Simulate Log-Normal Spot Path

```
S(t_{j}) = S(t_{j-1}) × exp(drift[j-1] + stdev[j-1] × Z_j)
```

where `Z_j ~ N(0,1)` drawn from the path generator.

This is the **Euler-Maruyama** discretization of GBM:

```
dS = (r - q) S dt + simVol S dW
```

#### 4b — Initialize Hedge at t = 0

```
initPrice, initDelta = priceAndDelta(S₀, 0, fwdratesToT_[0])
cash = -initPrice + initDelta × S₀
```

From the **buyer's perspective**:
- Pay the option premium: `-initPrice`
- Short `initDelta` shares and receive proceeds: `+initDelta × S₀`
- Net cash position at inception = `-initPrice + initDelta × S₀`

#### 4c — Rehedge at Interior Nodes j = 1, ..., N-1

At each rehedge time `t_j`:

1. **Accrue cash** at the risk-free rate over `[t_{j-1}, t_j]`:
   ```
   cash *= exp(fwdrates_[j-1] × dt)
   ```

2. **Recompute delta** using BS formula at current spot and time:
   ```
   priceAndDelta(S_j, t_j, fwdratesToT_[j], priceNow, deltaNow)
   ```

3. **Rebalance** the delta position:
   ```
   cash += (deltaNow - deltaPrev) × S_j
   ```
   Buy/sell `(Δ_new - Δ_old)` shares at current price `S_j`. The cash account
   reflects the cost of this trade.

#### 4d — Final Settlement at t = T

1. **Compute payoff** from the terminal spot:
   ```
   payoff = max(S_T - K, 0)   [call]
          = max(K - S_T, 0)   [put]
   ```
   via `prod_->eval(0, spots, 0.0)`

2. **Accrue final cash leg**:
   ```
   cash *= exp(fwdrates_[N-1] × (T - t_{N-1}))
   ```

3. **Unwind delta position** at S_T:
   ```
   cash -= deltaPrev × S_T
   ```

4. **P&L**:
   ```
   pnl = payoff + cash
   ```
   If the hedge were perfect (continuous rebalancing, markVol = simVol),
   `pnl` would be zero on every path. Discrete hedging and vol misspecification
   cause the P&L to deviate.

---

### Step 5 — Collect Statistics

```cpp
for (auto* sc : statsCalcs)
    sc->addSample(&pnl, &pnl+1);
```

Each `StatisticsCalculator` (e.g. `MeanVarCalculator`, `HistogramCalculator`)
receives the P&L of the current path. After all paths, they return:

- **Mean P&L** — should be near zero if markVol ≈ simVol
- **Std Dev of P&L** — measures hedging error; decreases as hedgeFreq increases
- **Histogram** — full P&L distribution shape

---

## P&L Formula Summary

$$
\text{P\&L} = \underbrace{\text{Payoff}(S_T)}_{\text{option received}}
+ \underbrace{(-V_0 + \Delta_0 S_0)}_{\text{initial position}}
\times e^{r \cdot T}
+ \sum_{j=1}^{N-1} \underbrace{(\Delta_j - \Delta_{j-1}) S_j}_{\text{rebalance cost}}
\times e^{r(T - t_j)}
- \Delta_{N-1} S_T
$$

---

## Key Design Decisions

| Decision | Reason |
|---|---|
| `simVol` ≠ `markVol` | Allows testing P&L sensitivity to vol misspecification |
| Buyer perspective | Cash account starts negative (paid premium) |
| Forward rates precomputed | Avoids ~10M redundant yield curve calls |
| `unique_ptr<PathGenerator>` | Exclusive ownership, zero overhead vs raw pointer |
| Iterator-range `addSample` | Uniform interface for scalar and vector samples |
