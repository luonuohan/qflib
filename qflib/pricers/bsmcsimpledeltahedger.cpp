#include "bsmcsimpledeltahedger.hpp"
#include <qflib/math/random/rng.hpp>
#include <qflib/methods/montecarlo/eulerpathgenerator.hpp>
#include <cmath>
#include <vector>                   // for std::vector<>
#include <qflib/math/stats/meanvarcalculator.hpp>   // or wherever your MeanVarCalculator lives
#include <qflib/math/stats/histogramcalculator.hpp>// ditto
#include <cmath>
#include <armadillo>

BEGIN_NAMESPACE(qf)

BsMcSimpleDeltaHedger::BsMcSimpleDeltaHedger(SPtrProduct    prod,
                                             SPtrYieldCurve discountCurve,
                                             double         divYield,
                                             double         simVol,
                                             double         markVol,
                                             double         spot,
                                             McParams       mcparams,
                                             OptionType     optType,
                                             int            payoffType,
                                             double         strike,
                                             double         maturity )
: prod_(prod), discyc_(discountCurve),
  divyld_(divYield), simVol_(simVol),
  markVol_(markVol), spot0_(spot),
  mcparams_(mcparams), optType_(optType),
  payoffType_(payoffType), K_(strike), T_(maturity)
{
  // nothing else yet—time grid will be built in hedge()
}

void BsMcSimpleDeltaHedger::priceAndDelta(double S, double t, double r,
                                           double& price, double& delta) const
{
  double tau = prod_->payTimes().back() - t;  // time to maturity
  Vector greeks;
  if (optType_ == European) {
    greeks = europeanOptionBS(payoffType_, S, K_, tau, r, divyld_, markVol_ );
    price = greeks[0];  delta = greeks[1];
  } else {
    greeks = digitalOptionBS(payoffType_, S, K_, tau, r, divyld_, markVol_ );
    price = greeks[0];  delta = greeks[1];
  }
}

void BsMcSimpleDeltaHedger::hedge(const std::vector<StatisticsCalculator<double*>*>& statsCalcs,
                                  unsigned hedgeFreq,
                                  unsigned long nPaths)
{
  // 1) Build uniform time grid [0, T] with N = hedgeFreq * T
  double T = prod_->payTimes().back();
  unsigned N = unsigned(std::ceil(hedgeFreq * T)); // number of re-hedges = number of time steps in the grid
  times_.resize(N+1);
  for (unsigned i=0; i<=N; ++i) times_[i] = double(i)/hedgeFreq;

  // 2) Precompute drifts, stdevs, and forward rates
  drifts_.resize(N);
  stdevs_.resize(N);
  fwdrates_.resize(N);
  fwdratesToT_.resize(N);
  double t0 = 0.0;
  double vol2 = simVol_ * simVol_;
  for (unsigned i=0; i<N; ++i) {
    double t1 = times_[i+1];
    double dt = t1 - t0;
    double var = vol2 * dt;
    stdevs_[i] = std::sqrt(var);
    fwdrates_[i] = discyc_->fwdRate(t0, t1);
    fwdratesToT_[i] = discyc_->fwdRate(t0, T);
    drifts_[i] = (fwdrates_[i] - divyld_)*dt - 0.5*var;
    t0 = t1;
  }

  // 3) Create Euler path generator over our grid
  using EulerMt = EulerPathGenerator<NormalRngMt19937>;
  if (mcparams_.pathGenType == McParams::PathGenType::EULER) {
      pathGen_ = std::make_unique<EulerMt>(times_.begin(), times_.end(), 1);
  } else {QF_ASSERT(false, "Only EULER pathGen supported");}

  // 4) Monte Carlo loop
  qf::Matrix path(N+1, 1);
  std::vector<double> spotPath(N+1); // Store the spot path

  for (unsigned long i=0; i<nPaths; ++i) {
    pathGen_->next(path);
    // --- simulate spot path ---
    // double S = spot0_;
    // simulate *and store* the entire stock path
    spotPath[0] = spot0_;

    for (unsigned j=1; j<=N; ++j) {
      double z = path(j,0);
      //S = S * std::exp(drifts_[j-1] + stdevs_[j-1]*z);
      double increment = std::exp(drifts_[j-1] + stdevs_[j-1]*z);
      spotPath[j] = spotPath[j-1] * increment;
    }

    // --- initialize hedge at t=0 ---
    double cash, deltaPrev;
    double initPrice, initDelta;
    priceAndDelta(spot0_, 0.0, fwdratesToT_[0], initPrice, initDelta);
    // buyer:  pay premium, receive proceeds from shorting Δ shares
    cash       = -initPrice + initDelta*spotPath[0]; // for buyer of option -initPrice + initDelta * spotPath[0];  
    //cash       = initPrice - initDelta*spot0_; for seller of option
    deltaPrev  = initDelta;
    double tPrev = 0.0;

    // --- rehedge at each grid node ---
    for (unsigned j=1; j<N; ++j) {
      double tNow = times_[j];
      // accrue cash at risk-free
      cash *= std::exp(fwdrates_[j-1] * (tNow - tPrev));

      // compute new Δ
      double priceNow, deltaNow;
      priceAndDelta(spotPath[j], tNow, fwdratesToT_[j], priceNow, deltaNow);

      // trade underlying
      double trade = (deltaNow - deltaPrev)*spotPath[j];
      // buyer:  sell (short) Δ_new–Δ_old shares → receive
      cash   += trade;
      deltaPrev = deltaNow;
      tPrev = tNow;
    }

    // 1) extract final spot ST
    double ST = spotPath[N];

    // 2) compute payoff by calling the 1‐fixing overload
    qf::Vector spots(1);
    spots[0] = ST;
    prod_->eval(/* idx = */ 0, spots, /* contVal = */ 0.0);
    double payoff = prod_->payAmounts()[0];

    
    // single‐payoff
    // accrue cash final leg (tPrev = times_[N-1])
    cash *= std::exp(fwdrates_[N-1] * (T - tPrev));
    cash -= deltaPrev * spotPath[N]; 

    double pnl = payoff + cash;

    // 5) Feed P&L to stats
    for (auto* sc : statsCalcs)
      sc->addSample(&pnl, &pnl+1);
  }
}

END_NAMESPACE(qf)
