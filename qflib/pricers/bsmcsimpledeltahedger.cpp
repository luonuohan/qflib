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

void BsMcSimpleDeltaHedger::priceAndDelta(double S, double t,
                                           double& price, double& delta) const
{
  double tau = prod_->payTimes().back() - t;  // time to maturity
  //double tau = T_ - t; 
  double r   = discyc_->fwdRate(t, t+tau);                  // instantaneous forward rate
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
                                  unsigned                                      hedgeFreq,
                                  unsigned long                                 nPaths)
{
  // 1) Build uniform time grid [0, T] with N = hedgeFreq * T
  double T = prod_->payTimes().back();
  unsigned N = unsigned(std::ceil(hedgeFreq * T));
  times_.resize(N+1);
  for (unsigned i=0; i<=N; ++i) times_[i] = double(i)/hedgeFreq;

  // 2) Precompute drifts & stdevs for log-normal Euler using simVol_
  drifts_.resize(N);
  stdevs_.resize(N);
  double t0 = 0.0;
  for (unsigned i=0; i<N; ++i) {
    double t1 = times_[i+1];
    double dt = t1 - t0;
    double var = simVol_*simVol_ * dt;
    stdevs_[i] = std::sqrt(var);
    double fwd = discyc_->fwdRate(t0, t1);
    drifts_[i] = (fwd - divyld_)*dt - 0.5*var;
    t0 = t1;
  }

  // 3) Create Euler path generator over our grid
  if (mcparams_.pathGenType == McParams::PathGenType::EULER) {
    pathGen_.reset(new EulerPathGenerator<NormalRngMt19937>(
                    times_.begin(), times_.end(), 1));
  } else {
    QF_ASSERT(false, "Only EULER pathGen supported");
  }

  // 4) Monte Carlo loop
  arma::mat path(N+1, 1);
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
    priceAndDelta(spot0_, 0.0, initPrice, initDelta);
    // buyer:  pay premium, receive proceeds from shorting Δ shares
    cash       = -initPrice + initDelta*spotPath[0]; // for buyer of option -initPrice + initDelta * spotPath[0];  
    //cash       = initPrice - initDelta*spot0_; for seller of option
    deltaPrev  = initDelta;
    double tPrev = 0.0;

    // --- rehedge at each grid node ---
    for (unsigned j=1; j<N; ++j) {
      double tNow = times_[j];
      // accrue cash at risk-free
      double rrf = discyc_->fwdRate(tPrev, tNow);
      cash *= std::exp(rrf*(tNow - tPrev));

      // compute new Δ
      double priceNow, deltaNow;
      priceAndDelta(spotPath[j], tNow, priceNow, deltaNow);

      // trade underlying
      double trade = (deltaNow - deltaPrev)*spotPath[j];
      // buyer:  sell (short) Δ_new–Δ_old shares → receive
      cash   += trade;
      deltaPrev = deltaNow;
      tPrev = tNow;

      // advance S for next step (already done)
    }

    // --- at maturity: settle option + liquidate hedge ---
    // copy spotPath to path
    /*
    for (unsigned j = 0; j <= N; ++j) path(j, 0) = spotPath[j];
    prod_->eval(path);
    double payoff = prod_->payAmounts()[0]; 
    */
         
    

    // 1) extract final spot ST
    double ST = spotPath[N];

    // 2) compute payoff by calling the 1‐fixing overload
    qf::Vector spots(1);
    spots[0] = ST;
    prod_->eval(/* idx = */ 0, spots, /* contVal = */ 0.0);
    double payoff = prod_->payAmounts()[0];

    
    // single‐payoff
    // accrue cash final leg
    double rrf = discyc_->fwdRate(tPrev, T);
    cash *= std::exp(rrf*(T - tPrev));
    cash -= deltaPrev * spotPath[N]; // change here

    double pnl = payoff + cash; // change here

    // 5) Feed P&L to stats
    for (auto* sc : statsCalcs)
      sc->addSample(&pnl, &pnl+1);
  }
}

END_NAMESPACE(qf)
