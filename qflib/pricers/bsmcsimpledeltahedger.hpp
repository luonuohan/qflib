#pragma once

#include <qflib/products/product.hpp>                    // for SPtrProduct
#include <qflib/market/yieldcurve.hpp>                   // for SPtrYieldCurve
#include <qflib/methods/montecarlo/mcparams.hpp>         // for McParams
#include <qflib/methods/montecarlo/pathgenerator.hpp>    // for PathGenerator
#include <qflib/math/stats/statisticscalculator.hpp>     // for StatisticsCalculator
#include <qflib/pricers/simplepricers.hpp>               // for europeanOptionBS, digitalOptionBS

BEGIN_NAMESPACE(qf)

class BsMcSimpleDeltaHedger {
public:
  enum OptionType { European, Digital };

  /// ctor: takes a product (European or digital), market curves, volatilities, spot, MC params
  BsMcSimpleDeltaHedger(SPtrProduct    prod,
                        SPtrYieldCurve discountCurve,
                        double         divYield,
                        double         simVol,
                        double         markVol,
                        double         spot,
                        McParams       mcparams,
                        OptionType     optType,
                        int            payoffType,
                        double         strike,
                        double         maturity );

  /// run the Δ-hedge: accumulate P&L into any number of stats calculators
  /// hedgeFreq = rehedges per year; nPaths = Monte Carlo paths
  /// in bsmcsimpledeltahedger.hpp
  void hedge(const std::vector<StatisticsCalculator<double*>*>& statsCalcs,
    unsigned                                 hedgeFreq,
    unsigned long                            nPaths);


private:
  SPtrProduct              prod_;
  SPtrYieldCurve           discyc_;
  double                   divyld_, simVol_, markVol_, spot0_;
  McParams                 mcparams_;
  OptionType optType_;
  int        payoffType_;  // +1 or –1
  double     K_, T_;


  std::unique_ptr<PathGenerator> pathGen_;
  std::vector<double>      times_, drifts_, stdevs_;

  /// compute analytic price & Δ (via simplepricers) at (S,t)
  void priceAndDelta(double S, double t, double& price, double& delta) const;
};

END_NAMESPACE(qf)
