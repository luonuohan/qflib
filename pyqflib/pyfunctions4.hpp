/**
@file  pyfunctions4.hpp
@brief Implementation of Python callable functions
*/
#include <pyqflib/pyutils.hpp>

#include <qflib/market/market.hpp>
#include <qflib/products/europeancallput.hpp>
#include <qflib/products/digitalcallput.hpp>
#include <qflib/products/americancallput.hpp>
#include <qflib/methods/pde/pde1dsolver.hpp>
#include <qflib/pricers/bsmcsimpledeltahedger.hpp>
#include <qflib/math/stats/meanvarcalculator.hpp>
#include <qflib/math/stats/histogramcalculator.hpp>

static
PyObject*  pyQfEuroBSPDE(PyObject* pyDummy, PyObject* pyArgs)
{
PY_BEGIN;

  PyObject* pyPayoffType(NULL);
  PyObject* pySpot(NULL);
  PyObject* pyStrike(NULL);
  PyObject* pyTimeToExp(NULL);
  PyObject* pyDiscountCrv(NULL);
  PyObject* pyDivYield(NULL);
  PyObject* pyVolatility(NULL);
  PyObject* pyPdeParams(NULL);
  PyObject* pyAllResults(NULL);

  if (!PyArg_ParseTuple(pyArgs, "OOOOOOOOO", &pyPayoffType, &pyStrike, &pyTimeToExp, 
    &pySpot, &pyDiscountCrv, &pyDivYield, &pyVolatility, &pyPdeParams, &pyAllResults))
    return NULL;

  int payoffType = asInt(pyPayoffType);
  double spot = asDouble(pySpot);
  double strike = asDouble(pyStrike);
  double timeToExp = asDouble(pyTimeToExp);

  std::string name = asString(pyDiscountCrv);
  qf::SPtrYieldCurve spyc = qf::market().yieldCurves().get(name);
  QF_ASSERT(spyc, "error: yield curve " + name + " not found");

  double divYield = asDouble(pyDivYield);
  double vol = asDouble(pyVolatility);

  // read the PDE parameters
  qf::PdeParams pdeparams = asPdeParams(pyPdeParams);
  // read the allresults flag
  bool allresults = asBool(pyAllResults);

  // create the product
  qf::SPtrProduct spprod(new qf::EuropeanCallPut(payoffType, strike, timeToExp));
  // create the PDE solver
  qf::Pde1DResults results;
  qf::Pde1DSolver solver(spprod, spyc, spot, divYield, vol, results);
  solver.solve(pdeparams);

  // write results
  PyObject* ret = PyDict_New();
  int ok = PyDict_SetItem(ret, asPyScalar("Price"), asPyScalar(results.prices[0]));

  if (allresults) {
    qf::Vector spots;
    results.getSpotAxis(0, spots);
    qf::Matrix values(results.times.size(), results.values.front().size());
    for (size_t i = 0; i < results.times.size(); ++i)
      for (size_t j = 0; j < results.values.front().size(); ++j)
        values(i, j) = results.values[i](j, 0);

    PyDict_SetItem(ret, asPyScalar("Times"), asNumpy(results.times));
    PyDict_SetItem(ret, asPyScalar("Spots"), asNumpy(spots));
    PyDict_SetItem(ret, asPyScalar("Values"), asNumpy(values));
  }
  return ret;

PY_END;
}

static
PyObject*  pyQfAmerBSPDE(PyObject* pyDummy, PyObject* pyArgs)
{
PY_BEGIN;

  PyObject* pyPayoffType(NULL);
  PyObject* pySpot(NULL);
  PyObject* pyStrike(NULL);
  PyObject* pyTimeToExp(NULL);
  PyObject* pyDiscountCrv(NULL);
  PyObject* pyDivYield(NULL);
  PyObject* pyVolatility(NULL);
  PyObject* pyPdeParams(NULL);
  PyObject* pyAllResults(NULL);

  if (!PyArg_ParseTuple(pyArgs, "OOOOOOOOO", &pyPayoffType, &pyStrike, &pyTimeToExp, 
    &pySpot, &pyDiscountCrv, &pyDivYield, &pyVolatility, &pyPdeParams, &pyAllResults))
    return NULL;

  int payoffType = asInt(pyPayoffType);
  double spot = asDouble(pySpot);
  double strike = asDouble(pyStrike);
  double timeToExp = asDouble(pyTimeToExp);

  std::string name = asString(pyDiscountCrv);
  qf::SPtrYieldCurve spyc = qf::market().yieldCurves().get(name);
  QF_ASSERT(spyc, "error: yield curve " + name + " not found");

  double divYield = asDouble(pyDivYield);
  double vol = asDouble(pyVolatility);

  // read the PDE parameters
  qf::PdeParams pdeparams = asPdeParams(pyPdeParams);
  // read the allresults flag
  bool allresults = asBool(pyAllResults);

  // create the product
  qf::SPtrProduct spprod(new qf::AmericanCallPut(payoffType, strike, timeToExp));
  // create the PDE solver
  qf::Pde1DResults results;
  bool storeAllResults = true;
  qf::Pde1DSolver solver(spprod, spyc, spot, divYield, vol, results, storeAllResults);
  solver.solve(pdeparams);

  // write results
  PyObject* ret = PyDict_New();
  int ok = PyDict_SetItem(ret, asPyScalar("Price"), asPyScalar(results.prices[0]));

  if (allresults) {
    qf::Vector spots;
    results.getSpotAxis(0, spots);
    qf::Matrix values(results.times.size(), results.values.front().size());
    for (size_t i = 0; i < results.times.size(); ++i)
      for (size_t j = 0; j < results.values.front().size(); ++j)
        values(i, j) = results.values[i](j, 0);

    PyDict_SetItem(ret, asPyScalar("Times"), asNumpy(results.times));
    PyDict_SetItem(ret, asPyScalar("Spots"), asNumpy(spots));
    PyDict_SetItem(ret, asPyScalar("Values"), asNumpy(values));
  }
  return ret;

PY_END;
}

// pyQfBSHedge for simulation hedging
static
PyObject* pyQfBSHedge(PyObject* pyDummy, PyObject* pyArgs)
{
  PY_BEGIN;

  // 16 args: optType, payoffType, strike, timeToExp,
  //           discName, divY, spot, simVol, markVol,
  //           pathGenType, urngType, hedgeFreq, nPaths,
  //           histMin, histMax, histNbins
  int    optType, payoffType;
  double strike, timeToExp;
  const char* discName;
  double divY, spot, simVol, markVol;
  int    pathGenType, urngType, hedgeFreq;
  unsigned long nPaths;
  double histMin, histMax;
  int    histNbins;

  if (!PyArg_ParseTuple(pyArgs,
        "iiddsddddiiikddi",
        &optType, &payoffType,
        &strike, &timeToExp,
        &discName,
        &divY, &spot, &simVol, &markVol,
        &pathGenType, &urngType,
        &hedgeFreq, &nPaths,
        &histMin, &histMax, &histNbins))
    return NULL;

  // 1) Get the yield curve
  std::string ycName(discName);
  qf::SPtrYieldCurve yc = qf::market().yieldCurves().get(ycName);
  QF_ASSERT(yc, "yield curve not found");

  // 2) Build the product
  qf::SPtrProduct prod;
  if (optType == 0)  // European
    //prod.reset(new qf::EuropeanCallPut(payoffType, strike, timeToExp));
    prod = std::make_shared<qf::EuropeanCallPut>(payoffType, strike, timeToExp);

  else               // Digital
    //prod.reset(new qf::DigitalCallPut(payoffType, strike, timeToExp));
    prod = std::make_shared<qf::DigitalCallPut>(payoffType, strike, timeToExp);

  // 3) MC params
  qf::McParams mc;
  mc.pathGenType = qf::McParams::PathGenType(pathGenType);
  mc.urngType    = qf::McParams::UrngType   (urngType);

  // 4) Instantiate the hedger
  qf::BsMcSimpleDeltaHedger hedger(
    prod, yc,
    divY, simVol, markVol,
    spot,
    mc,
    static_cast<qf::BsMcSimpleDeltaHedger::OptionType>(optType),
    payoffType,
    strike,
    timeToExp  ); 


  // 5) Prepare calculators
  qf::MeanVarCalculator<double*>   mv(1);
  // build binEdges vector
  std::vector<double> edges(histNbins+1);
  for (int i = 0; i <= histNbins; ++i){
    edges[i] = histMin + (histMax - histMin) * i / histNbins};
  qf::HistogramCalculator<double*> hist(edges, 1);

  // 6) Run the Δ-hedge
  hedger.hedge({ &mv, &hist }, hedgeFreq, nPaths);

  // 7) Collect results
  auto M  = mv.results();                 // 2×1: [ mean ; variance ]
  double mean = M(0,0);
  double var  = M(1,0);
  double std  = std::sqrt(var);

  auto H  = hist.results();               // histNbins×1

  // 8) Pack into Python dict
  PyObject* ret = PyDict_New();
  PyDict_SetItem(ret, asPyScalar("mean"),       asPyScalar(mean));
  PyDict_SetItem(ret, asPyScalar("std"),        asPyScalar(std));
  PyDict_SetItem(ret, asPyScalar("hist_counts"),asNumpy(H));
  // convert std::vector<double> → qf::Vector so asNumpy() can pick the right overload
  qf::Vector edgeVec(edges.size());
  for (size_t i = 0; i < edges.size(); ++i)
    edgeVec(i) = edges[i];
  PyDict_SetItem(ret, asPyScalar("hist_edges"), asNumpy(edgeVec));
  return ret;

  PY_END;
}

