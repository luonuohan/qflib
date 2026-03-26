#pragma once

#include <qflib/math/stats/statisticscalculator.hpp>
#include <vector>
#include <algorithm>

BEGIN_NAMESPACE(qf)

/**
 * Computes a frequency histogram over fixed bin intervals,
 * for N variables per sample.
 *
 * Template parameter ITER is the iterator type for your sample data,
 * e.g. double* if you call addSample(&pnl, &pnl+1).
 */
template <typename ITER>
class HistogramCalculator : public StatisticsCalculator<ITER> {
public:
  /**
   * @param binEdges  Vector of length B+1 giving the B bin boundaries in increasing order.
   *                  The i-th bin is [binEdges[i], binEdges[i+1]).
   * @param nvars     Number of variables per sample (columns).
   */
  HistogramCalculator(const std::vector<double>& binEdges, size_t nvars)
    : StatisticsCalculator<ITER>(nvars, binEdges.size()-1)
    , binEdges_(binEdges), nBins_(binEdges.size()-1)
  {
    QF_ASSERT(binEdges_.size() >= 2,
      "Need at least two edges to form one bin");
    // sanity: edges must be strictly increasing
    for (size_t i = 1; i < binEdges_.size(); ++i)
      QF_ASSERT(binEdges_[i] > binEdges_[i-1],
        "binEdges must be strictly increasing");
  }

  /// Clear all counts back to zero
  virtual void reset() override {
    StatisticsCalculator<ITER>::reset();
    this->results_.zeros(nBins_, this->nVariables());
  }

  /// For each variable j, find its bin and increment that count
  virtual void addSample(ITER begin, ITER end) override {
    QF_ASSERT(end - begin == this->nVariables(),
      "HistogramCalculator: sample size mismatch");
    for (size_t j = 0; j < this->nVariables(); ++j) {
      double x = *(begin + j);
      // find bin index b so that binEdges[b] <= x < binEdges[b+1]
      auto it = std::upper_bound(binEdges_.begin(),
                                 binEdges_.end(), x);
      if (it == binEdges_.begin() || it == binEdges_.end()) {
        // x is below first edge or >= last edge: skip or put in under/overflow?
        // here we choose to ignore values outside [min,max)
        continue;
      }
      size_t b = size_t((it - binEdges_.begin()) - 1);
      // results_ is inherited from StatisticsCalculator:
      // rows = bins, cols = variables
      this->results_(b, j) += 1.0;
    }
    ++this->nsamples_;
  }

  /// Return the counts matrix: rows = bins, cols = variables
  virtual Matrix const & results() override {
    return this->results_;
  }

private:
  std::vector<double> binEdges_;
  size_t nBins_;
};

END_NAMESPACE(qf)
