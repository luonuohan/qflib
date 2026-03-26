/**
 * @file digitalcallput.hpp
 * @brief Cash‑or‑nothing European digital call / put.
 */
#ifndef QF_DIGITALCALLPUT_HPP
#define QF_DIGITALCALLPUT_HPP

#include <qflib/defines.hpp>
#include <qflib/products/product.hpp>

BEGIN_NAMESPACE(qf)

/* ---------------------------------------------------------------
   Digital (cash‑or‑nothing) option:  +1 = call, –1 = put
   --------------------------------------------------------------- */
class DigitalCallPut : public Product
{
public:
    DigitalCallPut(int payoffType, double strike, double timeToExp);

    virtual size_t nAssets() const override { return 1; }

    virtual void eval(size_t       eventIdx,
                      Vector const& spots,
                      double       contValue) override;

    virtual void eval(Matrix const& pricePath) override;

private:
    int    payoffType_;        // +1 call, –1 put
    double strike_;
    double timeToExp_;
};

/* -------- inline implementation -------------------------------- */

inline
DigitalCallPut::DigitalCallPut(int payoffType,
                               double strike,
                               double timeToExp)
: payoffType_(payoffType),
  strike_(strike),
  timeToExp_(timeToExp)
{
    QF_ASSERT((payoffType == 1 || payoffType == -1),
              "DigitalCallPut: payoff type must be 1 (call) or -1 (put)");
    QF_ASSERT(strike    > 0.0, "DigitalCallPut: strike must be positive");
    QF_ASSERT(timeToExp > 0.0, "DigitalCallPut: time to exp must be positive");

    // only one fixing time, the expiration
    fixTimes_.resize(1);
    fixTimes_[0] = timeToExp_;

    // assume that it will settle (pay) at expiration
    // TODO allow payment time later than expiration
    payTimes_.resize(1);
    payTimes_[0] = timeToExp_;

    // this product generates only one payment
    payAmounts_.resize(1);
}

/* ---- PDE / tree evaluation ---- */
inline void
DigitalCallPut::eval(size_t       /*eventIdx*/,
                     Vector const& spots,
                     double        /*contValue*/)
{
    double S = spots[0];
    payAmounts_[0] =
        (payoffType_ == 1) ? ((S >= strike_) ? 1.0 : 0.0)
                           : ((S <  strike_) ? 1.0 : 0.0);
}

/* ---- Monte‑Carlo wrapper ---- */
inline void
DigitalCallPut::eval(Matrix const& pricePath)
{
    Vector spotVec(1);
    spotVec[0] = pricePath(0, 0);
    eval(0, spotVec, 0.0);
}

END_NAMESPACE(qf)

#endif  /* QF_DIGITALCALLPUT_HPP */
