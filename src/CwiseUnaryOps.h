#ifndef CWISEUNARYOPS_H
#define CWISEUNARYOPS_H

#include <cmath>

template<typename Scalar>
struct CwiseMaxOp {
  CwiseMaxOp(const Scalar& other) : m_other(other){}
  const Scalar operator()(const Scalar& x) const { return x < m_other ? m_other : x; }
  Scalar m_other;
};

template<typename Scalar>
struct CwiseCeilOp {
  CwiseCeilOp(){}
  const Scalar operator()(const Scalar& x) const { return ceil(x); }
};


template<typename Scalar>
struct CwiseAngleOp {
  CwiseAngleOp(){}
  const Scalar operator()(const Scalar& x) const { return atan2(x.imag(), x.real()); }
};


template<typename Scalar>
struct CwiseExpOp {
  CwiseExpOp(const Scalar& other) : m_other(other){}
  const Scalar operator()(const Scalar& x) const { return pow(m_other, x); }
  Scalar m_other;
};

#endif // CWISEUNARYOPS_H
