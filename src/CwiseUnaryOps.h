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
struct CwiseAngleOp {
  CwiseAngleOp(){}
  const Scalar operator()(const Scalar& x) const { return atan2(x.imag(), x.real()); }
};

#endif // CWISEUNARYOPS_H
