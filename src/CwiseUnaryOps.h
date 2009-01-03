#ifndef CWISEUNARYOPS_H
#define CWISEUNARYOPS_H

template<typename Scalar>
struct CwiseMaxOp {
  CwiseMaxOp(const Scalar& other) : m_other(other){}
  const Scalar operator()(const Scalar& x) const { return x < m_other ? m_other : x; }
  Scalar m_other;
};

#endif // CWISEUNARYOPS_H
