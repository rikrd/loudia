#ifndef FUNCTORSADDONS_H
#define FUNCTORSADDONS_H

template<typename Scalar> struct ei_scalar_ceil_op EIGEN_EMPTY_STRUCT {
  inline const Scalar operator() (const Scalar& a) const { return std::ceil(a); }
};

template<typename Scalar>
struct ei_functor_traits<ei_scalar_ceil_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = false }; };

#endif // FUNCTORSADDONS_H
