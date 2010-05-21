#ifndef CWISEADDONS_H
#define CWISEADDONS_H

EIGEN_STRONG_INLINE const CwiseUnaryOp<ei_scalar_angle_op<Scalar>, Derived>
angle() const
{
  return derived();
}

EIGEN_STRONG_INLINE const CwiseUnaryOp<ei_scalar_sgn_op<Scalar>, Derived>
sgn() const
{
  return derived();
}

EIGEN_STRONG_INLINE const CwiseUnaryOp<ei_scalar_ceil_op<Scalar>, Derived>
ceil() const
{
  return derived();
}

EIGEN_STRONG_INLINE const CwiseUnaryOp<ei_scalar_floor_op<Scalar>, Derived>
floor() const
{
  return derived();
}

EIGEN_STRONG_INLINE const CwiseUnaryOp<ei_scalar_isnan_op<Scalar>, Derived>
isnan() const
{
  return derived();
}

EIGEN_STRONG_INLINE const CwiseUnaryOp<ei_scalar_mod_n_op<Scalar>, Derived>
modN(const Scalar& divisor) const
{
  return CwiseUnaryOp<ei_scalar_mod_n_op<Scalar>,Derived>
          (derived(), ei_scalar_mod_n_op<Scalar>(divisor));
}


EIGEN_STRONG_INLINE const CwiseUnaryOp<ei_scalar_exp_n_op<Scalar>, Derived>
expN(const Scalar& base) const
{
  return CwiseUnaryOp<ei_scalar_exp_n_op<Scalar>,Derived>
          (derived(), ei_scalar_exp_n_op<Scalar>(base));
}


EIGEN_STRONG_INLINE const CwiseUnaryOp<ei_scalar_log_n_op<Scalar>, Derived>
logN(const Scalar& base) const
{
  return CwiseUnaryOp<ei_scalar_log_n_op<Scalar>,Derived>
          (derived(), ei_scalar_log_n_op<Scalar>(base));
}


EIGEN_STRONG_INLINE const CwiseUnaryOp<ei_scalar_clip_under_op<Scalar>, Derived>
clipUnder(const Scalar& under = 0) const
{
  return CwiseUnaryOp<ei_scalar_clip_under_op<Scalar>,Derived>
          (derived(), ei_scalar_clip_under_op<Scalar>(under));
}

EIGEN_STRONG_INLINE const CwiseUnaryOp<ei_scalar_clip_op<Scalar>, Derived>
clip(const Scalar& under, const Scalar& over) const
{
  return CwiseUnaryOp<ei_scalar_clip_op<Scalar>,Derived>
          (derived(), ei_scalar_clip_op<Scalar>(under, over));
}

#endif // CWISEADDONS_H
