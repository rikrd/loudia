#ifndef CWISEADDONS_H
#define CWISEADDONS_H

inline const EIGEN_CWISE_UNOP_RETURN_TYPE(ei_scalar_angle_op)
angle() const
{ 
  return _expression(); 
}

inline const EIGEN_CWISE_UNOP_RETURN_TYPE(ei_scalar_ceil_op)
ceil() const
{ 
  return _expression(); 
}

inline const EIGEN_CWISE_UNOP_RETURN_TYPE(ei_scalar_floor_op)
floor() const
{ 
  return _expression(); 
}

inline const EIGEN_CWISE_UNOP_RETURN_TYPE(ei_scalar_isnan_op)
isnan() const
{ 
  return _expression(); 
}

inline const EIGEN_CWISE_UNOP_RETURN_TYPE(ei_scalar_exp_n_op)
expN(const Scalar& base) const
{
  return EIGEN_CWISE_UNOP_RETURN_TYPE(ei_scalar_exp_n_op)(_expression(), ei_scalar_exp_n_op<Scalar>(base));
}


inline const EIGEN_CWISE_UNOP_RETURN_TYPE(ei_scalar_log_n_op)
logN(const Scalar& base) const
{
  return EIGEN_CWISE_UNOP_RETURN_TYPE(ei_scalar_log_n_op)(_expression(), ei_scalar_log_n_op<Scalar>(base));
}


inline const EIGEN_CWISE_UNOP_RETURN_TYPE(ei_scalar_clip_under_op)
clipUnder(const Scalar& under) const
{
  return EIGEN_CWISE_UNOP_RETURN_TYPE(ei_scalar_clip_under_op)(_expression(), ei_scalar_clip_under_op<Scalar>(under));
}

inline const EIGEN_CWISE_UNOP_RETURN_TYPE(ei_scalar_clip_op)
clip(const Scalar& under, const Scalar& over) const
{
  return EIGEN_CWISE_UNOP_RETURN_TYPE(ei_scalar_clip_op)(_expression(), ei_scalar_clip_op<Scalar>(under, over));
}

#endif // CWISEADDONS_H
