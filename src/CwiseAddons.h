#ifndef CWISEADDONS_H
#define CWISEADDONS_H

inline const EIGEN_CWISE_UNOP_RETURN_TYPE(ei_scalar_ceil_op)
ceil() const
{ 
  return _expression(); 
}

#endif // CWISEADDONS_H
