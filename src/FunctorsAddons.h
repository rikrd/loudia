#ifndef FUNCTORSADDONS_H
#define FUNCTORSADDONS_H


/**
 *
 * angle() operator
 *
 */
template<typename Scalar> struct ei_scalar_angle_op EIGEN_EMPTY_STRUCT {
  inline const Scalar operator() (const Scalar& a) const { return atan2(a.imag(), a.real()); }
};

template<typename Scalar>
struct ei_functor_traits<ei_scalar_angle_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = false }; };

/**
 *
 * ceil() operator
 *
 */
template<typename Scalar> struct ei_scalar_ceil_op EIGEN_EMPTY_STRUCT {
  inline const Scalar operator() (const Scalar& a) const { return std::ceil(a); }
};

template<typename Scalar>
struct ei_functor_traits<ei_scalar_ceil_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = false }; };

/**
 *
 * floor() operator
 *
 */
template<typename Scalar> struct ei_scalar_floor_op EIGEN_EMPTY_STRUCT {
  inline const Scalar operator() (const Scalar& a) const { return std::floor(a); }
};

template<typename Scalar>
struct ei_functor_traits<ei_scalar_floor_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = false }; };


/**
 *
 * expN(Scalar under) operator used for clipping
 *
 */
template<typename Scalar> struct ei_scalar_exp_n_op {
  // FIXME default copy constructors seems bugged with std::complex<>
  inline ei_scalar_exp_n_op(const ei_scalar_exp_n_op& other) : m_base(other.m_base) { }
  inline ei_scalar_exp_n_op(const Scalar& base) : m_base(base) {}
  inline Scalar operator() (const Scalar& a) const { 
    return ei_pow(m_base, a); 
  }
  const Scalar m_base;
};
template<typename Scalar>
struct ei_functor_traits<ei_scalar_exp_n_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = false }; };


/**
 *
 * logN(Scalar under) operator used for clipping
 *
 */
template<typename Scalar> struct ei_scalar_log_n_op {
  // FIXME default copy constructors seems bugged with std::complex<>
  inline ei_scalar_log_n_op(const ei_scalar_log_n_op& other) : m_log_base(other.m_log_base) { }
  inline ei_scalar_log_n_op(const Scalar& base) : m_log_base(ei_log(base)) {}
  inline Scalar operator() (const Scalar& a) const { 
    return ei_log(a) / m_log_base; 
  }
  const Scalar m_log_base;
};
template<typename Scalar>
struct ei_functor_traits<ei_scalar_log_n_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = false }; };


/**
 *
 * clipUnder(Scalar under) operator used for clipping
 *
 */
template<typename Scalar> struct ei_scalar_clip_under_op {
  // FIXME default copy constructors seems bugged with std::complex<>
  inline ei_scalar_clip_under_op(const ei_scalar_clip_under_op& other) : m_under(other.m_under) { }
  inline ei_scalar_clip_under_op(const Scalar& under) : m_under(under) {}
  inline Scalar operator() (const Scalar& a) const { 
    return a < m_under ? m_under : a; 
  }
  const Scalar m_under;
};
template<typename Scalar>
struct ei_functor_traits<ei_scalar_clip_under_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = false }; };


/**
 *
 * clipOver(Scalar over) operator used for clipping
 *
 */
template<typename Scalar> struct ei_scalar_clip_over_op {
  // FIXME default copy constructors seems bugged with std::complex<>
  inline ei_scalar_clip_over_op(const ei_scalar_clip_over_op& other) : m_over(other.m_over) { }
  inline ei_scalar_clip_over_op(const Scalar& over) : m_over(over) {}
  inline Scalar operator() (const Scalar& a) const { 
    return a > m_over ? m_over : a; 
  }
  const Scalar m_over;
};
template<typename Scalar>
struct ei_functor_traits<ei_scalar_clip_over_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = false }; };


/**
 *
 * clip(Scalar under, Scalar over) operator used for clipping
 *
 */
template<typename Scalar> struct ei_scalar_clip_op {
  // FIXME default copy constructors seems bugged with std::complex<>
  inline ei_scalar_clip_op(const ei_scalar_clip_op& other) : m_under(other.m_under), m_over(other.m_over) { }
  inline ei_scalar_clip_op(const Scalar& under, const Scalar& over) : m_under(under), m_over(over) {}
  inline Scalar operator() (const Scalar& a) const { 
    return (a < m_under ? m_under : (a > m_over ? m_over : a)); 
  }
  const Scalar m_under;
  const Scalar m_over;
};
template<typename Scalar>
struct ei_functor_traits<ei_scalar_clip_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = false }; };


#endif // FUNCTORSADDONS_H
