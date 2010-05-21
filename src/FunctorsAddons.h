#ifndef FUNCTORSADDONS_H
#define FUNCTORSADDONS_H


/*
 *
 * angle() operator
 *
 */
template<typename Scalar> struct ei_scalar_angle_op {
    EIGEN_EMPTY_STRUCT_CTOR(ei_scalar_angle_op)
    typedef typename NumTraits<Scalar>::Real result_type;
    EIGEN_STRONG_INLINE const result_type operator() (const Scalar& a) const
    {
        return atan2(a.imag(), a.real());
    }
};

template<typename Scalar>
struct ei_functor_traits<ei_scalar_angle_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = 0 }; };

/*
 *
 * sgn() operator
 *
 */
template<typename Scalar> struct ei_scalar_sgn_op {
    EIGEN_EMPTY_STRUCT_CTOR(ei_scalar_sgn_op)
    EIGEN_STRONG_INLINE const Scalar operator() (const Scalar& a) const
    {
        return (a>0 ? 1 : (a<0 ? -1 : 0));
    }
};

template<typename Scalar>
struct ei_functor_traits<ei_scalar_sgn_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = 0 }; };

/*
 *
 * ceil() operator
 *
 */
template<typename Scalar> struct ei_scalar_ceil_op {
    EIGEN_EMPTY_STRUCT_CTOR(ei_scalar_ceil_op)
    EIGEN_STRONG_INLINE const Scalar operator() (const Scalar& a) const
    {
        return std::ceil(a);
    }
};

template<typename Scalar>
struct ei_functor_traits<ei_scalar_ceil_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = 0 }; };

/*
 *
 * floor() operator
 *
 */
template<typename Scalar> struct ei_scalar_floor_op {
    EIGEN_EMPTY_STRUCT_CTOR(ei_scalar_floor_op)
    EIGEN_STRONG_INLINE const Scalar operator() (const Scalar& a) const
    {
        return std::floor(a);
    }
};

template<typename Scalar>
struct ei_functor_traits<ei_scalar_floor_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = 0 }; };


/*
 *
 * isnan() operator
 *
 */
template<typename Scalar> struct ei_scalar_isnan_op {
    EIGEN_EMPTY_STRUCT_CTOR(ei_scalar_isnan_op)
    EIGEN_STRONG_INLINE const Scalar operator() (const Scalar& a) const
    {
        return std::isnan(a);
    }
};

template<typename Scalar>
struct ei_functor_traits<ei_scalar_isnan_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = 0 }; };


/*
 *
 * modN(Scalar divisor) operator used for getting the remainder
 *
 */
template<typename Scalar> struct ei_scalar_mod_n_op {
    // FIXME default copy constructors seems bugged with std::complex<>
    EIGEN_STRONG_INLINE ei_scalar_mod_n_op(const ei_scalar_mod_n_op& other) : m_divisor(other.m_divisor) { }
    EIGEN_STRONG_INLINE ei_scalar_mod_n_op(const Scalar& divisor) : m_divisor(divisor) {}
    EIGEN_STRONG_INLINE Scalar operator() (const Scalar& a) const {
        Scalar div = a/m_divisor;
        return (div - floor(div)) * m_divisor;
    }
    const Scalar m_divisor;
};
template<typename Scalar>
struct ei_functor_traits<ei_scalar_mod_n_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = 0 }; };


/*
 *
 * expN(Scalar under) operator used for clipping
 *
 */
template<typename Scalar> struct ei_scalar_exp_n_op {
    // FIXME default copy constructors seems bugged with std::complex<>
    EIGEN_STRONG_INLINE ei_scalar_exp_n_op(const ei_scalar_exp_n_op& other) : m_base(other.m_base) { }
    EIGEN_STRONG_INLINE ei_scalar_exp_n_op(const Scalar& base) : m_base(base) {}
    EIGEN_STRONG_INLINE Scalar operator() (const Scalar& a) const {
        return ei_pow(m_base, a);
    }
    const Scalar m_base;
};
template<typename Scalar>
struct ei_functor_traits<ei_scalar_exp_n_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = 0 }; };


/*
 *
 * logN(Scalar under) operator used for clipping
 *
 */
template<typename Scalar> struct ei_scalar_log_n_op {
    // FIXME default copy constructors seems bugged with std::complex<>
    EIGEN_STRONG_INLINE ei_scalar_log_n_op(const ei_scalar_log_n_op& other) : m_log_base(other.m_log_base) { }
    EIGEN_STRONG_INLINE ei_scalar_log_n_op(const Scalar& base) : m_log_base(ei_log(base)) {}
    EIGEN_STRONG_INLINE Scalar operator() (const Scalar& a) const {
        return ei_log(a) / m_log_base;
    }
    const Scalar m_log_base;
};
template<typename Scalar>
struct ei_functor_traits<ei_scalar_log_n_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = 0 }; };


/*
 *
 * clipUnder(Scalar under) operator used for clipping
 *
 */
template<typename Scalar> struct ei_scalar_clip_under_op {
    // FIXME default copy constructors seems bugged with std::complex<>
    EIGEN_STRONG_INLINE ei_scalar_clip_under_op(const ei_scalar_clip_under_op& other) : m_under(other.m_under) { }
    EIGEN_STRONG_INLINE ei_scalar_clip_under_op(const Scalar& under) : m_under(under) {}
    EIGEN_STRONG_INLINE Scalar operator() (const Scalar& a) const {
        return a >= m_under ? a : m_under;
    }
    const Scalar m_under;
};
template<typename Scalar>
struct ei_functor_traits<ei_scalar_clip_under_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = 0 }; };


/*
 *
 * clipOver(Scalar over) operator used for clipping
 *
 */
template<typename Scalar> struct ei_scalar_clip_over_op {
    // FIXME default copy constructors seems bugged with std::complex<>
    EIGEN_STRONG_INLINE ei_scalar_clip_over_op(const ei_scalar_clip_over_op& other) : m_over(other.m_over) { }
    EIGEN_STRONG_INLINE ei_scalar_clip_over_op(const Scalar& over) : m_over(over) {}
    EIGEN_STRONG_INLINE Scalar operator() (const Scalar& a) const {
        return a <= m_over ? a : m_over;
    }
    const Scalar m_over;
};
template<typename Scalar>
struct ei_functor_traits<ei_scalar_clip_over_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = 0 }; };


/*
 *
 * clip(Scalar under, Scalar over) operator used for clipping
 *
 */
template<typename Scalar> struct ei_scalar_clip_op {
    // FIXME default copy constructors seems bugged with std::complex<>
    EIGEN_STRONG_INLINE ei_scalar_clip_op(const ei_scalar_clip_op& other) : m_under(other.m_under), m_over(other.m_over) { }
    EIGEN_STRONG_INLINE ei_scalar_clip_op(const Scalar& under, const Scalar& over) : m_under(under), m_over(over) {}
    EIGEN_STRONG_INLINE Scalar operator() (const Scalar& a) const {
        return (a < m_under ? m_under : (a > m_over ? m_over : a));
    }
    const Scalar m_under;
    const Scalar m_over;
};
template<typename Scalar>
struct ei_functor_traits<ei_scalar_clip_op<Scalar> >
{ enum { Cost = 5 * NumTraits<Scalar>::MulCost, PacketAccess = 0 }; };


#endif // FUNCTORSADDONS_H
