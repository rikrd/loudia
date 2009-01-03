#ifndef MATRIXBASEADDONS_H
#define MATRIXBASEADDONS_H

/*

inline Scalar at(uint i, uint j) const { return this->operator()(i,j); }
inline Scalar& at(uint i, uint j) { return this->operator()(i,j); }
inline Scalar at(uint i) const { return this->operator[](i); }
inline Scalar& at(uint i) { return this->operator[](i); }

inline RealScalar squaredLength() const { return squaredNorm(); }
inline RealScalar length() const { return norm(); }
inline RealScalar invLength(void) const { return fast_inv_sqrt(squaredNorm()); }

template<typename OtherDerived>
inline Scalar squaredDistanceTo(const MatrixBase<OtherDerived>& other) const
{ return (derived() - other.derived()).squaredNorm(); }

template<typename OtherDerived>
inline RealScalar distanceTo(const MatrixBase<OtherDerived>& other) const
{ return ei_sqrt(derived().squaredDistanceTo(other)); }

inline void scaleTo(RealScalar l) { RealScalar vl = norm(); if (vl>1e-9) derived() *= (l/vl); }

inline Transpose<Derived> transposed() {return transpose();}
inline const Transpose<Derived> transposed() const {return transpose();}

inline uint minComponentId(void) const  { int i; minCoeff(&i); return i; }
inline uint maxComponentId(void) const  { int i; maxCoeff(&i); return i; }

template<typename OtherDerived>
void makeFloor(const MatrixBase<OtherDerived>& other) { derived() = derived().cwise().min(other.derived()); }
template<typename OtherDerived>
void makeCeil(const MatrixBase<OtherDerived>& other) { derived() = derived().cwise().max(other.derived()); }

const typename Cwise<Derived>::ScalarAddReturnType
operator+(const Scalar& scalar) const { return cwise() + scalar }

friend const typename Cwise<Derived>::ScalarAddReturnType
operator+(const Scalar& scalar, const MatrixBase<Derived>& mat) { return mat + scalar; }

*/

#include "CwiseUnaryOps.h"
const CwiseUnaryOp< CwiseMaxOp< Scalar >, Derived > max(const Scalar s) { return derived().unaryExpr( CwiseMaxOp< Scalar >( s ) ); }



#endif // MATRIXBASEADDONS_H
