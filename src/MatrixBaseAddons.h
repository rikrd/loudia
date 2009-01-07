#ifndef MATRIXBASEADDONS_H
#define MATRIXBASEADDONS_H

#include "CwiseUnaryOps.h"
const CwiseUnaryOp< CwiseMaxOp< Scalar >, Derived > max(const Scalar s) const
{ return derived().unaryExpr( CwiseMaxOp< Scalar >( s ) ); }

const CwiseUnaryOp< CwiseAngleOp< Scalar >, Derived > angle() const
{ return derived().unaryExpr( CwiseAngleOp< Scalar >( ) ); }

#endif // MATRIXBASEADDONS_H
