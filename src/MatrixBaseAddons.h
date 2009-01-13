#ifndef MATRIXBASEADDONS_H
#define MATRIXBASEADDONS_H

#include "CwiseUnaryOps.h"

const CwiseUnaryOp< CwiseExpOp< Scalar >, Derived > exp(const Scalar s) const
{ return derived().unaryExpr( CwiseExpOp< Scalar >( s ) ); }

/*
const CwiseUnaryOp< CwiseAngleOp< Scalar >, Derived > angle() const
{ return derived().unaryExpr( CwiseAngleOp< Scalar >( ) ); }
*/

#endif // MATRIXBASEADDONS_H
