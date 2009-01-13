#ifndef MATRIXBASEADDONS_H
#define MATRIXBASEADDONS_H

#include "CwiseUnaryOps.h"

const CwiseUnaryOp< CwiseMaxOp< Scalar >, Derived > max(const Scalar s) const
{ return derived().unaryExpr( CwiseMaxOp< Scalar >( s ) ); }

const CwiseUnaryOp< CwiseExpOp< Scalar >, Derived > exp(const Scalar s) const
{ return derived().unaryExpr( CwiseExpOp< Scalar >( s ) ); }

const CwiseUnaryOp< CwiseAngleOp< Scalar >, Derived > angle() const
{ return derived().unaryExpr( CwiseAngleOp< Scalar >( ) ); }

const CwiseUnaryOp< CwiseCeilOp< Scalar >, Derived > ceil() const
{ return derived().unaryExpr( CwiseCeilOp< Scalar >( ) ); }

#endif // MATRIXBASEADDONS_H
