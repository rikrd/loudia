/*                                                         
** Copyright (C) 2008 Ricard Marxer <email@ricardmarxer.com>
**                                                                  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 2 of the License, or   
** (at your option) any later version.                                 
**                                                                     
** This program is distributed in the hope that it will be useful,     
** but WITHOUT ANY WARRANTY; without even the implied warranty of      
** MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the       
** GNU General Public License for more details.                        
**                                                                     
** You should have received a copy of the GNU General Public License   
** along with this program; if not, write to the Free Software         
** Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307, USA.
*/                                                                          

#include "typedefs.h"
#include "debug.h"

#include "utils.h"


/**
 * Given a matrix of polynomes (one per column)
 * returns a matrix of roots (a vector of roots per column)
 */
using namespace std;

void roots(const MatrixXR& poly, MatrixXC* result) {
  const int coeffs = poly.cols();
  
  if ( coeffs <= 1 ) {
    // Throw error about wrong input length
  }
  
  // Prepare the output
  (*result).resize(coeffs - 1, 1);

  // Build companion matrix and find its eigenvalues (the roots)
  MatrixXR companion(coeffs - 1, coeffs - 1);
  companion.corner( Eigen::TopRight, coeffs - 2, 1 ).setZero();
  companion.corner( Eigen::BottomLeft, coeffs - 2, coeffs - 2).setIdentity();
  companion.row(0) = -poly.corner( Eigen::TopRight, 1, coeffs - 1 ) / poly(0, 0);
  
  // Get the eigen values
  //(*result) = Eigen::EigenSolver<MatrixXR>(companion).eigenvalues();
}

/**
 * Given a matrix of roots (a vector of roots per column)
 * returns a matrix of polynomes (a polynome per vector of roots)
 */
template<typename InMatrixType>
void poly(const InMatrixType& roots, InMatrixType* result) {
  const int nroots = roots.cols();

  // Prepare the output
  (*result).resize(1, 1);
  (*result).setZero();
  (*result)(0, 0) = 1.0;

  //InMatrixType newA;
  //InMatrixType a = InMatrixType::Zero(1, nroots + 1);
  InMatrixType b(1, 2);

  for ( int i = 0; i < nroots; i++) {
    b << 1 , -roots(0, i);
    convolve(*result, b, result);
    //a.set(newA);
  }

  //(*result) = a;
}

void poly(const MatrixXC& roots, MatrixXC* result) {
  return poly<MatrixXC>(roots, result);
}

/**
 * Given two row matrices 
 * returns the convolution of both
 */
template<typename InMatrixType>
void convolve(const InMatrixType& a, const InMatrixType& b, InMatrixType* c) {
  const int asize = a.cols();
  const int bsize = b.cols();

  const int csize = asize + bsize - 1;
  
  const int rows = a.rows();

  if ( b.rows() != rows ) {
    // Throw an error the two inputs must have the same number of rows
    DEBUG("ERROR: the two inputs must have the same number of rows");
    return;
  }

  // Prepare the output
  (*c).resize( rows, csize );
  
  typename InMatrixType::Scalar s;

  int acol = 0;
  int bcol = 0;
  for ( int row = 0 ; row < rows ; row++ ) {

    for ( int ccol = 0 ; ccol < csize ; ccol++ ) {
      s = 0.0;
      
      int n_lo = 0 > (ccol - bsize + 1) ? 0 : ccol - bsize + 1;
      
      int n_hi = asize - 1 < ccol ? asize - 1 : ccol;
      
      acol = n_lo;
      
      bcol = ccol - n_lo;
      
      for ( int n = n_lo ; n <= n_hi ; n++ ) {
        s += a(row, acol) * b(row, bcol);
        
        acol++;
        
        bcol--;
        
      }
      
      (*c)(row, ccol) = s;
    }
  }
}

void convolve(const MatrixXC& a, const MatrixXC& b, MatrixXC* c) {
  return convolve<MatrixXC>(a, b, c);
}

void convolve(const MatrixXR& a, const MatrixXR& b, MatrixXR* c) {
  return convolve<MatrixXR>(a, b, c);
}


/**
 * Given two row matrices 
 * returns the correlation of both
 */
template<typename InMatrixType>
void correlate(const InMatrixType& _a, const InMatrixType& _b, InMatrixType* c) {
  // TODO: allow to calculate only one part of the correlation
  //       like in the case of autocorrelation where only half is needed
  // a must be the shortest and b the longuest
  const InMatrixType& a(_a.cols() > _b.cols() ? _b : _a);
  const InMatrixType& b(_a.cols() > _b.cols() ? _a : _b);
  
  const int asize = a.cols();
  const int bsize = b.cols();

  const int csize = asize + bsize - 1;
  
  const int rows = a.rows();

  if ( b.rows() != rows ) {
    // Throw an error the two inputs must have the same number of rows
    DEBUG("ERROR: the two inputs must have the same number of rows");
    return;
  }

  // Prepare the output
  (*c).resize( rows, csize );
  (*c).setZero();
  
  int minsize = min(asize, bsize);

  for (int lag = - bsize + 1; lag < asize; lag++ ) {    
    int astart = max(lag, 0);
    int bstart = max(-lag, 0);
    
    int len = min(minsize - astart, lag + bsize - astart);
    
    if (len != 0){
      (*c).col( csize - 1 - (lag + bsize - 1) ) = (a.block(0, astart, rows, len).cwise() * b.block(0, bstart, rows, len)).rowwise().sum();
    }
  }
}

void correlate(const MatrixXC& a, const MatrixXC& b, MatrixXC* c) {
  return correlate<MatrixXC>(a, b, c);
}

void correlate(const MatrixXR& a, const MatrixXR& b, MatrixXR* c) {
  return correlate<MatrixXR>(a, b, c);
}

/**
 * Reverse in place the order of the columns
 */
void reverseCols(MatrixXC* in) { 
  const int cols = (*in).cols();
  
  for(int i = 0; i < cols / 2; i++ ){
    (*in).col(i).swap((*in).col(cols - i - 1));
  }
}

void reverseCols(MatrixXR* in) {
  const int cols = (*in).cols();
  
  for(int i = 0; i < cols / 2; i++ ){
    (*in).col(i).swap((*in).col(cols - i - 1));
  }
}

/**
 * Convert from the b and a coefficients of an IIR filter to the
 * zeros, poles and gain of the filter
 */
void rowCumsum(MatrixXR* in) { 
  const int rows = (*in).rows();
  
  for(int i = 1; i < rows; i++ ){
    (*in).row(i) += (*in).row(i-1);
  }
}

void colCumsum(MatrixXR* in) { 
  const int cols = (*in).cols();
  
  for(int i = 1; i < cols; i++ ){
    (*in).col(i) += (*in).col(i-1);
  }
}

void range(Real start, Real end, int steps, MatrixXR* in){
  const Real step = (end - start) / steps;
  
  in->resize(1, steps);
  
  for (int i = 0; i<steps; i++) {
    (*in)(0, i) = i*step + start;
  }
}

void polar(const MatrixXR&  mag, const MatrixXR&  phase, MatrixXC* complex) {
  if ((mag.rows() != phase.rows()) || (mag.cols() != phase.cols())) {
    // Throw an error
  }

  (*complex).resize(mag.rows(), mag.cols());

  for(int i = 0; i < mag.rows(); i++){
    for(int j = 0; j < mag.cols(); j++){
      (*complex)(i, j) = polar(mag(i, j), phase(i, j));
    }
  }
}

void coeffsToZpk(const MatrixXR&  b, const MatrixXR&  a, MatrixXC* zeros, MatrixXC* poles, Real* gain){
  (*gain) = b(0, 0);
  MatrixXR bTemp = b;
  bTemp /= b(0, 0);
  roots(bTemp, zeros);
  roots(a, poles);
}

Real asinc(int M, Real omega) {
  return omega == 0 ? M : sin(M * omega / 2.0) / sin(omega / 2.0) ;
}


void raisedCosTransform(Real position, Real magnitude, 
                        int windowSize, int fftSize,
                        Real alpha, Real beta, 
                        MatrixXR* spectrum, int* begin, int* end, int bandwidth) {
  
  (*begin) = max((position - bandwidth / 2.0), 0.0);
  (*end) = min(ceil(position + bandwidth / 2.0 + 1), fftSize/2.0);
  
  if ( end <= begin ) {
    DEBUG("ERROR: end must be higher than begin");
    // Throw a ValueError end must be higher than begin
  }
  
  spectrum->resize(1, (*end) - (*begin));

  const Real omegaM = 2.0 * M_PI / fftSize;
  
  for ( int i = (*begin); i < (*end); i++ ) {
    
    Real omega = 2.0 * M_PI * (i - position) / fftSize;
    
    (*spectrum)(0, i - (*begin)) = magnitude * abs(alpha * asinc(windowSize, omega) + beta * (asinc(windowSize, omega - omegaM) + asinc(windowSize, omega + omegaM)));
    
  }
}

void raisedCosTransform(Real position, Real magnitude, 
                        int windowSize, int fftSize,
                        Real alpha, Real beta, 
                        MatrixXR* spectrum, int bandwidth) {
  int begin, end;

  return raisedCosTransform(position, magnitude, 
                            windowSize, fftSize,
                            alpha, beta, 
                            spectrum, &begin, &end, bandwidth);

}

void hannTransform(Real position, Real magnitude, 
                   int windowSize, int fftSize,
                   MatrixXR* spectrum, int* begin, int* end, int bandwidth) {

  return hannTransform(position, magnitude,
                       windowSize, fftSize,
                       spectrum, begin, end, bandwidth);
  
}

void hannTransform(Real position, Real magnitude, 
                   int windowSize, int fftSize,
                   MatrixXR* spectrum, int bandwidth) {

  int begin, end;

  return hannTransform(position, magnitude, 
                       windowSize, fftSize, 
                       spectrum, &begin, &end, bandwidth);
  
}


void hammingTransform(Real position, Real magnitude, 
                      int windowSize, int fftSize,
                      MatrixXR* spectrum, int* begin, int* end, int bandwidth) {

  return raisedCosTransform(position, magnitude,
                            windowSize, fftSize,
                            0.53836, 0.46164,
                            spectrum, begin, end, bandwidth);

}

void hammingTransform(Real position, Real magnitude, 
                      int windowSize, int fftSize,
                      MatrixXR* spectrum, int bandwidth) {
  
  int begin, end;
  
  return hammingTransform(position, magnitude,
                          windowSize, fftSize,
                          spectrum, &begin, &end, bandwidth);

}


void dbToMag(const MatrixXR& db, MatrixXR* mag) {
  mag->resize(db.rows(), db.cols());
  (*mag) = (db / 20.0).cwise().expN(10.0);
}

void magToDb(const MatrixXR& mag, MatrixXR* db, Real minMag) {
  db->resize(mag.rows(), mag.cols());
  (*db) = 20.0 * mag.cwise().clipUnder( minMag ).cwise().logN( 10.0 );
}

void unwrap(const MatrixXR& phases, MatrixXR* unwrapped) {
  const int nrows = phases.rows();
  const int ncols = phases.cols();
  
  unwrapped->resize(nrows, ncols);
  
  MatrixXR diff(nrows, ncols);
  diff << MatrixXR::Zero(nrows, 1), phases.block(0, 0, nrows, ncols - 1) - phases.block(0, 1, nrows, ncols - 1);

  MatrixXR upsteps = (diff.cwise() > M_PI).cast<Real>();
  MatrixXR downsteps = (diff.cwise() < -M_PI).cast<Real>();
  
  colCumsum(&upsteps);
  colCumsum(&downsteps);
  
  (*unwrapped) = phases + (2.0 * M_PI * (upsteps - downsteps));
  return;
}

void freqz(const MatrixXR& b, const MatrixXR& a, const MatrixXR& w, MatrixXC* resp) {
  const int coeffRows = max(b.rows(), a.rows());
  const int coeffCols = b.cols();
  
  if(coeffCols != a.cols()) {
    // Throw ValueError, b must have the same cols as a
  }
  
  const int nPoints = w.cols();
  
  MatrixXR k;
  range(0, coeffRows, coeffRows, &k);
 
  MatrixXC complexW(nPoints, coeffRows);

  resp->resize(nPoints, coeffCols);

  DEBUG("FREQZ: resp.shape: (" << resp->rows() << "," << resp->cols() << ")");
  DEBUG("FREQZ: complexW.shape: (" << complexW.rows() << "," << complexW.cols() << ")");

  complexW = (w.cast<Complex>().transpose() * Complex(0,-1)  * k).cwise().exp();
  
  for(int coeffCol = 0; coeffCol < coeffCols; coeffCol++ ) {
    resp->col(coeffCol) = ( complexW.block(0, 0, nPoints, b.rows()) * b.col(coeffCol) ).cwise() / \
      ( complexW.block(0, 0, nPoints, a.rows()) * a.col(coeffCol) );
  }  
}

void freqz(const MatrixXR& b, const MatrixXR& w, MatrixXC* resp) {
  MatrixXR a = MatrixXR::Ones(1,1);
  freqz(b, a, w, resp);
}

void derivate(const MatrixXR& a, MatrixXR* b) {
  const int rows = a.rows();
  const int cols = a.cols();
  
  (*b).resize(rows, cols);
  (*b).col(0).setZero();

  (*b).block(0, 1, rows, cols - 1) = a.block(0, 1, rows, cols - 1) - a.block(0, 0, rows, cols - 1);
}
