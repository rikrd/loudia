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

#include "Typedefs.h"
#include "Debug.h"

#include "Utils.h"


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
void poly(const InMatrixType& _roots, InMatrixType* result) {
  
  const int rows = _roots.rows();
  const int nroots = _roots.cols();
  
  // Prepare the output
  (*result).resize(1, 1);
  (*result).setZero();
  (*result).col(0).setOnes();
  
  InMatrixType b(rows, 2);

  for ( int i = 0; i < nroots; i++) {
    b.col(0).setOnes();
    b.col(1) = -_roots.col( i );
    
    InMatrixType temp = (*result);
    
    convolve(temp, b, result);
  }
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
void correlate(const InMatrixType& _a, const InMatrixType& _b, InMatrixType* c, int _minlag, int _maxlag) {
  /*
  // TODO: allow to calculate only one part of the correlation
  //       like in the case of autocorrelation where only half is needed
  // a must be the shortest and b the longuest
  const InMatrixType& a(_a.cols() > _b.cols() ? _b : _a);
  const InMatrixType& b(_a.cols() > _b.cols() ? _a : _b);
  
  const int asize = a.cols();
  const int bsize = b.cols();

  const int minlag = max(-bsize + 1, _minlag);
  const int maxlag = min(asize, _maxlag);

  const int csize = maxlag - minlag;
  
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

  for (int lag = minlag; lag < maxlag; lag++ ) {    
    int astart = max(lag, 0);
    int bstart = max(-lag, 0);
    
    int len = min(minsize - astart, lag + bsize - astart);
    
    if (len != 0){
      (*c).col( (lag - minlag) ) = (a.block(0, astart, rows, len).cwise() * b.block(0, bstart, rows, len)).rowwise().sum();
    }
  }
  */
  const int rows = _a.rows();
  InMatrixType temp;
  convolve<InMatrixType>(_a.rowwise().reverse(), _b, &temp);
  
  // In some cases users might want some zero coeffs of the autocorrelation
  // in those cases temp will not be big enough and then we only copy a part of it
  // (from colstart we copy 'colscopy' columns)
  const int colstart = max(_a.cols(), _b.cols()) - 1 + _minlag;
  const int colscopy = min(_maxlag - _minlag, temp.cols() - colstart);

  (*c) = InMatrixType::Zero(rows, _maxlag - _minlag);
  (*c).block(0, 0, rows, colscopy) = temp.block(0, colstart, rows, colscopy);
}

void correlate(const MatrixXC& a, const MatrixXC& b, MatrixXC* c, int _minlag, int _maxlag) {
  return correlate<MatrixXC>(a, b, c, _minlag, _maxlag);
}

void correlate(const MatrixXR& a, const MatrixXR& b, MatrixXR* c, int _minlag, int _maxlag) {
  return correlate<MatrixXR>(a, b, c, _minlag, _maxlag);
}

void autocorrelate(const MatrixXR& a, MatrixXR* c,  int _minlag, int _maxlag) {
  return correlate<MatrixXR>(a, a, c, _minlag, _maxlag);
}

void autocorrelate(const MatrixXC& a, MatrixXC* c,  int _minlag, int _maxlag) {
  return correlate<MatrixXC>(a, a, c, _minlag, _maxlag);
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

void rowShift(MatrixXR* in, int num) { 
  const int rows = (*in).rows();
  
  for(int i = -num; i < rows; i++ ){
    if(i >= 0){
      (*in).row( i ).swap( (*in).row( i + num ) );
    }
  }
}

void colShift(MatrixXR* in, int num) { 
  const int cols = (*in).cols();
  
  for(int i = -num; i < cols; i++ ){
    if(i >= 0){
      (*in).col( i ).swap( (*in).col( i + num ) );
    }
  }
}

template<typename InMatrixType>
void range(Real start, Real end, int steps, int rows, InMatrixType* in){
  const Real step = (end - start) / steps;
  
  in->resize(rows, steps);
  
  for (int i = 0; i<steps; i++) {
    (*in).col(i).setConstant( (typename InMatrixType::Scalar)(i*step + start) );
  }
}

template<typename InMatrixType>
void range(Real start, Real end, int steps, InMatrixType* in){
  return range<InMatrixType>(start, end, steps, 1, in);
}

void range(Real start, Real end, int steps, int rows, MatrixXC* in){
  range<MatrixXC>(start, end, steps, rows, in);
}

void range(Real start, Real end, int steps, MatrixXC* in){
  return range(start, end, steps, 1, in);
}

void range(Real start, Real end, int steps, int rows, MatrixXR* in){
  range<MatrixXR>(start, end, steps, rows, in);
}

void range(Real start, Real end, int steps, MatrixXR* in){
  return range(start, end, steps, 1, in);
}

void range(Real start, Real end, int steps, int rows, MatrixXI* in){
  range<MatrixXI>(start, end, steps, rows, in);
}

void range(Real start, Real end, int steps, MatrixXI* in){
  return range(start, end, steps, 1, in);
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

void zpkToCoeffs(const MatrixXC& zeros, const MatrixXC& poles, Real gain, MatrixXC*  b, MatrixXC*  a){
  poly( zeros, b );
  (*b) *= gain;
  
  poly( poles, a );
}

void lowPassToLowPass(const MatrixXC& b, const MatrixXC& a, Real freq, MatrixXC*  bout, MatrixXC*  aout) {
  const int asize = a.cols();
  const int bsize = b.cols();

  const int rows = a.rows();

  const int maxsize = max(asize, bsize);    

  *aout = a;
  *bout = b;

  MatrixXR pwo;
  range(maxsize-1, -1, maxsize, rows, &pwo);

  pwo = pwo.cwise().expN( freq );
 
  int start1 = max(bsize - asize, 0);
  int start2 = max(asize - bsize, 0);
  
  for ( int i = 0; i < bsize; i++ ) {
    (*bout).col(i) *= pwo.col( start2 + i ).cwise().inverse() * pwo.col( start1 );
  }

  for ( int i = 0; i < asize; i++ ) {
    (*aout).col(i) *= pwo.col( start1 + i ).cwise().inverse() * pwo.col( start1 );
  }

}

int combination(int N, int k) {
  if ((k > N) || (N < 0) || (k < 0)) return 0;
  
  int val = 1;
  for (int i = 0; i < min(k, N-k); i++) {
    val = (int)floor((val * (N - i)) / (i + 1));
  }

  return val;
}

void bilinear(const MatrixXC& b, const MatrixXC& a, Real fs, MatrixXR*  bout, MatrixXR*  aout) {
  const int asize = a.cols();
  const int bsize = b.cols();
  const int maxsize = max(asize, bsize);
  
  const int rows = a.rows();
  
  (*aout).resize(rows, maxsize);
  (*bout).resize(rows, maxsize);
  
  MatrixXC val;
  
  for ( int j = 0; j < maxsize; j++ ) {
    val = MatrixXC::Zero(rows, 1);
    
    for ( int i = 0; i < bsize; i++ ) {
      for ( int k = 0; k < i + 1; k++ ) {
        for ( int l = 0; l < (maxsize - i); l++ ) {
          
          if((k + l) == j)
            val += combination(i, k) * combination(maxsize - i - 1, l) * b.col(bsize - 1 - i) * pow(2*fs, (Real)i) * pow((Real)-1, (Real)k);
        }
      }
    }
    
    (*bout).col(j) = val.real();
  }

  for ( int j = 0; j < maxsize; j++ ) {
    val = MatrixXC::Zero(rows, 1);
    
    for ( int i = 0; i < asize; i++ ) {
      for ( int k = 0; k < i + 1; k++ ) {
        for ( int l = 0; l < (maxsize - i); l++ ) {
          
          if((k + l) == j)
            val += combination(i, k) * combination(maxsize - i - 1, l) * a.col(asize - 1 - i) * pow((Real)2.0*fs, (Real)i) * pow((Real)-1, (Real)k);
        }
      }
    }
    
    (*aout).col(j) = val.real();
  }
  
}

void bilinear(const MatrixXC& b, const MatrixXC& a, MatrixXR*  bout, MatrixXR*  aout) {
  bilinear(b, a, 1.0, bout, aout);
}

Real asinc(int M, Real omega) {
  return omega == 0 ? 1 : sin(M * omega / 2.0) / sin(omega / 2.0) / M;
}


void raisedCosTransform(Real position, Real magnitude, 
                        int windowSize, int fftSize,
                        Real alpha, Real beta, 
                        MatrixXR* spectrum, int* begin, int* end, int bandwidth) {
  
  (*begin) = (int)max((position - bandwidth / 2.0), 0.0);
  (*end) = (int)min(ceil(position + bandwidth / 2.0 + 1), fftSize/2.0);
  
  if ( (*end) <= (*begin) ) {
    DEBUG("ERROR: end (" << (*end) << ") must be higher than begin (" << (*begin) << ")");
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

int nextPowerOf2(Real a, int factor){
  if (a == 0) return (int)pow(2.0, 1.0 + (Real)factor);
  return (int)pow(2.0, ceil(log2(a)) + (Real)factor);
}


Real gaussian(Real x, Real mu, Real fi){
  return exp(- pow((x - mu), 2.0) / (2.0 * pow(fi, 2.0))) / (fi * sqrt(2.0 * M_PI));
}

void gaussian(Real x, MatrixXR mu, MatrixXR fi, MatrixXR* result){
  (*result) = (((-mu).cwise() + x).cwise().square().cwise() / (2.0 * fi.cwise().square())).cwise().exp().cwise() / (sqrt(2.0 * M_PI) * fi);
}
