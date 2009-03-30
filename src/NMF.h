/*                                                         
** Copyright (C) 2008, 2009 Ricard Marxer <email@ricardmarxer.com>
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

#ifndef NMF_H
#define NMF_H

#include "Typedefs.h"
#include "Debug.h"

/**
  * @class NMF
  *
  * @brief Algorithm to calculate the Non-negative Matrix Factorization of a 
  * matrix of Real values.
  *
  * This class represents an object to perform a 
  * Non-negative Matrix Factorization (NMF) on a matrix 
  * of Real values.  Which is a useful technique for decomposing multichannel
  * signals into linear independent components and their respective positive gains.
  * 
  * The algorithm estimates a set of M temporal independent components, each of which
  * is represented as a single row in a matrix of the same number of columns as the
  * input matrix.
  * The algorithm also estimates a set of M time-varying positive gains, 
  * each of which is represented as a single column in a matrix with as many rows
  * as the input matrix.
  *
  * In order to implement the NMF an iterative multiplicative update rule is applied in order
  * to minimize the distance from the multiplication of the gains and components matrices
  * to the input matrix.
  *
  * The distance minimized in this implementation is the Kullback-Liebler divergence.
  *
  * @author Ricard Marxer
  *
  * @sa INMF, ICA
  */
class NMF {
protected:
  // Internal parameters
  int _inputSize;
  int _componentCount;

  int _iterationCount;
  
  Real _epsilon;
  
  // Internal variables
  MatrixXR _xOverWH;
  MatrixXR _norms;

public:
  /**
     Constructs an NMF object with the specified @a inputSize, 
     @a componentCount, @a iterationCount and @a epsilon settings.
     
     @param inputSize size of the inputs arrays,
     must be > 0.
     
     @param componentCount number of linearly independnet components to be estimated.
     @param iterationCount number of update iterations performed to solve the NMF.
     @param epsilon parameter of the solver to clip the minimum values at each update.
     
  */
  NMF(int inputSize = 1024, int componentCount = 3, int iterationCount = 10, Real epsilon = 1e-9);

  /**
     Destroys the algorithm and frees its resources.
  */
  ~NMF();

  void reset();
  void setup();

  /**
     Performs an NMF on @a frames which represents multichannel signal where the rows
     are the time axis and the columns are the channels.
     Puts the resulting NMF components in the rows of @a components and the
     instantaneous gains of each components in @a gains.
     
     @param frames matrix of Real values.  The number of columns of @a 
     frames must be equal to the input size specified using setInputSize().
     
     @param components pointer to a matrix of Real values for the components.
     The matrix should have the same number of rows as componentCount and inputSize columns.

     @param gains pointer to a matrix of Real positive values for the instantaneous gains 
     of each component.
     The matrix should have the same number of rows as frames and componentCount columns.

     Note that if the output matrices are not of the required sizes they will be resized, 
     reallocating a new memory space if necessary.
  */
  void process(const MatrixXR& frames, MatrixXR* gains, MatrixXR* components);

  /**
     Returns the size of the input arrays.
     The default is 1024.
     
     @sa setInputSize()
  */
  int inputSize() const;  

  /**
     Specifies the @a size of the input.
     The given @a size must be higher than 0.
     
     @sa inputSize()
  */
  void setInputSize( int size, bool callSetup = true );

  /**
     Returns the number of components to be calculated.
     The default is 3.
     
     @sa setComponentCount()
  */
  int componentCount() const;

  /**
     Specifies the @a count of components to be calculated.
     The given @a count must be greater than 0.
          
     @sa componentCount()
  */
  void setComponentCount( int count, bool callSetup = true );

  /**
     Returns the number of iterations of the solver.
     The default is 10.
     
     @sa setIterationCount()
  */
  int iterationCount() const;

  /**
     Specifies the @a count of iterations of the solver.
     The given @a count must be greater than 0.
          
     @sa iterationCount()
  */
  void setIterationCount( int count, bool callSetup = true );

  /**
     Returns the epsilon.
     The default is 1e-6.
     
     @sa setEpsilon()
  */
  Real epsilon() const;

  /**
     Specifies the @a epsilon for the NMF solving update rule.
     
     @sa epsilon()
  */
  void setEpsilon( Real epsilon, bool callSetup = true );

};

#endif  /* NMF_H */
