/*                                                         
** Copyright (C) 2008 Ricard Marxer <email@ricardmarxer.com.com>
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


#ifndef TYPEDEFS_H
#define TYPEDEFS_H

#define EIGEN_MATRIXBASE_PLUGIN "MatrixBaseAddons.h"

//#define EIGEN_DEFAULT_TO_ROW_MAJOR

#include "debug.h"

#include <Eigen/Core>
#include <Eigen/Array>
#include <cmath>

// Types for scalar values
typedef int Integer;
typedef float Real;
typedef std::complex< Real > Complex;

// Types for vector values
typedef Eigen::Matrix< Integer, 1, Eigen::Dynamic > RowXI;
typedef Eigen::Matrix< Real, 1, Eigen::Dynamic > RowXR;
typedef Eigen::Matrix< Complex, 1, Eigen::Dynamic > RowXC;

typedef Eigen::Matrix< Integer, Eigen::Dynamic, 1 > ColXI;
typedef Eigen::Matrix< Real, Eigen::Dynamic, 1 > ColXR;
typedef Eigen::Matrix< Complex, Eigen::Dynamic, 1 > ColXC;

// Types for matrix values
typedef Eigen::Matrix< Integer, Eigen::Dynamic, Eigen::Dynamic > MatrixXI;
typedef Eigen::Matrix< Real, Eigen::Dynamic, Eigen::Dynamic > MatrixXR;
typedef Eigen::Matrix< Complex, Eigen::Dynamic, Eigen::Dynamic > MatrixXC;

// Types for mapping Scipy matrices (these are RowMajor)
typedef Eigen::Matrix< Integer, Eigen::Dynamic, Eigen::Dynamic, Eigen::Matrix_RowMajor, Eigen::Dynamic, Eigen::Dynamic > MatrixXIscipy;
typedef Eigen::Matrix< Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::Matrix_RowMajor, Eigen::Dynamic, Eigen::Dynamic > MatrixXRscipy;
typedef Eigen::Matrix< Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::Matrix_RowMajor, Eigen::Dynamic, Eigen::Dynamic > MatrixXCscipy;

#endif // TYPEDEFS_H
