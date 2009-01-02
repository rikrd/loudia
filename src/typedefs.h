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

#include <Eigen/Core>
#include <Eigen/Array>
#include<cmath>

// Type for the Real values
typedef int Integer;
typedef float Real;
typedef std::complex< Real > Complex;

// Normally used matrices
typedef Eigen::Matrix< Integer, Eigen::Dynamic,Eigen::Dynamic > MatrixXI;
typedef Eigen::Matrix< Real, Eigen::Dynamic,Eigen::Dynamic > MatrixXR;
typedef Eigen::Matrix< Complex, Eigen::Dynamic,Eigen::Dynamic > MatrixXC;

// Type for the matrices that come from scipy (these are RowMajor)
typedef Eigen::Matrix< Real, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor, Eigen::Dynamic, Eigen::Dynamic > MatrixXRscipy;
typedef Eigen::Matrix< Complex, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor, Eigen::Dynamic, Eigen::Dynamic > MatrixXCscipy;

#endif // TYPEDEFS_H
