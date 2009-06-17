/*                                                         
** Copyright (C) 2008, 2009 Ricard Marxer <email@ricardmarxer.com>
**                                                                  
** This program is free software; you can redistribute it and/or modify
** it under the terms of the GNU General Public License as published by
** the Free Software Foundation; either version 3 of the License, or   
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

#ifndef DEBUG_H
#define DEBUG_H

#include <sstream>
#include <cassert>
#include <exception>
#include <string.h>

/**
 * Exception class for Loudia. It has a whole slew of different constructors
 * to make it as easy as possible to throw an exception with a descriptive
 * message.
 */
class LoudiaException : public std::exception {

 public:
  LoudiaException(const char* msg) : exception(), _msg(msg) {}
  LoudiaException(const std::string& msg) : exception(), _msg(msg) {}
  LoudiaException(const std::ostringstream& msg) : exception(), _msg(msg.str()) {}

  template <typename T, typename U>
  LoudiaException(const T& a, const U& b) : exception() {
    std::ostringstream oss; oss << a << b; _msg = oss.str();
  }

  template <typename T, typename U, typename V>
  LoudiaException(const T& a, const U& b, const V& c) : exception() {
    std::ostringstream oss; oss << a << b << c; _msg = oss.str();
  }

  virtual ~LoudiaException() throw() {}
  virtual const char* what() const throw() { return _msg.c_str(); }

 protected:
  std::string _msg;

};

#define LOUDIA_ERROR(msg) ostringstream loudiaErrorMessage__;  loudiaErrorMessage__ << msg; throw LoudiaException(loudiaErrorMessage__);
#define LOUDIA_WARNING(msg) std::cerr << msg << std::endl;


// If none of the NO_DEBUG are defined we enable debugging
#if (!defined(LOUDIA_NO_DEBUG) && !defined(NDEBUG))

#include <iostream>
#define LOUDIA_DEBUG(msg) std::cerr << msg << std::endl;

// else we do nothing
#else

#define LOUDIA_DEBUG(msg)

#endif


#endif  /* DEBUG_H */
