#pragma once
#ifndef ALGORITHMCHECKER_HPP
#define ALGORITHMCHECKER_HPP


/***********************************************************************************
 ***********************************************************************************
 **                                                                               **
 **  Copyright (C) 2011 Tom Shorrock <t.h.shorrock@gmail.com> 
 **                                                                               **
 **                                                                               **
 **  This program is free software; you can redistribute it and/or                **
 **  modify it under the terms of the GNU General Public License                  **
 **  as published by the Free Software Foundation; either version 2               **
 **  of the License, or (at your option) any later version.                       **
 **                                                                               **
 **  This program is distributed in the hope that it will be useful,              **
 **  but WITHOUT ANY WARRANTY; without even the implied warranty of               **
 **  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                **
 **  GNU General Public License for more details.                                 **
 **                                                                               **
 **  You should have received a copy of the GNU General Public License            **
 **  along with this program; if not, write to the Free Software                  **
 **  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  **
 **                                                                               **
 ***********************************************************************************
 ***********************************************************************************/


#include <numeric> 

namespace std{
  /** Check that not inputting integer type into accumulate
   *  This is considered an error in this program (where a double was expected
   *  @tparam InputIterator The iterator to accumulate
   *  @tparam T The type to accumulate - will fail if integer.
   *  @param first The first iterator to accumulate from.
   *  @param last the iterator to acculate to,
   *  @param init The initial value
   *  @return The accumulated value as evaluated by std::accumulate.
   */
  template<class InputIterator, class T>
  inline
  T
  accumulate_checked(InputIterator first, InputIterator last, T init )
  {
    return std::accumulate(first,last, init);
  }

  //Not implemented for integers (will not compile if called).
  template<class InputIterator>
  inline
  int
  accumulate_checked(InputIterator first, InputIterator last, int init );

}

#endif  // guard for ALGORITHMCHECKER_HPP
