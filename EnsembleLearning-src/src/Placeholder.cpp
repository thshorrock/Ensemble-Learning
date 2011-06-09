
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




#include "EnsembleLearning/calculation_tree/Placeholder.hpp"
#include <iostream>

using namespace ICR::EnsembleLearning;

template<class T>
std::pair<T,T>
Placeholder<T>::Invert(const T rhs, subcontext_parameter c) const
{
  /** returns:
   *  The rhs - the summed terms
   *  The multiplicative factor on the rest.
   */

  if  (begin()==end()){
    std::cout<<"SOMETHING WRONG IN PLACEHOLDER.HPP"<<std::endl;
    BOOST_ASSERT(1==2);//This shouldn't happen 
  }
	
  const_iterator itr = begin(); // the first function.
  //Climb down to the root
  while(itr!=end())
    ++itr;
	
  //evaluate the whole tree
  //itr.Previous because itr is now at end().
  T lhs = itr.Previous()->Evaluate(c);


  //if this placeholder is summed, then undo the summation and return
  if (begin()->GetFunctionType()==FunctionName::PLUS)
    {
      return std::pair<T,T>(rhs - (lhs - c.Lookup(this)), 1);
    }
  //Okay so we have some mulitplication.
  //Need to evaluate the product


  itr = begin();
  //Decend until we are no long multiplying
  while(itr!=end() && itr->GetFunctionType()!=FunctionName::TIMES )
    { 
      ++itr;
    }
  T prod = itr.Previous()->Evaluate(c);  //The product
  //subtract from total
  T other = lhs-prod;
  //Evaluate the divisor 
  T Factor  = prod / c.Lookup(this);
  //Return
  return std::pair<T,T>((rhs-other), Factor);
}

template<class T>
size_t ICR::EnsembleLearning::Placeholder<T>::s_count = 0;


//The types we can use
template class ICR::EnsembleLearning::Placeholder<double>;
template class ICR::EnsembleLearning::Placeholder<float>;

