

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




#include "EnsembleLearning/calculation_tree/Factory.hpp"

template<class T>
typename ICR::EnsembleLearning::ExpressionFactory<T>::placeholder_t
ICR::EnsembleLearning::ExpressionFactory<T>::placeholder()
{
  typedef boost::shared_ptr<Placeholder<T> > Ptr;
  //create
  Ptr v(new Placeholder<T>());
  //store
  m_Expr.push_back(v);
  //return
  return v.get();
}

template<class T>
typename ICR::EnsembleLearning::ExpressionFactory<T>::expression_t
ICR::EnsembleLearning::ExpressionFactory<T>::Add(expression_parameter a, 
				    expression_parameter b)
{
  typedef boost::shared_ptr<Expression<T> > ExpPtr;
  ExpPtr r(new Plus<T>(a,b));
  m_Expr.push_back(r);
  return r.get();
}

template<class T>
typename ICR::EnsembleLearning::ExpressionFactory<T>::expression_t
ICR::EnsembleLearning::ExpressionFactory<T>::Multiply(expression_parameter a, 
					 expression_parameter b)
{
  typedef boost::shared_ptr<Expression<T> > ExpPtr;
  ExpPtr r(new Times<T>(a,b));
  m_Expr.push_back(r);
  return r.get();
}


//The types we can use
template class ICR::EnsembleLearning::ExpressionFactory<double>;
template class ICR::EnsembleLearning::ExpressionFactory<float>;

