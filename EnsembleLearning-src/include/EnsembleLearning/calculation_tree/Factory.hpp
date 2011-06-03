

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



#pragma once
#ifndef FACTORY_HPP
#define FACTORY_HPP


//Need these (rather than forward declarations)
//  as need to know the inheritence relationship
#include "Placeholder.hpp"
#include "Functions.hpp"
#include "Expression.hpp"

#include <boost/smart_ptr.hpp>
#include <boost/call_traits.hpp>
#include<list>

namespace ICR{
  namespace EnsembleLearning{

    
    /** @defgroup Calculation Evaluate the Calculation Nodes.
     *  
     */

    /** A Factory to help contruct expressions.
     *  This factory should be the user's sole interface with placeholders and mathmatical operations.
     *
     *  Example:  
     *  @code
     *  //Create an expression X*Y + Z.
     *  ExpressionFactory<double> Factory;  // Create a Factory.
     *  
     *  Placeholder<double>* X = Factory.placeholder();
     *  Placeholder<double>* Y = Factory.placeholder();
     *  Placeholder<double>* Z = Factory.placeholder();
     *  Expression<double>* XY = Factory.Multiply(X,Y);
     *  Expression<double>* Expr = Factory.Add(XY,Z);
     *  @endcode
     *
     *  A context is then provided to this expression with the Context class.
     *  The expression and context is then passed to a Calculation Variable.
     *
     * @tparam T The data type.  This will be either float or double.
     *
     * @attention  The Factory looks after the memory management of the placeholders 
     *  and expression for you.  
     *  The expression will be destroyed when the factory object is destroyed.  
     *  Do not use the delete keyword.
     *  
     * 
     *  @ingroup UserInterface
     *  @ingroup Calculation
     */
    template<class T>
    class ExpressionFactory{
    public:
      
      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{
      typedef typename boost::call_traits<Expression<T>*>::param_type
      expression_parameter;
      
      typedef typename boost::call_traits<Expression<T>*>::value_type
      expression_t;
      
      typedef typename boost::call_traits<Placeholder<T>*>::value_type
      placeholder_t;
      
      ///@}

      /** Create a placeholder
       *  @return A pointer to a placeholder.
       */
      placeholder_t
      placeholder();
      
      /** Add two sub expressions together.
       *  @param a The first half of the expression.
       *  @param b The second half of the expression.
       *  @return The resultant sum.
       */
      expression_t
      Add(expression_parameter a, 
	  expression_parameter b);
      
      /** Multiply two sub expressions together.
       *  @param a The first half of the expression.
       *  @param b The second half of the expression.
       *  @return The resultant product.
       */
      expression_t
      Multiply(expression_parameter a, 
	       expression_parameter b);
      

    private:
      //References to the expression - memory management hanled by the shared pointers.
      std::list<boost::shared_ptr<Expression<T> > >  m_Expr;
    };
  }
}



#endif  // guard for FACTORY_HPP
