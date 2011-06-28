
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
#ifndef FUNCTIONS_HPP
#define FUNCTIONS_HPP


#include "Expression.hpp"
#include <boost/call_traits.hpp>

#include <boost/proto/proto.hpp>

namespace ICR{

  namespace EnsembleLearning{

struct Plus
{
    template <class T1, class T2>
    struct apply
    {
      typedef typename boost::proto::plus<T1,T2>::type type;
    };
};
struct Times
{
    template <class T1, class T2>
    struct apply
    {
       typedef typename boost::proto::multiplies<T1,T2>::type type;
    };
};
    //forward declare
    
    // template<class T>
    // class SubContext;

    // /** Add two functions together.
    //  *  @tparam T The data type, either float or double.
    //  */
    // template<class T>
    // class Plus : public Function<T>
    // {
    // public:
      
      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{

      // typedef typename boost::call_traits<Expression<T>*>::param_type
    //   expression_parameter;
      
    //   typedef typename boost::call_traits<Expression<T>*>::value_type
    //   expression_t;
      
    //   typedef typename boost::call_traits<T>::value_type
    //   data_t;
      
    //   typedef typename boost::call_traits<Function<T>*>::param_type
    //   function_parameter;
      
    //   typedef typename boost::call_traits<Function<T>*>::value_type
    //   function_t;

    //   typedef typename boost::call_traits<SubContext<T> >::param_type
    //   subcontext_parameter;
      
    //   ///@}

    //   /** Constructor.
    //    *  @param a The left hand side of the addition.
    //    *  @param b The right-hand-side of the addition.
    //    */
    //   Plus(expression_parameter  a,  
    // 	   expression_parameter b);
 
    //   /** Evaluate the expression for the given context.
    //    *  @param C The context for a view of the Moments provided by SubContext object.
    //    *  @return The result of the addition.
    //    */
    //   data_t Evaluate(subcontext_parameter C) const;
    // private:
    //   /* Set the parent expression.
    //    *  In order to be able to evaluate the (different) 
    //    *  expression for the message to the parent node, 
    //    *  we need a hook to be able to invert the expression.
    //    */
    //   void SetParent(function_parameter p)
    //   {m_parent = p;}
      
    //   /*  Get the parent expression.
    //    *  Need to be able to traverse down the expression tree
    //    *  to be able to evaluate the inverse expression.
    //    */
    //   function_t
    //   GetParent() const
    //   {return m_parent;}

      
    //   /*  Get the Function type 
    //    *  The inversion is different depending on the type of the function.
    //    *  This function identifies what we are doing.
    //    *  This is ugly because it couples the functions together -
    //    *  would be good to have the inversion at every step.
    //    */
    //   FunctionName::Value
    //   GetFunctionType() const
    //   {return FunctionName::PLUS;}
      
    //   //Pointers to the expressions (the children)
    //   expression_t m_a, m_b;
    //   //and the parents.
    //   function_t m_parent;
    // };
    
    // /** Multiply two functions together.
    //  *  @tparam T The data type, either float or double.
    //  */
    // template<class T>
    // class Times : public Function<T>
    // {
    // public:
      
    //   /** @name Useful typdefs for types that are exposed to the user.
    //    */
    //   ///@{

    //   typedef typename boost::call_traits<Expression<T>*>::param_type
    //   expression_parameter;
      
    //   typedef typename boost::call_traits<Expression<T>*>::value_type
    //   expression_t;
      
    //   typedef typename boost::call_traits<T>::value_type
    //   data_t;
      
    //   typedef typename boost::call_traits<Function<T>*>::param_type
    //   function_parameter;
      
    //   typedef typename boost::call_traits<Function<T>*>::value_type
    //   function_t;

    //   typedef typename boost::call_traits<SubContext<T> >::param_type
    //   subcontext_parameter;

    //   ///@}

      
    //   /** Constructor.
    //    *  @param a The left hand side of the product.
    //    *  @param b The right-hand-side of the product.
    //    */
    //   Times(expression_parameter  a, 
    // 	    expression_parameter b);
      

    //   /** Evaluate the expression for the given context.
    //    *  @param C The context for a view of the Moments provided by SubContext object.
    //    *  @return The result of the addition.
    //    */
    //   data_t 
    //   Evaluate(subcontext_parameter C) const;
    // private:


    //   /* Set the parent expression.
    //    *  In order to be able to evaluate the (different) 
    //    *  expression for the message to the parent node, 
    //    *  we need a hook to be able to invert the expression.
    //    */
    //   void
    //   SetParent(function_parameter p)
    //   {m_parent = p;}
      
    //   /*  Get the parent expression.
    //    *  Need to be able to traverse down the expression tree
    //    *  to be able to evaluate the inverse expression.
    //    */
    //   function_t
    //   GetParent() const
    //   {return m_parent;}
      
      
    //   /*  Get the Function type 
    //    *  The inversion is different depending on the type of the function.
    //    *  This function identifies what we are doing.
    //    *  This is ugly because it couples the functions together -
    //    *  would be good to have the inversion at every step.
    //    */
    //   FunctionName::Value
    //   GetFunctionType() const
    //   {return FunctionName::TIMES;}
      
      
    //   //Pointers to the expressions (the children)
    //   expression_t m_a, m_b;
    //   //and the parents.
    //   function_t m_parent;
    // };
    
  }
}

// template<class T>
// inline
// ICR::EnsembleLearning::Plus<T>::Plus(expression_parameter a,
// 			expression_parameter b)
//   : Function<T>(), 
//     m_a(a),
//     m_b(b),
//     m_parent(0)
// {
//   m_a->SetParent(this);
//   m_b->SetParent(this);
// }

// template<class T>
// inline
// typename  ICR::EnsembleLearning::Plus<T>::data_t
// ICR::EnsembleLearning::Plus<T>::Evaluate(subcontext_parameter c) const
// {
//   return  
//     m_a -> Evaluate(c) +
//     m_b -> Evaluate(c);
// }

// template<class T>
// inline
// ICR::EnsembleLearning::Times<T>::Times(expression_parameter a, 
// 			  expression_parameter b)
//   : Function<T>(), 
//     m_a(a),
//     m_b(b),
//     m_parent(0)
// {
//   m_a->SetParent(this);
//   m_b->SetParent(this);
// }

// template<class T>
// inline
// typename ICR::EnsembleLearning::Times<T>::data_t
// ICR::EnsembleLearning::Times<T>::Evaluate(subcontext_parameter c) const
// {
//   return  
//     m_a -> Evaluate(c) *
//     m_b -> Evaluate(c);
// }

#endif  // guard for FUNCTIONS_HPP
