

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
#ifndef EXPRESSION_HPP
#define EXPRESSION_HPP


#include <boost/proto/proto.hpp>
#include <boost/call_traits.hpp>

namespace ICR{

  namespace EnsembleLearning{
    
    // struct plus_f
    // {
    //   template <class T1, class T2>
    //   struct apply
    //   {
    // 	typedef typename boost::proto::plus<T1,T2>::type type;
    //   };
    // };
    // struct times_f
    // {
    //   template <class T1, class T2>
    //   struct apply
    //   {
    // 	typedef typename boost::proto::multiplies<T1,T2>::type type;
    //   };
    // };

    // //Forward declarations
    // template<class>  class SubContext;
    // namespace detail{
    //   template<class>  class FunctionIterator;
    // }
    // template<class T>
    // class Plus;
    
    // template<class T>
    // class Times;
    
    // template<class T>
    // class Placeholder;
    
    // template<class> class Function;
    
    // /** A struct containing the Function Names.
    //  *  These couple the functions together and so will 
    //  *  hopefully be removed in a future version.
    //  */
    // struct FunctionName{
    //   /** An enumeration indicating the type of function. */
    //   enum Value{
    // 	PLUS,  //<!-- Addition function specifier.
    // 	TIMES, //<!-- Mutliplication function specifier.
    //   };
    // };
    
    // /** An interface for every expression (and sub expression).
    //  *  @tparam T The data type: float or double.
    //  */
    // template<class T>
    // class Expression
    // {
    // public:
    //   /** @name Useful typdefs for types that are exposed to the user.
    //    */
    //   ///@{
    //   typedef typename boost::call_traits<Function<T>*>::param_type
    //   function_parameter;
      
    //   typedef typename boost::call_traits<Function<T>*>::value_type
    //   function_t;

    //   typedef typename boost::call_traits<SubContext<T> >::param_type
    //   subcontext_parameter;

    //   typedef typename boost::call_traits<T>::value_type
    //   data_t;
      
    //   ///@}

    //   /** A destructor */
    //   virtual ~Expression(){};
      
    //   /** Evaluate the expression for a given context.
    //    * @param C The context from which to evaluate the expression.
    //    * @return  The result of the expression.
    //    */
    //   virtual data_t Evaluate(subcontext_parameter C) const = 0;

    // protected:
    //   friend class Plus<T>;
    //   friend class Times<T>;
    //   friend class Placeholder<T>;
    //   template<class> friend  class detail::FunctionIterator;
      
    //   /** Set the parent to the expression.
    //    * @param fp The function that is the parent
    //    */
    //   virtual void SetParent(function_parameter fp) = 0;
      
    //   /** Get the parent to the expression.
    //    * @return The function that is the parent
    //    */
    //   virtual function_t GetParent() const = 0;
    // };

    // /** An interface to Function expressions.
    //  * Function expressions are operators - Plus, Times.
    //  * The interface returns the function type.
    //  * In a future version this will be done by the compiler directly
    //  * without the need for the function call.
    //  */
    // template<class T>
    // class Function : public Expression<T>
    // {
    // public:
    //   /** A destructor */
    //   virtual ~Function(){};
    
    //   /** Get the function type.
    //    *  @return A member of the enumeration FunctionName::Value 
    //    * that indicates the function type.
    //    */
    //   virtual
    //   FunctionName::Value
    //   GetFunctionType() const = 0;
    // };
    

  }
}

#endif  // guard for EXPRESSION_HPP
