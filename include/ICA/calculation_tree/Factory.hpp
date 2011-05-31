#pragma once
#include "Placeholder.hpp"
#include "Functions.hpp"
#include "Expression.hpp"

#include <boost/smart_ptr.hpp>
#include<list>

namespace ICR{
  namespace ICA{


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
      
      typedef typename boost::call_traits<Expression<T>*>::param_type
      expression_parameter;
      
      typedef typename boost::call_traits<Expression<T>*>::value_type
      expression_t;
      
      typedef typename boost::call_traits<Placeholder<T>*>::value_type
      placeholder_t;
      
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



/*****************************************************
 *****************************************************
 ************ IMPLEMENTATION *************************
 *****************************************************
 *****************************************************/

template<class T>
typename ICR::ICA::ExpressionFactory<T>::placeholder_t
ICR::ICA::ExpressionFactory<T>::placeholder()
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
typename ICR::ICA::ExpressionFactory<T>::expression_t
ICR::ICA::ExpressionFactory<T>::Add(expression_parameter a, 
				    expression_parameter b)
{
  typedef boost::shared_ptr<Expression<T> > ExpPtr;
  ExpPtr r(new Plus<T>(a,b));
  m_Expr.push_back(r);
  return r.get();
}

template<class T>
typename ICR::ICA::ExpressionFactory<T>::expression_t
ICR::ICA::ExpressionFactory<T>::Multiply(expression_parameter a, 
					 expression_parameter b)
{
  typedef boost::shared_ptr<Expression<T> > ExpPtr;
  ExpPtr r(new Times<T>(a,b));
  m_Expr.push_back(r);
  return r.get();
}
