#pragma once


#include "Functions.hpp"
#include "FunctionsIterator.hpp"
#include "Expression.hpp"
#include "Context.hpp"
#include<string>

namespace ICR{

  namespace ICA{
    

    /**  A placeholder to a VariableNode in a function.
     *  The purpose of the placeholder is to be able to form an equation that can be used in multiple contexts.
     *  A placeholder indicates where a particular node will be in the equation.
     *  It is then linked to the Variable by a Context.
     *
     *  @tparam T The type of the stored data.  
     *   This will be either double or float - with the latter used only for memory contrained models.
     *
     *  Example:
     *  @code 
     *
     *  //Create an expression X*Y + Z.
     *  ExpressionFactory<double> Factory;  // Create a Factory.
     *  
     *  Placeholder<double>* X = Factory.placeholder();
     *  Placeholder<double>* Y = Factory.placeholder();
     *  Placeholder<double>* Z = Factory.placeholder();
     *  Expression<double>* XY = Factory.Multiply(X,Y);
     *  Expression<double>* Expr = Factory.Add(XY,Z);
     *  
     *  //Create some variables
     *  Builder builder;
     *  Builder::GaussianNode x = Builder.Gaussian(0.01,0.01);
     *  Builder::GaussianNode y = Builder.Gaussian(0.01,0.01);
     *  Builder::GaussianNode z = Builder.Gaussian(0.01,0.01);
     *  
     *  //Assign this context to the expression.
     *  Context context;
     *  context.Assign(X,x);
     *  context.Assign(Y,y);
     *  context.Assign(Z,z);
     *  
     *  @endcode
     *
     */
    template<class T>
    class Placeholder : public Expression<T>
    {
      
      
    public:
      
      /** Constructor */
      Placeholder()
	: Expression<T>(),
	  m_parent(0)
      {}
      
    private:
      
      typedef  FunctionIterator<const Function<T> > const_iterator;
      typedef  FunctionIterator<Function<T> > iterator;
      
      
      typedef typename boost::call_traits<Function<T>*>::param_type
      function_parameter;
      
      typedef typename boost::call_traits<Function<T>*>::value_type
      function_t;

      typedef typename boost::call_traits<SubContext<T> >::param_type
      subcontext_parameter;

      typedef typename boost::call_traits<T>::value_type
      data_t;

      template<class> friend class Gaussian;
      template<class> friend class RectifiedGaussian;

      /** A forward iterator pointing to the first operator (one below this placeholder).
       */
      const_iterator
      begin() const
      {
	return const_iterator(m_parent);
      }
      /** The iterator just beyond the root operator
       */
      const_iterator
      end() const
      {
	return const_iterator();
      }
      
      /** A forward iterator pointing to the first operator (one below this placeholder).
       */
      iterator
      begin() 
      {
	return const_iterator(m_parent);
      }
      
      /** The iterator just beyond the root operator
       */
      iterator
      end() 
      {
	return const_iterator();
      }
      
      /** Evaluate the expression.
       *  @param c The subcontext (i.e <x> or <x^2> of the variables that you wish to build)
       *  @see SubContext
       */
      data_t Evaluate(subcontext_parameter c) const
      {
	return c.Lookup(this);
      }
      
      /** Set the parent function.  This is needed for inversion, so that we can iterate through tree.*/
      void
      SetParent(function_parameter p)
      {m_parent = p;}
      
      /** Get the parent function.  This is needed for inversion, so that we can iterate through tree.*/
      function_t
      GetParent() const
      {return m_parent;}
      
      /** Invert the expression.
       *  The expression provided by the user is sufficient only to calclulate the messages from ParentNodes to ChildNodes.
       *  To evaluate the Message from the Child to the Parent, the expression needs to be inverted.
       *  This function inverts the expression.
       *  @param rhs The right-hand-side of the expression to be inverted.
       *  @param c The Subcontext of the expression to be inverted.
       *  @return A pair of values.
       *   The first is the rhs minus the summed terms.
       *   The second is the multiplicative factor on the rest.
       *   Both these numbers are required by the Models to invert the expression.
       */
      std::pair<T,T>
      Invert(const T rhs, const SubContext<T> c) const
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

      function_t m_parent;
      friend class SubContext<T>;
    };



  }
}


