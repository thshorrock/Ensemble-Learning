#pragma once

//#include "symbolic/visitor/visitor.hpp"

#include "Functions.hpp"
#include "FunctionsIterator.hpp"
#include "Expression.hpp"
#include "Context.hpp"
#include<string>

namespace ICR{

  namespace ICA{
    
    //forward declaration;

    template<class T>
    class DetConstant : public Expression<T>
    {
    public:
      T   Evaluate(const SubContext<T>& c) const
      {
	return m_value;
      }

    private:
      const T m_value;
    };


    /**
     *
     */
    template<class T>
    class Placeholder : public Expression<T>
    {
      
      
    public:
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
      
      Placeholder()
	: Expression<T>(),
	  m_parent(0)
      {}
      
      const_iterator
      begin() const
      {
	return const_iterator(m_parent);
      }
      
      const_iterator
      end() const
      {
	return const_iterator();
      }

      iterator
      begin() 
      {
	return const_iterator(m_parent);
      }
      
      iterator
      end() 
      {
	return const_iterator();
      }
      
      data_t   Evaluate(subcontext_parameter c) const
      {
	return c.Lookup(this);
      }
      
      void
      SetParent(function_parameter p)
      {m_parent = p;}
      
      function_t
      GetParent() const
      {return m_parent;}
      
      std::pair<T,T>
      Invert(const T rhs, const SubContext<T> c) const
      {
	//returns:
	//The rhs - the summed terms
	//The multiplicative factor on the rest.
	//std::cout<<"SubContext= "<<c<<std::endl;

	// T RHSminusSumed;
	T Factor;

	//Climb down to the root
	Function<T> const* root;
	// if  (m_parent==0){
	//   std::cout<<"SOMETHING WRONG IN PLACEHOLDER.HPP"<<std::endl;

	//   return std::pair<T,T>(rhs, 1);//This shouldn't happen 
	// }
	if  (begin()==end()){
	  std::cout<<"SOMETHING WRONG IN PLACEHOLDER.HPP"<<std::endl;
	  BOOST_ASSERT(1==2);
	  return std::pair<T,T>(rhs, 1);//This shouldn't happen 
	}
	

	// root = m_parent;
	// while(root->GetParent()!=0)
	//   { 
	//     root = root->GetParent();
	//   }
	const_iterator itr = begin(); // the first function.
	// const_iterator peek = begin();
	// ++peek;
	while(itr!=end())
	  ++itr;
	  //	  ++itr;
	
	//evaluate the whole tree
	//T lhs = root->Evaluate(c);
	T lhs = itr.Previous()->Evaluate(c);
	//std::cout<<"lhs = "<<lhs<<std::endl;


	//if this placeholder is summed, then undo the summation and return
	// if (m_parent->GetFunctionType()==FunctionName::PLUS)
	//   {
	//     return std::pair<T,T>(rhs - (lhs - c.Lookup(this)), 1);
	//   }
	if (begin()->GetFunctionType()==FunctionName::PLUS)
	  {
	    return std::pair<T,T>(rhs - (lhs - c.Lookup(this)), 1);
	  }
	//Okay so we have some mulitplication.
	//Need to evaluate the product

	// Function<T> const* p= m_parent;
	// bool descend = true;
	// while(descend)
	//   { 
	//     std::cout<<"Descend"<<std::endl;
	    
	//     Function<T> const* peek_root = p->GetParent(); //gradparent
	//     if (peek_root==0)
	//       descend = false;
	//     else if (peek_root->GetFunctionType()!=FunctionName::TIMES) 
	//       descend = false;
	//     else {
	//       p = peek_root;
	//     }
	    
	//   }

	// T prod = p->Evaluate(c); //evaluate the subexpression. (the product)
	// std::cout<<"prod = "<<prod<<std::endl;

	itr = begin();
	// peek = begin();
	// peek++; //look at next
	while(itr!=end() && itr->GetFunctionType()!=FunctionName::TIMES )
	  { 
	    ++itr;
	  }
	T prod = itr.Previous()->Evaluate(c); //
	//subtract from total
	T other = lhs-prod;
	//Evaluate the divisor 
	Factor  = prod / c.Lookup(this);
	//Return
	//return std::pair<T,T>((rhs-other)/Factor, 1);
	return std::pair<T,T>((rhs-other), Factor);
	  //   //std::cout<<"Descend"<<std::endl;
	  //   const_iterator peek = itr+1;
	  //   if (peek == end())
	  //     descend = false;
	  //   else if (peek->GetFunctionType()!=FunctionName::TIMES) 
	  //     descend = false;
	  //   else {
	  //     itr = peek;
	  //   }
	    
	  // }
	// while(itr->GetFunctionType()!=FunctionName::TIMES && (itr++)!=end() )
	//   {}



	//std::cout<<"other = "<<other<<std::endl;
	
	//std::cout<<"Factor = "<<Factor<<std::endl;
      }

    private:
      function_t m_parent;
      friend class SubContext<T>;
    };



  }
}


