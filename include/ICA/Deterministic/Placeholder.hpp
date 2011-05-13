#pragma once

//#include "symbolic/visitor/visitor.hpp"

#include "Functions.hpp"
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

    template<class T>
    class Placeholder : public Expression<T>
    {
    public:

      Placeholder()
	: Expression<T>(),
	  m_parent(0)
      {}
      
      T   Evaluate(const SubContext<T>& c) const
      {
	return c.Lookup(this);
      }
      
      void
      SetParent(const Function<T>* p)
      {m_parent = p;}
      
      const Function<T>*
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
	if  (m_parent==0){
	  std::cout<<"SOMETHING WRONG IN PLACEHOLDER.HPP"<<std::endl;

	  return std::pair<T,T>(rhs, 1);//This shouldn't happen 
	}
	
	root = m_parent;
	while(root->GetParent()!=0)
	  { 
	    root = root->GetParent();
	  }
	//evaluate the whole tree
	T lhs = root->Evaluate(c);
	//std::cout<<"lhs = "<<lhs<<std::endl;


	//if this placeholder is summed, then undo the summation and return
	if (m_parent->GetFunctionType()==FunctionName::PLUS)
	  {
	    return std::pair<T,T>(rhs - (lhs - c.Lookup(this)), 1);
	  }
	//Okay so we have some mulitplication.
	//Need to evaluate the product

	
	Function<T> const* p= m_parent;
	bool descend = true;
	while(descend)
	  { 
	    //std::cout<<"Descend"<<std::endl;
	    
	    Function<T> const* peek_root = p->GetParent(); //gradparent
	    if (peek_root==0)
	      descend = false;
	    else if (peek_root->GetFunctionType()!=FunctionName::TIMES) 
	      descend = false;
	    else {
	      p = peek_root;
	    }
	    
	  }
	T prod = p->Evaluate(c); //evaluate the subexpression. (the product)
	//std::cout<<"prod = "<<prod<<std::endl;

	//subtract from total
	T other = lhs-prod;
	//std::cout<<"other = "<<other<<std::endl;
	//Evaluate the divisor 
	
	
	Factor  = prod / c.Lookup(this);
	
	//std::cout<<"Factor = "<<Factor<<std::endl;
	//Return
	//return std::pair<T,T>((rhs-other)/Factor, 1);
	return std::pair<T,T>((rhs-other), Factor);
      }

    private:
      Function<T> const *m_parent;
      friend class SubContext<T>;
    };



  }
}


