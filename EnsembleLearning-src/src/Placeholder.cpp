
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


//The types we can use
template class ICR::EnsembleLearning::Placeholder<double>;
template class ICR::EnsembleLearning::Placeholder<float>;

