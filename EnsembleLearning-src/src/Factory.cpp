
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

