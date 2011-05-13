#pragma once
#include "Placeholder.hpp"
#include "Functions.hpp"
#include "Expression.hpp"

#include <boost/smart_ptr.hpp>
#include<list>

namespace ICR{
  namespace ICA{

    template<class T>
    class ExpressionFactory{
    public:
      
      Placeholder<T>*
      placeholder()
      {
	typedef boost::shared_ptr<Placeholder<T> > Ptr;
	
	Ptr v(new Placeholder<T>());
	m_Expr.push_back(v);
	return v.get();
      }
      
      
      typedef boost::shared_ptr<Expression<T> > ExpPtr;
      Expression<T>* 
      Add(Expression<T>* a, Expression<T>* b)
      {
	ExpPtr r(new Plus<T>(a,b));
	m_Expr.push_back(r);
	return r.get();
      }

      Expression<T>* 
      Multiply(Expression<T>* a, Expression<T>* b)
      {
	ExpPtr r(new Times<T>(a,b));
	m_Expr.push_back(r);
	return r.get();
      }

    private:
      std::list<boost::shared_ptr<Expression<T> > >  m_Expr;
    };
  }
}
