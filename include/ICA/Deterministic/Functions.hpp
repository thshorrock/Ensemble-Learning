#pragma once
#include "Expression.hpp"
#include "Context.hpp"


namespace ICR{

  namespace ICA{


    template<class T>
    class Plus : public Function<T>
    {
    public:
      
      Plus( Expression<T>*  a,  Expression<T>* b);
      
      void SetParent(const Function<T>* p)
      {m_parent = p;}
      
      const Function<T>*
      GetParent() const
      {return m_parent;}

      
      FunctionName::Value
      GetFunctionType() const
      {return FunctionName::PLUS;}

      T Evaluate(const SubContext<T>&) const;
    private:
      Expression<T> *m_a, *m_b;
      const Function<T>	*m_parent;
    };
    
    template<class T>
    class Times : public Function<T>
    {
    public:
      
      Times( Expression<T>*  a,  Expression<T>* b);
      
      void
      SetParent(const Function<T>* p)
      {m_parent = p;}
      
      const Function<T>*
      GetParent() const
      {return m_parent;}
      
      FunctionName::Value
      GetFunctionType() const
      {return FunctionName::TIMES;}
      
      

      T Evaluate(const SubContext<T>&) const;
    private:
      Expression<T>  *m_a, *m_b;
      const Function<T>	*m_parent;
    };
    
  }
}

template<class T>
inline
ICR::ICA::Plus<T>::Plus(Expression<T>* a, Expression<T>* b)
  : Function<T>(), 
    m_a(a),
    m_b(b),
    m_parent(0)
{
  m_a->SetParent(this);
  m_b->SetParent(this);
}

template<class T>
inline
T
ICR::ICA::Plus<T>::Evaluate(const SubContext<T>& c) const
{
  return  
    m_a -> Evaluate(c) +
    m_b -> Evaluate(c);
}

template<class T>
inline
ICR::ICA::Times<T>::Times(Expression<T>* a, Expression<T>* b)
  : Function<T>(), 
    m_a(a),
    m_b(b),
    m_parent(0)
{
  m_a->SetParent(this);
  m_b->SetParent(this);
}

template<class T>
inline
T
ICR::ICA::Times<T>::Evaluate(const SubContext<T>& c) const
{
  return  
    m_a -> Evaluate(c) *
    m_b -> Evaluate(c);
}
