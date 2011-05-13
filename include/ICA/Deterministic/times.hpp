#pragma once

#include "symbolic/MathsExp.hpp"
#include "symbolic/context.hpp"

#include <string>

namespace ICR{

  namespace maths{
    
    class times : public MathsExp
    {
    public:
      typedef boost::shared_ptr<MathsExp> MathsPtr;
      
      times(const MathsPtr& a,const MathsPtr& b);
      
      double evaluate(const context&);

      MathsPtr replace(const std::string&, MathsExp&);
      MathsPtr copy() const ;     
    private:
      MathsPtr m_a, m_b;
    };
    
  }
}

inline
ICR::maths::times::times(const MathsPtr& a,const  MathsPtr& b)
  : MathsExp(), 
    m_a(a),
    m_b(b)
{
}

inline
double
ICR::maths::times::evaluate(const context& c)
{
  return  
    m_a -> evaluate(c) *
    m_b -> evaluate(c);
}

inline
boost::shared_ptr<ICR::maths::MathsExp>
ICR::maths::times::copy() const
{
  MathsPtr r ( new times(m_a -> copy(), m_b -> copy() ) );
  return r;
}

inline
boost::shared_ptr<ICR::maths::MathsExp>
ICR::maths::times::replace(const std::string& name, MathsExp& exp)
{
  MathsPtr r ( new times(m_a -> replace(name,exp), m_b -> replace(name,exp) ) );
  return r;
}


