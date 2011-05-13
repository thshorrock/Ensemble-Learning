#pragma once

//#include "symbolic/visitor/visitor.hpp"

#include "MathsExp.hpp"


#include<boost/smart_ptr.hpp>
#include<string>

namespace ICR{

  namespace maths{
    
    //forward declaration;
    class context;

    class VariableExp : public MathsExp
    {
    public:
      typedef boost::shared_ptr<MathsExp> MathsPtr;
      
      VariableExp(const std::string& name);
      
      ~VariableExp();
      
      //virtual void accept(visitor&) = 0;
      double   evaluate(const context&);
      MathsPtr replace(const std::string&, MathsExp&);
      MathsPtr copy() const ;     

    private:
      const std::string m_name;
      friend class context;
    };



  }
}


