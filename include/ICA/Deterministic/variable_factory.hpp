#pragma once
#include "VariableExp.hpp"
#include "ConstExp.hpp"

namespace ICR{
  namespace maths{

    class variable_factory{
    public:
      typedef boost::shared_ptr<VariableExp> VarPtr;
      typedef boost::shared_ptr<ConstExp> ConstPtr;

      VarPtr
      make(const std::string& name)
      {
	VarPtr v(new VariableExp(name));
	return v;
      }
      
      
      ConstPtr
      make(const double& d)
      {
	ConstPtr v(new ConstExp(d));
	return v;
      }

    };
  }
}
