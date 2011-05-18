#pragma once

#include "ICA/node/Node.hpp"
#include "ICA/message/NaturalParameters.hpp"
#include "ICA/message/Moments.hpp"
#include "ICA/calculation_tree/Expression.hpp"
#include "ICA/calculation_tree/Context.hpp"

// #include "ICA/Message.hpp"

//#include "ICA/piping/Piping.hpp"
#include "ICA/node/variable/Dirichlet.hpp"
#include "ICA/node/variable/Hidden.hpp"
#include "ICA/node/variable/Calculation.hpp"
#include "ICA/exponential_model/Gaussian.hpp"
#include "ICA/exponential_model/Gamma.hpp"
#include "ICA/exponential_model/Discrete.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/none.hpp>
#include <vector>
#include <boost/assert.hpp> 
#include "rng.hpp"


namespace ICR{
  namespace ICA{
    namespace Details{
      
      
      template<template<class> class Model,  class T>
      class CalcGaussianFactor : public FactorNode<double>
      {
      public:
      
      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;

	CalcGaussianFactor( Expression<T>* Expr,  Context<T>& context,
			    DeterministicNode<Model<T>,T>* Child)
	  : m_expr(Expr),
	    m_context(context),
	    m_child_node(Child)
	{

	  context.AddChildFactor(this);
	  Child->SetParentFactor(this);
	  
	};
      
	Moments<T>
	InitialiseMoments() const
	{
	  //deterministic
	  return Model<T>::CalcMoments(Model<T>::CalcNP2Deterministic(m_expr,m_context));
	}
      
	NaturalParameters<T>
	GetNaturalNot( variable_parameter v) const
	{
	  if (v == m_child_node) 
	    {
	      return Model<T>::CalcNP2Deterministic(m_expr,m_context);
	    }
	  else
	    {
	      //parent node
	      return Model<T>::CalcNP2Parent(v,m_child_node,m_context);
	    }
	}
	T
	CalcLogNorm() const {return 0;}
      private: 
	Expression<T>* m_expr;
	Context<T> m_context;
	mutable DeterministicNode<Model<T>,T> *m_child_node;
      
      };
    

    }
    
  }
}
