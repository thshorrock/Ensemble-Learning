#pragma once

#include "EnsembleLearning/node/Node.hpp"
#include "EnsembleLearning/message/NaturalParameters.hpp"
#include "EnsembleLearning/message/Moments.hpp"
#include "EnsembleLearning/calculation_tree/Expression.hpp"
#include "EnsembleLearning/calculation_tree/Context.hpp"

// #include "EnsembleLearning/Message.hpp"

//#include "EnsembleLearning/piping/Piping.hpp"
#include "EnsembleLearning/node/variable/Hidden.hpp"
#include "EnsembleLearning/node/variable/Calculation.hpp"
#include "EnsembleLearning/exponential_model/Gaussian.hpp"
#include "EnsembleLearning/exponential_model/Gamma.hpp"
#include "EnsembleLearning/exponential_model/Discrete.hpp"
#include "EnsembleLearning/exponential_model/Dirichlet.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/none.hpp>
#include <vector>
#include <boost/assert.hpp> 
//#include "rng.hpp"


namespace ICR{
  namespace EnsembleLearning{
    namespace Details{
      
      
    /** A Deterministic Factor.
     *  The Calculation Nodes pass existing moments through an expression.
     *  @tparam Model  The model to use for the data data.
     *  @tparam T The data type (float or double)
     */
      template<template<class> class Model,  class T>
      class CalcGaussianFactor : public FactorNode<T>
      {
      public:
      
      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;
	
	/** Create A Deterministic factor.
	 *  @param Expr A pointer to the Expression to use to calculate the stored moments.
	 *  @param context The Context (The actual parent nodes) for the expression in this node.
	 *  @param Child The DeterministicNode That is the child to this factor.
	 */
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
