#pragma once
#ifndef FACTOR_CALCULATION_HPP
#define FACTOR_CALCULATION_HPP


/***********************************************************************************
 ***********************************************************************************
 **                                                                               **
 **  Copyright (C) 2011 Tom Shorrock <t.h.shorrock@gmail.com> 
 **                                                                               **
 **                                                                               **
 **  This program is free software; you can redistribute it and/or                **
 **  modify it under the terms of the GNU General Public License                  **
 **  as published by the Free Software Foundation; either version 2               **
 **  of the License, or (at your option) any later version.                       **
 **                                                                               **
 **  This program is distributed in the hope that it will be useful,              **
 **  but WITHOUT ANY WARRANTY; without even the implied warranty of               **
 **  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                **
 **  GNU General Public License for more details.                                 **
 **                                                                               **
 **  You should have received a copy of the GNU General Public License            **
 **  along with this program; if not, write to the Free Software                  **
 **  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  **
 **                                                                               **
 ***********************************************************************************
 ***********************************************************************************/


#include "EnsembleLearning/node/Node.hpp"
#include "EnsembleLearning/calculation_tree/Context.hpp"

#include <boost/call_traits.hpp> 


namespace ICR{
  namespace EnsembleLearning{
    //forward declare
    template<class T, int array_size> class Moments;
    template<class T, int array_size> class NaturalParameters;
    template<class T> class Expression; //only used as a pointer here
    template <class Model, class T,class List> class DeterministicNode; //only used as a pointer here
    
    namespace detail{
      
      
      /** A Deterministic Factor.
       *  The Calculation Nodes pass existing moments through an expression.
       *  @tparam Model  The model to use for the data data.
       *  @tparam T The data type (float or double)
       */
      template<template<class> class Model,  class T>
      class Deterministic : public FactorNode_basic
      {
      public:
      
	/** @name Useful typdefs for types that are exposed to the user.
	 */
	///@{

	typedef typename boost::call_traits< VariableNode<T>* const>::param_type
	variable_parameter;
	typedef typename boost::call_traits< VariableNode<T>* const>::value_type
	variable_t;
	
	///@}

	/** Create A Deterministic factor.
	 *  @param Expr A pointer to the Expression to use to calculate the stored moments.
	 *  @param context The Context (The actual parent nodes) for the expression in this node.
	 *  @param Child The DeterministicNode That is the child to this factor.
	 */
	Deterministic( Expression<T>* Expr,  
		       Context<T>& context,
		       DeterministicNode<Model<T>,T>* Child)
	  : m_expr(Expr),
	    m_context(context),
	    m_child_node(Child)
	{

	  Child->SetParentFactor(this);
	  context.AddChildFactor(this);
	  
	};
      
	Moments<T>
	InitialiseMoments() const
	{
	  //deterministic
	  return Model<T>::CalcMoments(Model<T>::CalcNP2Deterministic(m_expr,m_context));
	}
      
      /** Obtain the natural parameter destined for the variable_parameter v.
       * @param v A pointer to the  VariableNode for which the message is destined.
       *  The message is calculated from the moments of every node adjacent to the factor withe exception of v.
       * @return The natural parameter calculated for v.
       */
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
	      return Model<T>::CalcNP2Parent<v::id::value>(m_expr,m_context,m_child_node);
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
#endif  // guard for FACTOR_CALCULATION_HPP
