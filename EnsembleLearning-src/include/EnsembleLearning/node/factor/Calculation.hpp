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
//#include "EnsembleLearning/node/variable/Calculation.hpp"
#include "EnsembleLearning/detail/MACRO_defaults.hpp"//for ENSEMBLE_LEARNING_COMPONENTS
#include "EnsembleLearning/calculation_tree/Context.hpp"
#include "EnsembleLearning/detail/TypeList.hpp"

#include <boost/call_traits.hpp> 

#include <boost/shared_ptr.hpp>
#include <boost/none.hpp>
#include <boost/mpl/inherit_linearly.hpp>
#include <boost/mpl/inherit.hpp>
#include <boost/mpl/size.hpp>
#include <boost/fusion/include/for_each.hpp>

namespace ICR{
  namespace EnsembleLearning{
    //forward declare
    template<class T, int array_size> class Moments;
    template<class T, int array_size> class NaturalParameters;
    template<class T> class Expression; //only used as a pointer here
    template <template<class> class Model, class T,class List,int array_size,class Enable> class DeterministicNode; //only used as a pointer here
    
    namespace detail{
    
      /** A Deterministic Factor.
       *  The Calculation Nodes pass existing moments through an expression.
       *  @tparam Model  The model to use for the data data.
       *  @tparam T The data type (float or double)
       */
      // template<template<class> class Model, class T,
      // 	       class child_t>
      template<template<class> class Model,  class T, class context_t, class Expr_t>
      class Deterministic : public FactorNode_basic,
			    public boost::mpl::inherit_linearly<typename context_t::type,
								typename boost::mpl::inherit<boost::mpl::_1,
											     FactorNode<T,typename boost::mpl::_2> 
											     >::type 
			    >::type,				
			    public ParentFactorNode<T,DeterministicNode<Model,T,detail::TypeList::zeros> >
    
	
      {
      public:
      
	/** @name Useful typdefs for types that are exposed to the user.
	 */
	///@{

	typedef typename boost::call_traits< VariableNode<T>* const>::param_type
	variable_parameter;
	typedef typename boost::call_traits< VariableNode<T>* const>::value_type
	variable_t;
	
	typedef typename boost::call_traits< DeterministicNode<Model,T,detail::TypeList::zeros>*  const>::param_type
	deterministic_parameter;
	typedef typename boost::call_traits< DeterministicNode<Model,T,detail::TypeList::zeros>*  const>::value_type
	deterministic_t;
	///@}

	/** Create A Deterministic factor.
	 *  @param Expr A pointer to the Expression to use to calculate the stored moments.
	 *  @param context The Context (The actual parent nodes) for the expression in this node.
	 *  @param Child The DeterministicNode That is the child to this factor.
	 */
	Deterministic( Expr_t Expr,  
		       context_t& context,
		       deterministic_parameter Child)
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
	GetNaturalNot(deterministic_parameter  v) const
	{
	  return Model<T>::CalcNP2Deterministic(m_expr,m_context);
	}

	/* Define GetNaturalNot for each of the parent nodes.
	 * (Iterate with the preprocessor)
	 * This is done with an ugly preprocessor hack since these functions
	 * are overloads the base FactorNode class.
	 * They therefore must really be "written down" rather than just templated,
	 * which is where the preprocessor comes in.
	 */
#       include <boost/preprocessor/iteration/iterate.hpp>
#       define BOOST_PP_ITERATION_LIMITS (0, ENSEMBLE_LEARNING_PLACEHOLDERS  - 1)
#       define BOOST_PP_FILENAME_1       "EnsembleLearning/node/factor/Calculation_GetNaturalNot.hpp"
#       include BOOST_PP_ITERATE()
	/* End of pre-processor hack */
	
	// NaturalParameters<T>
	// GetNaturalNot( variable_parameter v) const
	// {
	// }

	T
	CalcLogNorm() const {return 0;}
      private: 
	context_t m_context;
	Expr_t m_expr;
	mutable deterministic_t m_child_node;
      
      };
    

    }
    
  }
}
#endif  // guard for FACTOR_CALCULATION_HPP
