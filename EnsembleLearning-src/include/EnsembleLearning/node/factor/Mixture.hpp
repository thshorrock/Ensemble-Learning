

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

#pragma once
#ifndef FACTOR_MIXTURE_HPP
#define FACTOR_MIXTURE_HPP


#include "EnsembleLearning/node/Node.hpp"
#include "EnsembleLearning/message/NaturalParameters.hpp"
#include "EnsembleLearning/message/Moments.hpp"
#include "EnsembleLearning/detail/MixtureVector.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/none.hpp>
#include <boost/mpl/inherit_linearly.hpp>
#include <boost/mpl/inherit.hpp>
#include <boost/mpl/size.hpp>
#include <vector>



namespace ICR{
  namespace EnsembleLearning{

    //forward declarations
    
    template<class T>
    class Discrete;
    
    template <template<class> class Model, class T,class List, class Enable>
    class HiddenNode; //only used as pointer here.
    
    namespace detail{

      /******************************************************************************
       * Default Specialisation
       ******************************************************************************/
      /** Factor Linking Mixture nodes */
      template<class Model, class T,
	       class parent1_t, class parent2_t, class child_t
	       >
      class Mixture : public FactorNode_basic,
		      public boost::mpl::inherit_linearly<typename parent1_t::type,
							  typename boost::mpl::inherit<boost::mpl::_1,
										       FactorNode<T,boost::mpl::_2> 
										       >::type 
		      >::type,
		      public boost::mpl::inherit_linearly<typename parent2_t::type,
							    typename boost::mpl::inherit<boost::mpl::_1,
											 FactorNode<T,boost::mpl::_2> 
											 >::type 
											 >::type, 
											   public FactorNode<T,child_t>	,							   public FactorNode<T,HiddenNode<Discrete,T > >
	
      {
	//Non-copieable
	//Mixture(const Mixture<Model, T>& f) {};
      public:
      
	/** @name Useful typdefs for types that are exposed to the user.
	 */
	///@{

	typedef typename boost::call_traits< VariableNode<T>* const>::param_type
	variable_parameter;
	typedef typename boost::call_traits< VariableNode<T>* const>::value_type
	variable_t;

      
	typedef typename boost::call_traits< parent1_t>::param_type
	p1_parameter;
	typedef typename boost::call_traits< parent1_t>::value_type
	p1_t;
	typedef typename boost::call_traits< parent2_t>::param_type
	p2_parameter;
	typedef typename boost::call_traits< parent2_t>::value_type
	p2_t;
	typedef typename boost::call_traits< child_t* const>::param_type
	c_parameter;
	typedef typename boost::call_traits< child_t* const>::value_type
	c_t;

	typedef typename boost::call_traits< std::vector<VariableNode<T>*> >::param_type
	variable_vector_parameter;
	typedef typename boost::call_traits< std::vector<VariableNode<T>*> >::value_type
	variable_vector_t;
	typedef typename boost::call_traits< Moments<T> >::value_type
	moments_t;
	typedef typename boost::call_traits< Moments<T> >::const_reference
	moments_const_reference;

	typedef typename boost::call_traits<HiddenNode<Discrete,T >* const>::value_type
	discrete_t;
	typedef typename boost::call_traits<HiddenNode<Discrete,T >* const>::value_type
	discrete_parameter;

	///@}
      
	/** Construct a mixture node.
	 * @param Parent1 A vector of parent Nodes (for example Means)
	 * @param Parent2 A vector of parent Nodes (for example Precisions)
	 * @param Weights The Weights on each of the vectors 
	 *   (The size of the weights need to be the same as sthe size of the vectors)
	 * @param child The child node.
	 */
	Mixture(p1_parameter Parent1, 
		p2_parameter Parent2,  
		discrete_parameter Weights,
		c_parameter child
	      
		)
	  :  m_parent1_nodes(Parent1), 
	     m_parent2_nodes(Parent2),
	     m_weights_node(Weights),
	     m_child_node(child),
	     m_LogNorm(0)
	{
	  //Need to be as many parents to both.
	  // BOOST_MPL_ASSERT_RELATION( mpl::size<p1_parameter::result_t::type>::value, ==,mpl::size<p2_parameter::result_t::type>::value  );
	  //Parent1.size() == Parent2.size());
	
	  child->SetParentFactor(this);
	  //boost::fusion::for_each(p1, boost::bin
	    

	//Add this factor to each parent
	//(Iterate with the preprocessor)
#       include <boost/preprocessor/iteration/iterate.hpp>
#       define BOOST_PP_ITERATION_LIMITS (0, ENSEMBLE_LEARNING_COMPONENTS - 1  )
#       define BOOST_PP_FILENAME_1       "EnsembleLearning/node/factor/Mixture_AddChildFactor.hpp"
#       include BOOST_PP_ITERATE()


	    

	  Weights->AddChildFactor(this);
	};
      
      
	Moments<T>
	InitialiseMoments() const
	{

	  std::vector<moments_t > moments1(boost::mpl::size<typename parent1_t::type>::value);
	  std::vector<moments_t > moments2(boost::mpl::size<typename parent2_t::type>::value);

	//Get the moments of each parent
	//(Iterate with the preprocessor)
#       include <boost/preprocessor/iteration/iterate.hpp>
#       define BOOST_PP_ITERATION_LIMITS (0, ENSEMBLE_LEARNING_COMPONENTS - 1 )
#       define BOOST_PP_FILENAME_1       "EnsembleLearning/node/factor/Mixture_InitialiseMoments.hpp"
#       include BOOST_PP_ITERATE()
	  
	  m_weights_node->InitialiseMoments();
	  return Model::CalcSample(moments1,
	  			 moments2, 
	  			 m_weights_node ->GetMoments()
	  			 );
	}

	//Request From the child Node
	NaturalParameters<T>
	GetNaturalNot( c_parameter v) const
	{
	  NaturalParameters<T> NP2Child(std::vector<T>(2));
	  m_LogNorm = 0;
	  
	  const size_t size = boost::mpl::size<typename parent1_t::type>::value;
	  std::vector<moments_t > moments1(size);
	  std::vector<moments_t > moments2(size);

	//Get the moments of each parent
	//(Iterate with the preprocessor)
#       include <boost/preprocessor/iteration/iterate.hpp>
#       define BOOST_PP_ITERATION_LIMITS (0, ENSEMBLE_LEARNING_COMPONENTS - 1 )
#       define BOOST_PP_FILENAME_1       "EnsembleLearning/node/factor/Mixture_InitialiseMoments.hpp"
#       include BOOST_PP_ITERATE()

	  const Moments<double>& weights = m_weights_node->GetMoments();
	  for(size_t i=0;i<size;++i){
	    NP2Child
	      += Model::CalcNP2Data(moments1[i], moments2[i]) * weights[i];
	    
	    m_LogNorm 
	      += Model::CalcLogNorm(moments1[i],moments2[i]) * weights[i];
	  }
	  return NP2Child;
	}

	//Request From the weights Node
	NaturalParameters<T>
	GetNaturalNot( discrete_parameter v) const
	{
	  
	  const size_t size = boost::mpl::size<typename parent1_t::type>::value;
	  std::vector<moments_t > moments1(size);
	  std::vector<moments_t > moments2(size);
	  NaturalParameters<T> NP2Weights(size);

	//Get the moments of each parent
	//(Iterate with the preprocessor)
#       include <boost/preprocessor/iteration/iterate.hpp>
#       define BOOST_PP_ITERATION_LIMITS (0, ENSEMBLE_LEARNING_COMPONENTS - 1 )
#       define BOOST_PP_FILENAME_1       "EnsembleLearning/node/factor/Mixture_InitialiseMoments.hpp"
#       include BOOST_PP_ITERATE()

	  const Moments<T>& child = m_child_node->GetMoments();
	  for(size_t i=0;i<size;++i){
	    NP2Weights[i] = Model::CalcAvLog(moments1[i], moments2[i],child);
	  }
	  return NP2Weights;
	}

	//Define GetNaturalNot for each of the parent nodes.
	//(Iterate with the preprocessor)
#       include <boost/preprocessor/iteration/iterate.hpp>
#       define BOOST_PP_ITERATION_LIMITS (0, ENSEMBLE_LEARNING_COMPONENTS  - 1)
#       define BOOST_PP_FILENAME_1       "EnsembleLearning/node/factor/Mixture_GetNaturalNot.hpp"
#       include BOOST_PP_ITERATE()


	T
	CalcLogNorm() const 
	{
	  return m_LogNorm;
	}


	/** Obtain the natural parameter destined for the variable_parameter v.
	 * @param v A pointer to the  VariableNode for which the message is destined.
	 *  The message is calculated from the moments of every node adjacent to the factor withe exception of v.
	 * @return The natural parameter calculated for v.
	 */
	NaturalParameters<T>
	GetNaturalNot( variable_parameter v) const;
      
      private: 
      
	p1_t m_parent1_nodes;
	p2_t m_parent2_nodes;
	discrete_t m_weights_node;
	c_t  m_child_node;
      
	mutable T m_LogNorm;
      };
    
    }
  }
}
#endif  // guard for FACTOR_MIXTURE_HPP
  
