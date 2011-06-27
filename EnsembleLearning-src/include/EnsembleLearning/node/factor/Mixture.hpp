

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
#include <boost/fusion/include/for_each.hpp>
#include <vector>



namespace ICR{
  namespace EnsembleLearning{

    //forward declarations
    
    template<class T>
    class Discrete;
    
    template <template<class> class Model, class T,class List,int array_size, class Enable>
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
		      								       FactorNode<T,typename boost::mpl::_2> 
		      								       >::type 
		      					  >::type ,
		      public boost::mpl::inherit_linearly<typename parent2_t::type,
		      					    typename boost::mpl::inherit<boost::mpl::_1,
		      									 FactorNode<T,typename boost::mpl::_2> 
		      									 >::type 
		      >::type, 
		      	public ParentFactorNode<T,child_t>	,				
		      	public FactorNode<T,HiddenNode<Discrete,T,detail::TypeList::zeros,ENSEMBLE_LEARNING_COMPONENTS >, typename boost::mpl::int_<ENSEMBLE_LEARNING_COMPONENTS>::type >
	
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
	typedef typename boost::call_traits< Moments<T,ENSEMBLE_LEARNING_COMPONENTS> >::value_type
	weights_moments_t;
	typedef typename boost::call_traits< Moments<T,ENSEMBLE_LEARNING_COMPONENTS> >::const_reference
	weights_moments_const_reference;

	typedef typename boost::call_traits<HiddenNode<Discrete,T,detail::TypeList::zeros,ENSEMBLE_LEARNING_COMPONENTS >* const>::value_type
	discrete_t;
	typedef typename boost::call_traits<HiddenNode<Discrete,T,detail::TypeList::zeros,ENSEMBLE_LEARNING_COMPONENTS >* const>::value_type
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
	  //Attach to the child
	  child->SetParentFactor(this);
	  //Add this factor to each parent
	  boost::fusion::for_each(m_parent1_nodes.data(), AttachToChildFactor(this));
	  boost::fusion::for_each(m_parent2_nodes.data(), AttachToChildFactor(this));
	  //And the Weights too
	  Weights->AddChildFactor(this);
	};
      
      
	Moments<T>
	InitialiseMoments() const
	{
	  const size_t size = boost::mpl::size<typename parent2_t::type>::value;
	  std::vector< moments_t > moments1, moments2;
	  //reserve the space so that push_back is not expensive
	  moments1.reserve(size);
	  moments2.reserve(size);

	  //add the moments to the vector (pre-reserved) (moments1)
	  boost::fusion::for_each( m_parent1_nodes.data(), 
				   GetMomentsInit(moments1) );
	  //add the moments to the vector (pre-reserved) (moments2)
	  boost::fusion::for_each( m_parent2_nodes.data(), 
				   GetMomentsInit(moments2) );
	  
	  m_weights_node->InitialiseMoments();
	  return Model::CalcSample(moments1,
				   moments2, 
				   m_weights_node -> GetMoments()
				   );
	}

	//Request From the child Node
	NaturalParameters<T>
	GetNaturalNot( c_parameter v) const
	{
	  NaturalParameters<T> NP2Child(std::vector<T>(2));
	  m_LogNorm = 0;
	  
	  const size_t size = boost::mpl::size<typename parent2_t::type>::value;
	  std::vector<const moments_t* > moments1, moments2;
	  //reserve the space so that push_back is not expensive
	  moments1.reserve(size);
	  moments2.reserve(size);

	  //add the moments to the vector (pre-reserved) (moments1)
	  boost::fusion::for_each( m_parent1_nodes.data(), 
				   GetMoments(moments1) );
	  //add the moments to the vector (pre-reserved) (moments2)
	  boost::fusion::for_each( m_parent2_nodes.data(), 
				   GetMoments(moments2) );

	  const weights_moments_t* weights = m_weights_node->GetMoments();
	  for(size_t i=0;i<size;++i){
	    const moments_t* m1i = moments1[i];
	    const moments_t* m2i = moments2[i];
	    const T wi  = weights->operator[](i);
	    
	    NP2Child
	      += Model::CalcNP2Data(m1i, m2i) * wi;
	    
	    m_LogNorm 
	      += Model::CalcLogNorm(m1i,m2i) * wi;
	  }
	  return NP2Child;
	}

	//Request From the weights Node
	NaturalParameters<T,ENSEMBLE_LEARNING_COMPONENTS>
	GetNaturalNot( discrete_parameter v) const
	{
	  
	  const size_t size = boost::mpl::size<typename parent2_t::type>::value;
	  std::vector<const moments_t* > moments1, moments2;
	  //reserve the space so that push_back is not expensive
	  moments1.reserve(size);
	  moments2.reserve(size);
	  NaturalParameters<T,ENSEMBLE_LEARNING_COMPONENTS> NP2Weights;

	  //add the moments to the vector (pre-reserved) (moments1)
	  boost::fusion::for_each( m_parent1_nodes.data(), 
				   GetMoments(moments1) );
	  //add the moments to the vector (pre-reserved) (moments2)
	  boost::fusion::for_each( m_parent2_nodes.data(), 
				   GetMoments(moments2) );

	  const moments_t* child = m_child_node->GetMoments();
	  for(size_t i=0;i<size;++i){
	    NP2Weights[i] = Model::CalcAvLog(moments1[i], moments2[i],child);
	  }
	  return NP2Weights;
	}

	/* Define GetNaturalNot for each of the parent nodes.
	 * (Iterate with the preprocessor)
	 * This is done with an ugly preprocessor hack since these functions
	 * are overloads the base FactorNode class.
	 * They therefore must really be "written down" rather than just templated,
	 * which is where the preprocessor comes in.
	 */
#       include <boost/preprocessor/iteration/iterate.hpp>
#       define BOOST_PP_ITERATION_LIMITS (0, ENSEMBLE_LEARNING_COMPONENTS  - 1)
#       define BOOST_PP_FILENAME_1       "EnsembleLearning/node/factor/Mixture_GetNaturalNot.hpp"
#       include BOOST_PP_ITERATE()
	/* End of pre-processor hack */

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
      
	struct GetMoments
	{
	  GetMoments(std::vector<const moments_t*>& m) : m_moments(m) {}
	  template <class Node>
	  void operator()(Node& node) const 
	  {
	    m_moments.push_back(node->GetMoments());
	  }
	private:
	  std::vector<const moments_t*>& m_moments;
	};
	
	//store as type rather than as pointer.
	struct GetMomentsInit
	{
	  GetMomentsInit(std::vector< moments_t>& m) : m_moments(m) {}
	  template <class Node>
	  void operator()(Node& node) const 
	  {
	    m_moments.push_back(*(node->GetMoments()));
	  }
	private:
	  std::vector<moments_t>& m_moments;
	};

	struct AttachToChildFactor
	{
	  AttachToChildFactor(Mixture* p) : m_ptr(p) {}
	  template <class Node>
	  void operator()(Node& node) const 
	  {
	    node->AddChildFactor(m_ptr);
	  }
	  Mixture* m_ptr;
	};

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
  
