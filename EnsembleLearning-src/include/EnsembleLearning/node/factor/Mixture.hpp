#pragma once
#ifndef FACTOR_MIXTURE_HPP
#define FACTOR_MIXTURE_HPP


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
#include "EnsembleLearning/node/variable/Hidden.hpp"
#include "EnsembleLearning/message/NaturalParameters.hpp"
#include "EnsembleLearning/message/Moments.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/none.hpp>
#include <vector>

namespace ICR{
  namespace EnsembleLearning{

    //forward declarations
    
    template<class T>
    class Discrete;
    
    
    /******************************************************************************
     * Default Specialisation
     ******************************************************************************/
    /** Factor Linking Mixture nodes */
    template<class Model, class T>
    class Mixture : public FactorNode<T>
    {
      //Non-copieable
      Mixture(const Mixture<Model, T>& f) {};
    public:
      
      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{

      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;
      typedef typename boost::call_traits< std::vector<VariableNode<T>*> >::param_type
      variable_vector_parameter;
      typedef typename boost::call_traits< std::vector<VariableNode<T>*> >::value_type
      variable_vector_t;

      typedef typename boost::call_traits<HiddenNode<Discrete<T>,T >*>::value_type
      discrete_t;
      typedef typename boost::call_traits<HiddenNode<Discrete<T>,T >*>::value_type
      discrete_parameter;

      ///@}
      
      /** Construct a mixture node.
       * @param Parent1 A vector of parent Nodes (for example Means)
       * @param Parent2 A vector of parent Nodes (for example Precisions)
       * @param Weights The Weights on each of the vectors 
       *   (The size of the weights need to be the same as sthe size of the vectors)
       * @param child The child node.
       */
      Mixture(variable_vector_parameter Parent1, 
	      variable_vector_parameter Parent2,  
	      discrete_parameter Weights,
	      variable_parameter child
	      
	       )
	:  m_parent1_nodes(Parent1), 
	   m_parent2_nodes(Parent2),
	   m_weights_node(Weights),
	   m_child_node(child),
	   m_LogNorm(0)
      {
	//Need to be as many parents to both.
	BOOST_ASSERT(Parent1.size() == Parent2.size());
	
	for(size_t i=0;i<Parent1.size();++i){
	  m_parent1_nodes[i]->AddChildFactor(this);
	  m_parent2_nodes[i]->AddChildFactor(this);
	}
	Weights->AddChildFactor(this);
    	child->SetParentFactor(this);
      };
      
      
      Moments<T>
      InitialiseMoments() const
      {
	//Initialise up the tree first
	//std::cout<<"p1size = "<<m_parent1_nodes.size()<<std::endl;

	for(size_t i=0;i<m_parent1_nodes.size();++i){
	  m_parent1_nodes[i]->InitialiseMoments();
	  m_parent2_nodes[i]->InitialiseMoments();
	}

	// std::for_each(m_parent1_nodes.begin(), m_parent1_nodes.end(),
	// 		 boost::bind(&VariableNode<T>::InitialiseMoments, _1)
	// 		 );
	// std::for_each(m_parent2_nodes.begin(), m_parent2_nodes.end(),
	// 		 boost::bind(&VariableNode<T>::InitialiseMoments, _1)
	// 		 );
	m_weights_node->InitialiseMoments();
	return Model::CalcSample(m_parent1_nodes,m_parent2_nodes, m_weights_node );
      }

      T
      CalcLogNorm() const 
      {
	return m_LogNorm;
      }


      NaturalParameters<T>
      GetNaturalNot( variable_parameter v) const;
      
    private: 

      variable_vector_t m_parent1_nodes, m_parent2_nodes;
      discrete_t m_weights_node;
      variable_t  m_child_node;
      
      mutable T m_LogNorm;
    };
    
  
    /**************************************************************************************
     **************************************************************************************
     **************************************************************************************
     *
     *  Default
     *
     **************************************************************************************
     **************************************************************************************
     **************************************************************************************/
    template<class Model, class T>
    inline
    NaturalParameters<T>
    Mixture< Model ,T>::GetNaturalNot( variable_parameter v) const
    {
      if (v == m_child_node) 
	{

	  NaturalParameters<T> NP2Child(std::vector<T>(2));
	  m_LogNorm = 0;
	  
	  const Moments<T>& weights = m_weights_node->GetMoments();
	  for(size_t i=0;i<m_parent1_nodes.size();++i){
	    const Moments<T>& parent1 = m_parent1_nodes[i]->GetMoments();
	    const Moments<T>& parent2 = m_parent2_nodes[i]->GetMoments();
	    NP2Child
	      += Model::CalcNP2Data( parent1, parent2) * weights[i];
	    
	    m_LogNorm 
	      += Model::CalcLogNorm(parent1,parent2) * weights[i];
	  }
	  return NP2Child;
	}
      else if (v == m_weights_node) 
	{
	  NaturalParameters<T> NP2Weights(m_parent1_nodes.size());
	  

	  const Moments<T>& child = m_child_node->GetMoments();
	  for(size_t i=0;i<m_parent1_nodes.size();++i){
	    const Moments<T>& parent1 = m_parent1_nodes[i]->GetMoments();
	    const Moments<T>& parent2 = m_parent2_nodes[i]->GetMoments();
	    NP2Weights[i] = Model::CalcAvLog(parent1, parent2,child);
	  }
	  return NP2Weights;
	}
      else 
	{
	  typename std::vector<VariableNode<T> *>::const_iterator it;
	  it = PARALLEL_FIND(m_parent1_nodes.begin(), m_parent1_nodes.end(), v);
	  if (it!=m_parent1_nodes.end()) 
	    {
	      size_t i = it-m_parent1_nodes.begin();
	      const Moments<T>& weights = m_weights_node->GetMoments();
	      const Moments<T>& parent2 = m_parent2_nodes[i]->GetMoments();
	      const Moments<T>& child = m_child_node->GetMoments();
	      return Model::CalcNP2Parent1(parent2,child)  * weights[i];
	    }
	  else
	    {
	      it = PARALLEL_FIND(m_parent2_nodes.begin(), m_parent2_nodes.end(), v);
	      size_t i = it-m_parent2_nodes.begin();
	      const Moments<T>& weights = m_weights_node->GetMoments();
	      const Moments<T>& parent1 = m_parent1_nodes[i]->GetMoments();
	      const Moments<T>& child = m_child_node->GetMoments();
	      return Model::CalcNP2Parent2(parent1,child) * weights[i];
	    }
	}
      throw ("Unknown Node in GetNaturalNot");
    }

  
  }
}

#endif  // guard for FACTOR_MIXTURE_HPP
