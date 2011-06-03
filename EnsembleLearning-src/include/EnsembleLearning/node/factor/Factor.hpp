#pragma once
#ifndef FACTOR_HPP
#define FACTOR_HPP


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
#include "EnsembleLearning/message/NaturalParameters.hpp"
#include "EnsembleLearning/message/Moments.hpp"
// #include "EnsembleLearning/Message.hpp"

//#include "EnsembleLearning/piping/Piping.hpp"
#include "EnsembleLearning/node/variable/Hidden.hpp"
#include "EnsembleLearning/exponential_model/Gaussian.hpp"
#include "EnsembleLearning/exponential_model/RectifiedGaussian.hpp"
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
    
   
    /******************************************************************************
     * Default Specialisation - Gaussian, RectifiedGaussian, Gamma
     ******************************************************************************/
    /** A Factor that joins various VariableNode.
     * @tparam Model The model used by the Factor.
     *  This general templated class is suitable for Gaussian, RectifiedGaussian and Gamma models.  The Dirichlet and Discrete models have partial specialisations.
     * @tparam T The data type used - either float or double.
     */
    template<class Model, class T>
    class Factor : public FactorNode<T>
    {
      //Non-copiable
      Factor(const Factor<Model,T>& f) {};
    public:
      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{

      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;

      typedef typename boost::call_traits< NaturalParameters<T> >::param_type
      NP_parameter;
      typedef typename boost::call_traits< NaturalParameters<T> >::value_type
      NP_t;
      
      typedef typename boost::call_traits<T>::value_type
      data_t;
      typedef typename boost::call_traits<T>::param_type
      data_parameter;
      

      typedef typename boost::call_traits< Moments<T> >::value_type
      moments_t;

      ///@}
      /** Constructor.
       *  @param Parent1 The first VariableNode attached to the Factor:
       *     This is the mean for the Gaussian/RectifiedGaussian or the Shape for the Gamma distribution.
       *  @param Parent2 The second VariableNode attached to the Factor:
       *     This is the precision for the Gaussian/RectifiedGaussian or the inverse scale for the Gamma distribution.
       *  @param Child The child parameter of the factor.
       */
      Factor( variable_parameter Parent1,  
	      variable_parameter Parent2,  
	      variable_parameter Child)
	: m_parent1_node(Parent1),
	  m_parent2_node(Parent2),
	  m_child_node(Child)
      {
    	Parent1->AddChildFactor(this);
    	Parent2->AddChildFactor(this);
    	Child->SetParentFactor(this);
      };
      
      /** Initialise the Moments for the child moment.
       *  This function should be called before the iterations start so that a random 
       *  Moments can be calculated for the child (based on the parent priors).
       *  @return The Child Moments.
       */
      moments_t
      InitialiseMoments() const
      {
	//Initialise up the tree first
	m_parent1_node->InitialiseMoments();
	m_parent2_node->InitialiseMoments();
	return  Model::CalcSample(m_parent1_node,m_parent2_node );
      }
      
      /** Calculate the Log Normalisation.
       *  This is called only from the child node.
       *  This is used to evaluate the cost.
       *  @return The log normalisation.
       */
      //Only called from the child node - can be no collision, no need for mutex.
      data_t
      CalcLogNorm() const 
      {
	return m_LogNorm;
      }

      /** Calculate the Natural Paramter to one of the connected VariableNodes.
       *  @param v The VaiableNode where the NaturalParameter is sent.
       *  @return The NaturalParameter.
       */
      NP_t
      GetNaturalNot(variable_parameter v) const
      {
	if (v==m_parent1_node)
	  {
	    const moments_t parent2 = m_parent2_node->GetMoments();
	    const moments_t child = m_child_node->GetMoments();
	    return Model::CalcNP2Parent1(parent2,child);
	  }
	else if (v==m_parent2_node)
	  {
	    const moments_t parent1 = m_parent1_node->GetMoments();
	    const moments_t child = m_child_node->GetMoments();
	    return Model::CalcNP2Parent2(parent1,child);
	  }
	else 
	  {
	    BOOST_ASSERT(v == m_child_node);
	    const moments_t parent1 = m_parent1_node->GetMoments();
	    const moments_t parent2 = m_parent2_node->GetMoments();
	    m_LogNorm = Model::CalcLogNorm(parent1,parent2);
	    return Model::CalcNP2Data(parent1,parent2);
	  }
      }

    private: 
      variable_t m_parent1_node, m_parent2_node, m_child_node;
      mutable data_t m_LogNorm;  //no need for mutex to protect this variable.
    };
   

    
    /******************************************************************************
     * DirichletModel Specialisation 
     ******************************************************************************/
    /** A partial specialisation to a Factor node using the Dirichlet Model.
     *  @tparam T The data type, double or float
     */
    template<class T>
    class Factor<Dirichlet<T>,T > : public FactorNode<T>
    {
      Factor(const Factor<Dirichlet<T>, T>& f) {};
    public:
      
      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{
      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;

      typedef typename boost::call_traits< NaturalParameters<T> >::param_type
      NP_parameter;
      typedef typename boost::call_traits< NaturalParameters<T> >::value_type
      NP_t;

      typedef typename boost::call_traits<T>::value_type
      data_t;
      typedef typename boost::call_traits<T>::param_type
      data_parameter;
      
      typedef typename boost::call_traits< Moments<T> >::value_type
      moments_t;

      ///@}

      /** Constructor.
       *  @param Prior The Parent VariableNode.
       *  @param Child The Child VariableNode.
       */
      Factor( variable_parameter Prior,  
	      variable_parameter Child)
	: m_prior_node(Prior),
	  m_child_node(Child)
      {
    	Prior->AddChildFactor(this);
    	Child->SetParentFactor(this);

      };
      
      /** Initialise the Moments for the child moment.
       *  This function should be called becore the iterations start so that a random 
       *  Moments can be calculated for the child (based on the parent priors).
       *  @return The Child Moments.
       */
      moments_t
      InitialiseMoments() const
      {
	//Initialise up the tree first
	m_prior_node->InitialiseMoments();
	return Dirichlet<T>::CalcSample(m_prior_node);
      }


      /** Calculate the Log Normalisation.
       *  This is called only from the child node.
       *  This is used to evaluate the cost.
       *  @return The log normalisation.
       */
      //Only called from the child node - can be no collision, no need for mutex.
      T
      CalcLogNorm() const 
      {
	return m_LogNorm;
      }


      /** Calculate the Natural Paramter to one of the connected VariableNodes.
       *  @param v The VaiableNode where the NaturalParameter is sent.
       *  @return The NaturalParameter.
       */
      NP_t
      GetNaturalNot( variable_parameter v) const
      {
	if (v==m_child_node)
	  {
	    moments_t prior = m_prior_node->GetMoments();
	    m_LogNorm = Dirichlet<T>::CalcLogNorm(prior);

	    return Dirichlet<T>::CalcNP2Data(prior);
	  }
       	throw ("Unknown Node in GetNaturalNot");
      }
      
    private: 
      variable_t m_prior_node,  m_child_node;
      
      mutable T m_LogNorm; //no need for mutex to protect this variable.
    };

    
    /******************************************************************************
     * Discrete Specialisation
     ******************************************************************************/
    /** A partial specialisation to a Factor node using the Discrete model.
     *  @tparam T The data type, double or float
     */
    template<class T>
    class Factor<Discrete<T>, T> : public FactorNode<T>
    {
    public:
      
      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{

      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;

      typedef typename boost::call_traits< NaturalParameters<T> >::param_type
      NP_parameter;
      typedef typename boost::call_traits< NaturalParameters<T> >::value_type
      NP_t;

      typedef typename boost::call_traits<T>::value_type
      data_t;
      typedef typename boost::call_traits<T>::param_type
      data_parameter;
  
      typedef typename boost::call_traits< Moments<T> >::param_type
      moments_parameter;
      typedef typename boost::call_traits< Moments<T> >::const_reference
      moments_const_reference;
      typedef typename boost::call_traits< Moments<T> >::value_type
      moments_t;
          
      ///@}

      /** Constructor.
       *  @param Prior The first VariableNode attached to the Factor:
       *     This is a Hidden Node for the Dirichlet distribution.
       *  @param Child The child parameter of the factor.
       */
      Factor( variable_parameter Prior,  variable_parameter Child)
	: m_prior_node(Prior),
	  m_child_node(Child)
      {
    	Prior->AddChildFactor(this);
    	Child->SetParentFactor(this);
      };
      
      /** Initialise the Moments for the child moment.
       *  This function should be called before the iterations start so that a random 
       *  Moments can be calculated for the child (based on the parent priors).
       *  @return The Child Moments.
       */
      Moments<T>
      InitialiseMoments() const
      {
	return Discrete<T>::CalcSample(m_prior_node);
      }



      /** Calculate the Log Normalisation.
       *  This is called only from the child node.
       *  This is used to evaluate the cost.
       *  @return The log normalisation.
       */
      //Only called from the child node - can be no collision, no need for mutex.
      T
      CalcLogNorm() const 
      {
	return m_LogNorm;
      }


      /** Calculate the Natural Paramter to one of the connected VariableNodes.
       *  @param v The VaiableNode where the NaturalParameter is sent.
       *  @return The NaturalParameter.
       */
      NP_t
      GetNaturalNot( variable_parameter v) const
      {
	if (v==m_prior_node)
	  {
	    const moments_t child  = m_child_node->GetMoments();
	    return Discrete<T>::CalcNP2Prior(child);/// = child;
	  }
	else 
	  {
	    BOOST_ASSERT(v==m_child_node);
	    const moments_t prior  = m_prior_node->GetMoments(); 
	    m_LogNorm = Discrete<T>::CalcLogNorm(prior);
	    return  Discrete<T>::CalcNP2Data(prior);/// = child;
	  }
      }
      
    private: 
      variable_t m_prior_node, m_child_node;
      mutable T  m_LogNorm; //no need for mutex to protect this variable.
    };
    
  
  }
}

#endif  // guard for FACTOR_HPP
