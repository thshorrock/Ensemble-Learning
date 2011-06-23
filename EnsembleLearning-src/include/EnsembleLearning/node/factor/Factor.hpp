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

#include "EnsembleLearning/message/Moments.hpp"
#include "EnsembleLearning/node/Node.hpp"
#include "EnsembleLearning/exponential_model/Discrete.hpp"
#include "EnsembleLearning/exponential_model/Dirichlet.hpp"

#include <boost/shared_ptr.hpp>
#include <boost/none.hpp>
#include <boost/assert.hpp> 
#include <boost/mpl/assert.hpp> 


#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/int.hpp>
#include <boost/mpl/bool.hpp>

#include <vector>


namespace ICR{
  namespace EnsembleLearning{
   
    //forward declaration
    
    template<class T>
    class NaturalParameters;
   

    /*  Need a NoSecondParent type for Discrete and Dirichlet Models that do not have 
     *  a second parent.  Calls to this parent won't be compiled and so 
     *  this struct can be empty.
     */
    struct
    NoSecondParent{};

    
    namespace detail{

      
      
    /******************************************************************************
     * Default Specialisation - Gaussian, RectifiedGaussian, Gamma
     ******************************************************************************/
    /** A Factor that joins various VariableNode.
     * @tparam Model The model used by the Factor.
     *  This general templated class is suitable for Gaussian, RectifiedGaussian and Gamma models.  The Dirichlet and Discrete models have partial specialisations.
     * @tparam T The data type used - either float or double.
     */
      template<template<class> class Model, class T,
	       class parent1_t, class parent2_t, class child_t
	       >
    class Factor : 
	public FactorNode<T,parent1_t>,
	public FactorNode<T,parent2_t>,
	public FactorNode<T,child_t>,
	public FactorNode_basic
    {
      struct enabler {};
      //Non-copiable
      //Factor(const Factor<Model,T>& f) {};
    public:
      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{

      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;

      typedef typename boost::call_traits< parent1_t* const>::param_type
      p1_parameter;
      typedef typename boost::call_traits< parent1_t* const>::value_type
      p1_t;
      typedef typename boost::call_traits< parent2_t* const>::param_type
      p2_parameter;
      typedef typename boost::call_traits< parent2_t* const>::value_type
      p2_t;
      typedef typename boost::call_traits< child_t* const>::param_type
      c_parameter;
      typedef typename boost::call_traits< child_t* const>::value_type
      c_t;

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
      Factor( p1_parameter Parent1,  
	      p2_parameter Parent2,  
	      c_parameter Child)
	: m_parent1_node(Parent1),
	  m_parent2_node(Parent2),
	  m_child_node(Child),
	  m_LogNorm(0)
      {
    	Child->SetParentFactor(this);
    	Parent1->AddChildFactor(this);
    	Parent2->AddChildFactor(this);
      };
       
      //For Discrete and Dirichlet Models that have 1 parent.
      Factor(parent1_t* Prior,  child_t* Child)
	: m_parent1_node(Prior),
	  m_parent2_node(0),
	  m_child_node(Child)
      {
    	m_parent1_node->AddChildFactor(this);
    	Child->SetParentFactor(this);

      };
      
      /** Initialise the Moments for the child moment.
       *  This function should be called before the iterations start so that a random 
       *  Moments can be calculated for the child (based on the parent priors).
       *  @return The Child Moments.
       */
      // moments_t
      // InitialiseMoments() const
      // {
      // 	//Initialise up the tree first
      // 	//m_parent1_node->InitialiseMoments();
      // 	//m_parent2_node->InitialiseMoments();
      // 	return  Model<T>::CalcSample(m_parent1_node->GetMoments(),
      // 				     m_parent2_node->GetMoments());
      // }
      
      template<class NoSecondParent_t>
      Moments<T>  
      InitialiseMoments(typename boost::disable_if<boost::is_same<parent2_t, NoSecondParent_t>, enabler >::type = enabler()) const
      {
	return Model<T>::CalcSample(m_parent1_node->GetMoments(),
				    m_parent2_node->GetMoments());
      }
      
      template<class NoSecondParent_t>
      Moments<T>  
      InitialiseMoments(typename boost::enable_if<boost::is_same<parent2_t, NoSecondParent_t>, enabler >::type = enabler()) const
      {
	return Model<T>::CalcSample(m_parent1_node->GetMoments());
      }

      Moments<T>  //The return type
      InitialiseMoments() const
      {
	return InitialiseMoments<NoSecondParent>();
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

      //Request From the Parent1 Node
      //parent2_t != NoSecondParent
      template<class NoSecondParent_t>
      NaturalParameters<T>
      GetNaturalNot(parent1_t* const,
		    typename boost::disable_if<boost::is_same<parent2_t, NoSecondParent_t>, enabler >::type = enabler()) const
      {
	Moments<T> parent2 = m_parent2_node->GetMoments();
	Moments<T> child = m_child_node->GetMoments();
	return Model<T>::CalcNP2Parent1(parent2,child);
      }

      //Request From the Parent1 Node
      template<class  NoSecondParent_t>
      NaturalParameters<T>
      GetNaturalNot(parent1_t* const,
		    typename boost::enable_if<boost::is_same<parent2_t, NoSecondParent_t>, enabler >::type = enabler()
		    ) const
      {
	Moments<T> child = m_child_node->GetMoments();
	return Model<T>::CalcNP2Prior(child);
      }

      //Request From the Parent1 Node
      NaturalParameters<T>
      GetNaturalNot( parent1_t* const v) const
      {
	return GetNaturalNot<NoSecondParent>(v);
      }

      /** Request From the Parent2 Node
       */
      template<class NoSecondParent_t>
      NaturalParameters<T>
      GetNaturalNot( parent2_t* const v,
		    typename boost::disable_if<boost::is_same<parent2_t, NoSecondParent_t>, enabler >::type = enabler()
		    ) const
      {
	Moments<T> parent1 = m_parent1_node->GetMoments();
	Moments<T> child = m_child_node->GetMoments();
	return Model<T>::CalcNP2Parent2(parent1,child);
      }

      /** Request From the Parent2 Node
       */
      template<class NoSecondParent_t>
      NaturalParameters<T>
      GetNaturalNot( parent2_t* const v,
		    typename boost::enable_if<boost::is_same<parent2_t, NoSecondParent_t>, enabler >::type = enabler()
		    ) const
      {
	BOOST_ASSERT(1==2);// This should never be called, but is required for interface
	return NaturalParameters<T>();
      }
      
      NaturalParameters<T>
      GetNaturalNot( parent2_t* const  v) const
      {
	return GetNaturalNot<NoSecondParent>(v);
      }
      
      //Request From the Child Node
      //parent2_t != NoSecondParent
      template<class NoSecondParent_t>
      NaturalParameters<T>
      GetNaturalNot( typename boost::call_traits< child_t* const>::param_type const v,
		    typename boost::disable_if<boost::is_same<parent2_t, NoSecondParent_t>, enabler >::type = enabler()
		    ) const
      {
	Moments<T> parent1 = m_parent1_node->GetMoments();
	Moments<T> parent2 = m_parent2_node->GetMoments();
	m_LogNorm = Model<T>::CalcLogNorm(parent1,parent2);
	return Model<T>::CalcNP2Data(parent1,parent2);
      }

      //Request From the Child Node
      //parent2_t == NoSecondParent
      template<class NoSecondParent_t>
      NaturalParameters<T>
      GetNaturalNot( typename boost::call_traits< child_t* const>::param_type v,
		    typename boost::enable_if<boost::is_same<parent2_t, NoSecondParent_t>, enabler >::type = enabler()
		    ) const
      {
	Moments<T> parent1 = m_parent1_node->GetMoments();
	m_LogNorm = Model<T>::CalcLogNorm(parent1);
	return Model<T>::CalcNP2Data(parent1);
      }

      //Request From the Child Node
      NaturalParameters<T>
      GetNaturalNot( typename boost::call_traits< child_t* const>::param_type const v) const
      {
	return GetNaturalNot<NoSecondParent>(v);
      }

      /** Calculate the Natural Paramter to one of the connected VariableNodes.
       *  @param v The VaiableNode where the NaturalParameter is sent.
       *  @return The NaturalParameter.
       */
      // template<class type>
      // NP_t
      // GetNaturalNot( typename boost::call_traits< type* const>::param_type v) const
      // {
      // 	if (boost::is_same<T,parent1_t>::value)
      // 	  {
      // 	    const moments_t parent2 = m_parent2_node->GetMoments();
      // 	    const moments_t child = m_child_node->GetMoments();
      // 	    return Model<T>::CalcNP2Parent1(parent2,child);
      // 	  }
      // 	else if (boost::is_same<T,parent2_t>::value)
      // 	  {
      // 	    const moments_t parent1 = m_parent1_node->GetMoments();
      // 	    const moments_t child = m_child_node->GetMoments();
      // 	    return Model<T>::CalcNP2Parent2(parent1,child);
      // 	  }
      // 	else 
      // 	  {
      // 	    BOOST_MPL_ASSERT((boost::is_same<T,child_t>));
      // 	    const moments_t parent1 = m_parent1_node->GetMoments();
      // 	    const moments_t parent2 = m_parent2_node->GetMoments();
      // 	    m_LogNorm = Model<T>::CalcLogNorm(parent1,parent2);
      // 	    return Model<T>::CalcNP2Data(parent1,parent2);
      // 	  }
      // }

    private: 
      p1_t m_parent1_node;
      p2_t m_parent2_node;
      c_t  m_child_node;
      mutable data_t m_LogNorm;  //no need for mutex to protect this variable.
    };
   

    
    // /******************************************************************************
    //  * DirichletModel Specialisation 
    //  ******************************************************************************/
    // /** A partial specialisation to a Factor node using the Dirichlet Model.
    //  *  @tparam T The data type, double or float
    //  */
    // template<class T>
    // class Factor<Dirichlet,T > : public FactorNode<T>
    // {
    //   Factor(const Factor<Dirichlet, T>& f) {};
    // public:
      
    //   /** @name Useful typdefs for types that are exposed to the user.
    //    */
    //   ///@{
    //   typedef typename boost::call_traits< VariableNode<T>* const>::param_type
    //   variable_parameter;
    //   typedef typename boost::call_traits< VariableNode<T>* const>::value_type
    //   variable_t;

    //   typedef typename boost::call_traits< NaturalParameters<T> >::param_type
    //   NP_parameter;
    //   typedef typename boost::call_traits< NaturalParameters<T> >::value_type
    //   NP_t;

    //   typedef typename boost::call_traits<T>::value_type
    //   data_t;
    //   typedef typename boost::call_traits<T>::param_type
    //   data_parameter;
      
    //   typedef typename boost::call_traits< Moments<T> >::value_type
    //   moments_t;

    //   ///@}

    //   /** Constructor.
    //    *  @param Prior The Parent VariableNode.
    //    *  @param Child The Child VariableNode.
    //    */
    //   Factor( variable_parameter Prior,  
    // 	      variable_parameter Child)
    // 	: m_prior_node(Prior),
    // 	  m_child_node(Child),
    // 	  m_LogNorm(0)
    //   {
    // 	Child->SetParentFactor(this);
    // 	Prior->AddChildFactor(this);

    //   };
      
    //   /** Initialise the Moments for the child moment.
    //    *  This function should be called becore the iterations start so that a random 
    //    *  Moments can be calculated for the child (based on the parent priors).
    //    *  @return The Child Moments.
    //    */
    //   moments_t
    //   InitialiseMoments() const
    //   {
    // 	//Initialise up the tree first
    // 	m_prior_node->InitialiseMoments();
    // 	return Dirichlet<T>::CalcSample(m_prior_node->GetMoments());
    //   }


    //   /** Calculate the Log Normalisation.
    //    *  This is called only from the child node.
    //    *  This is used to evaluate the cost.
    //    *  @return The log normalisation.
    //    */
    //   //Only called from the child node - can be no collision, no need for mutex.
    //   T
    //   CalcLogNorm() const 
    //   {
    // 	return m_LogNorm;
    //   }


    //   /** Calculate the Natural Paramter to one of the connected VariableNodes.
    //    *  @param v The VaiableNode where the NaturalParameter is sent.
    //    *  @return The NaturalParameter.
    //    */
    //   NP_t
    //   GetNaturalNot( variable_parameter v) const
    //   {
    // 	BOOST_ASSERT(v==m_child_node);
    // 	moments_t prior = m_prior_node->GetMoments();
    // 	m_LogNorm = Dirichlet<T>::CalcLogNorm(prior);

    // 	return Dirichlet<T>::CalcNP2Data(prior);
    //   }
      
    // private: 
    //   variable_t m_prior_node,  m_child_node;
      
    //   mutable T m_LogNorm; //no need for mutex to protect this variable.
    // };

    
    // /******************************************************************************
    //  * Discrete Specialisation
    //  ******************************************************************************/
    // /** A partial specialisation to a Factor node using the Discrete model.
    //  *  @tparam T The data type, double or float
    //  */
    // template<class T>
    // class Factor<Discrete, T> : public FactorNode<T>
    // {
    // public:
      
    //   /** @name Useful typdefs for types that are exposed to the user.
    //    */
    //   ///@{

    //   typedef typename boost::call_traits< VariableNode<T>* const>::param_type
    //   variable_parameter;
    //   typedef typename boost::call_traits< VariableNode<T>* const>::value_type
    //   variable_t;

    //   typedef typename boost::call_traits< NaturalParameters<T> >::param_type
    //   NP_parameter;
    //   typedef typename boost::call_traits< NaturalParameters<T> >::value_type
    //   NP_t;

    //   typedef typename boost::call_traits<T>::value_type
    //   data_t;
    //   typedef typename boost::call_traits<T>::param_type
    //   data_parameter;
  
    //   typedef typename boost::call_traits< Moments<T> >::param_type
    //   moments_parameter;
    //   typedef typename boost::call_traits< Moments<T> >::const_reference
    //   moments_const_reference;
    //   typedef typename boost::call_traits< Moments<T> >::value_type
    //   moments_t;
          
    //   ///@}

    //   /** Constructor.
    //    *  @param Prior The first VariableNode attached to the Factor:
    //    *     This is a Hidden Node for the Dirichlet distribution.
    //    *  @param Child The child parameter of the factor.
    //    */
    //   Factor( variable_parameter Prior,  variable_parameter Child)
    // 	: m_prior_node(Prior),
    // 	  m_child_node(Child),
    // 	  m_LogNorm(0)
    //   {
    // 	Child->SetParentFactor(this);
    // 	Prior->AddChildFactor(this);
    //   };
      
    //   /** Initialise the Moments for the child moment.
    //    *  This function should be called before the iterations start so that a random 
    //    *  Moments can be calculated for the child (based on the parent priors).
    //    *  @return The Child Moments.
    //    */
    //   Moments<T>
    //   InitialiseMoments() const
    //   {
    // 	return Discrete<T>::CalcSample(m_prior_node->GetMoments());
    //   }



    //   /** Calculate the Log Normalisation.
    //    *  This is called only from the child node.
    //    *  This is used to evaluate the cost.
    //    *  @return The log normalisation.
    //    */
    //   //Only called from the child node - can be no collision, no need for mutex.
    //   T
    //   CalcLogNorm() const 
    //   {
    // 	return m_LogNorm;
    //   }


    //   /** Calculate the Natural Paramter to one of the connected VariableNodes.
    //    *  @param v The VaiableNode where the NaturalParameter is sent.
    //    *  @return The NaturalParameter.
    //    */
    //   NP_t
    //   GetNaturalNot( variable_parameter v) const
    //   {
    // 	if (v==m_prior_node)
    // 	  {
    // 	    const moments_t child  = m_child_node->GetMoments();
    // 	    return Discrete<T>::CalcNP2Prior(child);/// = child;
    // 	  }
    // 	else 
    // 	  {
    // 	    BOOST_ASSERT(v==m_child_node);
    // 	    const moments_t prior  = m_prior_node->GetMoments(); 
    // 	    m_LogNorm = Discrete<T>::CalcLogNorm(prior);
    // 	    return  Discrete<T>::CalcNP2Data(prior);/// = child;
    // 	  }
    //   }
      
    // private: 
    //   variable_t m_prior_node, m_child_node;
    //   mutable T  m_LogNorm; //no need for mutex to protect this variable.
    // };
    
    }
  }
}

#endif  // guard for FACTOR_HPP
