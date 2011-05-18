#pragma once

#include "ICA/node/Node.hpp"
#include "ICA/message/NaturalParameters.hpp"
#include "ICA/message/Moments.hpp"
// #include "ICA/Message.hpp"

//#include "ICA/piping/Piping.hpp"
#include "ICA/node/variable/Dirichlet.hpp"
#include "ICA/node/variable/Hidden.hpp"
#include "ICA/exponential_model/Gaussian.hpp"
#include "ICA/exponential_model/RectifiedGaussian.hpp"
#include "ICA/exponential_model/Gamma.hpp"
#include "ICA/exponential_model/Discrete.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/none.hpp>
#include <vector>
#include <boost/assert.hpp> 
#include "rng.hpp"


namespace ICR{
  namespace ICA{
    

    // template<class Model, class T = double>
    // class Factor :  public FactorNode<T>
    // {
    //   //this will not compile if you try to call it directly.
    //   // You need to use one of the specializations.
    // };
   
    /******************************************************************************
     * Default Specialisation Double - Gaussian, RectifiedGaussian, Gamma
     ******************************************************************************/
    template<class Model, class T = double>
    class Factor : public FactorNode<T>
    {
      Factor(const Factor<Model>& f) {};
    public:
      typedef typename boost::call_traits< VariableNode<T>* const>::param_type
      variable_parameter;
      typedef typename boost::call_traits< VariableNode<T>* const>::value_type
      variable_t;

      
      Factor( variable_parameter Parent1,  variable_parameter Parent2,  variable_parameter Child)
	: m_parent1_node(Parent1),
	  m_parent2_node(Parent2),
	  m_child_node(Child)
      {
    	Parent1->AddChildFactor(this);
    	Parent2->AddChildFactor(this);
    	Child->SetParentFactor(this);
      };
      
      Moments<T>
      InitialiseMoments() const
      {
	//Initialise up the tree first
	m_parent1_node->InitialiseMoments();
	m_parent2_node->InitialiseMoments();
	return  Model::CalcSample(m_parent1_node,m_parent2_node );
      }
      
      T
      CalcLogNorm() const 
      {
	return m_LogNorm;
      }


      NaturalParameters<T>
      GetNaturalNot( variable_parameter const v) const
      {

	if (v==m_parent1_node)
	  {
	    Moments<T> parent2 = m_parent2_node->GetMoments();
	    Moments<T> child = m_child_node->GetMoments();
	    return Model::CalcNP2Parent1(parent2,child);
	  }
	else if (v==m_parent2_node)
	  {
	    Moments<T> parent1 = m_parent1_node->GetMoments();
	    Moments<T> child = m_child_node->GetMoments();
	    return Model::CalcNP2Parent2(parent1,child);
	  }
	else if (v == m_child_node) 
	  {
	    Moments<T> parent1 = m_parent1_node->GetMoments();
	    Moments<T> parent2 = m_parent2_node->GetMoments();
	    m_LogNorm = Model::CalcLogNorm(parent1,parent2);
	    return Model::CalcNP2Data(parent1,parent2);
	  }
	else{
	  throw ("Unknown Node in GetNaturalNot");
	}
      }

    private: 
      VariableNode<T> *m_parent1_node, *m_parent2_node, *m_child_node;
      mutable T m_LogNorm;
    private:
      mutable boost::mutex m_mutex;
    };
   
    /******************************************************************************
     * DirichletModel Specialisation Double
     ******************************************************************************/
    template<>
    class Factor<Dirichlet<double> > : public FactorNode<double>
    {
      Factor(const Factor<Dirichlet<double> >& f) {};
    public:
      
      typedef  boost::call_traits< VariableNode<double>* const>::param_type
      variable_parameter;
      typedef  boost::call_traits< VariableNode<double>* const>::value_type
      variable_t;


      Factor( variable_parameter Prior,  variable_parameter Child)
	: m_prior_node(Prior),
	  m_child_node(Child)
      {
    	Prior->AddChildFactor(this);
    	Child->SetParentFactor(this);

      };
      
      Moments<double>
      InitialiseMoments() const
      {
	//Initialise up the tree first
	m_prior_node->InitialiseMoments();
	return Dirichlet<double>::CalcSample(m_prior_node);
      }


      double
      CalcLogNorm() const 
      {
	return m_LogNorm;
      }


      NaturalParameters<double>
      GetNaturalNot( variable_parameter v) const;
      

      // void 
      // Iterate();
    private: 
       VariableNode<double> *m_prior_node,  *m_child_node;
      
      mutable double m_LogNorm;
    private:
      mutable boost::mutex m_mutex;
    };
    /******************************************************************************
     * Discrete Specialisation
     ******************************************************************************/
    template<>
    class Factor<Discrete<double> > : public FactorNode<double>
    {
    public:
      

      typedef  boost::call_traits< VariableNode<double>* const>::param_type
      variable_parameter;
      typedef  boost::call_traits< VariableNode<double>* const>::value_type
      variable_t;

      Factor( variable_parameter Prior,  variable_parameter Child)
	: m_prior_node(Prior),
	  m_child_node(Child)
      {
    	Prior->AddChildFactor(this);
    	Child->SetParentFactor(this);

      };
      
      Moments<double>
      InitialiseMoments() const
      {
	return Discrete<double>::CalcSample(m_prior_node);
      }




      double
      CalcLogNorm() const 
      {
	return m_LogNorm;
      }


      NaturalParameters<double>
      GetNaturalNot( variable_parameter v) const;
      
      // void
      // Reset();

      void 
      Iterate();
    private: 
      VariableNode<double> *m_prior_node, *m_child_node;
      mutable double  m_LogNorm;
    private:
      mutable boost::mutex m_mutex;
    };
    
  

   


    /**************************************************************************************
     **************************************************************************************
     **************************************************************************************
     *
     *  Dirichlet
     *
     **************************************************************************************
     **************************************************************************************
    //  **************************************************************************************/
    //template<>
    inline
    NaturalParameters<double>
    Factor< Dirichlet<double> ,double>::GetNaturalNot( variable_parameter v) const
    {
      if (v==m_child_node)
	{
	  Moments<double> prior = m_prior_node->GetMoments();
	  m_LogNorm = Dirichlet<double>::CalcLogNorm(prior);

	  return Dirichlet<double>::CalcNP2Data(prior);

	  //std::cout<<"returning child"<<m_NP2Child<<" from Dirichlet"<<std::endl;
	  //return m_NP2Child;
	}
       	throw ("Unknown Node in GetNaturalNot");
    }
    /**************************************************************************************
     **************************************************************************************
     **************************************************************************************
     *
     *  Discrete
     *
     **************************************************************************************
     **************************************************************************************
     **************************************************************************************/
  

    inline
    NaturalParameters<double>
    Factor<Discrete<double>,double>::GetNaturalNot( variable_parameter v) const
    {
      boost::mutex::scoped_lock lock(m_mutex);
      if (v==m_prior_node)
       	{
      Moments<double> child  = m_child_node->GetMoments();
       return Discrete<double>::CalcNP2Prior(child);/// = child;

      	}
      else 
       	{
	  Moments<double> prior  = m_prior_node->GetMoments(); 
	  m_LogNorm = Discrete<double>::CalcLogNorm(prior);

	  return  Discrete<double>::CalcNP2Data(prior);/// = child;
	  // return m_NP2Child;
	}
    }

  }
}
