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
      
      Factor( VariableNode<T>* Parent1,  VariableNode<T>* Parent2,  VariableNode<T>* Child)
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
      GetNaturalNot(const VariableNode<T>* v) const
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
    class Factor<DirichletModel<double> > : public FactorNode<double>
    {
      Factor(const Factor<DirichletModel<double> >& f) {};
    public:
      

      Factor( VariableNode<double>* Prior,  VariableNode<double>* Child)
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
	return DirichletModel<double>::CalcSample(m_prior_node);
      }


      double
      CalcLogNorm() const 
      {
	return m_LogNorm;
      }


      NaturalParameters<double>
      GetNaturalNot(const VariableNode<double>* v) const;
      

      // void 
      // Iterate();
    private: 
       VariableNode<double> *m_prior_node,  *m_child_node;
      
      mutable double m_LogNorm;
    private:
      mutable boost::mutex m_mutex;
    };
    /******************************************************************************
     * DiscreteModel Specialisation
     ******************************************************************************/
    template<>
    class Factor<DiscreteModel<double> > : public FactorNode<double>
    {
    public:
      

      Factor( VariableNode<double>* Prior,  VariableNode<double>* Child)
	: m_prior_node(Prior),
	  m_child_node(Child)
      {
    	Prior->AddChildFactor(this);
    	Child->SetParentFactor(this);

      };
      
      Moments<double>
      InitialiseMoments() const
      {
	return DiscreteModel<double>::CalcSample(m_prior_node);
      }




      double
      CalcLogNorm() const 
      {
	return m_LogNorm;
      }


      NaturalParameters<double>
      GetNaturalNot(const VariableNode<double>* v) const;
      
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
    Factor< DirichletModel<double> ,double>::GetNaturalNot(const VariableNode<double>* v) const
    {
      if (v==m_child_node)
	{
	  Moments<double> prior = m_prior_node->GetMoments();
	  m_LogNorm = DirichletModel<double>::CalcLogNorm(prior);

	  return DirichletModel<double>::CalcNP2Data(prior);

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
    Factor<DiscreteModel<double>,double>::GetNaturalNot(const VariableNode<double>* v) const
    {
      boost::mutex::scoped_lock lock(m_mutex);
      if (v==m_prior_node)
       	{
      Moments<double> child  = m_child_node->GetMoments();
       return DiscreteModel<double>::CalcNP2Prior(child);/// = child;

      	}
      else 
       	{
	  Moments<double> prior  = m_prior_node->GetMoments(); 
	  m_LogNorm = DiscreteModel<double>::CalcLogNorm(prior);

	  return  DiscreteModel<double>::CalcNP2Data(prior);/// = child;
	  // return m_NP2Child;
	}
    }

  }
}
