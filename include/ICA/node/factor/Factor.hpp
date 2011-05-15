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
    

    template<class Model, class T = double>
    class Factor :  public FactorNode<T>
    {
      //this will not compile if you try to call it directly.
      // You need to use one of the specializations.
    };
   
    /******************************************************************************
     * GaussianModel Specialisation Double
     ******************************************************************************/
    template<>
    class Factor<GaussianModel<double> > : public FactorNode<double>
    {
      Factor(const Factor<GaussianModel<double> >& f) {};
    public:
      
      Factor( VariableNode<double>* Parent1,  VariableNode<double>* Parent2,  VariableNode<double>* Child)
	: m_parent1_node(Parent1),
	  m_parent2_node(Parent2),
	  m_child_node(Child)
      {
    	Parent1->AddChildFactor(this);
    	Parent2->AddChildFactor(this);
    	Child->SetParentFactor(this);

      };
      
      Moments<double>
      InitialiseMoments() const
      {
	return  GaussianModel<double>::CalcSample(m_parent1_node,m_parent2_node );
      }
      
      double
      CalcLogNorm() const 
      {
	return m_LogNorm;
      }


      NaturalParameters<double>
      GetNaturalNot(const VariableNode<double>* v) const;
      

    private: 
       VariableNode<double> *m_parent1_node, *m_parent2_node, *m_child_node;
      mutable double m_LogNorm;
    private:
      mutable boost::mutex m_mutex;
    };
     /******************************************************************************
     * RectifiedGaussianModel Specialisation Double
     ******************************************************************************/
    template<>
    class Factor<RectifiedGaussianModel<double> > : public FactorNode<double>
    {
      Factor(const Factor<RectifiedGaussianModel<double> >& f) {};
    public:
      
      Factor( VariableNode<double>* Parent1,  VariableNode<double>* Parent2,  VariableNode<double>* Child)
	: m_parent1_node(Parent1),
	  m_parent2_node(Parent2),
	  m_child_node(Child)
      {
    	Parent1->AddChildFactor(this);
    	Parent2->AddChildFactor(this);
    	Child->SetParentFactor(this);

	//const Moments<double> M =RectifiedGaussianModel<double>::CalcSample(Parent1, Parent2);
      };
      
      
      Moments<double>
      InitialiseMoments() const
      {
	return  RectifiedGaussianModel<double>::CalcSample(m_parent1_node,m_parent2_node );
      }
      double
      CalcLogNorm() const 
      {
	return m_LogNorm;
      }


      NaturalParameters<double>
      GetNaturalNot(const VariableNode<double>* v) const;
      

    private: 
       VariableNode<double> *m_parent1_node, *m_parent2_node, *m_child_node;
      
      mutable double m_LogNorm;
    private:
      mutable boost::mutex m_mutex;
    };
    
    /******************************************************************************
     * GammaModel Specialisation
     ******************************************************************************/
    template<>
    class Factor<GammaModel<double> > : public FactorNode<double>
    {
      
    public:
      
      Factor( VariableNode<double>* Parent1,  VariableNode<double>* Parent2,  VariableNode<double>* Child)
	: m_parent1_node(Parent1),
	  m_parent2_node(Parent2),
	  m_child_node(Child)
      {
    	Parent1->AddChildFactor(this);
    	Parent2->AddChildFactor(this);
    	Child->SetParentFactor(this);
      };
      
      Moments<double>
      InitialiseMoments() const
      {
	
	//std::cout<<"here0"<<std::endl;
	return  GammaModel<double>::CalcSample(m_parent1_node, m_parent2_node);
      }

      double
      CalcLogNorm() const 
      {
	//boost::mutex::scoped_lock lock(m_mutex);
	return m_LogNorm;//m_Model.CalcLogNorm();
      }


      NaturalParameters<double>
      GetNaturalNot(const VariableNode<double>* v) const;
      
    private: 
      VariableNode<double> *m_parent1_node, *m_parent2_node, *m_child_node;
     mutable  double  m_LogNorm;
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
     *  GAUSSIAN
     *
     **************************************************************************************
     **************************************************************************************
    //  **************************************************************************************/
    //template<>
    inline
    NaturalParameters<double>
    Factor< GaussianModel<double> ,double>::GetNaturalNot(const VariableNode<double>* v) const
    {
      if (v==m_parent1_node)
	{
	  Moments<double> prec = m_parent2_node->GetMoments();
	  Moments<double> child = m_child_node->GetMoments();
	  return GaussianModel<double>::CalcNP2Parent1(prec,child);
	}
      else if (v==m_parent2_node)
	{
	  Moments<double> parent1 = m_parent1_node->GetMoments();
	  Moments<double> child = m_child_node->GetMoments();
	  return GaussianModel<double>::CalcNP2Parent2(parent1,child);
	}
      else if (v == m_child_node) 
	{
	  Moments<double> parent1 = m_parent1_node->GetMoments();
	  Moments<double> prec = m_parent2_node->GetMoments();
	  m_LogNorm = GaussianModel<double>::CalcLogNorm(parent1,prec);
	  return GaussianModel<double>::CalcNP2Data(parent1,prec);
	}
      else{
	throw ("Unknown Node in GetNaturalNot");
      }
    }

    /**************************************************************************************
     **************************************************************************************
     **************************************************************************************
     *
     *  Rectified GAUSSIAN
     *
     **************************************************************************************
     **************************************************************************************
    //  **************************************************************************************/
    //template<>
    inline
    NaturalParameters<double>
    Factor< RectifiedGaussianModel<double> ,double>::GetNaturalNot(const VariableNode<double>* v) const
    {
      if (v==m_parent1_node)
	{
	  Moments<double> prec = m_parent2_node->GetMoments();
	  Moments<double> child = m_child_node->GetMoments();
	  return RectifiedGaussianModel<double>::CalcNP2Parent1(prec,child);
	  //return m_NP2Parent1;
	}
      else if (v==m_parent2_node)
	{
	  Moments<double> parent1 = m_parent1_node->GetMoments();
	  Moments<double> child = m_child_node->GetMoments();
	  return RectifiedGaussianModel<double>::CalcNP2Parent2(parent1,child);
	}
      else if (v == m_child_node) 
	{
	  Moments<double> parent1 = m_parent1_node->GetMoments();
	  Moments<double> prec = m_parent2_node->GetMoments();
	  m_LogNorm = RectifiedGaussianModel<double>::CalcLogNorm(parent1,prec);
	  return RectifiedGaussianModel<double>::CalcNP2Data(parent1,prec);
	  
	  // std::cout<<"returning child"<<m_NP2Child<<" from gaussian"<<std::endl;
	  //return m_NP2Child;
	}
      else{
	throw ("Unknown Node in GetNaturalNot");
      }
    }

    /**************************************************************************************
     **************************************************************************************
     **************************************************************************************
     *
     *  GAMMA
     *
     **************************************************************************************
     **************************************************************************************
     **************************************************************************************/
  
  

    inline
    NaturalParameters<double>
    Factor<GammaModel<double>,double >::GetNaturalNot(const VariableNode<double>* v) const
    {
      boost::mutex::scoped_lock lock(m_mutex);
      if (v==m_parent1_node)
	{
	  std::cerr<<"Can't pass message to parent1 variable  of Gamma Distribution"<<v<<std::endl;
	  throw ("Can't pass message to parent1 variable  of Gamma Distribution");
	}
      else if (v==m_parent2_node)
	{
	  Moments<double> parent1  = m_parent1_node->GetMoments();
	  //Moments<double> parent2 = m_parent2_node->GetMoments();
	  Moments<double> child  = m_child_node->GetMoments();
	  return GammaModel<double>::CalcNP2Parent2(parent1,child);
	}
      else 
	{
	  Moments<double> parent1  = m_parent1_node->GetMoments();
	  Moments<double> parent2 = m_parent2_node->GetMoments();
	  //Moments<double> child  = m_child_node->GetMoments();
	  BOOST_ASSERT(parent2.size() == 2);
	  BOOST_ASSERT(parent2[0] > 0);
	  if (parent1[0] == 0) {
	    std::cout<<"parent1 "<<parent1<<std::endl;
	    std::cout<<"parent2"<<parent2<<std::endl;
	  
	  }
	  BOOST_ASSERT(parent1[0] != 0);
	  m_LogNorm = GammaModel<double>::CalcLogNorm(parent1,parent2);
	  return GammaModel<double>::CalcNP2Data(parent1, parent2);

	}
    }

    

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
