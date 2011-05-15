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
      
      Factor( VariableNode<double>* Mean,  VariableNode<double>* Precision,  VariableNode<double>* Child)
	: m_mean_node(Mean),
	  m_precision_node(Precision),
	  m_child_node(Child)
      {
    	Mean->AddChildFactor(this);
    	Precision->AddChildFactor(this);
    	Child->SetParentFactor(this);

      };
      
      Moments<double>
      InitialiseMoments() const
      {
	return  GaussianModel<double>::CalcSample(m_mean_node,m_precision_node );
      }
      
      double
      CalcLogNorm() const 
      {
	return m_LogNorm;
      }


      NaturalParameters<double>
      GetNaturalNot(const VariableNode<double>* v) const;
      

    private: 
       VariableNode<double> *m_mean_node, *m_precision_node, *m_child_node;
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
      
      Factor( VariableNode<double>* Mean,  VariableNode<double>* Precision,  VariableNode<double>* Child)
	: m_mean_node(Mean),
	  m_precision_node(Precision),
	  m_child_node(Child)
      {
    	Mean->AddChildFactor(this);
    	Precision->AddChildFactor(this);
    	Child->SetParentFactor(this);

	//const Moments<double> M =RectifiedGaussianModel<double>::CalcSample(Mean, Precision);
      };
      
      
      Moments<double>
      InitialiseMoments() const
      {
	return  RectifiedGaussianModel<double>::CalcSample(m_mean_node,m_precision_node );
      }
      double
      CalcLogNorm() const 
      {
	return m_LogNorm;
      }


      NaturalParameters<double>
      GetNaturalNot(const VariableNode<double>* v) const;
      

    private: 
       VariableNode<double> *m_mean_node, *m_precision_node, *m_child_node;
      
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
      
      Factor( VariableNode<double>* Shape,  VariableNode<double>* IScale,  VariableNode<double>* Child)
	: m_shape_node(Shape),
	  m_iscale_node(IScale),
	  m_child_node(Child)
      {
    	Shape->AddChildFactor(this);
    	IScale->AddChildFactor(this);
    	Child->SetParentFactor(this);
      };
      
      Moments<double>
      InitialiseMoments() const
      {
	
	//std::cout<<"here0"<<std::endl;
	return  GammaModel<double>::CalcSample(m_shape_node, m_child_node);
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
      VariableNode<double> *m_shape_node, *m_iscale_node, *m_child_node;
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
      if (v==m_mean_node)
	{
	  Moments<double> prec = m_precision_node->GetMoments();
	  Moments<double> child = m_child_node->GetMoments();
	  return GaussianModel<double>::CalcNP2Mean(prec,child);
	}
      else if (v==m_precision_node)
	{
	  Moments<double> mean = m_mean_node->GetMoments();
	  Moments<double> child = m_child_node->GetMoments();
	  return GaussianModel<double>::CalcNP2Precision(mean,child);
	}
      else if (v == m_child_node) 
	{
	  Moments<double> mean = m_mean_node->GetMoments();
	  Moments<double> prec = m_precision_node->GetMoments();
	  m_LogNorm = GaussianModel<double>::CalcLogNorm(mean,prec);
	  return GaussianModel<double>::CalcNP2Data(mean,prec);
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
      if (v==m_mean_node)
	{
	  Moments<double> prec = m_precision_node->GetMoments();
	  Moments<double> child = m_child_node->GetMoments();
	  return RectifiedGaussianModel<double>::CalcNP2Mean(prec,child);
	  //return m_NP2Mean;
	}
      else if (v==m_precision_node)
	{
	  Moments<double> mean = m_mean_node->GetMoments();
	  Moments<double> child = m_child_node->GetMoments();
	  return RectifiedGaussianModel<double>::CalcNP2Precision(mean,child);
	}
      else if (v == m_child_node) 
	{
	  Moments<double> mean = m_mean_node->GetMoments();
	  Moments<double> prec = m_precision_node->GetMoments();
	  m_LogNorm = RectifiedGaussianModel<double>::CalcLogNorm(mean,prec);
	  return RectifiedGaussianModel<double>::CalcNP2Data(mean,prec);
	  
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
      if (v==m_shape_node)
	{
	  std::cerr<<"Can't pass message to shape variable  of Gamma Distribution"<<v<<std::endl;
	  throw ("Can't pass message to shape variable  of Gamma Distribution");
	}
      else if (v==m_iscale_node)
	{
	  Moments<double> shape  = m_shape_node->GetMoments();
	  //Moments<double> iscale = m_iscale_node->GetMoments();
	  Moments<double> child  = m_child_node->GetMoments();
	  return GammaModel<double>::CalcNP2IScale(shape,child);
	}
      else 
	{
	  Moments<double> shape  = m_shape_node->GetMoments();
	  Moments<double> iscale = m_iscale_node->GetMoments();
	  //Moments<double> child  = m_child_node->GetMoments();
	  BOOST_ASSERT(iscale.size() == 2);
	  BOOST_ASSERT(iscale[0] > 0);
	  if (shape[0] == 0) {
	    std::cout<<"shape "<<shape<<std::endl;
	    std::cout<<"iscale"<<iscale<<std::endl;
	  
	  }
	  BOOST_ASSERT(shape[0] != 0);
	  m_LogNorm = GammaModel<double>::CalcLogNorm(shape,iscale);
	  return GammaModel<double>::CalcNP2Data(shape, iscale);

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
