#pragma once

#include "ICA/Node.hpp"
#include "ICA/NaturalParameters.hpp"
#include "ICA/Moments.hpp"
// #include "ICA/Message.hpp"

//#include "ICA/piping/Piping.hpp"
#include "ICA/variable/DirichletNode.hpp"
#include "ICA/variable/HiddenNode.hpp"
#include "ICA/exponential_models/GaussianModel.hpp"
#include "ICA/exponential_models/RectifiedGaussianModel.hpp"
#include "ICA/exponential_models/GammaModel.hpp"
#include "ICA/exponential_models/DiscreteModel.hpp"
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
      //doublehe model
      //GaussianModel m_Model;
      //doublehe variable nodes that are attached to this function node
       VariableNode<double> *m_mean_node, *m_precision_node, *m_child_node;
      
      
      //doublehe NaturalParameters that serve as messages to be passed around
      mutable double m_LogNorm;
     // NaturalParameters<double> m_NP2Child;
     //  NaturalParameters<double> m_NP2Mean;
     //  NaturalParameters<double> m_NP2Precision;
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

	const Moments<double> M =RectifiedGaussianModel<double>::CalcSample(Mean, Precision);
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
      //doublehe model
      //RectifiedGaussianModel m_Model;
      //doublehe variable nodes that are attached to this function node
       VariableNode<double> *m_mean_node, *m_precision_node, *m_child_node;
      
      
      //doublehe NaturalParameters that serve as messages to be passed around
      mutable double m_LogNorm;
     // NaturalParameters<double> m_NP2Child;
     //  NaturalParameters<double> m_NP2Mean;
     //  NaturalParameters<double> m_NP2Precision;
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
      
      // void
      // ReceiveMoments(const Moments& m);

      // void
      // Reset();

      // void 
      // Iterate();
    private: 
      //doublehe model
      // GammaModel m_Model;
      //doublehe variable nodes that are attached to this function node
      VariableNode<double> *m_shape_node, *m_iscale_node, *m_child_node;
      //Failed messages that will need to be resent when have more information
      
      //doublehe NaturalParameters that serve as messages to be passed around
     mutable  double  m_LogNorm;
      // NaturalParameters<double> m_NP2Child;
      // NaturalParameters<double> m_NP2IScale;
      //No message can be sent to Shape Node (no Conjutage Variable).
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
	//boost::mutex::scoped_lock lock(m_mutex);
	//return m_Model.CalcLogNorm();
      }


      NaturalParameters<double>
      GetNaturalNot(const VariableNode<double>* v) const;
      

      // void 
      // Iterate();
    private: 
      //doublehe model
      //DirichletModel m_Model;
      //doublehe variable nodes that are attached to this function node
       VariableNode<double> *m_prior_node,  *m_child_node;
      
      
      //doublehe NaturalParameters that serve as messages to be passed around
      mutable double m_LogNorm;
      // NaturalParameters<double> m_NP2Child;
      //NaturalParameters<double> m_NPprior;
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
      //doublehe model
      //DiscreteModel m_Model;
      //doublehe variable nodes that are attached to this function node
      // VariableNode *m_mean_node, *m_precision_node, *m_child_node;
      VariableNode<double> *m_prior_node, *m_child_node;
      mutable double  m_LogNorm;
      //doublehe NaturalParameters that serve as messages to be passed around
      // NaturalParameters<double> m_NP2Child;
      // NaturalParameters<double> m_NP2Prior;
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
	  //std::cout<<"returning mean"<<m_NP2Mean<<" from gaussian"<<std::endl;
	  
	  Moments<double> prec = m_precision_node->GetMoments();
	  Moments<double> child = m_child_node->GetMoments();
	  return GaussianModel<double>::CalcNP2Mean(prec,child);
	  //return m_NP2Mean;
	}
      else if (v==m_precision_node)
	{
	  Moments<double> mean = m_mean_node->GetMoments();
	  Moments<double> child = m_child_node->GetMoments();
	  return GaussianModel<double>::CalcNP2Precision(mean,child);
	  // std::cout<<"returning prec"<<m_NP2Precision<<" from gaussian"<<std::endl;
	  //return m_NP2Precision;
	}
      else if (v == m_child_node) 
	{
	  Moments<double> mean = m_mean_node->GetMoments();
	  Moments<double> prec = m_precision_node->GetMoments();
	  m_LogNorm = GaussianModel<double>::CalcLogNorm(mean,prec);
	  return GaussianModel<double>::CalcNP2Data(mean,prec);
	  
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
	  //std::cout<<"returning mean"<<m_NP2Mean<<" from gaussian"<<std::endl;
	  
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
	  // std::cout<<"returning prec"<<m_NP2Precision<<" from gaussian"<<std::endl;
	  //return m_NP2Precision;
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

    // inline
    // void 
    // Factor<GaussianModel<double>,double >::Iterate()
    // {
    //   //      GaussianModel<double> model;
    //   // Moments<double> mean = m_mean_node->GetMoments();
    //   // Moments<double> prec = m_precision_node->GetMoments();
    //   // Moments<double> child = m_child_node->GetMoments();
      
    //   // std::cout<<"GMean = "<<mean<<std::endl;
    //   // std::cout<<"GPrec = "<<prec<<std::endl;
    //   // std::cout<<"GChild = "<<child<<std::endl;


    //   // m_NP2Mean       = GaussianModel<double>::CalcNP2Mean(prec,child);
    //   // m_NP2Precision  = GaussianModel<double>::CalcNP2Precision(mean,child);
    //   // m_NP2Child      = GaussianModel<double>::CalcNP2Data(mean,prec);
      
    //   // m_LogNorm = GaussianModel<double>::CalcLogNorm(mean,prec);


    //   // std::cout<<"NP2Gaussian = "<<m_NP2Child<<std::endl;
    //   // m_mean_node -> ReceiveNaturalParameters(m_NP2Mean);
    //   // m_precision_node -> ReceiveNaturalParameters(m_NP2Precision);
    //   // m_child_node -> ReceiveNaturalParameters(m_NP2Child);
    // }

    /**************************************************************************************
     **************************************************************************************
     **************************************************************************************
     *
     *  GAMMA
     *
     **************************************************************************************
     **************************************************************************************
     **************************************************************************************/
  
  
    // inline
    // void 
    // Factor<GammaModel>::Reset()
    // {
    //   // m_Model.Reset();
    // }

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
	  // std::cout<<"returning scale"<<m_NP2IScale<<" from gamma"<<std::endl;
	  //return m_NP2IScale;
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

	  //std::cout<<"returning child"<<m_NP2Child<<" from gamma"<<std::endl;
	  //return m_NP2Child;
	}
    }

    
    // inline
    // void 
    // Factor<GammaModel<double>,double >::Iterate()
    // {
  
    //   // Moments<double> shape  = m_shape_node->GetMoments();
    //   // Moments<double> iscale = m_iscale_node->GetMoments();
    //   // Moments<double> child  = m_child_node->GetMoments();
    //   // //      GammaModel<double> model;
      
    //   // std::cout<<"Gshape = "<<shape<<std::endl;
    //   // std::cout<<"Giscale = "<<iscale<<std::endl;
    //   // std::cout<<"Gchild = "<<child<<std::endl;
      
    //   // m_LogNorm = GammaModel<double>::CalcLogNorm(shape,iscale);

    //   // m_NP2Child = GammaModel<double>::CalcNP2Data(shape, iscale);
    //   // m_NP2IScale = GammaModel<double>::CalcNP2IScale(shape,child);

    //   // std::cout<<"NP2Gamma = "<<m_NP2Child<<std::endl;

    //   // m_iscale_node -> ReceiveNaturalParameters(m_NP2IScale);
    //   // m_child_node -> ReceiveNaturalParameters(m_NP2Child);
    //   // m_messages->CreateFactorMessage(this, m_iscale_node);
    //   // m_messages->CreateFactorMessage(this, m_child_node);
      

    // }

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
      // else if (v==m_precision_node)
      // 	{
      // 	  return m_NP2Precision;
      // 	}
      // else if (v == m_child_node) 
      // 	{
      // 	  return m_NP2Child;
      // 	}
       	throw ("Unknown Node in GetNaturalNot");
    }

    // inline
    // void 
    // Factor<DirichletModel<double>,double >::Iterate()
    // {
    //   // std::cout<<"Dirichlet Iterate Factor "<<this<<std::endl;
    //   // boost::thread th1(boost::bind(&DirichletModel::CalcNP2Mean, boost::ref(m_Model), m_NP2Mean));
    //   // boost::thread th2(boost::bind(&DirichletModel::CalcNP2Precision, boost::ref(m_Model), m_NP2Precision));
    //   // boost::thread th3(boost::bind(&DirichletModel::CalcNP2Data, boost::ref(m_Model), m_NP2Child));
      
    //   // th1.join();
    //   // th2.join();
    //   // th3.join();

    //   //      DirichletModel<double> model;
    //   Moments<double> prior = m_prior_node->GetMoments();
    //   // Moments<double> child = m_child_node->GetMoments();

    //   // model.SetValues(m_prior_node->GetMoments());

    //   m_LogNorm = DirichletModel<double>::CalcLogNorm(prior);

    //   m_NP2Child = DirichletModel<double>::CalcNP2Data(prior);

    //   // std::cout<<"Dirichlet Child NP = "<<m_NP2Child<<std::endl;
    //   // std::cout<<"Dirichlet LogNorm = "<<m_LogNorm<<std::endl;
    //   // std::cout<<"Dirichlet Child NP = "<<m_NP2Child<<std::endl;
    //   // m_mean_node -> ReceiveNaturalParameters(m_NP2Mean);
    //   // m_precision_node -> ReceiveNaturalParameters(m_NP2Precision);
    //   // m_child_node -> ReceiveNaturalParameters(m_NP2Child);
    // }

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

	  //std::cout<<"returning prior "<<m_NP2Prior<<" from Discrete"<<std::endl;
       // return m_NP2Prior;
      	}
      else 
       	{
	  //std::cout<<"returning child "<<m_NP2Child<<" from Discrete"<<std::endl;

	  Moments<double> prior  = m_prior_node->GetMoments(); 
	  m_LogNorm = DiscreteModel<double>::CalcLogNorm(prior);

	  return  DiscreteModel<double>::CalcNP2Data(prior);/// = child;
	  // return m_NP2Child;
	}
    }

    // inline
    // void 
    // Factor<DiscreteModel<double>,double>::Iterate()
    // {
    //   // std::cout<<"Discrete Iterate Factor "<<this<<std::endl;


    //   Moments<double> prior  = m_prior_node->GetMoments(); 
    //   //std::cout<<"prior size = "<<prior.size()<<std::endl;

    //   Moments<double> child  = m_child_node->GetMoments();

    //   // std::cout<<"PRIOR = "<<prior<<std::endl;
    //   // std::cout<<"CHILD = "<<child<<std::endl;


    //   // m_NP2Child = prior;
    //   // m_NP2Prior = child;


      
    //    m_NP2Child  = DiscreteModel<double>::CalcNP2Data(prior);/// = child;
    //    m_NP2Prior  = DiscreteModel<double>::CalcNP2Prior(child);/// = child;
    //    // model.SetData(child);
      
    //    // std::cout<<"Discrete Prior NP = "<<m_NP2Prior<<std::endl;
    //    //std::cout<<"Discrete Child NP = "<<m_NP2Child<<std::endl;

    //   m_LogNorm = DiscreteModel<double>::CalcLogNorm(prior);

    //   // std::cout<<"Discrete LogNorm = "<<m_LogNorm<<std::endl;
    //   // model.CalcNP2Data(m_NP2Child);
    //   // m_Model.CalcNP2Prior(m_NP2Prior); 

    // }

  }
}

// template<class Model>
// void
// ICR::ICA::Mixture<Model>::Reset()
// {}
