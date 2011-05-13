#pragma once

#include "ICA/Node.hpp"
#include "ICA/NaturalParameters.hpp"
#include "ICA/Moments.hpp"
// #include "ICA/Message.hpp"

//#include "ICA/piping/Piping.hpp"
#include "ICA/variable/DirichletNode.hpp"
#include "ICA/variable/HiddenNode.hpp"
#include "ICA/exponential_models/GaussianModel.hpp"
#include "ICA/exponential_models/GammaModel.hpp"
#include "ICA/exponential_models/DiscreteModel.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/none.hpp>
#include <vector>

namespace ICR{
  namespace ICA{
    

    template<class Model, class T = double>
    class Mixture :  public FactorNode<T>
    {
      //this will not compile if you try to call it directly.
      // You need to use one of the specializations.
    };
   
    /******************************************************************************
     * GaussianModel Specialisation Double
     ******************************************************************************/
    template<>
    class Mixture<GaussianModel<double> > : public FactorNode<double>
    {
      Mixture(const Mixture<GaussianModel<double> >& f) {};
    public:
      
      Mixture( std::vector<VariableNode<double>*>& Mean, 
	       std::vector<VariableNode<double>*>& Precision,  
	       HiddenNode<DiscreteModel<double> >* Weights,
	       VariableNode<double>* child
	       
	       )
	:  m_mean_nodes(Mean), 
	   m_precision_nodes(Precision),
	   m_weights_node(Weights),
	   m_child_node(child),
	   m_LogNorm(0)
	   
	   // m_NP2Child(0.0,0.0), 
	   // m_NP2Weights(std::vector<double>(Mean.size())),
	   
	   // m_NP2Mean(Mean.size()), m_NP2Precision(Mean.size())
      {
	
	BOOST_ASSERT(Mean.size() == Precision.size());
	
	for(size_t i=0;i<Mean.size();++i){
	  Mean[i]->AddChildFactor(this);
	  m_precision_nodes[i]->AddChildFactor(this);
	}
	Weights->AddChildFactor(this);
    	child->SetParentFactor(this);
	
	
      };
      
      
      Moments<double>
      InitialiseMoments() const
      {
	return GaussianModel<double>::CalcSample(m_mean_nodes,m_precision_nodes, m_weights_node );
      }

      // double
      // GetFFunction() const
      // {
      // 	return m_FFunction;
      // }

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

      std::vector<VariableNode<double> *> m_mean_nodes, m_precision_nodes;
      HiddenNode<DiscreteModel<double> > *m_weights_node;
      VariableNode<double>  *m_child_node;
      
      
      //he NaturalParameters that serve as messages to be passed around
      mutable double m_LogNorm;
      // NaturalParameters<double> m_NP2Child;
      // NaturalParameters<double> m_NP2Weights;
      // std::vector<NaturalParameters<double> > m_NP2Mean;
      // std::vector<NaturalParameters<double> > m_NP2Precision;
      mutable boost::mutex m_mutex;
    };
    
    /******************************************************************************
     * RectifiedGaussianModel Specialisation Double
     ******************************************************************************/
    template<>
    class Mixture<RectifiedGaussianModel<double> > : public FactorNode<double>
    {
      Mixture(const Mixture<RectifiedGaussianModel<double> >& f) {};
    public:
      
      Mixture( std::vector<VariableNode<double>*>& Mean, 
	       std::vector<VariableNode<double>*>& Precision,  
	       HiddenNode<DiscreteModel<double> >* Weights,
	       VariableNode<double>* child
	       
	       )
	:  m_mean_nodes(Mean), 
	   m_precision_nodes(Precision),
	   m_weights_node(Weights),
	   m_child_node(child),
	   m_LogNorm(0)
	   
	   // m_NP2Child(0.0,0.0), 
	   // m_NP2Weights(std::vector<double>(Mean.size())),
	   
	   // m_NP2Mean(Mean.size()), m_NP2Precision(Mean.size())
      {
	
	BOOST_ASSERT(Mean.size() == Precision.size());
	
	for(size_t i=0;i<Mean.size();++i){
	  Mean[i]->AddChildFactor(this);
	  m_precision_nodes[i]->AddChildFactor(this);
	}
	Weights->AddChildFactor(this);
    	child->SetParentFactor(this);
      }
      
      Moments<double>
      InitialiseMoments() const
      {
	return GaussianModel<double>::CalcSample(m_mean_nodes,m_precision_nodes, m_weights_node );
      }

      // double
      // GetFFunction() const
      // {
      // 	return m_FFunction;
      // }

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

      std::vector<VariableNode<double> *> m_mean_nodes, m_precision_nodes;
      HiddenNode<DiscreteModel<double> > *m_weights_node;
      VariableNode<double>  *m_child_node;
      
      
      //he NaturalParameters that serve as messages to be passed around
      mutable double m_LogNorm;
      // NaturalParameters<double> m_NP2Child;
      // NaturalParameters<double> m_NP2Weights;
      // std::vector<NaturalParameters<double> > m_NP2Mean;
      // std::vector<NaturalParameters<double> > m_NP2Precision;
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
    Mixture< GaussianModel<double> ,double>::GetNaturalNot(const VariableNode<double>* v) const
    {
      if (v == m_child_node) 
	{

	  NaturalParameters<double> NP2Child(std::vector<double>(2));
	  m_LogNorm = 0;
	  
	  const Moments<double>& weights = m_weights_node->GetMoments();
	  for(size_t i=0;i<m_mean_nodes.size();++i){
	    const Moments<double>& mean = m_mean_nodes[i]->GetMoments();
	    const Moments<double>& prec = m_precision_nodes[i]->GetMoments();
	    NP2Child
	      += GaussianModel<double>::CalcNP2Data( mean, prec) * weights[i];
	    
	    m_LogNorm 
	      += GaussianModel<double>::CalcLogNorm(mean,prec) * weights[i];
	  }
	  return NP2Child;
	  // // std::cout<<"returning child"<<m_NP2Child<<" from mixture gaussian"<<std::endl;
	  // return m_NP2Child;
	}
      else if (v == m_weights_node) 
	{
	  NaturalParameters<double> NP2Weights(m_mean_nodes.size());
	  

	  const Moments<double>& child = m_child_node->GetMoments();
	  for(size_t i=0;i<m_mean_nodes.size();++i){
	    const Moments<double>& mean = m_mean_nodes[i]->GetMoments();
	    const Moments<double>& prec = m_precision_nodes[i]->GetMoments();
	    NP2Weights[i] = GaussianModel<double>::CalcAvLog(mean, prec,child);
	  }
	  return NP2Weights;
	  // std::cout<<"returning weight"<<m_NP2Weights<<" from mixture gaussian"<<std::endl;
	  //std::cout<<"RETURNING m_weights = "<<m_NP2Weights<<std::endl;

	  // return m_NP2Weights;
	}
      else 
	{
	  std::vector<VariableNode<double> *>::const_iterator it;
	  it = PARALLEL_FIND(m_mean_nodes.begin(), m_mean_nodes.end(), v);
	  if (it!=m_mean_nodes.end()) 
	    {
	      size_t i = it-m_mean_nodes.begin();
	      const Moments<double>& weights = m_weights_node->GetMoments();
	      const Moments<double>& prec = m_precision_nodes[i]->GetMoments();
	      const Moments<double>& child = m_child_node->GetMoments();
	      return GaussianModel<double>::CalcNP2Mean(prec,child)  * weights[i];
	      //	std::cout<<"returning mean["<<i<<"] "<<m_NP2Mean[i]<<" from mixture gaussian"<<std::endl;
	      // return m_NP2Mean[i];
	    }
	  else
	    {
	      it = PARALLEL_FIND(m_precision_nodes.begin(), m_precision_nodes.end(), v);
	      size_t i = it-m_precision_nodes.begin();
	      const Moments<double>& weights = m_weights_node->GetMoments();
	      const Moments<double>& mean = m_mean_nodes[i]->GetMoments();
	      const Moments<double>& child = m_child_node->GetMoments();
	      return GaussianModel<double>::CalcNP2Precision(mean,child) * weights[i];
	      //	std::cout<<"returning prec["<<i<<"] "<<m_NP2Precision[i]<<" from mixture gaussian"<<std::endl;
	      // return m_NP2Precision[i];
	    }
	}
      throw ("Unknown Node in GetNaturalNot");
    }

  
    /**************************************************************************************
     **************************************************************************************
     **************************************************************************************
     *
     *  RECTIFIED GAUSSIAN
     *
     **************************************************************************************
     **************************************************************************************
     //  **************************************************************************************/
    //template<>
    inline
    NaturalParameters<double>
    Mixture< RectifiedGaussianModel<double> ,double>::GetNaturalNot(const VariableNode<double>* v) const
    {
      if (v == m_child_node) 
	{

	  NaturalParameters<double> NP2Child(std::vector<double>(2));
	  m_LogNorm = 0;
	  
	  const Moments<double>& weights = m_weights_node->GetMoments();
	  for(size_t i=0;i<m_mean_nodes.size();++i){
	    const Moments<double>& mean = m_mean_nodes[i]->GetMoments();
	    const Moments<double>& prec = m_precision_nodes[i]->GetMoments();
	    NP2Child
	      += RectifiedGaussianModel<double>::CalcNP2Data( mean, prec) * weights[i];
	    
	    m_LogNorm 
	      += RectifiedGaussianModel<double>::CalcLogNorm(mean,prec) * weights[i];
	  }
	  return NP2Child;
	  // // std::cout<<"returning child"<<m_NP2Child<<" from mixture gaussian"<<std::endl;
	  // return m_NP2Child;
	}
      else if (v == m_weights_node) 
	{
	  NaturalParameters<double> NP2Weights(m_mean_nodes.size());
	  

	  const Moments<double>& child = m_child_node->GetMoments();
	  for(size_t i=0;i<m_mean_nodes.size();++i){
	    const Moments<double>& mean = m_mean_nodes[i]->GetMoments();
	    const Moments<double>& prec = m_precision_nodes[i]->GetMoments();
	    NP2Weights[i] = RectifiedGaussianModel<double>::CalcAvLog(mean, prec,child);
	  }
	  return NP2Weights;
	  // std::cout<<"returning weight"<<m_NP2Weights<<" from mixture gaussian"<<std::endl;
	  //std::cout<<"RETURNING m_weights = "<<m_NP2Weights<<std::endl;

	  // return m_NP2Weights;
	}
      else 
	{
	  std::vector<VariableNode<double> *>::const_iterator it;
	  it = PARALLEL_FIND(m_mean_nodes.begin(), m_mean_nodes.end(), v);
	  if (it!=m_mean_nodes.end()) 
	    {
	      size_t i = it-m_mean_nodes.begin();
	      const Moments<double>& weights = m_weights_node->GetMoments();
	      const Moments<double>& prec = m_precision_nodes[i]->GetMoments();
	      const Moments<double>& child = m_child_node->GetMoments();
	      return RectifiedGaussianModel<double>::CalcNP2Mean(prec,child)  * weights[i];
	      //	std::cout<<"returning mean["<<i<<"] "<<m_NP2Mean[i]<<" from mixture gaussian"<<std::endl;
	      // return m_NP2Mean[i];
	    }
	  else
	    {
	      it = PARALLEL_FIND(m_precision_nodes.begin(), m_precision_nodes.end(), v);
	      size_t i = it-m_precision_nodes.begin();
	      const Moments<double>& weights = m_weights_node->GetMoments();
	      const Moments<double>& mean = m_mean_nodes[i]->GetMoments();
	      const Moments<double>& child = m_child_node->GetMoments();
	      return RectifiedGaussianModel<double>::CalcNP2Precision(mean,child) * weights[i];
	      //	std::cout<<"returning prec["<<i<<"] "<<m_NP2Precision[i]<<" from mixture gaussian"<<std::endl;
	      // return m_NP2Precision[i];
	    }
	}
      throw ("Unknown Node in GetNaturalNot");
    }

    // inline
    // void 
    // Mixture<GaussianModel<double>,double >::Iterate()
    // {

    //   Moments<double> weights = m_weights_node->GetMoments();
      
    //   // Moments<double> av_mean;
    //   // Moments<double> av_prec;
    //    Moments<double> child = m_child_node->GetMoments();

    //    // std::cout<<"Weights  = "<<weights<<std::endl;
    //    // std::cout<<"mean size = "<<m_mean_nodes.size()<<std::endl;
      

    //   for(size_t i=0;i<m_mean_nodes.size();++i){
	
    // 	// GaussianModel<double> model;

    // 	Moments<double> mean = m_mean_nodes[i]->GetMoments();
    // 	Moments<double> prec = m_precision_nodes[i]->GetMoments();

    // 	// model.SetAvMean(mean[0]);
    // 	// model.SetAvMeanSquared(mean[1]);
    // 	// model.SetAvPrecision(prec[0]);
    // 	// model.SetData(child[0]);
	
    // 	m_NP2Mean[i]
    // 	  = GaussianModel<double>::CalcNP2Mean(prec,child)  * weights[i];
	
    // 	m_NP2Precision[i]
    // 	  = GaussianModel<double>::CalcNP2Precision(mean,child) * weights[i];
      
    // 	// Moments<double> mean = mean->GetMoments() * weights[i];
    // 	// Moments<double> prec = prec->GetMoments() * weights[i];

    // 	if (i == 0)
    // 	  {
    // 	    m_NP2Child 
    // 	      = GaussianModel<double>::CalcNP2Data( mean, prec) * weights[i];
	    
    // 	    m_LogNorm 
    // 	      = GaussianModel<double>::CalcLogNorm(mean, prec) * weights[i];
    // 	  } 
    // 	else
    // 	  {
    // 	    m_NP2Child
    // 	      += GaussianModel<double>::CalcNP2Data( mean, prec)  * weights[i];
	    
    // 	    m_LogNorm 
    // 	      += GaussianModel<double>::CalcLogNorm(mean,prec) * weights[i];
    // 	  }

    // 	m_NP2Weights[i] = GaussianModel<double>::CalcAvLog(mean, prec,child);

    // 	//std::cout<<"m_NP2Weights (b4 norm) ["<<i<<"] = "<<(m_NP2Weights[i])<<std::endl;

    // 	// 	m_LogNorm += model.CalcLogNorm();
    // 	// 	m_FFunction += model.GetFFunction();
    //   } 
    // }
      
    //Normalise weights
    // double norm = 0;
    // for(size_t i=0;i<m_NP2Weights.size();++i){
    // 	norm+=std::exp(m_NP2Weights[i]);
    // }
    // double lognorm = std::log(norm);
    // std::cout<<"log norm = "<<lognorm<<std::endl;

    // for(size_t i=0;i<m_NP2Weights.size();++i){
    // 	m_NP2Weights[i] -=lognorm;
    // }

    // {
    // 	GaussianModel<double> model;
	
    // 	model.SetAvMean(av_mean[0]);
    // 	model.SetAvMeanSquared(av_mean[1]);
    // 	model.SetAvPrecision(av_prec[0]);
    // 	model.SetData(child[0]);

    // 	model.CalcNP2Data(m_NP2Child);

      // 	m_LogNorm = model.CalcLogNorm();
      // 	m_FFunction = model.GetFFunction();
      // }
      
      // {


      // 	// DiscreteModel<double> model;
      // 	// model.SetValues(weights);

      // 	// model.CalcNP2Data(m_NP2Weights);

      // 	//std::cout<<"m_NP2Weights = "<<m_NP2Weights<<std::endl;

      // 	//	std::cout<<"m_NP2Wieghts size = "<<m_NP2Weights.size()<<std::endl;
	
      // }
      // GaussianModel<double> model;

      // 	Moments<double> mean = m_mean_nodes[i]->GetMoments();
      // 	Moments<double> prec = m_precision_node[i]->GetMoments();
      // 	Moments<double> weights = m_weights_node->GetMoments();
      // 	Moments<double> child = m_child_node->GetMoments();

      // 	model.SetAvMean(mean[0]);
      // 	model.SetAvMeanSquared(mean[1]);
      // 	model.SetAvPrecision(prec[0]);
      // 	model.SetData(child[0]);

      
      // 	m_LogNorm = model.CalcLogNorm();
      // 	m_FFunction = model.GetFFunction();

      // 	model.CalcNP2Mean(m_NP2Mean[i]);
      // 	model.CalcNP2Precision(m_NP2Precision[i]);
      // 	model.CalcNP2Data(m_NP2Child);

      // }
      // m_mean_node -> ReceiveNaturalParameters(m_NP2Mean);
      // m_precision_node -> ReceiveNaturalParameters(m_NP2Precision);
      // m_child_node -> ReceiveNaturalParameters(m_NP2Child);
  
  }
}

// template<class Model>
// void
// ICR::ICA::Mixture<Model>::Reset()
// {}
