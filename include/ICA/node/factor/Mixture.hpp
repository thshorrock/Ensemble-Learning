#pragma once

#include "ICA/node/Node.hpp"
#include "ICA/message/NaturalParameters.hpp"
#include "ICA/message/Moments.hpp"
// #include "ICA/Message.hpp"

//#include "ICA/piping/Piping.hpp"
#include "ICA/node/variable/Dirichlet.hpp"
#include "ICA/node/variable/Hidden.hpp"
#include "ICA/exponential_model/Gaussian.hpp"
#include "ICA/exponential_model/Gamma.hpp"
#include "ICA/exponential_model/Discrete.hpp"
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
      
      Mixture( std::vector<VariableNode<double>*>& Parent1, 
	       std::vector<VariableNode<double>*>& Parent2,  
	       HiddenNode<DiscreteModel<double> >* Weights,
	       VariableNode<double>* child
	       
	       )
	:  m_parent1_nodes(Parent1), 
	   m_parent2_nodes(Parent2),
	   m_weights_node(Weights),
	   m_child_node(child),
	   m_LogNorm(0)
	   
	   // m_NP2Child(0.0,0.0), 
	   // m_NP2Weights(std::vector<double>(Parent1.size())),
	   
	   // m_NP2Parent1(Parent1.size()), m_NP2Parent2(Parent1.size())
      {
	
	BOOST_ASSERT(Parent1.size() == Parent2.size());
	
	for(size_t i=0;i<Parent1.size();++i){
	  Parent1[i]->AddChildFactor(this);
	  m_parent2_nodes[i]->AddChildFactor(this);
	}
	Weights->AddChildFactor(this);
    	child->SetParentFactor(this);
	
	
      };
      
      
      Moments<double>
      InitialiseMoments() const
      {
	return GaussianModel<double>::CalcSample(m_parent1_nodes,m_parent2_nodes, m_weights_node );
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

      std::vector<VariableNode<double> *> m_parent1_nodes, m_parent2_nodes;
      HiddenNode<DiscreteModel<double> > *m_weights_node;
      VariableNode<double>  *m_child_node;
      
      
      //he NaturalParameters that serve as messages to be passed around
      mutable double m_LogNorm;
      // NaturalParameters<double> m_NP2Child;
      // NaturalParameters<double> m_NP2Weights;
      // std::vector<NaturalParameters<double> > m_NP2Parent1;
      // std::vector<NaturalParameters<double> > m_NP2Parent2;
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
      
      Mixture( std::vector<VariableNode<double>*>& Parent1, 
	       std::vector<VariableNode<double>*>& Parent2,  
	       HiddenNode<DiscreteModel<double> >* Weights,
	       VariableNode<double>* child
	       
	       )
	:  m_parent1_nodes(Parent1), 
	   m_parent2_nodes(Parent2),
	   m_weights_node(Weights),
	   m_child_node(child),
	   m_LogNorm(0)
	   
	   // m_NP2Child(0.0,0.0), 
	   // m_NP2Weights(std::vector<double>(Parent1.size())),
	   
	   // m_NP2Parent1(Parent1.size()), m_NP2Parent2(Parent1.size())
      {
	
	BOOST_ASSERT(Parent1.size() == Parent2.size());
	
	for(size_t i=0;i<Parent1.size();++i){
	  Parent1[i]->AddChildFactor(this);
	  m_parent2_nodes[i]->AddChildFactor(this);
	}
	Weights->AddChildFactor(this);
    	child->SetParentFactor(this);
      }
      
      Moments<double>
      InitialiseMoments() const
      {
	return RectifiedGaussianModel<double>::CalcSample(m_parent1_nodes,m_parent2_nodes, m_weights_node );
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

      std::vector<VariableNode<double> *> m_parent1_nodes, m_parent2_nodes;
      HiddenNode<DiscreteModel<double> > *m_weights_node;
      VariableNode<double>  *m_child_node;
      
      
      //he NaturalParameters that serve as messages to be passed around
      mutable double m_LogNorm;
      // NaturalParameters<double> m_NP2Child;
      // NaturalParameters<double> m_NP2Weights;
      // std::vector<NaturalParameters<double> > m_NP2Parent1;
      // std::vector<NaturalParameters<double> > m_NP2Parent2;
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
	  for(size_t i=0;i<m_parent1_nodes.size();++i){
	    const Moments<double>& parent1 = m_parent1_nodes[i]->GetMoments();
	    const Moments<double>& prec = m_parent2_nodes[i]->GetMoments();
	    NP2Child
	      += GaussianModel<double>::CalcNP2Data( parent1, prec) * weights[i];
	    
	    m_LogNorm 
	      += GaussianModel<double>::CalcLogNorm(parent1,prec) * weights[i];
	  }
	  return NP2Child;
	  // // std::cout<<"returning child"<<m_NP2Child<<" from mixture gaussian"<<std::endl;
	  // return m_NP2Child;
	}
      else if (v == m_weights_node) 
	{
	  NaturalParameters<double> NP2Weights(m_parent1_nodes.size());
	  

	  const Moments<double>& child = m_child_node->GetMoments();
	  for(size_t i=0;i<m_parent1_nodes.size();++i){
	    const Moments<double>& parent1 = m_parent1_nodes[i]->GetMoments();
	    const Moments<double>& prec = m_parent2_nodes[i]->GetMoments();
	    NP2Weights[i] = GaussianModel<double>::CalcAvLog(parent1, prec,child);
	  }
	  return NP2Weights;
	  // std::cout<<"returning weight"<<m_NP2Weights<<" from mixture gaussian"<<std::endl;
	  //std::cout<<"RETURNING m_weights = "<<m_NP2Weights<<std::endl;

	  // return m_NP2Weights;
	}
      else 
	{
	  std::vector<VariableNode<double> *>::const_iterator it;
	  it = PARALLEL_FIND(m_parent1_nodes.begin(), m_parent1_nodes.end(), v);
	  if (it!=m_parent1_nodes.end()) 
	    {
	      size_t i = it-m_parent1_nodes.begin();
	      const Moments<double>& weights = m_weights_node->GetMoments();
	      const Moments<double>& prec = m_parent2_nodes[i]->GetMoments();
	      const Moments<double>& child = m_child_node->GetMoments();
	      return GaussianModel<double>::CalcNP2Parent1(prec,child)  * weights[i];
	      //	std::cout<<"returning parent1["<<i<<"] "<<m_NP2Parent1[i]<<" from mixture gaussian"<<std::endl;
	      // return m_NP2Parent1[i];
	    }
	  else
	    {
	      it = PARALLEL_FIND(m_parent2_nodes.begin(), m_parent2_nodes.end(), v);
	      size_t i = it-m_parent2_nodes.begin();
	      const Moments<double>& weights = m_weights_node->GetMoments();
	      const Moments<double>& parent1 = m_parent1_nodes[i]->GetMoments();
	      const Moments<double>& child = m_child_node->GetMoments();
	      return GaussianModel<double>::CalcNP2Parent2(parent1,child) * weights[i];
	      //	std::cout<<"returning prec["<<i<<"] "<<m_NP2Parent2[i]<<" from mixture gaussian"<<std::endl;
	      // return m_NP2Parent2[i];
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
	  for(size_t i=0;i<m_parent1_nodes.size();++i){
	    const Moments<double>& parent1 = m_parent1_nodes[i]->GetMoments();
	    const Moments<double>& prec = m_parent2_nodes[i]->GetMoments();
	    NP2Child
	      += RectifiedGaussianModel<double>::CalcNP2Data( parent1, prec) * weights[i];
	    
	    m_LogNorm 
	      += RectifiedGaussianModel<double>::CalcLogNorm(parent1,prec) * weights[i];
	  }
	  return NP2Child;
	  // // std::cout<<"returning child"<<m_NP2Child<<" from mixture gaussian"<<std::endl;
	  // return m_NP2Child;
	}
      else if (v == m_weights_node) 
	{
	  NaturalParameters<double> NP2Weights(m_parent1_nodes.size());
	  

	  const Moments<double>& child = m_child_node->GetMoments();
	  for(size_t i=0;i<m_parent1_nodes.size();++i){
	    const Moments<double>& parent1 = m_parent1_nodes[i]->GetMoments();
	    const Moments<double>& prec = m_parent2_nodes[i]->GetMoments();
	    NP2Weights[i] = RectifiedGaussianModel<double>::CalcAvLog(parent1, prec,child);
	  }
	  return NP2Weights;
	  // std::cout<<"returning weight"<<m_NP2Weights<<" from mixture gaussian"<<std::endl;
	  //std::cout<<"RETURNING m_weights = "<<m_NP2Weights<<std::endl;

	  // return m_NP2Weights;
	}
      else 
	{
	  std::vector<VariableNode<double> *>::const_iterator it;
	  it = PARALLEL_FIND(m_parent1_nodes.begin(), m_parent1_nodes.end(), v);
	  if (it!=m_parent1_nodes.end()) 
	    {
	      size_t i = it-m_parent1_nodes.begin();
	      const Moments<double>& weights = m_weights_node->GetMoments();
	      const Moments<double>& prec = m_parent2_nodes[i]->GetMoments();
	      const Moments<double>& child = m_child_node->GetMoments();
	      return RectifiedGaussianModel<double>::CalcNP2Parent1(prec,child)  * weights[i];
	      //	std::cout<<"returning parent1["<<i<<"] "<<m_NP2Parent1[i]<<" from mixture gaussian"<<std::endl;
	      // return m_NP2Parent1[i];
	    }
	  else
	    {
	      it = PARALLEL_FIND(m_parent2_nodes.begin(), m_parent2_nodes.end(), v);
	      size_t i = it-m_parent2_nodes.begin();
	      const Moments<double>& weights = m_weights_node->GetMoments();
	      const Moments<double>& parent1 = m_parent1_nodes[i]->GetMoments();
	      const Moments<double>& child = m_child_node->GetMoments();
	      return RectifiedGaussianModel<double>::CalcNP2Parent2(parent1,child) * weights[i];
	      //	std::cout<<"returning prec["<<i<<"] "<<m_NP2Parent2[i]<<" from mixture gaussian"<<std::endl;
	      // return m_NP2Parent2[i];
	    }
	}
      throw ("Unknown Node in GetNaturalNot");
    }

    // inline
    // void 
    // Mixture<GaussianModel<double>,double >::Iterate()
    // {

    //   Moments<double> weights = m_weights_node->GetMoments();
      
    //   // Moments<double> av_parent1;
    //   // Moments<double> av_prec;
    //    Moments<double> child = m_child_node->GetMoments();

    //    // std::cout<<"Weights  = "<<weights<<std::endl;
    //    // std::cout<<"parent1 size = "<<m_parent1_nodes.size()<<std::endl;
      

    //   for(size_t i=0;i<m_parent1_nodes.size();++i){
	
    // 	// GaussianModel<double> model;

    // 	Moments<double> parent1 = m_parent1_nodes[i]->GetMoments();
    // 	Moments<double> prec = m_parent2_nodes[i]->GetMoments();

    // 	// model.SetAvParent1(parent1[0]);
    // 	// model.SetAvParent1Squared(parent1[1]);
    // 	// model.SetAvParent2(prec[0]);
    // 	// model.SetData(child[0]);
	
    // 	m_NP2Parent1[i]
    // 	  = GaussianModel<double>::CalcNP2Parent1(prec,child)  * weights[i];
	
    // 	m_NP2Parent2[i]
    // 	  = GaussianModel<double>::CalcNP2Parent2(parent1,child) * weights[i];
      
    // 	// Moments<double> parent1 = parent1->GetMoments() * weights[i];
    // 	// Moments<double> prec = prec->GetMoments() * weights[i];

    // 	if (i == 0)
    // 	  {
    // 	    m_NP2Child 
    // 	      = GaussianModel<double>::CalcNP2Data( parent1, prec) * weights[i];
	    
    // 	    m_LogNorm 
    // 	      = GaussianModel<double>::CalcLogNorm(parent1, prec) * weights[i];
    // 	  } 
    // 	else
    // 	  {
    // 	    m_NP2Child
    // 	      += GaussianModel<double>::CalcNP2Data( parent1, prec)  * weights[i];
	    
    // 	    m_LogNorm 
    // 	      += GaussianModel<double>::CalcLogNorm(parent1,prec) * weights[i];
    // 	  }

    // 	m_NP2Weights[i] = GaussianModel<double>::CalcAvLog(parent1, prec,child);

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
	
    // 	model.SetAvParent1(av_parent1[0]);
    // 	model.SetAvParent1Squared(av_parent1[1]);
    // 	model.SetAvParent2(av_prec[0]);
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

      // 	Moments<double> parent1 = m_parent1_nodes[i]->GetMoments();
      // 	Moments<double> prec = m_parent2_node[i]->GetMoments();
      // 	Moments<double> weights = m_weights_node->GetMoments();
      // 	Moments<double> child = m_child_node->GetMoments();

      // 	model.SetAvParent1(parent1[0]);
      // 	model.SetAvParent1Squared(parent1[1]);
      // 	model.SetAvParent2(prec[0]);
      // 	model.SetData(child[0]);

      
      // 	m_LogNorm = model.CalcLogNorm();
      // 	m_FFunction = model.GetFFunction();

      // 	model.CalcNP2Parent1(m_NP2Parent1[i]);
      // 	model.CalcNP2Parent2(m_NP2Parent2[i]);
      // 	model.CalcNP2Data(m_NP2Child);

      // }
      // m_parent1_node -> ReceiveNaturalParameters(m_NP2Parent1);
      // m_parent2_node -> ReceiveNaturalParameters(m_NP2Parent2 );
      // m_child_node -> ReceiveNaturalParameters(m_NP2Child);
  
  }
}

// template<class Model>
// void
// ICR::ICA::Mixture<Model>::Reset()
// {}
