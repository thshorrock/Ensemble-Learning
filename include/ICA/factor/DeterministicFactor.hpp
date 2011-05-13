#pragma once

#include "ICA/Node.hpp"
#include "ICA/NaturalParameters.hpp"
#include "ICA/Moments.hpp"
// #include "ICA/Message.hpp"

//#include "ICA/piping/Piping.hpp"
#include "ICA/variable/DirichletNode.hpp"
#include "ICA/variable/HiddenNode.hpp"
#include "ICA/variable/DeterministicNode.hpp"
#include "ICA/exponential_models/GaussianModel.hpp"
#include "ICA/exponential_models/GammaModel.hpp"
#include "ICA/exponential_models/DiscreteModel.hpp"
#include <boost/shared_ptr.hpp>
#include <boost/none.hpp>
#include <vector>
#include <boost/assert.hpp> 
#include "rng.hpp"


namespace ICR{
  namespace ICA{
    

    template<template<class> class Model, template<class> class Operation, class T>
    class DeterministicFactor : public FactorNode<double>
    {
    public:
      
      DeterministicFactor(VariableNode<T>* ParentA,  
			   VariableNode<T>* ParentB,  
			   DeterministicNode<Model<T>,T>* Child)
	: m_parentA_node(ParentA),
	  m_parentB_node(ParentB),
	  m_child_node(Child)
      {
	// std::cout<<"DFC A = "  <<m_parentA_node<<std::endl;
	// std::cout<<"DFC A = "  <<ParentA<<std::endl;
	ParentA->AddChildFactor(this);
	ParentB->AddChildFactor(this);

    	Child->SetParentFactor(this);

	Moments<double> a = ParentA->GetMoments();
	Moments<double> b = ParentB->GetMoments();

	//This code shouldn't be here need to refactor it somewhere else...
	
	T prec = 1.0/( Operation<T>::Forward( a[1], b[1]) - Operation<T>::Forward(a[0]*a[0], b[0]*b[0]) );
	T mean = Operation<T>::Forward( a[0], b[0]);
	
	Moments<double> M(mean,mean*mean+1/prec);
	Child->InitialiseMoments(M);

      };
      
      NaturalParameters<T>
      GetNaturalNot(const VariableNode<T>* v) const
      {
	if (v==m_parentA_node)
	  {
	    return Model<T>::template CalcNP2Parent<Operation >(m_parentB_node,m_child_node);
	  }
	else if (v==m_parentB_node)
	  {
	    return Model<T>::template CalcNP2Parent<Operation >(m_parentA_node,m_child_node);
	  }
	else if (v == m_child_node) {
	  // std::cout<<"DF A = "  <<m_parentA_node<<std::endl;

	  return Model<T>::template CalcNP2Deterministic<Operation >(m_parentA_node,m_parentB_node);
	  
	}
	else{
	  throw ("Unknown Node in GetNaturalNot");
	}
      }
      // void 
      // Iterate();
      
      Moments<T>
      InitialiseMoments() const
      {
	return Deterministic<T>::CalcSample(m_prior_node);
      }

      T
      CalcLogNorm() const {return 0;}
    private: 
      //The model
      //GaussianModel m_Model;
      //doublehe variable nodes that are attached to this function node
      VariableNode<T> *m_parentA_node, *m_parentB_node ;
      DeterministicNode<Model<T>,T> *m_child_node;
      
      
      //The NaturalParameters that serve as messages to be passed around
     // NaturalParameters<T> m_NP2Child;
     //  NaturalParameters<T> m_NP2Mean;
     //  NaturalParameters<T> m_NP2Precision;
    private:
      mutable boost::mutex m_mutex;
    };
          
  }
}
