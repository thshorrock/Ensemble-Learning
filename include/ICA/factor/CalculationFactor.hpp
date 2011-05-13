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
#include <boost/assert.hpp> 
#include "rng.hpp"


namespace ICR{
  namespace ICA{
    
    /******************************************************************************
     * CalculationFactor
     ******************************************************************************/
    template<class Calculation, class T = double>
    class CalculationFactor<Calculation> : public FactorNode<double>
    {
    public:
      
      CalculationFactor( VariableNode<double>* First,  VariableNode<double>* Second,  VariableNode<double>* Child)
	: m_first_node(First),
	  m_second_node(Second),
	  m_child_node(Child)
      {
    	First->AddChildFactor(this);
    	Second->AddChildFactor(this);
    	Child->SetParentFactor(this);
      };

      NaturalParameters<double>
      GetNaturalNot(const VariableNode<double>* v) const;
      
    private: 
       VariableNode<double> *m_first_node, *m_second_node, *m_child_node;
    };
    
    
  
    /**************************************************************************************
     **************************************************************************************
     **************************************************************************************
     *
     *  Implementaion
     *
     **************************************************************************************
     **************************************************************************************
     **************************************************************************************/
    template<class Calculation, class T>
    inline
    NaturalParameters<T>
    CalculationFactor<Calculation,T>::GetNaturalNot(const VariableNode<double>* v) const
    {
      if (v==m_first_node)
	{
	  const Moments<double>& second = m_second_node->GetMoments();
	  const Moments<double>& child  = m_child_node->GetMoments();
	  return Calculation::CalcNP2First(second,child);
	}
      else if (v==m_second_node)
	{
	  const Moments<double>& frist  = m_first_node->GetMoments();
	  const Moments<double>& child  = m_child_node->GetMoments();
	  return Calculation::CalcNP2Second(first,child);
	}
      else if (v == m_child_node) 
	{
	  const Moments<double>& frist  = m_first_node->GetMoments();
	  const Moments<double>& second = m_second_node->GetMoments();
	  return Calculation::CalcNP2Data(first,second);
	}
      else{
	throw ("Unknown Node in GetNaturalNot");
      }
    }
  }
}
