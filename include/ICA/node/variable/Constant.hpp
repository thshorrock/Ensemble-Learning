#pragma once
#include "ICA/message/Moments.hpp"
#include "ICA/node/Node.hpp"
#include "ICA/node/variable/ConstantModels.hpp"

namespace ICR{
  namespace ICA{
    
    template<class Model, class T=double>
    class ConstantNode : public VariableNode<T>
    {
    public:
      
      ConstantNode( const T& value)
	: m_Model( value), m_children()
      {}

      ConstantNode(const size_t& elements, const T& value)
	: m_Model(elements, value), m_children()
      {}
      
      
      void
      ReceiveNaturalParameters(const  NaturalParameters<T>& np) {}
      
      //No parent node possible
      void
      SetParentFactor(FactorNode<T>* f){}

      void
      AddChildFactor(FactorNode<T>* f);
      
      void
      InitialiseMoments(){}

      const Moments<T>&
      GetMoments() const;
      
      
      //No contribution to cost, nothing to update
      void 
      Iterate(Coster& C){};
      
    private:
      Model m_Model;
      std::list<FactorNode<T>*> m_children;
      mutable boost::mutex m_mutex;
    };


  }
}

template<class model,class T>
inline
const ICR::ICA::Moments<T>&
ICR::ICA::ConstantNode<model,T>::GetMoments() const
{
  //Model is thread safe
  return m_Model.GetMoments();
}

template<class model,class T>
inline
void
ICR::ICA::ConstantNode<model,T>::AddChildFactor(FactorNode<T>* f)
{
  //Could be many factors, potentially added with many threads
  boost::lock_guard<boost::mutex> lock(m_mutex);
  m_children.push_back(f);
}
      
