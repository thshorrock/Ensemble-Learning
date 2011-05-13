#pragma once
#include "ICA/Node.hpp"
#include "ICA/Moments.hpp"

namespace ICR{
  namespace ICA{
    
    template<class Model, class T=double>
    class ForwardingNode : public VariableNode<T>
    {
    public:
      
      ForwardingNode(const T& value)
	: m_Model(value), m_children()
      {}
      
      
      void
      ReceiveNaturalParameters(const  NaturalParameters<T>& np) {}
      
      //No parent node possible
      void
      SetParentFactor(FactorNode<T>* f)

      void
      AddChildFactor(FactorNode<T>* f);
      
      Moments<T>
      GetMoments() const;
      
      
      //No contribution to cost, nothing to update
      void 
      Iterate(Coster& C){};
      
    private:
      Model m_Model;
      FactorNode<T>* m_parent;
      std::vector<FactorNode<T>*> m_children;
      mutable boost::mutex m_mutex;
    };


  }
}

template<class model,class T>
inline
ICR::ICA::Moments<T>
ICR::ICA::ForwardingNode<model,T>::GetMoments() const
{
  //Model is thread safe
  return m_Moments;
}

template<class model,class T>
inline
void
ICR::ICA::ForwardingNode<model,T>::SetParentFactor(FactorNode<T>* f)
{
  //Could be many factors, potentially added with many threads
  boost::lock_guard<boost::mutex> lock(m_mutex);
  m_parent = f;
}
      


template<class model,class T>
inline
void
ICR::ICA::ForwardingNode<model,T>::AddChildFactor(FactorNode<T>* f)
{
  //Could be many factors, potentially added with many threads
  boost::lock_guard<boost::mutex> lock(m_mutex);
  m_children.push_back(f);
}
      
