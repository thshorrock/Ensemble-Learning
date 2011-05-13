#pragma once
#include "ICA/Node.hpp"
#include "ICA/Moments.hpp"
#include "ICA/exponential_models/DirichletModel.hpp"

namespace ICR{
  namespace ICA{
    
    template<class T=double>
    class DirichletNode : public VariableNode<T>
    {
    public:
      
      DirichletNode(const size_t elements);
      
      Moments<T>
      GetMoments() const;

      void
      SetParentFactor(FactorNode<T>* f){}

      void
      AddChildFactor(FactorNode<T>* f);
      
      
      void
      InitialiseMoments(const Moments<T>& m)
      {
      }

      //No contribution to cost, nothing to update
      void 
      Iterate(Coster& C){};
      
    private:
      const Moments<T> m_moments;
      std::list<FactorNode<T>*> m_children;
      mutable boost::mutex m_mutex;
    };


  }
}

namespace {
  template<class T>
  ICR::ICA::Moments<T> make_moments(const size_t elements) 
  {

    // ICR::ICA::DirichletModel<T> model;
    // std::vector<T> values(elements,1.0);
    // model.SetValues(values);
    
    //    T v = std::log(1.0/(T(elements)));
    std::vector<T> the_moments(elements, 1.0);
    return ICR::ICA::Moments<T>(the_moments);
  }
  
}

template<class T>
inline
ICR::ICA::DirichletNode<T>::DirichletNode(const size_t elements)
  :  m_moments( make_moments<T>(elements) ),  m_children()
{
}
  
template<class T>
inline    
ICR::ICA::Moments<T>
ICR::ICA::DirichletNode<T>::GetMoments() const
{
  //m_moments is const, threadsafe
  return m_moments;
}
      
template<class T>
inline
void
ICR::ICA::DirichletNode<T>::AddChildFactor(FactorNode<T>* f)
{
  //Could be many factors, potentially added with many threads
  boost::lock_guard<boost::mutex> lock(m_mutex);
  m_children.push_back(f);
}
      
