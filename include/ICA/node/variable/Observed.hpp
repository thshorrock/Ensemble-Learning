#pragma once
#include "ICA/message/Moments.hpp"
#include "ICA/node/Node.hpp"
#include "ICA/node/variable/ConstantModels.hpp"

namespace ICR{
  namespace ICA{
    
    template<class Model, class T=double>
    class ObservedNode : public VariableNode<T>
    {
    public:
      
      ObservedNode( const T& value)
	: m_Moments(2, make_Moments(value), Model() ), 
	  m_parent(0),
	  m_children()
      {}

      ObservedNode(const size_t& elements, const T& value)
	: m_Model(make_Moments(elements,value)), 
	  m_parent(0), 
	  m_children()
      {}
      
      
      //No parent node possible
      void
      SetParentFactor(FactorNode<T>* f);

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
      //default
      Moments<T>
      make_Moments(const size_t s,const T& d, const Gaussian)
      {
	return Moments<T>(d, d*d);
      }
      Moments<T>
      make_Moments(const size_t s,const T& d, const RectifiedGaussian)
      {
	return Moments<T>(d, d*d);
      }
      Moments<T>
      make_Moments(const size_t s,const T& d, const Gamma)
      {
	return Moments<T>(d, std::log(d));
      }
      Moments<T>
      make_Moments(const size_t s,const T& d, const Dirichlet)
      {
	return Moments<T>(std::vector<T>(s,d));
      }
      const Moments<T> m_Moments;
      FactorNode<T>* m_parent;
      std::list<FactorNode<T>*> m_children;
    };
    
  }
}

template<class model,class T>
inline
const ICR::ICA::Moments<T>&
ICR::ICA::ObservedNode<model,T>::GetMoments() const
{
  //Model is thread safe
  return m_Model.GetMoments();
}

template<class model,class T>
inline
void
ICR::ICA::ObservedNode<model,T>::AddChildFactor(FactorNode<T>* f)
{
  //Could be many factors, potentially added with many threads
  boost::lock_guard<boost::mutex> lock(m_mutex);
  m_children.push_back(f);
}
    

template<class Model, class T>
inline
void
ICR::ICA::ObvservedNode<Model,T>::SetParentFactor(FactorNode<T>* f)
{
  //Assuming only one factor: this should be called once per iteration 
  //  and therefore be thread safe.  
  //To be conservative (but very wastefull) uncomment the following
  //boost::lock_guard<boost::mutex> lock(m_mutex);
  m_parent = f;
}
  

template<class Model, class T>
inline
void
ICR::ICA::DataNode<Model,T>::Iterate(Coster& Total)
{
  //To be conservative (but very wastefull) uncomment the following
  //boost::lock_guard<boost::mutex> lock(m_mutex);

  //Get the natural parameters from the parent factor for this node.
  // (Node, the Parent Factor must have been called before this node)
  // (this is easily guareenteed by calling ALL the factor nodes prior to ANY of the variableNodes.
  // (This is done by the Builder
  if (m_parent!=0)
    {
      NaturalParameters<T> ParentNP =m_parent->GetNaturalNot(this); 
      //See page 41 of Winn's thesis for this formula
      const T Cost = ParentNP*GetMoments() +  m_parent->CalcLogNorm();// + Model::GetDataPenalty(GetMoments()[0]);
      Total += Cost;
    }
  // std::cout<<"DataCost = "<<Cost<<std::endl;
  //  std::cout<<"ParentNP = "<<ParentNP<<std::endl;
  //  std::cout<<"Moments = "<<GetMoments()<<std::endl;
  //  std::cout<<"ParentNorm = "<<m_parent->CalcLogNorm()<<std::endl;


}
