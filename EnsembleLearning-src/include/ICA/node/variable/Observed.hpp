#pragma once
#include "ICA/message/Moments.hpp"
#include "ICA/node/Node.hpp"
#include "ICA/exponential_model/Gaussian.hpp"
#include "ICA/exponential_model/RectifiedGaussian.hpp"
#include "ICA/exponential_model/Gamma.hpp"
#include "ICA/exponential_model/Discrete.hpp"
#include "ICA/exponential_model/Dirichlet.hpp"

namespace ICR{
  namespace ICA{
    
    /** An observed node.
     *  A node that contains data that is a known constant.
     *  Examples include data nodes (where the data is experimentally determined)
     *   and prior knowledge - the mean of the data is 4.0, say).
     *  @tparam Model The model catagory into which the data falls.
     *    For example, Data might be modelled by a Gaussian distribution.
     *  @tparam T The data type used in calculations - either float or double.
     */
    template<class Model, class T>
    class ObservedNode : public VariableNode<T>
    {
    public:
      /** A Constructor.
       *  @param value The observed value of the node.
       */
      ObservedNode( const T& value)
	: m_Moments(make_Moments(2, value, Model() ) ), 
	  m_parent(0),
	  m_children()
      {}

      /** A Constructor.
       * @param  elements The number of elements in the observed node.
       * @param  value The value of each of the elements 
       */
      ObservedNode(const size_t& elements, const T& value)
	: m_Moments(make_Moments(elements,value, Model() ) ), 
	  m_parent(0), 
	  m_children()
      {}
      
      
      void
      SetParentFactor(FactorNode<T>* f);

      void
      AddChildFactor(FactorNode<T>* f);
      
      void
      InitialiseMoments(){};

      const Moments<T>&
      GetMoments() const;
      
      const std::vector<T>
      GetMean() ;
      
      const std::vector<T>
      GetVariance() ;

      void 
      Iterate(Coster& C);
      
    private:
      Moments<T>
      make_Moments(const size_t s,const T& d, const Gaussian<T> )
      {
	return Moments<T>(d, d*d);
      }
      Moments<T>
      make_Moments(const size_t s,const T& d, const RectifiedGaussian<T>)
      {
	return Moments<T>(d, d*d);
      }
      Moments<T>
      make_Moments(const size_t s,const T& d, const Gamma<T>)
      {
	return Moments<T>(d, std::log(d));
      }
      Moments<T>
      make_Moments(const size_t s,const T& d, const Dirichlet<T>)
      {
	return Moments<T>(std::vector<T>(s,d));
      }
      const Moments<T> m_Moments;
      FactorNode<T>* m_parent;
      std::vector<FactorNode<T>*> m_children;
    };
    
  }
}

template<class model,class T>
inline
const ICR::ICA::Moments<T>&
ICR::ICA::ObservedNode<model,T>::GetMoments() const
{
  //Obvserved moments are not modified and so this is thead safe.
  return m_Moments;
}

template<class model,class T>
inline
void
ICR::ICA::ObservedNode<model,T>::AddChildFactor(FactorNode<T>* f)
{
  //Could be many factors, potentially added with many threads,
  //therefore make the following critical
#pragma omp critical
  {
    m_children.push_back(f);
  }
}


template<class Model, class T>
inline
void
ICR::ICA::ObservedNode<Model,T>::SetParentFactor(FactorNode<T>* f)
{
  //Only one factor: this should be called once (therefore thread safe)
  m_parent = f;
}
  

template<class Model,class T>
inline
const std::vector<T>
ICR::ICA::ObservedNode<Model,T>::GetMean() 
{
  std::vector<T> Mean(1, m_Moments[0]);
  return Mean;
}
   
template<class Model,class T>
inline
const std::vector<T>
ICR::ICA::ObservedNode<Model,T>::GetVariance() 
{
  std::vector<T> Var(1,0.0);
  return Var;
}

template<class Model, class T>
inline
void
ICR::ICA::ObservedNode<Model,T>::Iterate(Coster& Total)
{
  //Constant nodes have no parents (and contribute nothing to the cost)
  if (m_parent!=0)
    {
      //This is a data node...
      //Assume thead-safety of other nodes so do not need to worry here.
      const NaturalParameters<T> ParentNP =m_parent->GetNaturalNot(this); 
      //See page 41 of Winn's thesis for this formula
      const T Cost = ParentNP*GetMoments() +  m_parent->CalcLogNorm();
      Total += Cost;
    }
}
