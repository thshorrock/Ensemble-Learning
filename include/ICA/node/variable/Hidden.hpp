#pragma once
#include "ICA/node/Node.hpp"
#include "ICA/message/Moments.hpp"
#include "ICA/message/NaturalParameters.hpp"
#include "ICA/detail/parallel_algorithms.hpp"
#include "ICA/detail/Mutex.hpp"

#include <boost/bind.hpp>
#include <vector>



namespace ICR{
  namespace ICA{
    
    /** A Hidden Variable Node.
     *  HiddenNode's store moments that are inferred by Ensemble Learning.
     *  The model used (for example Gaussian/Gamma) is passed by template parameter.
     *  @tparam Model  The model to use for the inferred data.
     *  @tparam T The data type (float or double)
     */
    template <class Model, class T>
    class HiddenNode : public VariableNode<T>
    {
    public:
      /** A constructor.
       *  @param moment_size The number of elements in the stored moments.
       *  This is usually two but varies for discrete nodes.
       */
      HiddenNode(const size_t moment_size = 2); //mostly two, but variable for Discrete model

      void
      SetParentFactor(FactorNode<T>* f);
      
      void
      AddChildFactor(FactorNode<T>* f);

      void 
      Iterate(Coster& C);

      void
      InitialiseMoments()
      {
	m_Moments = m_parent->InitialiseMoments();
      }


      const Moments<T>&
      GetMoments() const;

      const std::vector<T>
      GetMean() ;
      
      const std::vector<T>
      GetVariance() ;

      /** The number of elements in the stored Moments */
      size_t 
      size() const {return m_Moments.size();}
      
    private:
      
      const NaturalParameters<T>
      GetNP();

      FactorNode<T>* m_parent;
      std::vector<FactorNode<T>*> m_children;
      Moments<T> m_Moments;
      mutable Mutex m_mutex;
    };

  }
}


template<class Model,class T>
ICR::ICA::HiddenNode<Model,T>::HiddenNode(const size_t moment_size) 
  :   m_parent(0), m_children(), m_Moments(moment_size)
{}


template<class Model,class T>
inline 
void
ICR::ICA::HiddenNode<Model,T>::SetParentFactor(FactorNode<T>* f)
{
  //This should only be called once, so should get no collisions here
  m_parent=f;
  // Can now initialise
  InitialiseMoments();
  
}



   
template<class Model,class T>
inline
void
ICR::ICA::HiddenNode<Model,T>::AddChildFactor(FactorNode<T>* f)
{ 
  //This could be called simultaneously by different threads, so call it critical
#pragma omp critical
  {
    m_children.push_back(f);
  }
}

template<class Model,class T>
inline
const ICR::ICA::Moments<T>&
ICR::ICA::HiddenNode<Model,T>::GetMoments() const
{
  /*This value is update in Iterate and called to evaluate other Hidden Nodes
   * (also in iterate mode).  It therefore needs to be protected by a mutex.
   */
  Lock lock(m_mutex);
  return m_Moments;
}
   
template<class Model,class T>
inline
const ICR::ICA::NaturalParameters<T>
ICR::ICA::HiddenNode<Model,T>::GetNP()
{
  //EVALUATE THE COST

  //first get the NP from the parent
  const NaturalParameters<T> ParentNP = (m_parent->GetNaturalNot(this));
  //The initial value in the total
  NaturalParameters<T> NP = ParentNP; //total

  //get the Natural parameters from all the children.
  std::vector<NaturalParameters<T> > ChildrenNP(m_children.size());
  //and put them in the vector
  PARALLEL_TRANSFORM(m_children.begin(), m_children.end(), 
		     ChildrenNP.begin(), boost::bind(&FactorNode<T>::GetNaturalNot, _1, this)
		     ); //from children

  //Add them up
  //(Note: Bug 48750 in gcc in parallel accumulate causes crash at next line: 
  //  fixed in 4.6.1)
  NP = PARALLEL_ACCUMULATE(ChildrenNP.begin(), ChildrenNP.end(), NP ); 
  return NP;
}

template<class Model,class T>
inline
const std::vector<T>
ICR::ICA::HiddenNode<Model,T>::GetMean() 
{
  return Model::CalcMean(GetNP());
}
   
template<class Model,class T>
inline
const std::vector<T>
ICR::ICA::HiddenNode<Model,T>::GetVariance() 
{
  std::vector<T> prec = Model::CalcPrecision(GetNP());
  std::vector<T> var(prec.size());
  for(size_t i=0;i<prec.size();++i){
    var[i] = 1.0/prec[i];
  }
  return var;
}

template<class Model,class T>
inline
void 
ICR::ICA::HiddenNode<Model,T>::Iterate(Coster& C)
{
  const NaturalParameters<T> NP = GetNP();
  //Get the moments and update the model
  const T LogNorm = Model::CalcLogNorm(NP);
  {
    Lock lock(m_mutex);
    m_Moments = Model::CalcMoments(NP);  //update the moments and the model
  }
  //first get the NP from the parent
  const NaturalParameters<T> ParentNP = (m_parent->GetNaturalNot(this));
  C +=  (ParentNP - NP)*m_Moments +m_parent->CalcLogNorm() -  LogNorm;


}
