#pragma once
#include "EnsembleLearning/node/Node.hpp"
#include "EnsembleLearning/message/Moments.hpp"
#include "EnsembleLearning/message/NaturalParameters.hpp"
#include "EnsembleLearning/detail/Mutex.hpp"
#include "EnsembleLearning/detail/parallel_algorithms.hpp"

#include <boost/bind.hpp>
#include <vector>

namespace ICR{
  namespace EnsembleLearning{
    
    /** A Deterministic Node.
     *  DeterministicNode's are not inferred,
     *  rather they are staging posts to pass moments that are calculated 
     *  from other (hidden) nodes.
     *  The model used (for example Gaussian/Gamma) is passed by template parameter.
     *  @tparam Model  The model to use for the inferred data.
     *  @tparam T The data type (float or double)
     */
    template <class Model, class T>
    class DeterministicNode : public VariableNode<T>
    {
    public:
      //mostly two, but variable for Discrete model
      /** A constructor.
       *  @param moment_size The number of elements in the stored moments.
       *  This is usually two but varies for discrete nodes.
       */
      DeterministicNode(const size_t moment_size = 2); 

      void
      SetParentFactor(FactorNode<T>* f);
      
      void
      AddChildFactor(FactorNode<T>* f);

      void 
      Iterate(Coster& C);

      const Moments<T>&
      GetMoments() const;

      //should make a const version of this 
      const Moments<T>
      GetForwardedMoments() ;

      /** The number of elements in the stored Moments */
      size_t 
      size() const {return m_Moments.size();}
      
      const std::vector<T>
      GetMean() ;
      
      const std::vector<T>
      GetVariance() ;

      void
      InitialiseMoments()
      {
	m_Moments = m_parent->InitialiseMoments();
      }


      
    private:
      
      FactorNode<T>* m_parent;
      std::vector<FactorNode<T>*> m_children;
      Moments<T> m_Moments;
      mutable Mutex m_mutex;
    };

  }
}


template<class Model,class T>
ICR::EnsembleLearning::DeterministicNode<Model,T>::DeterministicNode(const size_t moment_size) 
  :   m_parent(0), m_children(), m_Moments(moment_size) //, m_ForwardedMoments(moment_size)
{
}


template<class Model,class T>
inline 
void
ICR::EnsembleLearning::DeterministicNode<Model,T>::SetParentFactor(FactorNode<T>* f)
{
  //This should only be called once, so should get no collisions here
  m_parent=f;
  // Can now initialise
  InitialiseMoments();
}
   
template<class Model,class T>
inline
void
ICR::EnsembleLearning::DeterministicNode<Model,T>::AddChildFactor(FactorNode<T>* f)
{ 
  //Could be many factors, potentially added with many threads,
  //therefore make the following critical
#pragma omp critical
  {
    m_children.push_back(f);
  }
}

template<class Model,class T>
inline
const ICR::EnsembleLearning::Moments<T>&
ICR::EnsembleLearning::DeterministicNode<Model,T>::GetMoments() const
{

  /*This value is update in Iterate and called to evaluate other Hidden Nodes
   * (also in iterate mode).  It therefore needs to be protected by a mutex.
   */
  Lock lock(m_mutex);
  return m_Moments;
}
   

template<class Model,class T>
inline
const ICR::EnsembleLearning::Moments<T>
ICR::EnsembleLearning::DeterministicNode<Model,T>::GetForwardedMoments() 
{
  /* No stored value updated here so thread safe.
   */

  //Need to collect this fresh, the moments from other parts 
  //of the graph may have been updated since the last call.

  std::vector<NaturalParameters<T> > vChildrenNP(m_children.size());

  PARALLEL_TRANSFORM(m_children.begin(), m_children.end(), 
  		     vChildrenNP.begin(), 
		     boost::bind(&FactorNode<T>::GetNaturalNot, _1, this)
  		     ); 

  BOOST_ASSERT(vChildrenNP.size() >0);
  NaturalParameters<T> ChildrenNP = PARALLEL_ACCUMULATE(vChildrenNP.begin(), vChildrenNP.end(), NaturalParameters<T>(vChildrenNP[0].size() )); 
  const Moments<T> ForwardedMoments = Model::CalcMoments(ChildrenNP);// ;  //update the moments and the model

  return ForwardedMoments;
}
   
 
  

template<class Model,class T>
inline
const std::vector<T>
ICR::EnsembleLearning::DeterministicNode<Model,T>::GetMean() 
{
  std::vector<T> Mean(1, m_Moments[0]);
  return Mean;
}
   
template<class Model,class T>
inline
const std::vector<T>
ICR::EnsembleLearning::DeterministicNode<Model,T>::GetVariance() 
{
  std::vector<T> Var(1,0.0);
  return Var;
}
   
template<class Model,class T>
void 
ICR::EnsembleLearning::DeterministicNode<Model,T>::Iterate(Coster& C)
{
  //EVALUATE THE COST

  //first get the NP from the parent
  NaturalParameters<T> ParentNP = (m_parent->GetNaturalNot(this));
  //Calcualate the moments
  {
    Lock lock(m_mutex);
    m_Moments = Model::CalcMoments(ParentNP);  
  }
}
