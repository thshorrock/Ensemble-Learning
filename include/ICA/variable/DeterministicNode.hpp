#pragma once
#include "ICA/Node.hpp"
#include "ICA/Moments.hpp"
#include "ICA/NaturalParameters.hpp"
#include <boost/thread/locks.hpp>
#include <numeric>
#include "ICA/parallel_algorithms.hpp"

namespace ICR{
  namespace ICA{
    
    template <class Model, class T = double>
    class DeterministicNode : public VariableNode<T>
    {
    public:
      
      DeterministicNode(const size_t moment_size = 2); //mostly two, but variable for Discrete model

      void
      SetParentFactor(FactorNode<T>* f);
      
      void
      AddChildFactor(FactorNode<T>* f);

      void 
      Iterate(Coster& C);

      const Moments<T>&
      GetMoments() const;

      const Moments<T>
      GetForwardedMoments() const;

      size_t 
      size() const {return m_Moments.size();}
      
      void
      InitialiseMoments()
      {
	m_Moments = m_parent->InitialiseMoments();
      }


      
    private:
      
      //Model<T> m_Model; 
      FactorNode<T>* m_parent;
      std::vector<FactorNode<T>*> m_children;
      Moments<T> m_Moments;//, m_ForwardedMoments;
      mutable boost::mutex m_mutex;
    };

  }
}


template<class Model,class T>
ICR::ICA::DeterministicNode<Model,T>::DeterministicNode(const size_t moment_size) 
  :   m_parent(0), m_children(), m_Moments(moment_size) //, m_ForwardedMoments(moment_size)
{
}


template<class Model,class T>
inline 
void
ICR::ICA::DeterministicNode<Model,T>::SetParentFactor(FactorNode<T>* f)
{
  //This should only be called once, so should get no collisions here
  m_parent=f;
  // Can now initialise
  InitialiseMoments();
}
   
template<class Model,class T>
inline
void
ICR::ICA::DeterministicNode<Model,T>::AddChildFactor(FactorNode<T>* f)
{ 
  //This could be called by different threads, so lock
  boost::lock_guard<boost::mutex> lock(m_mutex);
  m_children.push_back(f);
}

template<class Model,class T>
inline
const ICR::ICA::Moments<T>&
ICR::ICA::DeterministicNode<Model,T>::GetMoments() const
{
  /*This value is update once in update stage,
   * and called in a Factor update stage.
   * These never overlap and so this is thead safe
   */


  return m_Moments;
}
   

template<class Model,class T>
inline
const ICR::ICA::Moments<T>
ICR::ICA::DeterministicNode<Model,T>::GetForwardedMoments() const
{
  /*This value is update once in update stage,
   * and called in a Factor update stage.
   * These never overlap and so this is thead safe
   */

  //Need to collect this fresh!

  std::vector<NaturalParameters<T> > vChildrenNP(m_children.size());
  //Get the NPs and put them in the vector
  PARALLEL_TRANSFORM(m_children.begin(), m_children.end(), 
		     vChildrenNP.begin(), boost::bind(&FactorNode<T>::GetNaturalNot, _1, this)
		     ); 

  BOOST_ASSERT(vChildrenNP.size() >0);
  NaturalParameters<T> ChildrenNP = PARALLEL_ACCUMULATE(vChildrenNP.begin(), vChildrenNP.end(), NaturalParameters<T>(vChildrenNP[0].size() )); 
  Moments<T> ForwardedMoments = Model::CalcMoments(ChildrenNP);// ;  //update the moments and the model
  //std::cout<<"chilren = "<<vChildrenNP.size()<<" FMoment = "<<ForwardedMoments<<std::endl;

  return ForwardedMoments;
}
   

   
template<class Model,class T>
void 
ICR::ICA::DeterministicNode<Model,T>::Iterate(Coster& C)
{
  //EVALUATE THE COST

  //first get the Natural parameters from all the children.
  //This assumes that the Children (factors) have already run.
  //This is guarenteed by the Builder class.
  NaturalParameters<T> ParentNP = (m_parent->GetNaturalNot(this));
  //Add it to the total
  
  // std::vector<NaturalParameters<T> > vChildrenNP(m_children.size());
  // //Get the NPs and put them in the vector
  // PARALLEL_TRANSFORM(m_children.begin(), m_children.end(), 
  // 		     vChildrenNP.begin(), boost::bind(&FactorNode<T>::GetNaturalNot, _1, this)
  // 		     ); //from children

  // //Add them up
  // NaturalParameters<T> ChildrenNP = PARALLEL_ACCUMULATE(vChildrenNP.begin(), vChildrenNP.end(), NaturalParameters<T>(ParentNP.size() )); 
  // //seems to be a bug here, do this one serially
  // //NaturalParameters<T> NP = std::accumulate(ChildrenNP.begin(), ChildrenNP.end(), init ); //

  //Get the moments and update the model
  m_Moments = Model::CalcMoments(ParentNP);  //update the moments and the model
  //m_ForwardedMoments = Model<t>::CalcMoments(ChildrenNP);// ;  //update the moments and the model
  
  //std::cout<<"m_Moments = "<<m_Moments<<std::endl;
  // std::cout<<"m_ForwdM  = "<<m_ForwardedMoments<<std::endl;



}
