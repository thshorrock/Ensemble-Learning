#pragma once
#include "ICA/Node.hpp"
#include "ICA/Moments.hpp"
#include "ICA/NaturalParameters.hpp"
#include <boost/thread/locks.hpp>
#include <numeric>
#include "parallel_algorithms.hpp"

namespace ICR{
  namespace ICA{
    
    template <class Model, class T = double>
    class MixtureNode : public VariableNode<T>, 
		       public LearningNode<T>,
		       public CanMakeMixture<T>
    {
    public:
      
      MixtureNode(const size_t moment_size = 2); //mostly two, but variable for Discrete model

      void
      SetParentFactor(FactorNode<T>* f);
      
      void
      AddChildFactor(FactorNode<T>* f);

      void 
      Iterate(Coster& C);

      Moments<T>
      GetMoments() const;

      
      /** The Av Log of exponential distribution.
       *
       * This is needed when calculating the log proability distribution
       * in a Mixture model (passed to Discrete Node).
       * (See page 46 of Winn's thesis).
       */
      double
      GetAvLog() const;
      
      size_t 
      size() const {return m_Moments.size();}
      
    private:
      
      //Model<T> m_Model; 
      FactorNode<T>* m_parent;
      std::vector<FactorNode<T>*> m_children;
      Moments<T> m_Moments;
      double m_AvLog;
      mutable boost::mutex m_mutex;
    };

    template <class T>
    inline
    std::ostream&
    operator<<( std::ostream& out, const NaturalParameters<T>& m)
    {
      for(size_t i=0;i<m.size();++i)
	{
	  out<<m[i]<<" ";
	}
      return out;
    }
  }
}


template<class Model,class T>
ICR::ICA::MixtureNode<Model,T>::MixtureNode(const size_t moment_size) 
  :   m_parent(), m_children(), m_Moments(moment_size), m_AvLog()
{}


template<class Model,class T>
inline 
void
ICR::ICA::MixtureNode<Model,T>::SetParentFactor(FactorNode<T>* f)
{
  //This should only be called once, so should get no collisions here
  m_parent=f;
}
   
template<class Model,class T>
inline
void
ICR::ICA::MixtureNode<Model,T>::AddChildFactor(FactorNode<T>* f)
{ 
  //This could be called by different threads, so lock
  boost::lock_guard<boost::mutex> lock(m_mutex);
  m_children.push_back(f);
}

template<class Model,class T>
inline
ICR::ICA::Moments<T>
ICR::ICA::MixtureNode<Model,T>::GetMoments() const
{
  /*This value is update once in update stage,
   * and called in a Factor update stage.
   * These never overlap and so this is thead safe
   */
  return m_Moments;
}
   
template<class Model,class T> 
inline  
double
ICR::ICA::MixtureNode<Model,T>::GetAvLog() const 
{
  /*This value is update once in update stage,
   * and called in a Factor update stage.
   * These never overlap and so this is thead safe
   */
  return m_AvLog;
}
      
// T
// ICR::ICA::MixtureNode<Model,T>::GetNormalisation() const 
// {
  
//   return m_Model.GetNormalisation();
// }



   
template<class Model,class T>
void 
ICR::ICA::MixtureNode<Model,T>::Iterate(Coster& C)
{
  //EVALUATE THE COST

  //first get the Natural parameters from all the children.
  //This assumes that the Children (factors) have already run.
  //This is guarenteed by the Builder class.

  std::vector<NaturalParameters<T> > ChildrenNP(m_children.size());
  //Get the NPs and put them in the vector
  PARALLEL_TRANSFORM(m_children.begin(), m_children.end(), 
		     ChildrenNP.begin(), boost::bind(&FactorNode<T>::GetNaturalNot, _1, this)
		     ); //from children

  //Add them up
  NaturalParameters<T> init(0.0,0.0);
  NaturalParameters<T> NP = PARALLEL_ACCUMULATE(ChildrenNP.begin(), ChildrenNP.end(), init ); 
  //seems to be a bug here, do this one serially
  //NaturalParameters<T> NP = std::accumulate(ChildrenNP.begin(), ChildrenNP.end(), init ); 

  // typename std::vector<NaturalParameters<T> >::iterator it = ChildrenNP.begin();
  // while ( it!= ChildrenNP.end()){
  //   init = init + *it++;  // or: init=binary_op(init,*first++) for the binary_op version
    
  // }
  // NaturalParameters<T> NP = init;

  //std::cout<<"init = "<<init<<std::endl;
  //NP from the parent factor
  
  //Do the same for the parents
  
  std::vector<NaturalParameters<T> > ParentsNP(m_parents.size());
  //Get the NPs and put them in the vector
  PARALLEL_TRANSFORM(m_parents.begin(), m_parents.end(), 
		     ParentsNP.begin(), boost::bind(&FactorNode<T>::GetNaturalNot, _1, this)
		     ); //from parents

  //Add them up to the total
  // NaturalParameters<T> initP(0.0,0.0);
  //  NaturalParameters<T> 
  //NP = PARALLEL_ACCUMULATE(ParentsNP.begin(), ParentsNP.end(),NP ); 


  NaturalParameters<T> ParentNP = (m_parent->GetNaturalNot(this));
  //Add it to the total
  NP += ParentNP; //total

  //Get the moments and update the model
  Model m_Model;
  m_Moments = m_Model.CalculateMoments(NP);  //update the moments and the model
  
  //The common part to AvLog and Cost
  double val1 = ParentNP*m_Moments+ m_parent->GetNormalisation();  

  //The AvLog
  m_AvLog = val1 +  m_parent->GetFFunction();
  //The cost -- see page 40 (eqn 2.34 of Winn's thesis)
  C += val1 - NP*m_Moments -  m_Model.GetNormalisation();


}
