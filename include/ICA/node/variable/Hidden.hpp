#pragma once
#include "ICA/node/Node.hpp"
#include "ICA/message/Moments.hpp"
#include "ICA/message/NaturalParameters.hpp"
#include <boost/thread/locks.hpp>
#include "ICA/detail/parallel_algorithms.hpp"

namespace ICR{
  namespace ICA{
    
    template <class Model, class T = double>
    class HiddenNode : public VariableNode<T>// , 
		       // public LearningNode<T>// ,
		       // public CanMakeMixture<T>
    {
    public:
      
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

      
      /** The Av Log of exponential distribution.
       *
       * This is needed when calculating the log proability distribution
       * in a Mixture model (passed to Discrete Node).
       * (See page 46 of Winn's thesis).
       */
      // double
      // GetAvLog() const;
      
      size_t 
      size() const {return m_Moments.size();}
      
    private:
      
      //Model<T> m_Model; 
      FactorNode<T>* m_parent;
      std::vector<FactorNode<T>*> m_children;
      Moments<T> m_Moments;
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
  //This could be called by different threads, so lock
  boost::lock_guard<boost::mutex> lock(m_mutex);
  m_children.push_back(f);
}

template<class Model,class T>
inline
const ICR::ICA::Moments<T>&
ICR::ICA::HiddenNode<Model,T>::GetMoments() const
{
  /*This value is update once in update stage,
   * and called in a Factor update stage.
   * These never overlap and so this is thead safe
   */
  //std::cout<<"Moments from "<<this<<" = "<<m_Moments<<std::endl;

  return m_Moments;
}
   
// template<class Model,class T> 
// inline  
// double
// ICR::ICA::HiddenNode<Model,T>::GetAvLog() const 
// {
//   /*This value is update once in update stage,
//    * and called in a Factor update stage.
//    * These never overlap and so this is thead safe
//    */
//   return m_AvLog;
// }
      
// T
// ICR::ICA::HiddenNode<Model,T>::CalcLogNorm() const 
// {
  
//   return m_Model.CalcLogNorm();
// }



   
template<class Model,class T>
inline
void 
ICR::ICA::HiddenNode<Model,T>::Iterate(Coster& C)
{
  //EVALUATE THE COST

  //first get the Natural parameters from all the children.
  //This assumes that the Children (factors) have already run.
  //This is guarenteed by the Builder class.

  // std::cout<<"THIS = "<<this<<std::endl;
  NaturalParameters<T> ParentNP = (m_parent->GetNaturalNot(this));
  //Add it to the total
  NaturalParameters<T> NP = ParentNP; //total

  // std::cout<<"P NP"<<NP<<std::endl;


  std::vector<NaturalParameters<T> > ChildrenNP(m_children.size());
  //  NaturalParameters<T> TotalChildrenNP(ParentNP.size())
  //Get the NPs and put them in the vector
  PARALLEL_TRANSFORM(m_children.begin(), m_children.end(), 
		     ChildrenNP.begin(), boost::bind(&FactorNode<T>::GetNaturalNot, _1, this)
		     ); //from children

  // //Add them up
  // //std::cout<<"MM size = "<<m_Moments.size()<<std::endl;
  // NaturalParameters<T> init(m_Moments.size());
  // //std::cout<<"init size= "<<init.size()<<std::endl;

  NP = PARALLEL_ACCUMULATE(ChildrenNP.begin(), ChildrenNP.end(), NP ); 
  //seems to be a bug here, do this one serially
  //NaturalParameters<T> NP = std::accumulate(ChildrenNP.begin(), ChildrenNP.end(), init ); 
  // std::cout<<"total NP"<<NP<<std::endl;


  // typename std::vector<NaturalParameters<T> >::iterator it = ChildrenNP.begin();
  // while ( it!= ChildrenNP.end()){
  //   init = init + *it++;  // or: init=binary_op(init,*first++) for the binary_op version
    
  // }
  // NaturalParameters<T> NP = init;

  //std::cout<<"init = "<<init<<std::endl;
  //NP from the parent factor
  

  //Get the moments and update the model
  //std::cout<<"NP size = "<<NP.size()<<std::endl;

  T LogNorm = Model::CalcLogNorm(NP);

  m_Moments = Model::CalcMoments(NP);  //update the moments and the model

    // std::cout<<"Cost Contr = "<<(ParentNP - NP)*m_Moments +m_parent->CalcLogNorm() -  LogNorm<<std::endl;


  C +=  (ParentNP - NP)*m_Moments +m_parent->CalcLogNorm() -  LogNorm;


   // std::cout<<"Parent Norm = "<<m_parent->CalcLogNorm()<<std::endl;
   // std::cout<<"Model Norm  = "<<LogNorm<<std::endl;
   // std::cout<<"m_Moments  = "<<m_Moments<<std::endl;
   // std::cout<<"Parent NP = "<<ParentNP<<std::endl;
   // std::cout<<"NP  = "<<NP<<std::endl;
  


  //printf("NODE: %p \t",  this);
  //std::cout<<"Moments: "<<this<<std::endl;


  //The common part to AvLog and Cost
  // double val1 = ParentNP*m_Moments+ m_parent->CalcLogNorm();  
  

  //The AvLog
  // m_AvLog = val1 +  m_parent->GetFFunction();
  //The cost -- see page 40 (eqn 2.34 of Winn's thesis)
  // std::cout<<"Cost Contr = "<<val1 - NP*m_Moments -  Model::CalcLogNorm()<<std::endl;

  

}
