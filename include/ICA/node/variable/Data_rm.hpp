#pragma once
#include "ICA/node/Node.hpp"
#include "ICA/message/Moments.hpp"
#include "ICA/node/factor/Factor.hpp"
// #include "ICA/Message.hpp"

#include <set>

namespace ICR{
  namespace ICA{
    
    /** A Data Node.
     *
     *  @tpram Model The model.
     */
    template<class Model, class T = double>
    class DataNode : public VariableNode<T>
    {
    public:
      /** Constructor.
       *
       * @Param value: The value of the data.
       */
      DataNode(const T& value);
      

      void
      InitialiseMoments(){}

      /** Get the moments from the data.*/
      const Moments<T>&
      GetMoments() const;

      /** Set the Parent Factor.*/
      void
      SetParentFactor(FactorNode<T>* f);

      //Cannot add child to data.
      void
      AddChildFactor(FactorNode<T>* f){}

      /**Iterate this node and calculate Cost.
       * @param Coster cost.  An object that calculates the global cost
       *  in a thread safe way.
       */
      void 
      Iterate(Coster&);


      private:
      Model m_Model;
      FactorNode<T>* m_parent;

      //Uncomment if want mutex (shouldn't need it).
      //mutable boost::mutex mutex;
    };


  }
}

template<class Model, class T>
inline
ICR::ICA::DataNode<Model,T>::DataNode(const T& value)
  : m_Model(value) , m_parent(0)
{}

template<class Model, class T>
inline
const ICR::ICA::Moments<T>&
ICR::ICA::DataNode<Model,T>::GetMoments() const
{
  /* Each Data Node will be connected to only 1 factor node.
   * Therefore there is no need for a mutex here.
   * (this function will be called only once per iteration)
   */
  //To be conservative (but very wastefull) uncomment the following
  //boost::lock_guard<boost::mutex> lock(m_mutex);
  return m_Model.GetMoments();
}

template<class Model, class T>
inline
void
ICR::ICA::DataNode<Model,T>::SetParentFactor(FactorNode<T>* f)
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
  NaturalParameters<T> ParentNP =m_parent->GetNaturalNot(this); 
  //See page 41 of Winn's thesis for this formula
  double Cost = ParentNP*GetMoments() +  m_parent->CalcLogNorm();// + Model::GetDataPenalty(GetMoments()[0]);
  // std::cout<<"DataCost = "<<Cost<<std::endl;
  //  std::cout<<"ParentNP = "<<ParentNP<<std::endl;
  //  std::cout<<"Moments = "<<GetMoments()<<std::endl;
  //  std::cout<<"ParentNorm = "<<m_parent->CalcLogNorm()<<std::endl;


  Total += Cost;
}
 
