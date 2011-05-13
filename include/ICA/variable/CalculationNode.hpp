


namespace ICR{
  namespace ICA{

    template <class Model, class T = double>
    class ForwardingNode : public VariableNode<T>
    {
    public:
      
      CalculationNode(const size_t moment_size = 2); //mostly two, but variable for Discrete model

      void
      SetParentFactor(FactorNode<T>* f);
      
      void
      AddChildFactor(FactorNode<T>* f);

      void 
      Iterate(Coster& C);

      Moments<T>
      GetMoments() const;

      size_t 
      size() const {return m_Moments.size();}
      
      
    private:
      
      //Model<T> m_Model; 
      FactorNode<T>* m_parent;
      std::vector<FactorNode<T>*> m_children;
      Moments<T> m_Moments;
      mutable boost::mutex m_mutex;
    };

    
  }
}




template<class Model,class T>
ICR::ICA::CalculationNode<Model,T>::CalculationNode(const size_t moment_size) 
  :   m_parent(0), m_children(), m_Moments(moment_size)
{}


template<class Model,class T>
inline 
void
ICR::ICA::CalculationNode<Model,T>::SetParentFactor(FactorNode<T>* f)
{
  //This should only be called once, so should get no collisions here
  m_parent=f;
}
   
template<class Model,class T>
inline
void
ICR::ICA::CalculationNode<Model,T>::AddChildFactor(FactorNode<T>* f)
{ 
  //This could be called by different threads, so lock
  boost::lock_guard<boost::mutex> lock(m_mutex);
  m_children.push_back(f);
}

template<class Model,class T>
inline
ICR::ICA::Moments<T>
ICR::ICA::CalculationNode<Model,T>::GetMoments() const
{
  /*This value is update once in update stage,
   * and called in a Factor update stage.
   * These never overlap and so this is thead safe
   */
  return m_Moments;
}
   

   
template<class Model,class T>
void 
ICR::ICA::CalculationNode<Model,T>::Iterate(Coster& C)
{
  //EVALUATE THE COST

  //first get the Natural parameters from all the children.
  //This assumes that the Children (factors) have already run.
  //This is guarenteed by the Builder class.
  NaturalParameters<T> NP = (m_parent->GetNaturalNot(this));
  //Add it to the total

  
  
  std::vector<NaturalParameters<T> > ChildrenNP(m_children.size());
  //Get the NPs and put them in the vector
  PARALLEL_TRANSFORM(m_children.begin(), m_children.end(), 
		     ChildrenNP.begin(), boost::bind(&FactorNode<T>::GetNaturalNot, _1, this)
		     ); //from children

  //Add them up
  NaturalParameters<T> NP = PARALLEL_ACCUMULATE(ChildrenNP.begin(), ChildrenNP.end(), NP ); 
  //seems to be a bug here, do this one serially
  //NaturalParameters<T> NP = std::accumulate(ChildrenNP.begin(), ChildrenNP.end(), init ); 

  //Get the moments and update the model
  Model m_Model;
  m_Moments = m_Model.CalculateMoments(NP);  //update the moments and the model
  


}
