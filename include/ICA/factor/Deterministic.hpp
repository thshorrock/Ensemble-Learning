
namespace ICR{
  namespace ICA{

    <template class T>
    class Expression
    {

    public:
      Expression(std::vector<VariableNode*> parents,
		 VariableNode* child)
      
      
    }
    
    template <class T>
    class Add : public FactorNode<T>
    {
      
    public:
      Add(GaussianNode G1, GaussianNode G2, ForwardingNode child);
      
      

    }

    template<class Expression, class T = double>
    class Determinisic : public FactorNode<T>
    {
    public:
      typedef HiddenNode<Model> ParentNode;
      typedef HiddenNode<DiscreteModel<T> > WeightsNode;
      typedef VariableNode<T> ChildNode;
      
      Determinisic(const std::vector<ParentNode*> parents,
	      WeightsNode* weights, 
	      ChildNode* child);
      
      T
      GetFFunction() const;
      
      T
      GetNormalisation() const;


      NaturalParameters<T>
      GetNaturalNot(const ChildNode* v) const;

      void 
      Iterate();
    private:
      std::vector<ParentNode*> m_parents;
      WeightsNode* m_weights;
      ChildNode* m_child;
      
      T m_FFunction, m_Normalisation;
      NaturalParameters<T> m_NP2Child;
      NaturalParameters<T> m_NP2Weight;
      std::vector<NaturalParameters<T> > m_NP2Parent;
    private:
      mutable boost::mutex m_mutex;
    };

  }
}
