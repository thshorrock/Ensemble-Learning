#include "ICA/factor/GaussianFactor.hpp"

template<class Model>
ICR::ICA::FactorNode::FactorNode( VariableNode& v1,  VariableNode& v2,  VariableNode& child)
  : Piping<Model>(&v1,&v2, &child);
{
  v1.AddChildFactor(*this);
  v2.AddChildFactor(*this);
  child.SetParentFactor(*this);
}


void
ICR::ICA::FactorNode::Iterate( boost::thread_group& messages)
{
  std::vector<VariableNode*> Neighbours = Piping<Model>::GetNeighbours();
  for(std::list<FactorNode*>::iterator it = Neighbours.begin();
      it != Neighbours.end();
      ++it)
    {
      messages.create_thread(FactorMessage(*this,**it));
    }
}

	
}
