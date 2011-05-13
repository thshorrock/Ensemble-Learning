#include "ICA/variable/HiddenNode.hpp"

void
ICR::ICA::HiddenNode::Iterate(MessageGroup& messages)
{

  //This should only be called once s should get no collisions here.

  for_each(m_children.begin(), m_children.end(), 
	   boost::bind(MessageGroup::CreateVariableMessage,
		       boost::ref(messages),
		       *this, **it)
	   );
  messages.CreateVariableMessage(*this, *m_parent);

}

void
ICR::ICA::HiddenNode::EvaluateCost(Coster& Total)
{
  //This should only be called once s should get no collisions here.
  
  NaturalParameters init(0,0);
  NaturalParameters NP = std::accumulate(m_NPs.begin(), m_NPs.end(), init);
  m_Moments = m_Model.CalculateMoments(NP,this);

  NaturalParameters ParentNP = m_parent->GetNaturalNot(*this);
  double Cost = (ParentNP - NP)*m_moments + m_parent->GetNormalisation() - GetNormalisation();
  Total += Cost;
  // std::cout<<"Hidden Node "<<this<<" Cost = "<<Cost<<std::endl;
  // std::cout<<"Of which:"
  // 	   <<"\n\t GetNormalisation() = "<<GetNormalisation()
  // 	   <<"\n\t P  Normalisation() = "<<m_parent->GetNormalisation()
  // 	   <<"\n\t ParentNP           = "<<ParentNP
  // 	   <<"\n\t CurrentNP          = "<<m_CurrentNP
  // 	   <<"\n\t Diff               = "<< (ParentNP - m_CurrentNP)
  // 	   <<"\n\t m_moments          = "<<m_moments
  // 	   <<"\n\t Prod =             = "<<(ParentNP - m_CurrentNP)*m_moments
  // 	   <<std::endl;
      

}

void
ICR::ICA::HiddenNode::Reset()
{
  //boost::mutex mutex;
  //boost::mutex::scoped_lock lock(mutex);
  if (m_has_iterated) {
    m_EvaluatedCost = false;
    m_has_iterated = false;
    m_NPs.clear();
    for_each(m_children.begin(), m_children.end(), boost::bind(&FactorNode::Reset, _1));
    m_parent->Reset();
  }
	
}
