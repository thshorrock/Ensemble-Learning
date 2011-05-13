#include "ICA/variable/ConstantNode.hpp"

void
ICR::ICA::ConstantNode::Iterate(boost::thread_group& messages)
{

  //boost::mutex::scoped_lock lock(mutex);
  if (!m_has_iterated)
    {
      //std::cout<<"Const It "<<this<<", size = "<<messages.size()<<std::endl;
      m_has_iterated = true;
      for_each(m_children.begin(), m_children.end(), boost::bind(&FactorNode::Iterate, _1, boost::ref(messages)));
      //boost::shared_ptr<ConstantNode> tmp(this);
	
       //std::cout<<"num children = "<<m_children.size()<<std::endl;

      for(std::list<FactorNode*>::iterator it = m_children.begin();
	  it != m_children.end();
	  ++it)
	{
	  messages.create_thread(VariableMessage( *this,**it));
	}
      //std::cout<<"Const It end "<<this<<", size = "<<messages.size()<<std::endl;
    }
}

void
ICR::ICA::ConstantNode::EvaluateCost(double& Total)
{

  //boost::mutex::scoped_lock lock(mutex);
  if (!m_EvaluatedCost)
    {
      m_EvaluatedCost = true;
      for_each(m_children.begin(), m_children.end(), boost::bind(&FactorNode::EvaluateCost, _1, boost::ref(Total)));
    }
}

void
ICR::ICA::ConstantNode::Reset()
{
  //boost::mutex mutex;
  //boost::mutex::scoped_lock lock(mutex);
  if (m_has_iterated) {
    m_has_iterated = false;
    for_each(m_children.begin(), m_children.end(), boost::bind(&FactorNode::Reset, _1));
  }
	
}
