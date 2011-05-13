#include "ICA/factor/GammaFactor.hpp"
#include <gsl/gsl_sf_gamma.h>

ICR::ICA::GammaFactor::GammaFactor(VariableNode& shape, VariableNode& scale, VariableNode& child)
  : m_shape_node(&shape), m_scale_node(&scale), m_child_node(&child),  m_has_iterated(false), m_EvaluatedCost(false)
{

  // boost::shared_ptr<GammaFactor> tmp(this);
  m_shape_node -> AddChildFactor(*this);
  m_scale_node -> AddChildFactor(*this);
  m_child_node -> SetParentFactor(*this);
	
}

ICR::ICA::NaturalParameters
ICR::ICA::GammaFactor::GetNaturalNot(const VariableNode& v) const
{
  // boost::mutex mutex;
  //  boost::condition cond;
  // std::cout<<"shapep = "<<m_shape_node<<std::endl;
  // std::cout<<"scalep = "<<m_scale_node<<std::endl;
  // std::cout<<"childp = "<<m_child_node<<std::endl;
      
  //std::cerr<<"GetNaturalNot "<<std::endl;
  if (&v==m_shape_node)
    {
      std::cerr<<"Can't pass message to shape variable  of Gamma Distribution"<<&v<<std::endl;
      throw ("Can't pass message to shape variable  of Gamma Distribution");
    }
  else if (&v==m_scale_node)
    {
      //std::cerr<<"scale node"<<std::endl;

      boost::mutex::scoped_lock lock(mutex);
      while (!m_NP2Scale) //lock the thread until parameter is set.
	cond.wait(lock);
      //std::cerr<<"Returning Scale NP"<<std::endl;
      return *m_NP2Scale;
    }
  else if (&v == m_child_node) 
    {
      //std::cerr<<"child node"<<std::endl;
      //boost::mutex::scoped_lock lock(mutex);
      boost::mutex::scoped_lock lock(mutex);
      while (!m_NP2Child) //lock the thread until parameter is set.
	cond.wait(lock);
      //std::cerr<<"Returning Child NP"<<std::endl;
      return *m_NP2Child;
    }
  else{
    //std::cout<<"v = "<<&v<<std::endl;
    //std::cout<<"shape = "<<m_shape_node<<std::endl;

    throw ("Unknown Node in GetNaturalNot");
  }
}

void
ICR::ICA::GammaFactor::ReceiveMoments(const Moments& m)
{
  //std::cout<<"Receive Moments"<<std::endl;
  //std::cout<<" source = "<<m.source()<<std::endl;

  if (m.source() == m_shape_node ) 
    { 
      boost::mutex::scoped_lock lock(mutex);
        //std::cerr<<"shape Moment"<<std::endl;
      while (m_shape) 
	cond.wait(lock); //already got this message. block thread till reset.
	    
      m_shape = m.first();
      cond.notify_all();
      //std::cerr<<"got shape Moment"<<std::endl;
    }
  else if (m.source() == m_scale_node) 
    {
      boost::mutex::scoped_lock lock(mutex);
      //std::cerr<<"scale Moment"<<std::endl;
      while (m_scale) 
	cond.wait(lock); //already got this message. block thread till reset.
      m_scale = m.first();
      cond.notify_all();
      //std::cerr<<"got scale Moment"<<std::endl;
	    
    }
  else if (m.source() == m_child_node) 
    {
      boost::mutex::scoped_lock lock(mutex);
      //std::cout<<"child Moment"<<std::endl;
      while (m_child) 
	cond.wait(lock); //already got this message. block thread till reset.
	    
      m_child = m.first();
      cond.notify_all();
      //std::cerr<<"got child Moment"<<std::endl;
    }
  else
    {
      //std::cout<<"source = "<<m.source()<<std::endl;
    }
  if (m_shape && m_scale) {
    {
      boost::mutex::scoped_lock lock(mutex);
      m_NP2Child = NaturalParameters(-*m_scale, *m_shape - 1);
    }
    cond.notify_all();
    //if (m_NP2Child) 
       //std::cerr<<"NP2Child set"<<std::endl;
	    
  }
  if (m_shape && m_child) {
    {
      boost::mutex::scoped_lock lock(mutex);
      m_NP2Scale = NaturalParameters(-*m_child, *m_shape);
    }
    cond.notify_all();
    //if (m_NP2Scale) 
       //std::cerr<<"NP2Scale set"<<std::endl;
	  
  }
}

void
ICR::ICA::GammaFactor::EvaluateCost(double& Total)
{
  if (!m_EvaluatedCost)
    {
      // m_EvaluatedCost=true;
      // const NaturalParameters NP = *m_NP2Child;
      // const double a = NP[1]+1;
      // const double b = -NP[0] ;
      // const Moments m = m_child_node->GetMoments();
      // const double CostFromThisIteration  = NP[0]*m[0] + NP[1]*m[1] + a*std::log(b) - std::log(gsl_sf_gamma (a));
      // if (m_child_node->IsData()) 
      // 	Total+= CostFromThisIteration;
      // else
      // 	Total+= CostFromThisIteration - m_CostFromLastIteration;
      // m_CostFromLastIteration = CostFromThisIteration;

      m_shape_node->EvaluateCost(Total);
      m_scale_node->EvaluateCost(Total);
      m_child_node->EvaluateCost(Total);
    }
}


void
ICR::ICA::GammaFactor::Iterate( boost::thread_group& messages)
{
  // boost::mutex mutex;
  //boost::mutex::scoped_lock lock(mutex);
  {
    if (!m_has_iterated)
      {
	//std::cout<<"Gamma It, size = "<<messages.size()<<std::endl;
	//std::cout<<"GAMMA iterate!"<<std::endl;

	m_has_iterated = true;
	//boost::shared_ptr<GammaFactor> tmp(this);
	messages.create_thread(FactorMessage(*this, *m_scale_node));
	messages.create_thread(FactorMessage(*this, *m_child_node));
	//std::cout<<"SHAPE IT"<<std::endl;

	m_shape_node->Iterate(messages);
	//std::cout<<"SCALE IT"<<std::endl;
	m_scale_node->Iterate(messages);
	//std::cout<<"CHILD IT"<<std::endl;
	m_child_node->Iterate(messages);
	//std::cout<<"Done IT"<<std::endl;
	    
	//std::cout<<"Gamma It end, size = "<<messages.size()<<std::endl;
      }
  }
}

void
ICR::ICA::GammaFactor::Reset()
{
  // boost::mutex mutex;
  //boost::mutex::scoped_lock lock(mutex);
  if (m_has_iterated) {
    m_has_iterated = false;
    m_EvaluatedCost = false;
    m_shape=boost::none;
    m_scale=boost::none;
    m_child=boost::none;
    m_NP2Child=boost::none;
    m_NP2Scale=boost::none;

    m_shape_node->Reset();
    m_scale_node->Reset();
    m_child_node->Reset();
	    

  }
	
}
