#include "ICA/factor/GaussianFactor.hpp"


ICR::ICA::GaussianFactor::GaussianFactor( VariableNode& mean,  VariableNode& precision,  VariableNode& child)
  : m_mean_node(&mean), m_precision_node(&precision), m_child_node(&child),  m_has_iterated(false), m_EvaluatedCost(false)
{

  // boost::shared_ptr<GaussianFactor> tmp(this);
  m_mean_node      -> AddChildFactor(*this);
  m_precision_node -> AddChildFactor(*this);
  m_child_node     -> SetParentFactor(*this);
	
}

ICR::ICA::NaturalParameters
ICR::ICA::GaussianFactor::GetNaturalNot(const VariableNode& v) const
{
  // boost::mutex mutex;
  //  boost::condition cond;
  // //std::cout<<"meanp = "<<m_mean_node<<std::endl;
  // std::cout<<"precisionp = "<<m_precision_node<<std::endl;
  // std::cout<<"childp = "<<m_child_node<<std::endl;
      
  //std::cerr<<"GetNaturalNot "<<std::endl;
  if (&v==m_mean_node)
    {
      
      //std::cerr<<"mean node"<<std::endl;

      boost::mutex::scoped_lock lock(mutex);
      while (!m_NP2Mean) //lock the thread until parameter is set.
	cond.wait(lock);
      //std::cerr<<"Returning Mean NP"<<std::endl;
      return *m_NP2Mean;
    }
  else if (&v==m_precision_node)
    {
      //std::cerr<<"precision node"<<std::endl;

      boost::mutex::scoped_lock lock(mutex);
      while (!m_NP2Precision) //lock the thread until parameter is set.
	cond.wait(lock);
      //std::cerr<<"Returning Precision NP"<<std::endl;
      return *m_NP2Precision;
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
    //std::cout<<"mean = "<<m_mean_node<<std::endl;

    throw ("Unknown Node in GetNaturalNot");
  }
}

void
ICR::ICA::GaussianFactor::ReceiveMoments(const Moments& m)
{
  //std::cout<<"Receive Moments"<<std::endl;
  //std::cout<<" source = "<<m.source()<<std::endl;

  if (m.source() == m_mean_node ) 
    { 
      boost::mutex::scoped_lock lock(mutex);
        //std::cerr<<"mean Moment"<<std::endl;
      while (m_mean) 
	cond.wait(lock); //already got this message. block thread till reset.
	    
      m_mean = m.first();
      m_mean_squared = m.second();
      cond.notify_all();
      //std::cerr<<"got mean Moment"<<std::endl;
    }
  else if (m.source() == m_precision_node) 
    {
      boost::mutex::scoped_lock lock(mutex);
      //std::cerr<<"precision Moment"<<std::endl;
      while (m_precision) 
	cond.wait(lock); //already got this message. block thread till reset.
      m_precision = m.first();
      cond.notify_all();
      //std::cerr<<"got precision Moment"<<std::endl;
	    
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
  if (m_mean && m_precision) {
    {
      boost::mutex::scoped_lock lock(mutex);
      double mean = *m_mean;
      double prec = *m_precision;
      m_NP2Child = NaturalParameters(mean*prec, -0.5*prec);
    }
    cond.notify_all();
    //if (m_NP2Child) 
       //std::cerr<<"NP2Child set"<<std::endl;
	    
  }
  if (m_mean && m_child) {
    {
      boost::mutex::scoped_lock lock(mutex);
      double mean = *m_mean; //average of mean
      double mean_squared = *m_mean_squared; //average of mean squared, so not mean*mean
      double child = *m_child;
      m_NP2Precision = NaturalParameters(-0.5*(child*child - 2*child*mean + mean_squared), 0.5 );
    }
    cond.notify_all();
    //if (m_NP2Precision) 
       //std::cerr<<"NP2Precision set"<<std::endl;
	  
  }
  if (m_precision && m_child) {
    {
      boost::mutex::scoped_lock lock(mutex);
      double child = *m_child;
      double prec  = *m_precision;
      m_NP2Mean = NaturalParameters(prec*child, -0.5*prec);
    }
    cond.notify_all();
    //if (m_NP2Mean) 
       //std::cerr<<"NP2Mean set"<<std::endl;
	  
  }
}

void
ICR::ICA::GaussianFactor::EvaluateCost(double& Total)
{
  if (!m_EvaluatedCost)
    {
      // m_EvaluatedCost=true;
      // const NaturalParameters NP = *m_NP2Child;
      // const double prec = -2.0* NP[1];
      // const double mean = NP[0]/prec ;
      // const Moments m = m_child_node->GetMoments();
      // const double CostFromThisIteration  = NP[0]*m[0] + NP[1]*m[1] + 0.5*(std::log(prec) - prec*mean*mean - std::log(2*M_PI) );
      // if (m_child_node->IsData()) 
      // 	Total+= CostFromThisIteration;
      // else
      // 	Total+= CostFromThisIteration - m_CostFromLastIteration;
      // m_CostFromLastIteration = CostFromThisIteration;

      m_mean_node->EvaluateCost(Total);
      m_precision_node->EvaluateCost(Total);
      m_child_node->EvaluateCost(Total);
    }
}


void
ICR::ICA::GaussianFactor::Iterate( boost::thread_group& messages)
{
  // boost::mutex mutex;
  //boost::mutex::scoped_lock lock(mutex);
  {
    if (!m_has_iterated)
      {
	//std::cout<<"Gaussian It , size = "<<messages.size()<<std::endl;
	//std::cout<<"GAUSSIAN iterate!"<<std::endl;

	m_has_iterated = true;
	//boost::shared_ptr<GaussianFactor> tmp(this);
	messages.create_thread(FactorMessage(*this, *m_mean_node));
	messages.create_thread(FactorMessage(*this, *m_precision_node));
	messages.create_thread(FactorMessage(*this, *m_child_node));
	//std::cout<<"MEAN IT"<<std::endl;

	m_mean_node->Iterate(messages);
	//std::cout<<"PRECISION IT"<<std::endl;
	m_precision_node->Iterate(messages);
	//std::cout<<"CHILD IT"<<std::endl;
	m_child_node->Iterate(messages);
	//std::cout<<"Done IT"<<std::endl;
	    
	//std::cout<<"Gaussian It end, size = "<<messages.size()<<std::endl;
      }
  }
}

void
ICR::ICA::GaussianFactor::Reset()
{
  // boost::mutex mutex;
  //boost::mutex::scoped_lock lock(mutex);
  if (m_has_iterated) {
    m_has_iterated = false;
    m_mean=boost::none;
    m_precision=boost::none;
    m_child=boost::none;
    m_NP2Mean=boost::none;
    m_NP2Child=boost::none;
    m_NP2Precision=boost::none;

    m_mean_node->Reset();
    m_precision_node->Reset();
    m_child_node->Reset();
	    

  }
	
}
