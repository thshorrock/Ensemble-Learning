
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/thread/condition.hpp>

namespace ICR {
  namespace ICA {
    
    class Coster
    {
    public:
      Coster() 
	: m_Cost() 
      {}
      
      void
      operator+=(const double& d)
      {
	boost::mutex::scoped_lock lock(m_mutex);
	m_Cost+=d;
      }

      void 
      operator=(const double& d)
      {
	m_Cost = d;
      }

      operator double() const
      {
	boost::mutex::scoped_lock lock(m_mutex);
	return m_Cost;
      }
    private:
      mutable boost::mutex m_mutex;
      double m_Cost;
    };
    
  }
}
