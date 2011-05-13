
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/call_traits.hpp>

namespace ICR {
  namespace ICA {
    
    /** Store the current (global) evidence bound in a thread safe way.
     *   The evidence is evaluated at each iteration of the algorithm.
     *   This class is passed to every VariableNode, potentially in parallel.
     *   Coster guaranees that the cost is handled in a thread safe way.
     */
    class Coster
    {
      typedef boost::call_traits<double>::param_type pDouble;
    public:
      /** Constructor. 
       *  @param cost The initial cost (default 0).
       */
      Coster(pDouble cost = 0.0) 
	: m_Cost(cost) 
      {}
      
      /** Copy Constructor. 
       */
      Coster(const Coster& other) 
	: m_Cost(other.m_Cost) 
      {}
      
      /** Add a double to the cost.
       * @param local The local cost on a variable node that is to be added to the global cost.
       */
      void
      operator+=(pDouble local)
      { 
	boost::lock_guard<boost::mutex> lock(m_mutex); 
	m_Cost+=local;
      }

      /** Assign the cost.
       * @param cost The new cost stored.
       */
      void 
      operator=(pDouble cost)
      {
	boost::lock_guard<boost::mutex> lock(m_mutex); 
	m_Cost = cost;
      }

      /** Implicitly convert the stored global evidence to a double.
       */
      operator double() const
      {
	boost::lock_guard<boost::mutex> lock(m_mutex); 
	return m_Cost;
      }
    private:
      mutable boost::mutex m_mutex;
      double m_Cost;
    };
    
  }
}
