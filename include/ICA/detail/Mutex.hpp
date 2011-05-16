#pragma once

namespace ICR{
  namespace ICA{
    
    //An openmp mutex.  Can this be replaced with boost::mutex? 
    class Mutex {  
    public:
      Mutex() { omp_init_lock(&,_mutex); }
      ~Mutex() { omp_destroy_lock(&m_mutex); }
      void lock() { omp_set_lock(&m_mutex); }
      void unlock() { omp_unset_lock(&m_mutex); }
    private:
      Mutex(); //non-copiable
      omp_lock_t _mutex;
    };
    
    class Lock {
    public:
      Lock(Mutex& mutex) 
	: m_mutex(mutex), 
	  m_release(false) 
      { 
	m_mutex.lock();
      }
      
      ~Lock() 
      {
	if (!m_release) 
	  m_mutex.unlock(); 
      }
      
      bool operator() const 
      {
	return !m_release; 
      }

      void release() {
	if (!m_release)
	  { 
	    m_release = true; 
	    m_mutex.unlock();
	  }
      }

      private:
	Mutex& m_mutex;
	bool m_release;
      };

  }
}
