#pragma once

namespace ICR{
  namespace ICA{
    
    //An openmp mutex.  Can this be replaced with boost::mutex? 
    class Mutex {  
    public:
      Mutex() { omp_init_lock(&m_mutex); }
      ~Mutex() { omp_destroy_lock(&m_mutex); }
      void lock() { omp_set_lock(&m_mutex); }
      void unlock() { omp_unset_lock(&m_mutex); }
    private:
      Mutex(const Mutex&); //non-copiable
      omp_lock_t m_mutex;
    };
    
    class Lock {
    public:
      Lock(Mutex& mutex) 
	: m_mutex(mutex)
      { 
	m_mutex.lock();
      }
      
      ~Lock() 
      {
	m_mutex.unlock(); 
      }
      
    private:
      Mutex& m_mutex;
    };

  }
}
