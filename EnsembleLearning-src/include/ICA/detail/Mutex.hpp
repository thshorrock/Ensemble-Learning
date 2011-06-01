#pragma once
#include <omp.h>

namespace ICR{
  namespace ICA{
    
    //An openmp mutex.  
    /** An opemmp mutex class.
     */
    class Mutex {  
    public:
      /** Contructor.
       *  Create the openmp mutex.
       */
      Mutex() { omp_init_lock(&m_mutex); }

      /** Destructor.
       *  Destroy the openmp mutex.
       */
      ~Mutex() { omp_destroy_lock(&m_mutex); }
      /** Lock the mutex */
      void lock() { omp_set_lock(&m_mutex); }
      /** Unlock the mutex */
      void unlock() { omp_unset_lock(&m_mutex); }
    private:
      Mutex(const Mutex&); //non-copiable
      omp_lock_t m_mutex;
    };
    
    /** Create an openmp mutex lock.
     */
    class Lock {
    public:
      /** Constructor.
       * Lock the openmp mutex on construction.
       * @param mutex The mutex to lock.
       */
      Lock(Mutex& mutex) 
	: m_mutex(mutex)
      { 
	m_mutex.lock();
      }

      /** Destructor
       * Unlock the openmp mutex on destruction.
       */
      ~Lock() 
      {
	m_mutex.unlock(); 
      }
      
    private:
      Mutex& m_mutex;
    };

  }
}
