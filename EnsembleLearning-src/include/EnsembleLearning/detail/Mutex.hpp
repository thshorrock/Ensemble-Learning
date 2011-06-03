#pragma once
#ifndef MUTEX_HPP
#define MUTEX_HPP


/***********************************************************************************
 ***********************************************************************************
 **                                                                               **
 **  Copyright (C) 2011 Tom Shorrock <t.h.shorrock@gmail.com> 
 **                                                                               **
 **                                                                               **
 **  This program is free software; you can redistribute it and/or                **
 **  modify it under the terms of the GNU General Public License                  **
 **  as published by the Free Software Foundation; either version 2               **
 **  of the License, or (at your option) any later version.                       **
 **                                                                               **
 **  This program is distributed in the hope that it will be useful,              **
 **  but WITHOUT ANY WARRANTY; without even the implied warranty of               **
 **  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the                **
 **  GNU General Public License for more details.                                 **
 **                                                                               **
 **  You should have received a copy of the GNU General Public License            **
 **  along with this program; if not, write to the Free Software                  **
 **  Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.  **
 **                                                                               **
 ***********************************************************************************
 ***********************************************************************************/




#include <omp.h>

namespace ICR{
  namespace EnsembleLearning{
    
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
#endif  // guard for MUTEX_HPP
