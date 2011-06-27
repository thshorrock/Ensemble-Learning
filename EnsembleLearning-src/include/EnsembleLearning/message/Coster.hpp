#pragma once
#ifndef COSTER_HPP
#define COSTER_HPP


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




#include "EnsembleLearning/detail/Mutex.hpp"
#include <boost/call_traits.hpp>
#include <omp.h>

namespace ICR {
  namespace EnsembleLearning {
    
    /** Store the current (global) evidence bound in a thread safe way.
     *   The evidence is evaluated at each iteration of the algorithm.
     *   This class is passed to every VariableNode, potentially in parallel.
     *   Coster guaranees that the cost is handled in a thread safe way.
     */
    class Coster
    {
      typedef boost::call_traits<double>::param_type double_p;
    public:
      /** Constructor. 
       *  @param cost The initial cost (default 0).
       */
      Coster(double_p cost = 0.0) 
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
      operator+=(double_p local)
      { 
	//Lock lock(m_mutex); 
	#pragma omp atomic
	m_Cost+=local;
      }

      /** Assign the cost.
       * @param cost The new cost stored.
       */
      void 
      operator=(double_p cost)
      {
	//Lock lock(m_mutex); 
	m_Cost = cost;
      }

      /** Implicitly convert the stored global evidence to a double.
       */
      operator double() const
      {
	//Lock lock(m_mutex); 
	return m_Cost;
      }
    private:
      //mutable Mutex m_mutex;
      double m_Cost;
    };
    
  }
}
#endif  // guard for COSTER_HPP
