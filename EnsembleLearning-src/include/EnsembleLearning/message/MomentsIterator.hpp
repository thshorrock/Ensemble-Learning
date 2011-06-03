#pragma once
#ifndef MOMENTSITERATOR_HPP
#define MOMENTSITERATOR_HPP



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

#include <boost/iterator/iterator_facade.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/call_traits.hpp>

namespace ICR{
  namespace EnsembleLearning{
    
    //forward declaration
    template<class> class Moments;
    
    /**Threadsafe iterator for the Moments class.
     * @tparam moments The moments class, intended to be either Moments<T> or Moments<T> const.
     * @tparam T The data type used in the iterator, intended to be either float or double.
     */ 
    template<class moments, class T>
    class  MomentsIterator
      : public boost::iterator_facade<MomentsIterator<moments,T>, T, boost::random_access_traversal_tag >
    {
      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{

    typedef  typename boost::iterator_facade<MomentsIterator<moments,T>, T, boost::random_access_traversal_tag > F;
      struct enabler {};
    public:
      
      typedef typename boost::call_traits<T>::param_type
      data_parameter;
      
      typedef typename boost::call_traits<T>::reference  
      data_reference;
      
      typedef typename boost::call_traits<T>::const_reference 
      data_const_reference;
      
      typedef typename boost::call_traits<T>::value_type
      data_type;

      ///@}

      /** A constructor */
      MomentsIterator()
	: m_Moments(0),
	  m_index(0),
	  m_mutex() 
      {}
      
      /** A constructor
       *  @param p A pointer to the Moments class that you want to iterate.
       *  @param index The index where you want to start the iterator.  
       *     For example, begin() would start at 0, end() would start at the size of the Moments. container
       */
      explicit MomentsIterator(moments* p,  typename F::difference_type  index)
	: m_Moments(p),
	  m_index(index),
	  m_mutex() 
      {}
      
      /** A copy constructor.
       *  @param other A Moments class of the same type (i.e. same level of constness) as the current iterator.
       */
      MomentsIterator(
      		      const MomentsIterator<moments,T>& other
      		      )
      : m_Moments(other.m_Moments),
       	m_index(other.m_index),
       	m_mutex() //don't copy the mutex (non-copiable)
      {}
      
      
      /** A copy constructor.
       *  @tparam Othermoments A Moments class of different constness.
       *  @tparam other_data_type The data type of different constness.
       *  @param other A Moments class of stricly different constness.  
       *  This function will fail if the constness is decreased, 
       *  (i.e. converting from const to non-const)
       *  with the  error  that 
       *  the Moments class (and not the MomentsIterator) cannot go from const to non-const.
       *  This error will be fixed in a future version.
       */
      //To fix the erronious compilation message when going from const to non-const
      //need to make this function vanish if decreasing constness.
      template <class Othermoments, class other_data_type>
      MomentsIterator(
      		      const MomentsIterator<Othermoments, other_data_type>& other,
      		      typename boost::disable_if<
      		      boost::is_same<Othermoments*,moments*>
      		      , enabler
      		      >::type = enabler()
      			)
      	 : m_Moments(other.m_Moments),
      	   m_index(other.m_index),
      	   m_mutex() //don't copy the mutex (non-copiable)
       {}


      
      /** Assignment.
       *  @tparam Othermoments A Moments class of different constness.
       *  @tparam other_data_type The data type of different constness.
       *  @param other A Moments class. 
       *  @return A reference to the current MomentsIterator.
       *  This function will fail if the constness is decreased, 
       *  (i.e. converting from const to non-const)
       *  with the  error  that 
       *  the Moments class (and not the MomentsIterator) cannot go from const to non-const.
       *  This error will be fixed in a future version.
       */
      //To fix the erronious compilation message when going from const to non-const
      //need to make this function vanish if decreasing constness.
      MomentsIterator<moments,T>&
      operator=(const MomentsIterator<moments,T>& other)
      {
      	if (this!=&other)
      	  {
      	    m_Moments = other.m_Moments;
      	    m_index   = other.m_index;
      	    //don't copy mutex;
      	  }
      	return *this;
      }


    private:
      /** Increment the iterator */
      void 
      increment() 
      {
	Lock lock(m_mutex); 
	++m_index;
      };
      /** Decrement the iterator */
      void 
      decrement() 
      { 
	Lock lock(m_mutex); 
	--m_index; 
      };
      
      /** Increment by a distance
       *  @param n The distance to increment
       */
      void 
      advance(const typename F::difference_type n)
      {
	Lock lock(m_mutex); 
	m_index+=n;
      };
      
      /** Find the distance to another iterator
       *  @param n The distance to increment
       */
      template < class Othermoments,class other_type>
      typename F::difference_type 
      distance_to(const MomentsIterator<Othermoments,other_type>& other   ) const 
      {
	Lock lock(m_mutex); 
	return other.m_index - m_index;
      };
      //Equality check.
      template <class Othermoments, class other_type>
      bool 
      equal(const MomentsIterator<Othermoments,other_type>& other) const
      {
      	Lock lock(m_mutex); 
      	return (other.m_index == m_index && m_Moments==other.m_Moments);
      }
      
      data_reference dereference() const
      { 
	Lock lock(m_mutex);
	return m_Moments->operator[](m_index);
      };

      data_reference dereference() 
      { 
	Lock lock(m_mutex);
	return m_Moments->operator[](m_index);
      };


      friend class boost::iterator_core_access;
      template <class, class> friend class MomentsIterator;
      moments* m_Moments;
      mutable  typename F::difference_type  m_index;
      mutable Mutex m_mutex;
    };
    
    
    
  }
}
#endif  // guard for MOMENTSITERATOR_HPP
