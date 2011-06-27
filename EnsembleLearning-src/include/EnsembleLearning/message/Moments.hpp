#pragma once
#ifndef MOMENTS_HPP
#define MOMENTS_HPP


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



//#include "MomentsIterator.hpp"
//#include "EnsembleLearning/detail/Mutex.hpp"
//#include "EnsembleLearning/detail/parallel_algorithms.hpp"

#include <boost/call_traits.hpp>
#include <boost/array.hpp>
#include <iostream>
#include <algorithm>
#include <vector>

namespace ICR{
  namespace EnsembleLearning{
    
    
    /** A threadsafe container for the Moments.
     *  @tparam T The datatype to be used for storing the moments (typically double or float).
     */
    template<class T=double, int array_size = 2>
    class 
    Moments
    {
    public:
      
      /** @name Useful typdefs for types that are exposed to the user.
       */
      ///@{
      typedef typename boost::call_traits<std::vector<T> >::param_type 
      vector_parameter;
      
      typedef typename boost::call_traits<T>::param_type
      data_parameter;
      
      typedef typename boost::call_traits<T>::reference  
      data_reference;
      
      typedef typename boost::call_traits<T>::const_reference 
      data_const_reference;
      
      typedef typename boost::call_traits<T>::value_type
      data_type;
      
      
      typedef typename boost::call_traits<Moments<T,array_size> >::param_type 
      parameter;
      
      typedef typename boost::call_traits<Moments<T,array_size> >::reference  
      reference;
      
      typedef typename boost::call_traits<Moments<T,array_size> >::const_reference  
      const_reference;
      
      typedef typename boost::call_traits<Moments<T,array_size> >::value_type  
      type;
      
      typedef typename boost::call_traits<size_t>::param_type  
      size_parameter;
      
      typedef typename boost::call_traits<size_t>::value_type  
      size_type;
      
      typedef typename boost::call_traits<size_t>::reference  
      size_reference;
      
      typedef typename boost::call_traits<size_t>::const_reference  
      size_const_reference;
      
      // typedef  MomentsIterator<type, data_type>  
      // iterator;
      
      // typedef  MomentsIterator<const type, const data_type>  
      // const_iterator;
       
      typedef typename boost::array<data_type,array_size>::iterator
       iterator;
      
      typedef typename boost::array<data_type,array_size>::const_iterator
       const_iterator;
      
      ///@}

      /** Constructor 
       * @param size The size of the set of moments to store 
       */
      Moments();

      /** Constructor.
       *  Construct a set of moments based upon the passed vector.
       *  @param data The set of moments to be copied.
       */
      Moments(vector_parameter data);
      
      /** Constructor.
       *  Construct a pair of moments.
       * @param d1 The zeroth moment in the set.
       * @param d2 The first moment in the set.
       */
      Moments(data_parameter  d1, data_parameter d2);

      /** Copy constructor.
       * @param other The other set of moments to be copied.
       */
      Moments(parameter other);
      
      /** Assingment operator.
       *  @param other The  Moments container to be copied.
       *  @return A reference to the current container.
       */
      reference operator=(parameter other);
      
      /** Obtain an iterator for the first moment.
       *  @return An iterator pointing to the first moment.
       *  The returned iterator is thread safe 
       */
      iterator
      begin();

      /** Obtain an const_iterator for the first moment.
       *  @return A const_iterator pointing to the first moment.
       *  The returned iterator is thread safe 
       */
      const_iterator
      begin() const;

      /** Obtain an iterator for the last+1 moment.
       *  @return An iterator pointing to the last+1 moment.
       *  The returned iterator is thread safe 
       */
      iterator
      end();

      /** Obtain a const_iterator for the last+1 moment.
       *  @return An const_iterator pointing to the last+1 moment.
       *  The returned iterator is thread safe 
       */
      const_iterator
      end() const;

      /** The number of Moments stored.
       *  @return The number of moments stored.
       */
      size_type
      size() const {return m_data.size();}
      
      
      /** Access the moments data.
       *  @param index The index of the moment to return.
       *    For example, the average will be zeroth element,
       *    and the average of the squares will be the first element.
       *  @return The requested moment.
       */
      data_const_reference 
      operator[](size_parameter index) const;
      
      /** Access the moments data.
       *  @param index The index of the moment to return.
       *    For example, the average will be zeroth element,
       *    and the average of the squares will be the first element.
       *  @return The requested moment.
       */
      data_reference 
      operator[](size_parameter index);


      /** Add another Moments to this.
       * @param other The Moments to add.
       * @return A reference to the current Moments.
       */
      reference
      operator+=(parameter other);
      

      /** Multiply the moments by a scalar.
       * @param other The scalar by which to multiply the moments.
       * @return A reference to the current Moments.
       */
      reference
      operator*=(const data_parameter other);
      
      /** Multiply the moments by another set of moments.
       * @param other The other moments by which to multiply the stored moments.
       * @return A reference to the current, stored, Moments.
       */
      reference
      operator*=(const parameter other);



    private:
      
      
      struct plus
      {
	typename Moments<T>::data_type
	operator()(data_parameter i, 
		   data_parameter j) 
	{
	  return i+j;
	}
      };

      struct times
      {
      	typename Moments<T>::data_type
      	operator()(data_parameter i, 
      		   data_parameter j) 
      	{
      	  return i*j;
      	}
      };

      struct times_by{
	times_by(data_parameter t) : m_t(t) {};
  
	data_type
	operator()(data_parameter d)
	{ 
	  return d*m_t; 
	}
      private:
	T m_t;
      };


      boost::array<data_type,array_size> m_data;
      //mutable Mutex m_mutex;

    };


    /** Output the contents of the Moments Container.
     * @param out The output stream.
     * @param m The moments to output.
     * @return A reference to the output stream.
     */
    template <class T>
    inline
    std::ostream&
    operator<<( std::ostream& out, const Moments<T>& m)
    {
      for(size_t i=0;i<m.size();++i)
    	{
    	  out<<m[i]<<" ";
    	}
      return out;
    }

    /** Add two moments together.
     *  @param a The first moments container to add.
     *  @param b The second moments container to add
     *  @return The sum of the two moments.
     */
    template<class T,int array_size>
    inline
    const Moments<T,array_size>
    operator+(const Moments<T,array_size>& a, const Moments<T,array_size>& b)
    {
      Moments<T,array_size> tmp = a;
      return tmp+=b;
    }  
    
    /** Multiply a Moments contatiner by a scaler.
     *  @param m The  moments container.
     *  @param d The scalar.
     *  @return The product of the moments and the scaler.
     */
    template<class T,int array_size>  
    inline 
    const Moments<T,array_size>
    operator*(const Moments<T,array_size>& m,  
	      const T d)
    {
      Moments<T,array_size> tmp = m;
      return tmp*=d;
    }
      

    /** Multiply a Moments contatiner by a scaler.
     *  @param m The  moments container.
     *  @param d The scalar.
     *  @return The product of the moments and the scaler.
     */
    template<class T,int array_size>  
    inline 
    const Moments<T,array_size>
    operator*(const T d, 
	      const Moments<T,array_size>& m)
    {
      Moments<T,array_size> tmp = m;
      return tmp*=d;
    }
  }
}


template<class T,int array_size>
inline
ICR::EnsembleLearning::Moments<T,array_size>::Moments( vector_parameter data)
  : m_data()// ,
    // m_mutex()
{
  std::copy(data.begin(),data.end(),m_data.begin());
}

template<class T,int array_size> 
inline   
ICR::EnsembleLearning::Moments<T,array_size>::Moments()
  : m_data()
{
}

template<class T,int array_size> 
inline   
ICR::EnsembleLearning::Moments<T,array_size>::Moments(data_parameter d1,
						      data_parameter d2)
  : m_data()
    // m_mutex() //non-copiable
{
  m_data[0] = d1;
  m_data[1] = d2;
}
      
template<class T,int array_size> 
inline   
ICR::EnsembleLearning::Moments<T,array_size>::Moments(parameter other)
  : m_data(other.m_data)
    // m_mutex() //non-copiable
{}
      
template<class T,int array_size> 
inline   
typename ICR::EnsembleLearning::Moments<T,array_size>::reference
ICR::EnsembleLearning::Moments<T,array_size>::operator=(parameter other)
{    
  if (this!= &other) {
    m_data = (other.m_data);
    // m_mutex non-copiable  
  }
  return *this;
}

template<class T,int array_size>  
inline  
typename ICR::EnsembleLearning::Moments<T,array_size>::iterator
ICR::EnsembleLearning::Moments<T,array_size>::begin()
{
  return m_data.begin(); //iterator(this,0);
}

template<class T,int array_size>
inline
typename ICR::EnsembleLearning::Moments<T,array_size>::const_iterator
ICR::EnsembleLearning::Moments<T,array_size>::begin() const
{
  return m_data.begin();
  //  return const_iterator (this,0);
}

template<class T,int array_size>
inline
typename ICR::EnsembleLearning::Moments<T,array_size>::iterator
ICR::EnsembleLearning::Moments<T,array_size>::end()
{ 
  return m_data.end();
  //return iterator(this,m_data.size());
}


template<class T,int array_size>
inline
typename ICR::EnsembleLearning::Moments<T,array_size>::const_iterator
ICR::EnsembleLearning::Moments<T,array_size>::end() const
{
  return m_data.end();
  //  return const_iterator(this,m_data.size());
}


template<class T,int array_size>
inline
typename ICR::EnsembleLearning::Moments<T,array_size>::data_const_reference
ICR::EnsembleLearning::Moments<T,array_size>::operator[](size_parameter index) const
{
  // Lock lock(m_mutex); 
  return m_data[index];
}
      
template<class T,int array_size>
inline
typename ICR::EnsembleLearning::Moments<T,array_size>::data_reference
ICR::EnsembleLearning::Moments<T,array_size>::operator[](size_parameter index)
{
  // Lock lock(m_mutex); 
  return m_data[index];
}
      

template<class T,int array_size>
typename ICR::EnsembleLearning::Moments<T,array_size>::reference  
ICR::EnsembleLearning::Moments<T,array_size>::operator+=(parameter other)
{
  //more efficient than the std::transform
  for(size_t i=0;i<array_size;++i){
    m_data[i]+=other[i];
  }
  //std::transform(m_data.begin(), m_data.end(), other.begin(), m_data.begin(), plus());
  return *this;
}
      
template<class T,int array_size>  
inline 
typename ICR::EnsembleLearning::Moments<T,array_size>::reference
ICR::EnsembleLearning::Moments<T,array_size>::operator*=(parameter other)
{
  //more efficient than the std::transform
  for(size_t i=0;i<array_size;++i){
    m_data[i]*=other[i];
  }
  //std::transform(m_data.begin(), m_data.end(), other.begin(), m_data.begin(), times() );
  return *this;
}

template<class T,int array_size>  
inline 
typename ICR::EnsembleLearning::Moments<T,array_size>::reference
ICR::EnsembleLearning::Moments<T,array_size>::operator*=(data_parameter d)
{
  //more efficient than the std::transform
  for(size_t i=0;i<array_size;++i){
    m_data[i]*=d;
  }
  //std::transform(m_data.begin(), m_data.end(),  m_data.begin(), times_by(d));
  return *this;
}
      

#endif  // guard for MOMENTS_HPP
