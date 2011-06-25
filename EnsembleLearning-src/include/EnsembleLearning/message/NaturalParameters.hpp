#pragma once
#ifndef NATURALPARAMETERS_HPP
#define NATURALPARAMETERS_HPP



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

#include "EnsembleLearning/detail/parallel_algorithms.hpp"
#include "EnsembleLearning/message/Moments.hpp"

#include <boost/call_traits.hpp>


#include <boost/array.hpp>
#include <iostream>
#include <algorithm>
//#include<vector>


namespace ICR{
  namespace EnsembleLearning{
    

    /** A container for the Natural Paramemeters.
     *  @tparam T The data type to be used.
     *  This is intended to be either float or double.
     *  @attention Natural Parameters is not thread safe - it is intended to be a temporary container for passing messages between nodes.
     */
    template<class T, int array_size = 2>
    class NaturalParameters
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

      typedef typename boost::call_traits<NaturalParameters<T,array_size> >::param_type 
      parameter;
      
      typedef typename boost::call_traits<NaturalParameters<T,array_size> >::reference  
      reference;
      
      typedef typename boost::call_traits<NaturalParameters<T,array_size> >::const_reference  
      const_reference;
      
      typedef typename boost::call_traits<NaturalParameters<T,array_size> >::value_type  
      type;
      
      typedef typename boost::call_traits<size_t>::param_type  
      size_parameter;
      
      typedef typename boost::call_traits<size_t>::value_type  
      size_type;
      
      typedef typename boost::array<data_type,array_size>::iterator
       iterator;
      
      typedef typename boost::array<data_type,array_size>::const_iterator
       const_iterator;
   
      // typedef typename std::vector<T>::iterator
      // iterator;
      
      // typedef typename std::vector<T>::const_iterator
      // const_iterator;

      ///@}

      /** Constructor.
       *  @param size The number of elements to be held in the container.
       */
      NaturalParameters();
      
      /** Constructor.
       *  @param NP A vector containing the the Natural Parameters to be copied into the container.
       */
      NaturalParameters(vector_parameter NP);

      /** Constructor.
       *  Construct a container holding two elements.
       *  @param d1 The zeroth element of the container.
       *  @param d2 The first element of the container.
       */
      NaturalParameters(data_parameter d1, 
			data_parameter d2);

      /** Obtain a constant reference to an element of the  Natural Parameter container.
       *  @param element The index of the container to return.
       *  @return The element of the container requested by element.
       *  @attention Natural Parameters is not thread safe - it is intended to be a temporary container for passing messages between nodes.
       */
      data_const_reference 
      operator[](size_parameter element) const;
      
      /** Obtain a reference to an element of the  Natural Parameter container.
       *  @param element The index of the container to return.
       *  @return The element of the container requested by element.
       *  @attention Natural Parameters is not thread safe - it is intended to be a temporary container for passing messages between nodes.
       */
      data_reference
      operator[](size_parameter element) ;

      /** Obtain an iterator for the zeroth natural parameter.
       *  @return An iterator pointing to the zeroth natural parameter.
       *  The returned iterator is not thread safe 
       */
      iterator
      begin();
      
      /** Obtain an const_iterator for the zeroth natural parameter.
       *  @return A const_iterator pointing to the  zeroth natural parameter.
       *  The returned iterator is not thread safe 
       */
      const_iterator
      begin() const;

      /** Obtain an iterator for the end of the container.
       *  @return An iterator pointing to the end of the container
       *  The returned iterator is not thread safe 
       */
      iterator
      end();

      /** Obtain an const_iterator for the end of the container.
       *  @return A const_iterator pointing to the end of the container
       *  The returned iterator is not thread safe 
       */
      const_iterator
      end() const;

      /** The number of Natural Parameters stored.
       *  @return The number of Natural Parameters stored.
       */
      size_type
      size() const {return m_data.size();}
      
      /** Add another set of Natural Parameters to this container.
       *  @param other The other Natural Parameters container.
       *  @return A reference to the current (and now altered) Natural Paramters container.
       */
      reference
      operator+=(parameter other);
      
      /** Subtract another set of Natural Parameters to this container.
       *  @param other The other Natural Parameters container.
       *  @return A reference to the current (and now altered) Natural Paramters container.
       */
      reference
      operator-=(parameter other);

      /** Multiply this container by a scalar.
       *  @param other The scalar
       *  @return A reference to the current (and now altered) Natural Paramters container.
       */
      reference
      operator*=(data_parameter other);
      
      /** Multiply a Natural Paramers container by a Moments Container.
       *  This is an inner-product like operation.
       * @param a The Natural Parameters container.
       * @param b The Moments container.
       * @return The scalar product from the two containers.
       */
      template<class U,int array_size2>
      friend
      U
      operator*(const NaturalParameters<U,array_size2>&  a, 
		const Moments<U,array_size2>& b);
    private:

      struct plus
      {
	data_type
	operator()(data_parameter i, 
		   data_parameter j) 
	{
	  return i+j;
	}
      };

      struct minus
      {
	data_type
	operator()(data_parameter i, 
		   data_parameter j) 
	{
	  return i-j;
	}
      };

      struct times
      {
      	data_type
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
	data_type m_t;
      };

      boost::array<data_type,array_size> m_data;
      

    };

    /** Output the contents of the Natural Parameters Container.
     * @param out The output stream.
     * @param m The Natural Parameters to output.
     * @return A reference to the output stream.
     */
    template<class T,int array_size>
    inline
    std::ostream&
    operator<<(std::ostream& out, const NaturalParameters<T,array_size>& NP)
    {
      for(size_t i=0;i<NP.size();++i)
	{
	  out<<NP[i]<<" ";
	}

      return out;
    }
    

    
    /** Add two Natural Parameters containers together.
     *  @param a The first Natural Parameters  container to add.
     *  @param b The second Natural Parameters container to add
     *  @return The sum of the two Natural Parameters.
     */
    template<class T,int array_size>
    inline
    NaturalParameters<T,array_size>
    operator+(const NaturalParameters<T,array_size>& a,
	      const NaturalParameters<T,array_size>& b)
    {
      NaturalParameters<T,array_size> tmp = a;
      return tmp+=b;
    }  
    
    /** Subtract two Natural Parameters containers.
     *  @param a The first Natural Parameters  container.
     *  @param b The second Natural Parameters container.
     *  @return The result of the subtraction.
     */
    template<class T,int array_size>
    inline
    NaturalParameters<T,array_size>
    operator-(const NaturalParameters<T,array_size>& a,
	      const NaturalParameters<T,array_size>& b)
    {
      ICR::EnsembleLearning::NaturalParameters<T,array_size> tmp = a;
      return tmp-=b;
    }  

    
 
      /** Multiply a Natural Paramers container by a Moments Container.
       *  This is an inner-product like operation.
       * @param a The Natural Parameters container.
       * @param b The Moments container.
       * @return The scalar product from the two containers.
       */

    template<class T,int array_size>
    inline
    T
    operator*(const NaturalParameters<T,array_size>&  a, 
	      const Moments<T,array_size>& b)
    { 
      return std::inner_product(a.begin(), a.end(), b.begin(), 0.0);
      // std::vector<T,array_size> prod(a.size());
      // PARALLEL_TRANSFORM(a.begin(), a.end(), b.begin(), prod.begin(), 
      // 			 typename NaturalParameters<T,array_size>::times());
      // T sum = 0.0;
      // return PARALLEL_ACCUMULATE(prod.begin(), prod.end(), sum);
    }  
      

    /** Multiply a Natural Paramers container by a Moments Container.
     *  This is an inner-product like operation.
     * @param p The Natural Parameters container.
     * @param m The Moments container.
     * @return The scalar product from the two containers.
     */
    template<class T,int array_size>
    inline
    T
    operator*(const Moments<T,array_size>& m,const NaturalParameters<T,array_size>&  p)
    {
      return p*m;
    }  

    /** Multiply a Natural Parameters container by a scalar.
     * @param n The Natural Parameters container.
     * @param d The scalar.
     * @return The product of the natural parameters and the scalar.
     */
    template<class T,int array_size>
    inline
    const NaturalParameters<T,array_size>
    operator*(const NaturalParameters<T,array_size>& n,const T  d)
    {
      NaturalParameters<T,array_size> tmp(n);
      return tmp*=d;
    }  



  }
}

template<class T,int array_size>
inline
ICR::EnsembleLearning::NaturalParameters<T,array_size>::NaturalParameters()
  : 
  m_data()
{
}

template<class T,int array_size>
inline
ICR::EnsembleLearning::NaturalParameters<T,array_size>::NaturalParameters(vector_parameter data)
  : //m_mutex_ptr(new boost::mutex), 
    m_data()
{
  std::copy(data.begin(),data.end(),m_data.begin());
}

template<class T,int array_size>
inline
ICR::EnsembleLearning::NaturalParameters<T,array_size>::NaturalParameters(data_parameter d1,
						  data_parameter d2 )
//  : 
  //  m_mutex_ptr(new boost::mutex)
  : m_data()
{
  //boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  m_data[0] = d1;
  m_data[1] = d2;
}

template<class T,int array_size>
typename ICR::EnsembleLearning::NaturalParameters<T,array_size>::iterator
ICR::EnsembleLearning::NaturalParameters<T,array_size>::begin()
{
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.begin();
}


template<class T,int array_size>
typename ICR::EnsembleLearning::NaturalParameters<T,array_size>::const_iterator
ICR::EnsembleLearning::NaturalParameters<T,array_size>::begin() const
{
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.begin();
}


template<class T,int array_size>
typename ICR::EnsembleLearning::NaturalParameters<T,array_size>::iterator
ICR::EnsembleLearning::NaturalParameters<T,array_size>::end()
{
  // boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.end();
}


template<class T,int array_size>
typename ICR::EnsembleLearning::NaturalParameters<T,array_size>::const_iterator
ICR::EnsembleLearning::NaturalParameters<T,array_size>::end() const
{
  // boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.end();
}


 

template<class T,int array_size>
inline
typename ICR::EnsembleLearning::NaturalParameters<T,array_size>::data_const_reference 
ICR::EnsembleLearning::NaturalParameters<T,array_size>::operator[](size_parameter i) const
{
  //boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data[i];
}
   
template<class T,int array_size>
inline
typename ICR::EnsembleLearning::NaturalParameters<T,array_size>::data_reference
ICR::EnsembleLearning::NaturalParameters<T,array_size>::operator[](size_parameter i) 
{
  //boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data[i];
}
   


template<class T,int array_size>
inline
typename ICR::EnsembleLearning::NaturalParameters<T,array_size>::reference
ICR::EnsembleLearning::NaturalParameters<T,array_size>::operator+=(parameter other)
{
  //This could be called by different threads, so lock
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  //DO KEEP THIS ONE!
  //  std::cout<<"size = "<<m_data.size()<<"other = "<<m_data.size<<std::endl;

  std::transform(m_data.begin(), m_data.end(), other.begin(), m_data.begin(), plus());
  // for(size_t i=0;i<other.size();++i){
  //   m_data[i]+=other[i];
  // }

  return *this;
}
   

template<class T,int array_size>  
inline 
typename ICR::EnsembleLearning::NaturalParameters<T,array_size>::reference
ICR::EnsembleLearning::NaturalParameters<T,array_size>::operator-=(parameter other)
{
  //This could be called by different threads, so lock
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  //DO KEEP THIS ONE!
  std::transform(m_data.begin(), m_data.end(), other.begin(), m_data.begin(), minus());
  return *this;
}
  
template<class T,int array_size>
inline
typename ICR::EnsembleLearning::NaturalParameters<T,array_size>::reference
ICR::EnsembleLearning::NaturalParameters<T,array_size>::operator*=(data_parameter other)
{
  
  std::transform(m_data.begin(), m_data.end(), m_data.begin(), times_by(other));
  return *this;
}    



#endif  // guard for NATURALPARAMETERS_HPP
