#pragma once
//#include "NaturalParameters.hpp"
#include "MomentsIterator.hpp"
#include "ICA/node/Node.hpp"
#include "ICA/detail/Mutex.hpp"
#include "ICA/detail/parallel_algorithms.hpp"
#include <boost/assign/list_of.hpp>
#include <boost/assign/std/vector.hpp>
#include <iostream>

namespace ICR{
  namespace ICA{
    
    //Forward declaration.
    template<class> class Moments;
    
    /** A threadsafe container for the Moments.
     *  @tparam T The datatype to be used for storing the moments (typically double or float).
     */
    template<class T=double>
    class 
    Moments
    {
    public:
      
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
      
      
      typedef typename boost::call_traits<Moments<T> >::param_type 
      parameter;
      
      typedef typename boost::call_traits<Moments<T> >::reference  
      reference;
      
      typedef typename boost::call_traits<Moments<T> >::const_reference  
      const_reference;
      
      typedef typename boost::call_traits<Moments<T> >::value_type  
      type;
      
      typedef typename boost::call_traits<size_t>::param_type  
      size_parameter;
      
      typedef typename boost::call_traits<size_t>::value_type  
      size_type;
      
      typedef typename boost::call_traits<size_t>::reference  
      size_reference;
      
      typedef typename boost::call_traits<size_t>::const_reference  
      size_const_reference;
      
      typedef  MomentsIterator<type, data_type>  
      iterator;
      
      typedef  MomentsIterator<const type, const data_type>  
      const_iterator;
      
      /** Constructor 
       * @param size The size of the set of moments to store 
       */
      Moments(size_parameter size = 0);

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


      std::vector<T> m_data;
      mutable Mutex m_mutex;

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
    template<class T>
    inline
    const Moments<T>
    operator+(const Moments<T>& a, const Moments<T>& b)
    {
      Moments<T> tmp = a;
      return tmp+=b;
    }  
    
    /** Multiply a Moments contatiner by a scaler.
     *  @param m The  moments container.
     *  @param d The scalar.
     *  @return The product of the moments and the scaler.
     */
    template<class T>  
    inline 
    const Moments<T>
    operator*(const Moments<T>& m,  
	      const T d)
    {
      Moments<T> tmp = m;
      return tmp*=d;
    }
      

    /** Multiply a Moments contatiner by a scaler.
     *  @param m The  moments container.
     *  @param d The scalar.
     *  @return The product of the moments and the scaler.
     */
    template<class T>  
    inline 
    const Moments<T>
    operator*(const T d, 
	      const Moments<T>& m)
    {
      Moments<T> tmp = m;
      return tmp*=d;
    }
  }
}


template<class T>
inline
ICR::ICA::Moments<T>::Moments( vector_parameter data)
  : m_data(data),
    m_mutex()
{}

template<class T> 
inline   
ICR::ICA::Moments<T>::Moments( size_parameter size)
  : m_data(std::vector<T>(size)),
    m_mutex() //non-copiable
{}

template<class T> 
inline   
ICR::ICA::Moments<T>::Moments(data_parameter d1,
			      data_parameter d2)
  : m_data(),
    m_mutex() //non-copiable
{
  using namespace boost::assign;
  m_data += d1,d2; // insert values at the end of the container
}
      
template<class T> 
inline   
ICR::ICA::Moments<T>::Moments(parameter other)
  : m_data(other.m_data),
    m_mutex() //non-copiable
{}
      
template<class T> 
inline   
typename ICR::ICA::Moments<T>::reference
ICR::ICA::Moments<T>::operator=(parameter other)
{    
  if (this!= &other) {
    m_data = (other.m_data);
    // m_mutex non-copiable  
  }
  return *this;
}

template<class T>  
inline  
typename ICR::ICA::Moments<T>::iterator
ICR::ICA::Moments<T>::begin()
{
  return iterator(this,0);
}

template<class T>
inline
typename ICR::ICA::Moments<T>::const_iterator
ICR::ICA::Moments<T>::begin() const
{
  return const_iterator (this,0);
}

template<class T>
inline
typename ICR::ICA::Moments<T>::iterator
ICR::ICA::Moments<T>::end()
{
  return iterator(this,m_data.size());
}


template<class T>
inline
typename ICR::ICA::Moments<T>::const_iterator
ICR::ICA::Moments<T>::end() const
{
  return const_iterator(this,m_data.size());
}


template<class T>
inline
typename ICR::ICA::Moments<T>::data_const_reference
ICR::ICA::Moments<T>::operator[](size_parameter index) const
{
  Lock lock(m_mutex); 
  return m_data[index];
}
      
template<class T>
inline
typename ICR::ICA::Moments<T>::data_reference
ICR::ICA::Moments<T>::operator[](size_parameter index)
{
  Lock lock(m_mutex); 
  return m_data[index];
}
      

template<class T>
typename ICR::ICA::Moments<T>::reference  
ICR::ICA::Moments<T>::operator+=(parameter other)
{
  //This could be called by different threads, so lock
  Lock lock(m_mutex);  
  PARALLEL_TRANSFORM(m_data.begin(), m_data.end(), other.begin(), m_data.begin(), plus());
  return *this;
}
      
template<class T>  
inline 
typename ICR::ICA::Moments<T>::reference
ICR::ICA::Moments<T>::operator*=(parameter other)
{
  //This could be called by different threads, so lock
  Lock lock(m_mutex); 
  PARALLEL_TRANSFORM(m_data.begin(), m_data.end(), other.begin(), m_data.begin(), times() );
  return *this;
}

template<class T>  
inline 
typename ICR::ICA::Moments<T>::reference
ICR::ICA::Moments<T>::operator*=(data_parameter d)
{
  //This could be called by different threads, so lock
  Lock lock(m_mutex); 
  PARALLEL_TRANSFORM(m_data.begin(), m_data.end(),  m_data.begin(), times_by(d));
  return *this;
}
      
   

