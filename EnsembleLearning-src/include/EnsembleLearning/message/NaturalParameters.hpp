#pragma once
#include <iostream>

#include "Moments.hpp"

#include<vector>


#include "EnsembleLearning/detail/parallel_algorithms.hpp"

namespace ICR{
  namespace EnsembleLearning{

    /** A container for the Natural Paramemeters.
     *  @tparam T The data type to be used.
     *  This is intended to be either float or double.
     *  @attention Natural Parameters is not thread safe - it is intended to be a temporary container for passing messages between nodes.
     */
    template<class T>
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

      typedef typename boost::call_traits<NaturalParameters<T> >::param_type 
      parameter;
      
      typedef typename boost::call_traits<NaturalParameters<T> >::reference  
      reference;
      
      typedef typename boost::call_traits<NaturalParameters<T> >::const_reference  
      const_reference;
      
      typedef typename boost::call_traits<NaturalParameters<T> >::value_type  
      type;
      
      typedef typename boost::call_traits<size_t>::param_type  
      size_parameter;
      
      typedef typename boost::call_traits<size_t>::value_type  
      size_type;
      
   
      typedef typename std::vector<T>::iterator
      iterator;
      
      typedef typename std::vector<T>::const_iterator
      const_iterator;

      ///@}

      /** Constructor.
       *  @param size The number of elements to be held in the container.
       */
      NaturalParameters(size_parameter size =0 );
      
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
      template<class U>
      friend
      U
      operator*(const NaturalParameters<U>&  a, 
		const Moments<U>& b);
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

      std::vector<data_type> m_data;
      

    };

    /** Output the contents of the Natural Parameters Container.
     * @param out The output stream.
     * @param m The Natural Parameters to output.
     * @return A reference to the output stream.
     */
    template<class T>
    inline
    std::ostream&
    operator<<(std::ostream& out, const NaturalParameters<T>& NP)
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
    template<class T>
    inline
    NaturalParameters<T>
    operator+(const NaturalParameters<T>& a,
	      const NaturalParameters<T>& b)
    {
      NaturalParameters<T> tmp = a;
      return tmp+=b;
    }  
    
    /** Subtract two Natural Parameters containers.
     *  @param a The first Natural Parameters  container.
     *  @param b The second Natural Parameters container.
     *  @return The result of the subtraction.
     */
    template<class T>
    inline
    NaturalParameters<T>
    operator-(const NaturalParameters<T>& a,
	      const NaturalParameters<T>& b)
    {
      ICR::EnsembleLearning::NaturalParameters<T> tmp = a;
      return tmp-=b;
    }  

    
 
      /** Multiply a Natural Paramers container by a Moments Container.
       *  This is an inner-product like operation.
       * @param a The Natural Parameters container.
       * @param b The Moments container.
       * @return The scalar product from the two containers.
       */

    template<class T>
    inline
    T
    operator*(const NaturalParameters<T>&  a, 
	      const Moments<T>& b)
    { 
      
      std::vector<T> prod(a.size());
      PARALLEL_TRANSFORM(a.begin(), a.end(), b.begin(), prod.begin(), 
			 typename NaturalParameters<T>::times());
      T sum = 0.0;
      return PARALLEL_ACCUMULATE(prod.begin(), prod.end(), sum);
    }  
      

    /** Multiply a Natural Paramers container by a Moments Container.
     *  This is an inner-product like operation.
     * @param p The Natural Parameters container.
     * @param m The Moments container.
     * @return The scalar product from the two containers.
     */
    template<class T>
    inline
    T
    operator*(const Moments<T>& m,const NaturalParameters<T>&  p)
    {
      return p*m;
    }  

    /** Multiply a Natural Parameters container by a scalar.
     * @param n The Natural Parameters container.
     * @param d The scalar.
     * @return The product of the natural parameters and the scalar.
     */
    template<class T>
    inline
    const NaturalParameters<T>
    operator*(const NaturalParameters<T>& n,const T  d)
    {
      NaturalParameters<T> tmp(n);
      return tmp*=d;
    }  



  }
}

template<class T>
inline
ICR::EnsembleLearning::NaturalParameters<T>::NaturalParameters(size_parameter size)
  : //m_mutex_ptr(new boost::mutex), 
  m_data(size)
{
}

template<class T>
inline
ICR::EnsembleLearning::NaturalParameters<T>::NaturalParameters(vector_parameter NP)
  : //m_mutex_ptr(new boost::mutex), 
    m_data(NP)
{
}

template<class T>
inline
ICR::EnsembleLearning::NaturalParameters<T>::NaturalParameters(data_parameter d1,
						  data_parameter d2 )
//  : 
  //  m_mutex_ptr(new boost::mutex)
  : m_data(2)
{
  //boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  m_data[0] = d1;
  m_data[1] = d2;
}

template<class T>
typename ICR::EnsembleLearning::NaturalParameters<T>::iterator
ICR::EnsembleLearning::NaturalParameters<T>::begin()
{
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.begin();
}


template<class T>
typename ICR::EnsembleLearning::NaturalParameters<T>::const_iterator
ICR::EnsembleLearning::NaturalParameters<T>::begin() const
{
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.begin();
}


template<class T>
typename ICR::EnsembleLearning::NaturalParameters<T>::iterator
ICR::EnsembleLearning::NaturalParameters<T>::end()
{
  // boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.end();
}


template<class T>
typename ICR::EnsembleLearning::NaturalParameters<T>::const_iterator
ICR::EnsembleLearning::NaturalParameters<T>::end() const
{
  // boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data.end();
}


 

template<class T>
inline
typename ICR::EnsembleLearning::NaturalParameters<T>::data_const_reference 
ICR::EnsembleLearning::NaturalParameters<T>::operator[](size_parameter i) const
{
  //boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data[i];
}
   
template<class T>
inline
typename ICR::EnsembleLearning::NaturalParameters<T>::data_reference
ICR::EnsembleLearning::NaturalParameters<T>::operator[](size_parameter i) 
{
  //boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  
  return m_data[i];
}
   


template<class T>
inline
typename ICR::EnsembleLearning::NaturalParameters<T>::reference
ICR::EnsembleLearning::NaturalParameters<T>::operator+=(parameter other)
{
  //This could be called by different threads, so lock
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  //DO KEEP THIS ONE!
  //  std::cout<<"size = "<<m_data.size()<<"other = "<<m_data.size<<std::endl;

  PARALLEL_TRANSFORM(m_data.begin(), m_data.end(), other.begin(), m_data.begin(), plus());
  // for(size_t i=0;i<other.size();++i){
  //   m_data[i]+=other[i];
  // }

  return *this;
}
   

template<class T>  
inline 
typename ICR::EnsembleLearning::NaturalParameters<T>::reference
ICR::EnsembleLearning::NaturalParameters<T>::operator-=(parameter other)
{
  //This could be called by different threads, so lock
  //  boost::lock_guard<boost::mutex> lock(*m_mutex_ptr);  //DO KEEP THIS ONE!
  PARALLEL_TRANSFORM(m_data.begin(), m_data.end(), other.begin(), m_data.begin(), minus());
  return *this;
}
  
template<class T>
inline
typename ICR::EnsembleLearning::NaturalParameters<T>::reference
ICR::EnsembleLearning::NaturalParameters<T>::operator*=(data_parameter other)
{
  
  PARALLEL_TRANSFORM(m_data.begin(), m_data.end(), m_data.begin(), times_by(other));
  return *this;
}    



