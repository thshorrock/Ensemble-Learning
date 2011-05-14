#pragma once
#include <boost/thread.hpp>
#include <boost/thread/mutex.hpp>
#include <boost/iterator/iterator_facade.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>

namespace ICR{
  namespace ICA{
    
    //forward declaration
    template<class> class Moments;
    
    /**Threadsafe iterator for the Moments class. */ 
    template<class moments, class T>
    class  MomentsIterator
      : public boost::iterator_facade<MomentsIterator<moments,T>, T, boost::random_access_traversal_tag >
    {
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

      MomentsIterator()
	: m_Moments(0),
	  m_index(0)
	  //	  m_mutex() 
      {}
      
      explicit MomentsIterator(moments* p,  typename F::difference_type  index)
	: m_Moments(p),
	  m_index(index)
	  // m_mutex() 
      {}
      
      MomentsIterator(
      		      const MomentsIterator<moments,T>& other
      		      )
      : m_Moments(other.m_Moments),
       	m_index(other.m_index),
       	m_mutex() //don't copy the mutex (non-copiable)
      {}
      

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

      // template <class Othermoments, class other_data_type>
      // typename boost::disable_if<
      // 	boost::is_same<Othermoments*,moments*>
      // 	, MomentsIterator<moments,T>&
      // 	>::type
      // operator=(
      // 		const MomentsIterator<Othermoments, other_data_type>& other
      // 		)
      // {
      // 	//certainly not same pointer (they are not even the same type!)
      // 	m_Moments
      // 	m_Moments = other.m_Moments;
      // 	m_index   = other.m_index;
      // 	//don't copy mutex;
      // 	return *this;
      // }

    private:
      /** Increment the iterator */
      void 
      increment() 
      {
	//	boost::lock_guard<boost::mutex> lock(m_mutex); 
	++m_index;
      };
      /** Decrement the iterator */
      void 
      decrement() 
      { 
	//boost::lock_guard<boost::mutex> lock(m_mutex); 
	--m_index; 
      };
      
      /** Increment by a distance
       *  @param n The distance to increment
       */
      void 
      advance(const typename F::difference_type n)
      {
	//boost::lock_guard<boost::mutex> lock(m_mutex); 
	m_index+=n;
      };
      
      /** Find the distance to another iterator
       *  @param n The distance to increment
       */
      template < class Othermoments,class other_type>
      typename F::difference_type 
      distance_to(const MomentsIterator<Othermoments,other_type>& other   ) const 
      {
	//boost::lock_guard<boost::mutex> lock(m_mutex); 
	return other.m_index - m_index;
      };
    
      template <class Othermoments, class other_type>
      bool 
      equal(const MomentsIterator<Othermoments,other_type>& other) const
      {
      	//boost::lock_guard<boost::mutex> lock(m_mutex); 
      	return (other->m_index - m_index && m_Moments==other.m_Moments);
      }
      
      data_reference dereference() const
      { 
	//boost::lock_guard<boost::mutex> lock(m_mutex);
	return m_Moments->operator[](m_index);
      };

      data_reference dereference() 
      { 
	//boost::lock_guard<boost::mutex> lock(m_mutex);
	return m_Moments->operator[](m_index);
      };


      friend class boost::iterator_core_access;
      template <class, class> friend class MomentsIterator;
      moments* m_Moments;
      mutable  typename F::difference_type  m_index;
      mutable boost::mutex m_mutex;
    };
    
    
    
  }
}
