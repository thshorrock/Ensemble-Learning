
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



#pragma once
#ifndef MIXTUREVECTOR_HPP
#define MIXTUREVECTOR_HPP

#include "EnsembleLearning/node/variable/Hidden.hpp"
#include "EnsembleLearning/detail/TypeList.hpp"

#include <boost/mpl/vector.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/at.hpp> 
#include <boost/mpl/transform.hpp> 

#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/adapted/mpl.hpp>
#include <boost/fusion/include/mpl.hpp>

#ifndef ENSEMBLE_LEARNING_COMPONENTS
#define   ENSEMBLE_LEARNING_COMPONENTS 5
#endif
#define ENSEMBLE_LEARNING_COMPONENTS_MINUS_1 BOOST_PP_SUB(ENSEMBLE_LEARNING_COMPONENTS,1)

namespace ICR{
  namespace EnsembleLearning{
      
    
    /** Store a vector of components for a Mixture model.
     *  The types of each of the component are different,
     *  (They have a different component id)
     *  and so the mixture models cannot be stored in a regular vector.
     *  This class resolves the problem, 
     *  enabling the components to be stored together,
     *  while maintaining their type information.
     *  
     *  To obtain a std::vector use the Mixture2Vector function.
     *  @see Mixture2Vector.
     */
    template<template<class> class Model, class T, int I, class D = detail::TypeList::component>
    class MixtureVector_impl
    {
      //Make all other MixtureVector_impl friends.
      template< template<class> class Model2, class T2,int I2, class D2>
      friend class  MixtureVector_impl;

      //create an mpl vector that increments the component number for each type from 1 to #components.
      typedef typename boost::mpl::push_back<typename MixtureVector_impl<Model,T,I-1,D>::base::type, HiddenNode<Model,T,typename detail::TypeList::incr_component<D,I>::type > >::type base;
      //Make the same vector but this time as pointers.
      typedef  boost::mpl::transform<typename base::type ,boost::add_pointer<boost::mpl::_1> > pointer_t;
      //The data is stored as a boost fusion vector. 
      //Obtain the fusion vector class from the mpl vector: base.
      typedef typename boost::fusion::result_of::as_vector<typename pointer_t::type> data_base;
      //Apply the fusion vector type to our data_t
      struct data_t : data_base::type
      {  };
      // The data.
      data_t m_data;
    public:

      //@Convenient typedefs
      //@{
      typedef  base type;
      typedef  pointer_t pointer_type;

      //@}
      const data_t& 
      data() const {return m_data;}
      data_t& 
      data() {return m_data;}

      //conversion operator
      operator const data_t &() const  { return m_data;}
      operator data_t& () {return m_data;}
      

      /** Get the type for the object stored at the given index.
       *  @tparam index The index of the object whose type you want.
       */
      template<int index>
      struct
      get_t 
      {
	typedef typename  boost::mpl::at<typename base::type,typename boost::mpl::int_<index>::type >::type type;
      };

      /** Get the object stored at the given index.
       *  @tparam index The index of the object whose type you want.
       *  @return The object stored at the given index.
       */
      template<int index>
      typename boost::add_pointer<typename get_t<index>::type>::type&
      get()
      {
	return boost::fusion::at_c<index>(m_data);
      }
      

      /** Get the object stored at the given index.
       *  @tparam index The index of the object whose type you want.
       *  @return The object stored at the given index.
       */
      template<int index>
      const typename  boost::add_pointer<typename get_t<index>::type>::type& 
      get() const
      {
	return boost::fusion::at_c<index>(m_data);
      }

    };

    /** Specialisation for when the vector is empty.
     *  This terminates the recursive generation of a vector.
     */
    template<template<class> class Model, class T,class D>
    class MixtureVector_impl<Model,T,0,D>
    {
      //Make all other MixtureVector_impl friends.
      template< template<class> class Model2, class T2,int I2, class D2>
      friend class  MixtureVector_impl;
      
      //create an mpl vector that increments the component number for each type from 1 to #components.
      typedef  boost::mpl::vector<HiddenNode<Model,T,typename detail::TypeList::incr_component<D,0>::type> > base;
      //Make the same vector but this time as pointers.
      typedef  boost::mpl::transform<typename base::type ,boost::add_pointer<boost::mpl::_1> > pointer_t;
      //The data is stored as a boost fusion vector. 
      //Obtain the fusion vector class from the mpl vector: base.
      typedef typename boost::fusion::result_of::as_vector<typename pointer_t::type> data_base;
      //Apply the fusion vector type to our data_t
      struct data_t : data_base::type 
      { };  
      
    public:
      //@Convenient typedefs
      //@{
      typedef  base type;
      typedef  pointer_t pointer_type;
      //@}

      /** Get the type for the object stored at the given index.
       *  @tparam index The index of the object whose type you want.
       */
      template<int index>
      struct
      get_t 
      {
	/** The type. */
	typedef typename   boost::mpl::at<typename base::type,typename boost::mpl::int_<index>::type >::type type;
      };
    };


    template<template<class> class Model, class T, int size = ENSEMBLE_LEARNING_COMPONENTS>
    struct MixtureVector : public MixtureVector_impl<Model, T, size-1 ,detail::TypeList::component>
{};

    /** Get the id for a particular class */
template<template<class> class Model, class T,class D>
int
get_id(HiddenNode<Model,T,D>* a)
{
  return detail::TypeList::get_id_t<D>::type::value;
}
 /** Get the component number for a particular class */
template<template<class> class Model, class T,class D>
int
get_component(HiddenNode<Model,T,D>* a)
{
  return detail::TypeList::get_component_t<D>::type::value;
}

 /** Get the positional number for a particular class */
template<template<class> class Model, class T,class D>
int
get_position(HiddenNode<Model,T,D>* a)
{
  return detail::TypeList::get_component_t<D>::type::value;
}
}
}
#endif  // guard for MIXTUREVECTOR_HPP
