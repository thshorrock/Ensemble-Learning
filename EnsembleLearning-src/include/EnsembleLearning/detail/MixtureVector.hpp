
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
#include <boost/mpl/placeholders.hpp> 

#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <boost/fusion/include/at.hpp>
#include <boost/fusion/include/as_vector.hpp>
#include <boost/fusion/adapted/mpl.hpp>
#include <boost/fusion/include/mpl.hpp>

//The maximum number of components is given by FUSION_MAX_VECTOR_SIZE (default is 20)
//If we want to increase the number of componets beyond this, then need to increase this size.
#ifdef ENSEMBLE_LEARNING_MAX_COMPONENTS
#  define FUSION_MAX_VECTOR_SIZE ENSEMBLE_LEARNING_MAX_COMPONENTS
#endif
#ifndef ENSEMBLE_LEARNING_COMPONENTS
#  define ENSEMBLE_LEARNING_COMPONENTS 5
#endif
#ifndef ENSEMBLE_LEARNING_PLACEHOLDERS
#  define ENSEMBLE_LEARNING_PLACEHOLDERS 5
#endif
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
    template<template<template<class> class,class,class,int,class> class ModelType,
	     template<class> class Model, 
	     class T, 
	     int I, 
	     template<class, int> class Op = detail::TypeList::Identity,
	     class D = detail::TypeList::zeros>
    class Vector_impl
    {
      //Make all other ector_impl friends.
      template<template<template<class> class,class,class,int,class> class ModelType2,
	     template<class> class Model2, 
	     class T2, 
	     int I2, 
	     template<class, int> class Operation2,
	       class D2>
    friend class Vector_impl;

      //create an mpl vector that increments the component number for each type from 1 to #components.
      typedef typename boost::mpl::push_back<typename Vector_impl<ModelType,Model,T,I-1,Op,D>::base::type, ModelType<Model,T,typename Op<D,I>::type,2,void > >::type base;
      //Make the same vector but this time as pointers.
      typedef  boost::mpl::transform<typename base::type ,boost::add_pointer<boost::mpl::_1> > pointer_t;
      //The data is stored as a boost fusion vector. 
      //Obtain the fusion vector class from the mpl vector: base.
      typedef typename boost::fusion::result_of::as_vector<typename pointer_t::type> data_base;
      //Apply the fusion vector type to our data_t
      // struct data_t : data_base::type
      // {  };
      // // The data.
      // data_t 
      typedef typename  data_base::type data_t;
      typename data_base::type m_data;
    public:

      //@Convenient typedefs
      //@{
      typedef  base type;
      typedef  pointer_t pointer_type;
      typedef  data_base data_type;

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
    
    template<template<template<class> class,class,class,int,class> class ModelType,
	     template<class> class Model, 
	     class T, 
	     template<class, int> class Op,
	     class D>
    class Vector_impl<ModelType,Model,T,0,Op,D>
    {
      //Make all other ector_impl friends.
      template<template<template<class> class,class,class,int,class> class ModelType2,
	     template<class> class Model2, 
	     class T2, 
	     int I2, 
	     template<class, int> class Operation2,
	     class D2>
    friend class Vector_impl;

      //create an mpl vector that increments the component number for each type from 1 to #components.
      typedef  boost::mpl::vector<ModelType<Model,T,typename Op<D,0>::type,2,void> > base;
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
    struct MixtureVector : public Vector_impl<HiddenNode,Model, T, size-1 ,detail::TypeList::incr_component,detail::TypeList::component>
    {};


    template<template<class> class Model, class T, int size = ENSEMBLE_LEARNING_PLACEHOLDERS>
    struct CalculationVector : public Vector_impl<HiddenNode,Model, T, size-1 ,detail::TypeList::incr_position,detail::TypeList::position>
    {};

    template<template<template<class> class,class,class,int,class> class ModelType,
	     template<class> class Model, 
	     class T, 
	     int size, 
	     template<class, int> class Op = detail::TypeList::Identity,
	     class D = detail::TypeList::zeros
	     >
    struct Vector : public Vector_impl<ModelType,Model, T, size-1 ,Op, D>
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
