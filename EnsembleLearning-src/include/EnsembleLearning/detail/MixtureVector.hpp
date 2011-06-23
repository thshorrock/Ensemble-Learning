
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


#include <iostream>
#include <typeinfo>


#include "EnsembleLearning/node/variable/Hidden.hpp"
#include "EnsembleLearning/detail/TypeList.hpp"

#include <boost/mpl/vector.hpp>
#include <boost/mpl/push_back.hpp>
#include <boost/mpl/back.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/type_traits.hpp>
#include <boost/mpl/at.hpp> 
#include <boost/mpl/transform.hpp> 

#include <boost/tuple/tuple.hpp> 
#include <boost/mpl/copy.hpp> 
#include <boost/type_traits/is_same.hpp> 
#include <boost/preprocessor/repetition/enum_params.hpp> 


#include <boost/fusion/adapted/boost_tuple.hpp>
#include <boost/fusion/include/boost_tuple.hpp>
#include <boost/fusion/container/vector.hpp>
#include <boost/fusion/include/vector.hpp>
#include <boost/fusion/container/generation/make_vector.hpp>
#include <boost/fusion/include/make_vector.hpp>
#include <boost/fusion/sequence/intrinsic/at.hpp>
#include <boost/fusion/include/at.hpp>

#ifndef ENSEMBLE_LEARNING_COMPONENTS
#define   ENSEMBLE_LEARNING_COMPONENTS 5
#endif
#define ENSEMBLE_LEARNING_COMPONENTS_MINUS_1 BOOST_PP_SUB(ENSEMBLE_LEARNING_COMPONENTS,1)

namespace ICR{
  namespace EnsembleLearning{
      
    template<class T, class Tuple> 
    struct tuple_push_front; 
    template<class T, BOOST_PP_ENUM_PARAMS(ENSEMBLE_LEARNING_COMPONENTS, class T)> 
    struct tuple_push_front<T, boost::fusion::vector<BOOST_PP_ENUM_PARAMS(ENSEMBLE_LEARNING_COMPONENTS, T)> > { 
      typedef boost::fusion::vector<T, BOOST_PP_ENUM_PARAMS(ENSEMBLE_LEARNING_COMPONENTS_MINUS_1, T)> type; 
    }; 
    template<class Sequence> 
    struct make_my_tuple { 
      typedef typename boost::mpl::reverse_copy< 
	Sequence, boost::mpl::inserter<boost::fusion::vector<>, 
				       tuple_push_front<boost::mpl::_2, boost::mpl::_1> > 
	>::type type; 
}; 
      
template<template<class> class Model, class T, int I, class D = detail::TypeList::component>
class MixtureVector_impl
{
  template< template<class> class Model2, class T2,int I2, class D2>
  friend class  MixtureVector_impl;

  typedef typename boost::mpl::push_back<typename MixtureVector_impl<Model,T,I-1,D>::result_t::type, 
				 HiddenNode<Model,T,typename detail::TypeList::incr_component<D,I>::type > >::type result_t;
typedef  boost::mpl::transform<typename result_t::type ,boost::add_pointer<boost::mpl::_1> > pointer_t;
typedef make_my_tuple<typename pointer_t::type> data_base;
  struct data_t : data_base::type
{  };
private:

    data_t m_data;
public:

  typedef  result_t type;
  typedef  pointer_t pointer_type;

  const data_t& 
  data() const {return m_data;}
   data_t& 
  data() {return m_data;}
  
template<int index>
struct
get_t 
{
  typedef typename  boost::mpl::at<typename result_t::type,typename boost::mpl::int_<index>::type >::type type;
};



  template<int index>
  typename boost::add_pointer<typename get_t<index>::type>::type&
  get()
  {
    return boost::fusion::at_c<index>(m_data);
  }
  template<int index>
  const typename  boost::add_pointer<typename get_t<index>::type>::type& 
  get() const
  {
    return boost::fusion::at_c<index>(m_data);
  }

};


template<template<class> class Model, class T,class D>
class MixtureVector_impl<Model,T,0,D>
{
public:
  typedef  boost::mpl::vector<HiddenNode<Model,T,typename detail::TypeList::incr_component<D,0>::type> > result_t;
  typedef  boost::mpl::transform<typename result_t::type ,boost::add_pointer<boost::mpl::_1> > pointer_t;

  typedef make_my_tuple<typename pointer_t::type> data_t;
  
  struct type : data_t::type 
  { };  
  template<int index>
  struct
  get_t 
  {
   typedef typename   boost::mpl::at<typename result_t::type,typename boost::mpl::int_<index>::type >::type type;
  };
};

template<template<class> class Model, class T>
class MixtureVector_impl<Model,T,0,detail::TypeList::component>
{
public:
  typedef  boost::mpl::vector<HiddenNode<Model,T,typename detail::TypeList::incr_component<detail::TypeList::component,0>::type> > result_t;
  typedef  boost::mpl::transform<typename result_t::type ,boost::add_pointer<boost::mpl::_1> > pointer_t;

  typedef make_my_tuple<typename pointer_t::type> data_t;
  
  struct type : data_t::type 
  { };  
  template<int index>
  struct
  get_t
  {
    typedef typename    boost::mpl::at<typename result_t::type,typename boost::mpl::int_<index>::type >::type type;
  };
};


// template<template<class> class Model, class T,int I>
// class MixtureVector_impl
// {
// public:
// 	typedef  boost::mpl::push_back<typename MixtureVector_impl<Model,T,I-1>::result_t::type, HiddenNode<Model,T,detail::MixtureList <typename boost::mpl::back<typename MixtureVector_impl<Model,T,I-1>::result_t::type>::type> > > result_t;
	
// 	typedef  boost::mpl::transform<typename result_t::type ,boost::add_pointer<boost::mpl::_1> > pointer_t;
	
// 	typedef make_my_tuple<typename pointer_t::type> data_t;

// 	struct type : data_t::type 
// 	{
// 	};

  
// 	template<int index>
// 	struct
// 	get_t : boost::mpl::at<typename result_t::type,typename boost::mpl::int_<index>::type >
// 	{};

// };


// template<template<class> class Model, class T>
// class MixtureVector_impl<Model,T,1>
// {
// public:
// 	typedef  boost::mpl::vector<HiddenNode<Model,T> > result_t;
// 	typedef  boost::mpl::transform<typename result_t::type ,boost::add_pointer<boost::mpl::_1> > pointer_t;
	
// 	typedef  make_my_tuple<typename pointer_t::type> data_t;
  
// 	struct type : data_t::type 
// 	{
// 	};  
// 	template<int index>
// 	struct
// 	get_t : boost::mpl::at<typename result_t::type,typename boost::mpl::int_<index>::type >::type
// 	{};
// };

template<template<class> class Model, class T>
struct MixtureVector : public MixtureVector_impl<Model, T, ENSEMBLE_LEARNING_COMPONENTS - 1,detail::TypeList::component>
{};

    
template<template<class> class Model, class T,class D>
int
get_id(HiddenNode<Model,T,D>* a)
{
  return detail::TypeList::get_id_t<D>::type::value;
}
 
template<template<class> class Model, class T,class D>
int
get_component(HiddenNode<Model,T,D>* a)
{
  return detail::TypeList::get_component_t<D>::type::value;
}

template<template<class> class Model, class T,class D>
int
get_position(HiddenNode<Model,T,D>* a)
{
  return detail::TypeList::get_component_t<D>::type::value;
}
}
}
#endif  // guard for MIXTUREVECTOR_HPP
