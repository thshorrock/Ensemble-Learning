
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
#ifndef TYPELIST_HPP
#define TYPELIST_HPP

#include <boost/mpl/vector.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/at.hpp> 
#include <boost/mpl/plus.hpp> 
#include <boost/mpl/times.hpp> 
#include <boost/mpl/apply_wrap.hpp>

namespace ICR{
  namespace EnsembleLearning{

    namespace detail{
      namespace TypeList {
	typedef boost::mpl::vector_c<int,0,0,0> zeros;
	typedef boost::mpl::vector_c<int,1,0,0> id;
	typedef boost::mpl::vector_c<int,0,1,0> component;
	typedef boost::mpl::vector_c<int,0,0,1> position;


	struct plus_f
	{
	  template <class T1, class T2>
	  struct apply 
	  { 
	    typedef boost::mpl::vector_c<int,
				  boost::mpl::plus<typename boost::mpl::at_c<T1,0>::type, typename boost::mpl::at_c<T2,0>::type>::type::value,
				  boost::mpl::plus<typename boost::mpl::at_c<T1,1>::type, typename boost::mpl::at_c<T2,1>::type>::type::value,
				  boost::mpl::plus<typename boost::mpl::at_c<T1,2>::type, typename boost::mpl::at_c<T2,2>::type>::type::value> type;
	};
      };


      struct times_f
      {
	template <class T1, class T2>
	struct apply 
	{ 
	  typedef boost::mpl::vector_c<int,
				boost::mpl::times<typename boost::mpl::at_c<T1,0>::type, typename boost::mpl::at_c<T2,0>::type>::type::value,
				boost::mpl::times<typename boost::mpl::at_c<T1,1>::type, typename boost::mpl::at_c<T2,1>::type>::type::value,
				boost::mpl::times<typename boost::mpl::at_c<T1,2>::type, typename boost::mpl::at_c<T2,2>::type>::type::value> type;
	};
      };


      struct scale
      {
	template <class D, class I>
	struct apply
	{
	  typedef boost::mpl::vector_c<int,I::value,I::value,I::value> sf;
	  typedef typename boost::mpl::apply_wrap2<times_f,D,sf>::type type;
	};
      };


      template<class D, int I=1>
      struct incr_id
      {
	typedef typename boost::mpl::apply_wrap2<plus_f,D,
					  typename boost::mpl::apply_wrap2<scale,id,typename boost::mpl::int_<I>::type>
					  ::type>::type type;
    };
    template<class D, int I=1>
    struct incr_position
    {
      typedef typename boost::mpl::apply_wrap2<plus_f,D,
					typename boost::mpl::apply_wrap2<scale,position,typename boost::mpl::int_<I>::type>
					::type>::type type;
  };
  template<class D, int I=1>
  struct incr_component
  {
    typedef typename boost::mpl::apply_wrap2<plus_f,D,
				      typename boost::mpl::apply_wrap2<scale,component,typename boost::mpl::int_<I>::type>
				      ::type>::type type;
  };
  template<class D, int I=1>
  struct Identity
  {
    typedef D type;
  };

  template <typename D>
  struct get_id_t
  {
    typedef typename boost::mpl::at<D,boost::mpl::int_<0>::type >::type type;
  };
	
  template <typename D>
  struct get_component_t
  {
    typedef typename boost::mpl::at<D,boost::mpl::int_<1>::type >::type type;
  };	
  template <typename D>
  struct get_position_t
  {
    typedef typename boost::mpl::at<D,boost::mpl::int_<2>::type >::type type;
  };
	
}
      // template <int id_t, int pos_t, int comp_t>
      // struct TypeList_Exact
      // {
      // 	enum
      // 	  {
      // 	    id        = id_t,
      // 	    position  = pos_t,
      // 	    component = comp_t,
      // 	  };
      // };
      
      // template <class Prev>
      // struct TypeList
      // {
      // 	enum
      // 	  {
      // 	    id        = (Prev::id) + 1,
      // 	    position  = (Prev::position),
      // 	    component = (Prev::component),
      // 	  };
      // };
      

      // template <>
      // struct TypeList<void>
      // {
      // 	enum
      // 	  {
      // 	    id        = 0,
      // 	    position  = 0,
      // 	    component = 0,
      // 	  };
      // };
      
      
      // template <class Prev>
      // struct MixtureList : TypeList<Prev>
      // {
      // public:
      // 	enum
      // 	  {
      // 	    id        = (Prev::id) + 1,
      // 	    position  = (Prev::position),
      // 	    component = (Prev::component) + 1,
      // 	  };
	
      // };
      
      // template <class Prev>
      // struct ContextList : TypeList<Prev>
      // {
      // public:
      // 	enum
      // 	  {
      // 	    id        = (Prev::id) + 1,
      // 	    position  = (Prev::position) + 1,
      // 	    component = (Prev::component),
      // 	  };
	
      // };

  }
  }
}
#endif  // guard for TYPELIST_HPP
