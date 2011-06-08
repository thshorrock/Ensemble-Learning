
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
#ifndef OBSERVED_HPP
#define OBSERVED_HPP


#include "EnsembleLearning/message/Moments.hpp"
#include "EnsembleLearning/node/Node.hpp"
#include "EnsembleLearning/exponential_model/Gaussian.hpp"
#include "EnsembleLearning/exponential_model/RectifiedGaussian.hpp"
#include "EnsembleLearning/exponential_model/Gamma.hpp"
#include "EnsembleLearning/exponential_model/Discrete.hpp"
#include "EnsembleLearning/exponential_model/Dirichlet.hpp"



#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/logical.hpp>
#include <boost/call_traits.hpp>

namespace ICR{
  namespace EnsembleLearning{

    

    namespace detail{
      
    template<template<class> class Model,class T>
    struct is_observable
      : boost::mpl::or_<boost::is_same<Model<T>,Dirichlet<T> >,
			boost::is_same<Model<T>,Gamma<T> >,
			boost::is_same<Model<T>,RectifiedGaussian<T> >,
			boost::is_same<Model<T>,Gaussian<T> >
			>
    {};
      
    //forward declare
    template<template<class> class Model, class T>
    struct GetMean_impl;

    template<template<class> class Model, class T>
    struct GetVariance_impl;
    
    }

    /** An observed node.
     *  A node that contains data that is a known constant.
     *  Examples include data nodes (where the data is experimentally determined)
     *   and prior knowledge - the mean of the data is 4.0, say).
     *  @tparam Model The model catagory into which the data falls.
     *    For example, Data might be modelled by a Gaussian distribution.
     *    The Model may be either Dirichlet, Gamma, RectifiedGaussian or Gaussian.
     *    Attempting to call any other class will not compile.
     *  @tparam T The data type used in calculations - either float or double.
     */
    template<template<class> class Model, class T, class Enable = void>
    class ObservedNode ; //uninitialised
    
    /** An observed node.
     *  A node that contains data that is a known constant.
     *  Examples include data nodes (where the data is experimentally determined)
     *   and prior knowledge - the mean of the data is 4.0, say).
     *  @tparam Model The model catagory into which the data falls.
     *    For example, Data might be modelled by a Gaussian distribution.
     *    The Model may be either Dirichlet, Gamma, RectifiedGaussian or Gaussian.
     *    Attempting to call any other class will not compile.
     *  @tparam T The data type used in calculations - either float or double.
     */
    template<template<class> class Model, class T>
    class ObservedNode<Model,T,
		       typename boost::enable_if<detail::is_observable<Model,T> >::type>
	: public VariableNode<T>
    {
    public:
      /** A Constructor.
       *  @param value The observed value of the node.
       *  This constructor is the default. 
       *  However, it is not appropriate, and therefore not available for  Discrete or Dirichlet models.
       */
      ObservedNode( const T& value 
	)
	: m_Moments(make_Moments(value, Model<T>() ) ), 
	  m_parent(0),
	  m_children()
      {}

      /** A Constructor.
       * @param  elements The number of elements in the observed node.
       * @param  value The value of each of the elements 
       *  This constructor is only available for Discrete and Dirichlet models.
       */
      ObservedNode(const size_t& elements, const T& value
		   )
	: m_Moments(make_Moments(elements,value, Model<T>() ) ), 
	  m_parent(0), 
	  m_children()
      {}
      
      
      void
      SetParentFactor(FactorNode<T>* f);

      void
      AddChildFactor(FactorNode<T>* f);
      
      void
      InitialiseMoments(){};

      const Moments<T>&
      GetMoments() ;
      
      const std::vector<T>
      GetMean() ;
      
      const std::vector<T>
      GetVariance() ;

      void 
      Iterate(Coster& C);
      
    private:
      friend struct detail::GetMean_impl<Model,T>;
      friend struct detail::GetVariance_impl< Model,T >;
      
      Moments<T>
      make_Moments(const T& d, const Gaussian<T> )
      {
	return Moments<T>(d, d*d);
      }
      Moments<T>
      make_Moments(const T& d, const RectifiedGaussian<T>)
      {
	return Moments<T>(d, d*d);
      }
      Moments<T>
      make_Moments(const T& d, const Gamma<T>)
      {
	return Moments<T>(d, std::log(d));
      }
      Moments<T>
      make_Moments(const size_t s,const T& d, const Dirichlet<T>)
      {
	return Moments<T>(std::vector<T>(s,d));
      }
      const Moments<T> m_Moments;
      FactorNode<T>* m_parent;
      std::vector<FactorNode<T>*> m_children;
    };
    
    namespace detail{
      
      /** Specialised structure to get the mean of the observed node.
       *  @tparam Model The model.
       *  @tparam T The data type (float or double)
       */
      template<template<class> class Model, class T>
      struct GetMean_impl
      {
	/** Get the mean.
	 *  @param t The node for which to obtain the mean
	 *  @return The vector of means (size of 1)
	 */
	static const std::vector<T>
	GetMean(ObservedNode<Model,T>& t){
	  return std::vector<T>(1, t.m_Moments[0]);
	}
      };

      /** Specialised structure to get the mean of Dirichlet model
       *  @tparam Model The model.
       *  @tparam T The data type (float or double)
       */
      template<class T>
      struct GetMean_impl<Dirichlet,T>
      {
	/** Get the mean.
	 *  @param t The Dirichlet node for which to obtain the mean
	 *  @return The vector of means (size of 1)
	 */
	static
	const std::vector<T>
	GetMean(ObservedNode<Dirichlet,T>& t){
	  return std::vector<T>(t.m_Moments.size(), t.m_Moments[0]);
	}
      };
      
      /** Specialised structure to get the variance of the observed node.
       *  @tparam Model The model.
       *  @tparam T The data type (float or double)
       */
      template<template<class> class Model, class T>
      struct GetVariance_impl
      {
	/** Get the variance.
	 *  @param t The  node for which to obtain the variance
	 *  @return The vector of means (size of 1)
	 */
	static  
	const std::vector<T>
	GetVariance(ObservedNode<Model,T>& t){
	  return std::vector<T>(1, 0);
	}
      };

      /** Specialised structure to get the variance of Dirichlet model
       *  @tparam Model The model.
       *  @tparam T The data type (float or double)
       */
      template<class T>
      struct GetVariance_impl<Dirichlet,T>
      {
	/** Get the variance.
	 *  @param t The Dirichlet node for which to obtain the variance
	 *  @return The vector of means (size of 1)
	 */
	static 
	const std::vector<T>
	GetVariance(ObservedNode<Dirichlet,T>& t){
	  return  std::vector<T>(t.m_Moments.size(), 0);
	}
      };
      
    }
    

  } 
}

template<template<class> class Model,class T>
inline
const ICR::EnsembleLearning::Moments<T>&
ICR::EnsembleLearning::ObservedNode<Model,T,
		       typename boost::enable_if<ICR::EnsembleLearning::detail::is_observable<Model,T> >::type>
::GetMoments() 
{
  //Obvserved moments are not modified and so this is thead safe.
  return m_Moments;
}

template<template<class> class Model,class T>
inline
void
ICR::EnsembleLearning::ObservedNode<Model,T,
		       typename boost::enable_if<ICR::EnsembleLearning::detail::is_observable<Model,T> >::type>
::AddChildFactor(FactorNode<T>* f)
{
  //Could be many factors, potentially added with many threads,
  //therefore make the following critical
#pragma omp critical
  {
    m_children.push_back(f);
  }
}


template<template<class> class Model, class T>
inline
void
ICR::EnsembleLearning::ObservedNode<Model,T,
		       typename boost::enable_if<ICR::EnsembleLearning::detail::is_observable<Model,T> >::type>
::SetParentFactor(FactorNode<T>* f)
{
  //Only one factor: this should be called once (therefore thread safe)
  m_parent = f;
}
  

template<template<class> class Model,class T>
inline
const std::vector<T>
ICR::EnsembleLearning::ObservedNode<Model,T,
		       typename boost::enable_if<ICR::EnsembleLearning::detail::is_observable<Model,T> >::type>
::GetMean() 
{
  return detail::GetMean_impl<Model,T>::GetMean(*this);
}
   
template<template<class> class Model,class T>
inline
const std::vector<T>
ICR::EnsembleLearning::ObservedNode<Model,T,
		       typename boost::enable_if<ICR::EnsembleLearning::detail::is_observable<Model,T> >::type>
::GetVariance() 
{
  return detail::GetVariance_impl<Model,T>::GetVariance(*this);
}

template<template<class> class Model, class T>
inline
void
ICR::EnsembleLearning::ObservedNode<Model,T,
		       typename boost::enable_if<ICR::EnsembleLearning::detail::is_observable<Model,T> >::type>
::Iterate(Coster& Total)
{
  //Constant nodes have no parents (and contribute nothing to the cost)
  if (m_parent!=0)
    {
      //This is a data node...
      //Assume thead-safety of other nodes so do not need to worry here.
      const NaturalParameters<T> ParentNP =m_parent->GetNaturalNot(this); 
      //See page 41 of Winn's thesis for this formula
      const T Cost = ParentNP*GetMoments() +  m_parent->CalcLogNorm();
      Total += Cost;
    }
}
#endif  // guard for OBSERVED_HPP
