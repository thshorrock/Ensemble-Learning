
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
#ifndef OUTPUT_HPP
#define OUTPUT_HPP


#include "EnsembleLearning/node/Node.hpp"
#include "EnsembleLearning/detail/MixtureVector.hpp"
#include <boost/fusion/include/for_each.hpp>


#include <cmath>  
namespace ICR{
  namespace EnsembleLearning{
    
    /** Evaluate the mean of a variable node. 
     * @param node A pointer to the VariableNode to evaluate.
     * @param index The index of the mean to return.
     *   Typically there is only one mean in the node (the index defaults to zero),
     *   however, Dirichlet models contain a mean for every component of the mixture.
     * @return The mean for the given node and index.
     * @ingroup UserInterface
     */
    template<class T,int array_size>
    T
    Mean(VariableNode<T,array_size>* node, size_t index = 0)
    {
      return node->GetMean()[index];
    }

    

    template<template<template<class> class,class,int> class VectorType,template<class> class Model, class T,int array_size>
    std::vector<VariableNode<T>* >
    to_std_vector(VectorType<Model,T,array_size> v)
    {
      std::vector<VariableNode<T>* > ret;
  
      boost::fusion::for_each(v.data(), 
			      boost::bind(&std::vector<VariableNode<T>* >::push_back,
					  boost::ref(ret),
					  (_1) 
					  ));
      return ret;
			  
    }
    

    template<template<template<class> class,class,int,int> class VectorType,template<class> class Model, class T,int from, int array_size>
    std::vector<VariableNode<T>* >
    to_std_vector(VectorType<Model,T,from,array_size> v)
    {
      std::vector<VariableNode<T>* > ret;
  
      boost::fusion::for_each(v.data(), 
			      boost::bind(&std::vector<VariableNode<T>* >::push_back,
					  boost::ref(ret),
					  (_1) 
					  ));
      return ret;
			  
    }
    /** Evaluate the standard deviation of a variable node. 
     * @param node A pointer to the VariableNode to evaluate.
     * @param index The index of the mean to return.
     *   Typically there is only one mean in the node (the index defaults to zero),
     *   however, Dirichlet models contain a mean for every component of the mixture.
     * @return The variance for the given node and index.
     * @ingroup UserInterface
     */
    template<class T,int array_size>
    T
    StandardDeviation(VariableNode<T,array_size>* node, size_t index = 0)
    {
      return std::sqrt(node->GetVariance()[index]);
    }
    
  }
}


#endif //OUTPUT_HPP guard
