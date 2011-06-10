

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
#ifndef INPUT_HPP
#define INPUT_HPP




#include "EnsembleLearning/node/variable/Hidden.hpp"
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
    template<template<class> class Model, class T>
    void
    SetMean(HiddenNode<Model,T>* node, double v)
    {
      node->SetMean(std::vector<T>(1, v));
    }
    /** Evaluate the standard deviation of a variable node. 
     * @param node A pointer to the VariableNode to evaluate.
     * @param index The index of the mean to return.
     *   Typically there is only one mean in the node (the index defaults to zero),
     *   however, Dirichlet models contain a mean for every component of the mixture.
     * @return The variance for the given node and index.
     * @ingroup UserInterface
     */
    template<template<class> class Model, class T>
    void
    SetStandardDeviation(HiddenNode<Model,T>* node,double v)
    {
      node->SetVariance(std::vector<T>(1, std::pow(v,2)));
    }
    
  }
}

#endif  // guard for INPUT_HPP
