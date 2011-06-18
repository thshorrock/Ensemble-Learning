
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
#ifndef PLACEHOLDER_HPP
#define PLACEHOLDER_HPP
#include <boost/proto/proto.hpp>
#include <iostream>
#include <vector>
#include <cmath>
#include <cassert>

#include <boost/type_traits.hpp>
#include <boost/mpl/vector.hpp>
#include <boost/mpl/range_c.hpp>
#include <boost/mpl/if.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/placeholders.hpp>
#include <boost/mpl/transform.hpp>
#include <boost/mpl/accumulate.hpp>
#include <boost/mpl/vector_c.hpp>
#include <boost/mpl/plus.hpp>
#include <boost/mpl/times.hpp>
#include <boost/bind.hpp>
#include <algorithm>

#include <boost/fusion/include/transform.hpp>
#include <boost/fusion/include/for_each.hpp>
// Create some namespace aliases
namespace mpl = boost::mpl;
namespace fusion = boost::fusion;
namespace proto = boost::proto;



#include "FunctionsIterator.hpp"
//#include "Context.hpp"
#include "Expression.hpp"

#include <boost/assert.hpp>
#include <boost/call_traits.hpp>

#include <boost/mpl/int.hpp>
namespace ICR{

  namespace EnsembleLearning{
    

    /**  A placeholder to a VariableNode in a function.
     *  The purpose of the placeholder is to be able to form an equation that can be used in multiple contexts.
     *  A placeholder indicates where a particular node will be in the equation.
     *  It is then linked to the Variable by a Context.
     *
     *  @tparam T The type of the stored data.  
     *   This will be either double or float - with the latter used only for memory contrained models.
     *
     *  Example:
     *  @code 
     *
     *  //Create an expression X*Y + Z.
     *  ExpressionFactory<double> Factory;  // Create a Factory.
     *  
     *  Placeholder<double>* X = Factory.placeholder();
     *  Placeholder<double>* Y = Factory.placeholder();
     *  Placeholder<double>* Z = Factory.placeholder();
     *  Expression<double>* XY = Factory.Multiply(X,Y);
     *  Expression<double>* Expr = Factory.Add(XY,Z);
     *  
     *  //Create some variables
     *  Builder builder;
     *  Builder::GaussianNode x = Builder.Gaussian(0.01,0.01);
     *  Builder::GaussianNode y = Builder.Gaussian(0.01,0.01);
     *  Builder::GaussianNode z = Builder.Gaussian(0.01,0.01);
     *  
     *  //Assign this context to the expression.
     *  Context context;
     *  context.Assign(X,x);
     *  context.Assign(Y,y);
     *  context.Assign(Z,z);
     *  
     *  @endcode
     *
     *  @attention Not thread safe in initialisation.
     */
    template<int I>
    struct placeholder :  mpl::int_<I>
    {};
      
  }
}


#endif  // guard for PLACEHOLDER_HPP
