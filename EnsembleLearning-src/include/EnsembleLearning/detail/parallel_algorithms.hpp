#pragma once
#ifndef PARALLEL_ALGORITHMS_HPP
#define PARALLEL_ALGORITHMS_HPP


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




//The msvc versions of these definitions have not been tested, should be easy enough to fix up though if your having problems.

//Node that you often need extra compile time directives .e.g. -openmp, to make these work

#ifdef __GNUC__
#include <parallel/algorithm> 
#include <algorithm> 
#elif _MSC_VER
#include <ppl.h>
#else  
#include <algorithm>
#endif

#ifdef __GNUC__
#include <parallel/numeric>
#include <numeric> 
#elif _MSC_VER
#include <ppl.h>
#else  
#include <numeric>
#endif

#include "AlgorithmChecker.hpp"

#ifdef __GNUC__
#define  PARALLEL_FOREACH __gnu_parallel::for_each  //std::for_each
//#define  PARALLEL_FOREACH  std::for_each
#elif _MSC_VER
#define  PARALLEL_FOREACH parallel_for_each //are these in std?
#else  
#define  PARALLEL_FOREACH  std::for_each
#endif

#ifdef __GNUC__
#define  PARALLEL_FIND  __gnu_parallel::find // std::find
//#define  PARALLEL_FIND  std::find
#elif _MSC_VER
#define  PARALLEL_FIND parallel_find //are these in std?
#else  
#define  PARALLEL_FIND  std::find
#endif

#ifdef __GNUC__
#define  PARALLEL_GENERATE __gnu_parallel::generate //std::generate
//#define  PARALLEL_GENERATE  std::generate
#elif _MSC_VER
#define  PARALLEL_GENERATE parallel_generate //are these in std?
#else  
#define  PARALLEL_GENERATE  std::generate
#endif


#ifdef __GNUC__
#define  PARALLEL_COPY  std::copy
#elif _MSC_VER
#define  PARALLEL_COPY parallel_copy //are these in std?
#else  
#define  PARALLEL_COPY  std::copy
#endif


#ifdef __GNUC__
#define  PARALLEL_ACCUMULATE __gnu_parallel::accumulate //std::accumulate_checked
//#define  PARALLEL_ACCUMULATE  std::accumulate_checked
#elif _MSC_VER
#define  PARALLEL_ACCUMULATE parallel_accumulate
#else  
#define  PARALLEL_ACCUMULATE  std::accumulate_checked
#endif

#ifdef __GNUC__
#define  PARALLEL_TRANSFORM __gnu_parallel::transform //std::transform
//#define  PARALLEL_TRANSFORM  std::transform
#elif _MSC_VER
#define  PARALLEL_TRANSFORM parallel_transform
#else  
#define  PARALLEL_TRANSFORM  std::transform
#endif

#ifdef __GNUC__
#define  PARALLEL_MAX __gnu_parallel::max_element  //std::max_element 
//#define  PARALLEL_MAX  std::max_element
#elif _MSC_VER
#define  PARALLEL_MAX parallel_max_element
#else  
#define  PARALLEL_MAX  std::max_element
#endif

#ifdef __GNUC__
#define  PARALLEL_INNERPRODUCT __gnu_parallel::inner_product  //std::max_element 
//#define  PARALLEL_INNERPRODUCT  std::inner_product
#elif _MSC_VER
#define  PARALLEL_INNERPRODUCT parallel_inner_product
#else  
#define  PARALLEL_INNERPRODUCT  std::inner_product
#endif


#endif  // guard for PARALLEL_ALGORITHMS_HPP
