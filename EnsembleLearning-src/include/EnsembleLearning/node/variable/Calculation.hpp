#pragma once
#ifndef VARIABLE_CALCULATION_HPP
#define VARIABLE_CALCULATION_HPP


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



#include "EnsembleLearning/node/Node.hpp"
#include "EnsembleLearning/message/Moments.hpp"
#include "EnsembleLearning/message/NaturalParameters.hpp"
#include "EnsembleLearning/detail/Mutex.hpp"
#include "EnsembleLearning/detail/parallel_algorithms.hpp"
#include "EnsembleLearning/detail/TypeList.hpp"

#include <boost/assert.hpp> 
#include <boost/bind.hpp>
#include <vector>
#include <set>

namespace ICR{
  namespace EnsembleLearning{
    
    /** A Deterministic Node.
     *  DeterministicNode's are not inferred,
     *  rather they are staging posts to pass moments that are calculated 
     *  from other (hidden) nodes.
     *  The model used (for example Gaussian/Gamma) is passed by template parameter.
     *  @tparam Model  The model to use for the inferred data.
     *  @tparam T The data type (float or double)
     */
    template <template<class> class Model, class T,class List=detail::TypeList::zeros
	      , int array_size =2,class Enable = void >
    class DeterministicNode : public VariableNode<T, boost::mpl::int_<array_size>::value>//, public List
    {
      DeterministicNode(const DeterministicNode<Model,T,List,array_size,Enable>&); //non-copyable

      typedef typename boost::call_traits<FactorNode<T, DeterministicNode, typename boost::mpl::int_<array_size>::type >* >::param_type 
      F_parameter;
      
      typedef typename boost::call_traits<ParentFactorNode<T, DeterministicNode, typename boost::mpl::int_<array_size>::type>* >::param_type 
      ParentF_parameter;
      
      typedef typename boost::call_traits<FactorNode<T, DeterministicNode, typename boost::mpl::int_<array_size>::type >* >::value_type 
      F_t;
      
      typedef typename boost::call_traits<ParentFactorNode<T, DeterministicNode, typename boost::mpl::int_<array_size>::type>* >::value_type 
      ParentF_t;
      
      typedef typename boost::call_traits<Moments<T,array_size> >::value_type 
      moments_t;
    public:
      //mostly two, but variable for Discrete model
      /** A constructor.
       *  @param moment_size The number of elements in the stored moments.
       *  This is usually two but varies for discrete nodes.
       */
      DeterministicNode(); 
      ~DeterministicNode(); 

      void
      SetParentFactor(ParentF_parameter f);
      
      void
      AddChildFactor(F_parameter f);

      void 
      Iterate(Coster& C);

      const moments_t*
      GetMoments() ;

      //should make a const version of this 
      /** Forward moments to the deterministic function.
       *  @return The moments from the child node.
       */
      const moments_t
      GetForwardedMoments() ;

      /** The number of elements in the stored Moments */
      size_t 
      size() const {return m_Moments->size();}
      
      const std::vector<T>
      GetMean() ;
      
      const std::vector<T>
      GetVariance() ;

      void
      InitialiseMoments()
      {
	//before iteration stage - no chance for GetMoments to be called by another thread,
	// therefore thread considerations don't apply here.
	*m_Moments = m_parent->InitialiseMoments();
      }


      
    private:
      
      ParentF_t m_parent;
      std::vector<F_t> m_children;
      mutable moments_t* m_Moments;
      //store the moments that might still being called by other threads after moments updated, deleted at the beginning of each iteration.
      //mutable moments_t* m_oldMoments; 
      // std::set<moments_t*> m_MomentsFromLastIteration;
      // std::set<moments_t*> m_MomentsFromCurrentIteration;
      
      // mutable Mutex m_mutex;
    };

  }
}


template<template<class> class Model,class T,class List,int array_size,class Enable>
ICR::EnsembleLearning::DeterministicNode<Model,T,List,array_size,Enable>::DeterministicNode() 
  :   m_parent(0), m_children(), m_Moments(0)//, m_MomentsFromLastIteration(),m_MomentsFromCurrentIteration()
{
   m_Moments = new moments_t();
  // m_oldMoments = new moments_t(); //dummy, deleted first iteration
}


template<template<class> class Model,class T,class List,int array_size,class Enable>
inline 
void
ICR::EnsembleLearning::DeterministicNode<Model,T,List,array_size,Enable>::SetParentFactor(ParentF_parameter f)
{
  //This should only be called once, so should get no collisions here
  m_parent=f;
  // Can now initialise
  InitialiseMoments();
}
template<template<class> class Model,class T,class List,int array_size,class Enable>
ICR::EnsembleLearning::DeterministicNode<Model,T,List,array_size,Enable>::~DeterministicNode() 
{
  
  //std::cout<<"~"<<std::endl;
  //make sure not deleted twice.
  // m_MomentsFromCurrentIteration.erase(m_Moments);
  delete m_Moments;
  // for(typename std::set<moments_t*>::iterator it = m_MomentsFromCurrentIteration.begin();
  //     it!=m_MomentsFromCurrentIteration.end();
  //     ++it)
  //   {
  //     // std::cout<<*it<<std::endl;
  //     delete *it;
  //   }
  //oldMoments deleted in iteration.
}

   
template<template<class> class Model,class T,class List,int array_size,class Enable>
inline
void
ICR::EnsembleLearning::DeterministicNode<Model,T,List,array_size,Enable>::AddChildFactor(F_parameter f)
{ 
  //Could be many factors, potentially added with many threads,
  //therefore make the following critical
#pragma omp critical
  {
    m_children.push_back(f);
  }
}

template<template<class> class Model,class T,class List,int array_size,class Enable>
inline
const typename ICR::EnsembleLearning::DeterministicNode<Model,T,List,array_size,Enable>::moments_t*
ICR::EnsembleLearning::DeterministicNode<Model,T,List,array_size,Enable>::GetMoments() 
{

  /* This node is called exactly once per iteration.
   * Therefore can atomically update the moments via their pointers
   *  with a simple deletion mechanism for old pointers.
   * (The old pointer from the previous iteration is for sure no longer being used).
   *  Also, since called exactly once, 
   *  don't need a CAS to check that right pointer is updated.
   */
  // Lock(m_mutex);
  # pragma omp critical
   {
     delete m_Moments;
     //first get the NP from the parent
     NaturalParameters<T> ParentNP = (m_parent->GetNaturalNot(this));
     //Calcualate the moments
     m_Moments =  new moments_t(Model<T>::CalcMoments(ParentNP));
   }
//   //first of all, delete the old moments.
//   //delete m_oldMoments;
//   //first get the NP from the parent
//   NaturalParameters<T> ParentNP = (m_parent->GetNaturalNot(this));
//   //Calcualate the moments
//   //Lock lock(m_mutex);
//   //atomically swap in new pointer
//   moments_t *newMoments, *oldMoments;
//   newMoments = new moments_t(Model<T>::CalcMoments(ParentNP));  
// # pragma omp critical
//   {
//     //store the old pointer so that it can be deleted next iteration.
//     m_MomentsFromCurrentIteration.push_back( newMoments ); 
//   }
//   do {
//     oldMoments = m_Moments;
//   } while (!CAS(&m_Moments,oldMoments ,newMoments ));

  return m_Moments;
}
   

template<template<class> class Model,class T,class List,int array_size,class Enable>
inline
const typename ICR::EnsembleLearning::DeterministicNode<Model,T,List,array_size,Enable>::moments_t
ICR::EnsembleLearning::DeterministicNode<Model,T,List,array_size,Enable>::GetForwardedMoments() 
{
  /* No stored data here so threadsafe
   */

  //Need to collect this fresh, the moments from other parts 
  //of the graph may have been updated since the last call.

  std::vector<NaturalParameters<T> > vChildrenNP(m_children.size());

  PARALLEL_TRANSFORM(m_children.begin(), m_children.end(), 
  		     vChildrenNP.begin(), 
		     boost::bind(&FactorNode<T,DeterministicNode>::GetNaturalNot, _1, this)
  		     ); 

  BOOST_ASSERT(vChildrenNP.size() >0);
  NaturalParameters<T> ChildrenNP = PARALLEL_ACCUMULATE(vChildrenNP.begin(), vChildrenNP.end(), NaturalParameters<T>()); 
  return Model<T>::CalcMoments(ChildrenNP) ;// ;  //update the moments and the model

}
   
 
  

template<template<class> class Model,class T,class List,int array_size,class Enable>
inline
const std::vector<T>
ICR::EnsembleLearning::DeterministicNode<Model,T,List,array_size,Enable>::GetMean() 
{
  std::vector<T> Mean(1, m_Moments->operator[](0));
  return Mean;
}
   
template<template<class> class Model,class T,class List,int array_size,class Enable>
inline
const std::vector<T>
ICR::EnsembleLearning::DeterministicNode<Model,T,List,array_size,Enable>::GetVariance() 
{
  std::vector<T> Var(1,0.0);
  return Var;
}
   
template<template<class> class Model,class T,class List,int array_size,class Enable>
void 
ICR::EnsembleLearning::DeterministicNode<Model,T,List,array_size,Enable>::Iterate(Coster& C)
{
  //std::cout<<"ITERATION"<<std::endl;

  /* This node is called exactly once per iteration.
   * Therefore can atomically update the moments via their pointers
   *  with a simple deletion mechanism for old pointers.
   * (The old pointer from the previous iteration is for sure no longer being used).
   *  Also, since called exactly once, 
   *  don't need a CAS to check that right pointer is updated.
   */
  // std::cout<<"Last size = "<<m_MomentsFromLastIteration.size()<<std::endl;
  // std::cout<<"Current size = "<<m_MomentsFromCurrentIteration.size()<<std::endl;

  // size_t size  = m_MomentsFromCurrentIteration.size();
  // for(size_t i=0;i<m_number_to_delete;++i){
  //   delete m_MomentsFromCurrentIteration.front();
  //   m_MomentsFromCurrentIteration.pop_front;
  // }
  // m_number_to_delete = size-m_number_to_delete;

  // for(typename std::set<moments_t*>::iterator it = m_MomentsFromLastIteration.begin();
  //     it!=m_MomentsFromLastIteration.end();
  //     ++it){
  //   //remove old pointers from current list
  //   m_MomentsFromCurrentIteration.erase(*it);
  //   //delete redundent pointers.
  //   delete *it;
  // }

  // //copy the list still being used to past moments (to be deleted next time)
  // m_MomentsFromLastIteration = m_MomentsFromCurrentIteration;
	   
    
  
  //m_MomentsFromCurrentIteration.insert( newMoments ); 
  //first of all, delete the old moments.
  // delete m_oldMoments;
  // //first get the NP from the parent
  // NaturalParameters<T> ParentNP = (m_parent->GetNaturalNot(this));
  // //Calcualate the moments
  // //Lock lock(m_mutex);
  // moments_t* newMoments = new moments_t(Model<T>::CalcMoments(ParentNP));  
  // //store the old pointer so that it can be deleted next iteration.
  // m_oldMoments = m_Moments; 
  // //atomically swap in new pointer
  // m_Moments = newMoments;
}

#endif  // guard for VARIABLE_CALCULATION_HPP
