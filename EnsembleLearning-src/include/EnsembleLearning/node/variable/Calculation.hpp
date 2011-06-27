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
    template <class Model, class T, class List=detail::TypeList::zeros >
    class DeterministicNode : public VariableNode<T>//, public List
    {
    public:
      //mostly two, but variable for Discrete model
      /** A constructor.
       *  @param moment_size The number of elements in the stored moments.
       *  This is usually two but varies for discrete nodes.
       */
      DeterministicNode(const size_t moment_size = 2); 

      void
      SetParentFactor(ParentFactorNode<T,DeterministicNode>* f);
      
      void
      AddChildFactor(FactorNode<T,DeterministicNode>* f);

      void 
      Iterate(Coster& C);

      const Moments<T>
      GetMoments() ;

      //should make a const version of this 
      /** Forward moments to the deterministic function.
       *  @return The moments from the child node.
       */
      const Moments<T>
      GetForwardedMoments() ;

      /** The number of elements in the stored Moments */
      size_t 
      size() const {return m_Moments.size();}
      
      const std::vector<T>
      GetMean() ;
      
      const std::vector<T>
      GetVariance() ;

      void
      InitialiseMoments()
      {
	m_Moments = m_parent->InitialiseMoments();
      }


      
    private:
      
      ParentFactorNode<T,DeterministicNode>* m_parent;
      std::vector<FactorNode<T,DeterministicNode>*> m_children;
      mutable Moments<T> m_Moments;
      mutable Mutex m_mutex;
    };

  }
}


template<class Model,class T,class List>
ICR::EnsembleLearning::DeterministicNode<Model,T,List>::DeterministicNode(const size_t moment_size) 
  :   m_parent(0), m_children(), m_Moments(moment_size) //, m_ForwardedMoments(moment_size)
{
}


template<class Model,class T,class List>
inline 
void
ICR::EnsembleLearning::DeterministicNode<Model,T,List>::SetParentFactor(ParentFactorNode<T,DeterministicNode>* f)
{
  //This should only be called once, so should get no collisions here
  m_parent=f;
  // Can now initialise
  InitialiseMoments();
}
   
template<class Model,class T,class List>
inline
void
ICR::EnsembleLearning::DeterministicNode<Model,T,List>::AddChildFactor(FactorNode<T,DeterministicNode>* f)
{ 
  //Could be many factors, potentially added with many threads,
  //therefore make the following critical
#pragma omp critical
  {
    m_children.push_back(f);
  }
}

template<class Model,class T,class List>
inline
const ICR::EnsembleLearning::Moments<T>
ICR::EnsembleLearning::DeterministicNode<Model,T,List>::GetMoments() 
{

  /*This value is update in Iterate and called to evaluate other Hidden Nodes
   * (also in iterate mode).  It therefore needs to be protected by a mutex.
   */
  //first get the NP from the parent
  NaturalParameters<T> ParentNP = (m_parent->GetNaturalNot(this));
  //Calcualate the moments
  Lock lock(m_mutex);
  m_Moments =  Model::CalcMoments(ParentNP);  
  return m_Moments;
}
   

template<class Model,class T,class List>
inline
const ICR::EnsembleLearning::Moments<T>
ICR::EnsembleLearning::DeterministicNode<Model,T,List>::GetForwardedMoments() 
{
  /* No stored value updated here so thread safe.
   */

  //Need to collect this fresh, the moments from other parts 
  //of the graph may have been updated since the last call.

  std::vector<NaturalParameters<T> > vChildrenNP(m_children.size());

  PARALLEL_TRANSFORM(m_children.begin(), m_children.end(), 
  		     vChildrenNP.begin(), 
		     boost::bind(&FactorNode<T,DeterministicNode>::GetNaturalNot, _1, this)
  		     ); 

  BOOST_ASSERT(vChildrenNP.size() >0);
  NaturalParameters<T> ChildrenNP = PARALLEL_ACCUMULATE(vChildrenNP.begin(), vChildrenNP.end(), NaturalParameters<T>(vChildrenNP[0].size() )); 
  const Moments<T> ForwardedMoments = Model::CalcMoments(ChildrenNP);// ;  //update the moments and the model

  return ForwardedMoments;
}
   
 
  

template<class Model,class T,class List>
inline
const std::vector<T>
ICR::EnsembleLearning::DeterministicNode<Model,T,List>::GetMean() 
{
  std::vector<T> Mean(1, m_Moments[0]);
  return Mean;
}
   
template<class Model,class T,class List>
inline
const std::vector<T>
ICR::EnsembleLearning::DeterministicNode<Model,T,List>::GetVariance() 
{
  std::vector<T> Var(1,0.0);
  return Var;
}
   
template<class Model,class T,class List>
void 
ICR::EnsembleLearning::DeterministicNode<Model,T,List>::Iterate(Coster& C)
{}

#endif  // guard for VARIABLE_CALCULATION_HPP
