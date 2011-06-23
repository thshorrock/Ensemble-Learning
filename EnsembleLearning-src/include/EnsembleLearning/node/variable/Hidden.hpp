#pragma once
#ifndef HIDDEN_HPP
#define HIDDEN_HPP


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
#include "EnsembleLearning/detail/parallel_algorithms.hpp"
#include "EnsembleLearning/detail/Mutex.hpp"
#include "EnsembleLearning/detail/TypeList.hpp"

#include <boost/assert.hpp> 
#include <boost/bind.hpp>
#include <vector>



namespace ICR{
  namespace EnsembleLearning{
    
    /** A Hidden Variable Node.
     *  HiddenNode's store moments that are inferred by Ensemble Learning.
     *  The model used (for example Gaussian/Gamma) is passed by template parameter.
     *  @tparam Model  The model to use for the inferred data.
     *  @tparam T The data type (float or double)
     */
    template <template<class> class Model, class T, class List=detail::TypeList::zeros
	     ,class Enable = void >
    class HiddenNode : public VariableNode<T> //, public List
    {
      HiddenNode(HiddenNode<Model,T, List>& other); //non-copyable

    public:
      /** A constructor.
       *  @param moment_size The number of elements in the stored moments.
       *  This is usually two but varies for discrete nodes.
       */
      HiddenNode(const size_t moment_size = 2); //mostly two, but variable for Discrete model

      void
      SetParentFactor(FactorNode<T, HiddenNode>* f);
      
      void
      AddChildFactor(FactorNode<T, HiddenNode>* f);

      void 
      Iterate(Coster& C);

      void
      InitialiseMoments()
      {
	m_Moments = m_parent->InitialiseMoments();
      }


      const Moments<T>&
      GetMoments() ;

      const std::vector<T>
      GetMean() ;
      
      const std::vector<T>
      GetVariance() ;

      void
      SetMean(const std::vector<T>& mean) ;
      
      void
      SetVariance(const std::vector<T>& v) ;

      

      /** The number of elements in the stored Moments */
      size_t 
      size() const {return m_Moments.size();}
      
    private:
      
      const NaturalParameters<T>
      GetNP();

      FactorNode<T, HiddenNode>* m_parent;
      std::vector<FactorNode<T, HiddenNode>*> m_children;
      Moments<T> m_Moments;
      mutable Mutex m_mutex;
    };

  }
}


template<template<class> class Model,class T,class List,class Enable>
ICR::EnsembleLearning::HiddenNode<Model,T,List,Enable>::HiddenNode(const size_t moment_size) 
  :   m_parent(0), m_children(), m_Moments(moment_size)
{}


template<template<class> class Model,class T,class List,class Enable>
inline 
void
ICR::EnsembleLearning::HiddenNode<Model,T,List,Enable>::SetParentFactor(FactorNode<T, HiddenNode>* f)
{
  //This should only be called once, so should get no collisions here
  m_parent=f;
  // Can now initialise
  InitialiseMoments();
  
}



   
template<template<class> class Model,class T,class List,class Enable>
inline
void
ICR::EnsembleLearning::HiddenNode<Model,T,List,Enable>::AddChildFactor(FactorNode<T, HiddenNode>* f)
{ 
  //This could be called simultaneously by different threads, so call it critical
#pragma omp critical
  {
    m_children.push_back(f);
  }
}

template<template<class> class Model,class T,class List,class Enable>
inline
const ICR::EnsembleLearning::Moments<T>&
ICR::EnsembleLearning::HiddenNode<Model,T,List,Enable>::GetMoments() 
{
  /*This value is update in Iterate and called to evaluate other Hidden Nodes
   * (also in iterate mode).  It therefore needs to be protected by a mutex.
   */
  Lock lock(m_mutex);
  return m_Moments;
}
   
template<template<class> class Model,class T,class List,class Enable>
inline
const ICR::EnsembleLearning::NaturalParameters<T>
ICR::EnsembleLearning::HiddenNode<Model,T,List,Enable>::GetNP()
{
  //EVALUATE THE COST
  BOOST_ASSERT(m_parent != 0); 
  //first get the NP from the parent
  const NaturalParameters<T> ParentNP = (m_parent->GetNaturalNot(this));
  //The initial value in the total
  NaturalParameters<T> NP = ParentNP; //total

  //get the Natural parameters from all the children.
  std::vector<NaturalParameters<T> > ChildrenNP(m_children.size());
  //and put them in the vector
  for(size_t i=0;i<m_children.size();++i){
    ChildrenNP[i] = m_children[i]->GetNaturalNot(this);
  }

  // PARALLEL_TRANSFORM(m_children.begin(), m_children.end(), 
  // 		     ChildrenNP.begin(), boost::bind(&FactorNode<T, HiddenNode>::GetNaturalNot, boost::ref(_1), this)
  // 		     ); //from children

  //Add them up
  //(Note: Bug 48750 in gcc in parallel accumulate causes crash at next line: 
  //  fixed in 4.6.1)
  NP = PARALLEL_ACCUMULATE(ChildrenNP.begin(), ChildrenNP.end(), NP ); 
  return NP;
}

template<template<class> class Model,class T,class List,class Enable>
inline
const std::vector<T>
ICR::EnsembleLearning::HiddenNode<Model,T,List,Enable>::GetMean() 
{
  return Model<T>::CalcMean(GetNP());
}
   
template<template<class> class Model,class T,class List,class Enable>
inline
const std::vector<T>
ICR::EnsembleLearning::HiddenNode<Model,T,List,Enable>::GetVariance() 
{
  std::vector<T> prec = Model<T>::CalcPrecision(GetNP());
  std::vector<T> var(prec.size());
  for(size_t i=0;i<prec.size();++i){
    var[i] = 1.0/prec[i];
  }
  return var;
}

template<template<class> class Model,class T,class List,class Enable>
inline
void
ICR::EnsembleLearning::HiddenNode<Model,T,List,Enable>::SetMean(const std::vector<T>& m) 
{
  m_Moments = Model<T>::CalcMoments(m,GetVariance());
}
   
template<template<class> class Model,class T,class List,class Enable>
inline
void
ICR::EnsembleLearning::HiddenNode<Model,T,List,Enable>::SetVariance(const std::vector<T>& v) 
{
  m_Moments = Model<T>::CalcMoments(GetMean(),v);
}

template<template<class> class Model,class T,class List,class Enable>
inline
void 
ICR::EnsembleLearning::HiddenNode<Model,T,List,Enable>::Iterate(Coster& C)
{
  const NaturalParameters<T> NP = GetNP();
  //Get the moments and update the model
  const T LogNorm = Model<T>::CalcLogNorm(NP);
  {
    Lock lock(m_mutex);
    m_Moments = Model<T>::CalcMoments(NP);  //update the moments and the model
  }
  //first get the NP from the parent
  const NaturalParameters<T> ParentNP = (m_parent->GetNaturalNot(this));
  C +=  (ParentNP - NP)*m_Moments +m_parent->CalcLogNorm() -  LogNorm;

}

#endif  // guard for HIDDEN_HPP
