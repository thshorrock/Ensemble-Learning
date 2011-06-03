
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
#ifndef FUNCTIONSITERATOR_HPP
#define FUNCTIONSITERATOR_HPP

#include "Expression.hpp"

#include <boost/iterator/iterator_facade.hpp>
#include <boost/type_traits/is_convertible.hpp>
#include <boost/type_traits/is_same.hpp>
#include <boost/utility/enable_if.hpp>


namespace ICR{
  namespace EnsembleLearning{

    
    /**Iterator to iterate down the function tree (towards the root).*/ 
    template <class function>
    class FunctionIterator
      : public boost::iterator_facade<
      FunctionIterator<function>
      , function
      , boost::forward_traversal_tag
      >
    {
    private:
      struct enabler {}; //stop class being abused.
    public:
      /** Constructor */
      FunctionIterator()
	: m_function(0), m_previous(0) {}

      /** Constructor
       * @param f The function at which to start the iterator.
       */
      //This function vanishes if cannot convert other_function into function.
      // The constructor won't be present if you try to go from const->non-const.
      template <class other_function>
      explicit FunctionIterator(other_function* f
				, typename boost::enable_if<
				boost::is_convertible<other_function*,function*>
				, enabler
				>::type = enabler()
				)
	: m_function(f),m_previous(0) {}

      /** Return back one iterator.
       *  Although FunctionIterator is strictly a forward iterator. 
       *   a reference to the previous location is retained 
       *   and may be accessed here.
       *   These calls cannot be strung together to iterate backwards.
       *   Doing so produces undefined behaviour.
       *
       *   @return An iterator pointing to the previous location if such a location exists,
       *    otherwise the current location.
       */
      FunctionIterator<function> Previous() const
      {
	if (m_previous!=0) 
	  return FunctionIterator<function>(m_previous);
	else
	  return *this;
      }
    private:
      friend class boost::iterator_core_access;

      //equality
      template <class other_function>
      bool equal(FunctionIterator<other_function> const& other) const
      {
        return this->m_function == other.m_function;
      }
      
      //desend the tree
      void increment()
      { 
	m_previous = m_function;
	m_function = m_function->GetParent(); 
      }

      //dereference
      function& dereference() const
      { return *m_function; }
      
      //pointer to a particular function.
      function* m_function;
      function* m_previous;
    };
    
  }
}

#endif  // guard for FUNCTIONSITERATOR_HPP
