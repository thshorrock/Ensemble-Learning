

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
#ifndef TERNARYOP_HPP
#define TERNARYOP_HPP


#include <boost/fusion/sequence/intrinsic/begin.hpp>
#include <boost/fusion/sequence/intrinsic/end.hpp>
#include <boost/fusion/iterator/equal_to.hpp>
#include <boost/fusion/iterator/next.hpp>
#include <boost/fusion/iterator/deref.hpp>
#include <boost/fusion/iterator/distance.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/fusion/support/category_of.hpp>



namespace ICR{
  namespace EnsembleLearning{

    
    namespace detail{
      
      
      template <typename First, typename Last,typename First2, typename First3, typename F>
      inline void
      TernaryOp_linear(First const&, Last const&,
		       First2 const&,
		       First3 const&, F const&, boost::mpl::true_)
      {
      }
      
      template <typename First, typename Last,typename First2, typename First3,typename F>
      inline void
      TernaryOp_linear(First const& first, 
		       Last const& last,
		       First2 const& first2,
		       First3 const& first3, F const& f, boost::mpl::false_)
      {
        f(*first,*first2,*first3);
        TernaryOp_linear(boost::fusion::next(first), last,
			 boost::fusion::next(first2),
			 boost::fusion::next(first3),
			 f, 
			 boost::fusion::result_of::equal_to<typename boost::fusion::result_of::next<First>::type, Last>());
      }

      template <typename Sequence1,typename Sequence2,typename Sequence3, typename F>
      inline void
      TernaryOp(Sequence1& seq1,Sequence2& seq2,Sequence3& seq3, F const& f)
      {
        TernaryOp_linear(
			 boost::fusion::begin(seq1)
			 , boost::fusion::end(seq1)
			 , boost::fusion::begin(seq2)
			 , boost::fusion::begin(seq3)
			 , f
			 , boost::fusion::result_of::equal_to<
			 typename boost::fusion::result_of::begin<Sequence1>::type
			 , typename boost::fusion::result_of::end<Sequence1>::type>());
    }



  }
    


}
}
#endif  // guard for TERNARYOP_HPP
  
