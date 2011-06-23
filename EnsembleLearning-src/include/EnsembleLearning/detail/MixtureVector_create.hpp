
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
#include <boost/type_traits/remove_reference.hpp>
#include <boost/type_traits/remove_pointer.hpp>
#define n BOOST_PP_ITERATION()
{
  
  // typedef ObservedNode<Model, T, ICR::EnsembleLearning::detail::TypeList_Exact<n,0,n> > p0_t;
  // typedef ObservedNode<Gamma, T, ICR::EnsembleLearning::detail::TypeList_Exact<n,0,n> >    p1_t;
  // boost::shared_ptr<p0_t> P0(new p0_t(v1[n]));
  // boost::shared_ptr<p1_t> P1(new p1_t(v2[n]));
  
  // m_Nodes.push_back(P0);
  // m_Nodes.push_back(P1);

  typedef typename boost::remove_pointer<typename boost::remove_reference<typename boost::fusion::result_of::at_c<p0_t,n>::type>::type>::type p0n_t;
  typedef typename boost::remove_pointer<typename boost::remove_reference<typename boost::fusion::result_of::at_c<p1_t,n>::type>::type>::type p1n_t;
 
  typedef  typename MV_t::template get_t<n>::type child_t;

  boost::shared_ptr<child_t> child( new child_t() );
  
  typedef detail::Factor<Model,T,p0n_t,p1n_t,child_t> Factor_t;
  boost::shared_ptr<Factor_t > F(new Factor_t(boost::fusion::at_c<n>(v0), 
					      boost::fusion::at_c<n>(v1), //v1.get(), 
					      child.get()));

  m_Factors.push_back(F);
  m_Nodes.push_back(child);
  MV.template get<n>()  = child.get();
}
//MV.template get<n>() = new typename MV_t::template get_t<n>(v1[n],v2[n]);
//m_Nodes.push_back(MV.template get<n>() );
#undef n
