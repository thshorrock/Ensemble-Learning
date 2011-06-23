
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

#include <boost/preprocessor/cat.hpp>

#define n BOOST_PP_ITERATION()

boost::shared_ptr<p0_t> BOOST_PP_CAT(shared_p0_,n)(new p0_t(v1[n]));
boost::shared_ptr<p1_t> BOOST_PP_CAT(shared_p1_,n)(new p1_t(v2[n]));
p0_t* BOOST_PP_CAT(p0_,n) = BOOST_PP_CAT(shared_p0_,n).get();
p1_t* BOOST_PP_CAT(p1_,n) = BOOST_PP_CAT(shared_p1_,n).get();

//MV.template get<n>() = new typename MV_t::template get_t<n>(v1[n],v2[n]);
//m_Nodes.push_back(MV.template get<n>() );
#undef n
