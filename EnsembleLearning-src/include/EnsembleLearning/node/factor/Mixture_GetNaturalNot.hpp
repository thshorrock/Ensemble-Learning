
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


#define n BOOST_PP_ITERATION()

      inline
      NaturalParameters<T>
      GetNaturalNot(typename boost::add_pointer<typename boost::mpl::at_c<typename parent1_t::type,n>::type>::type v) const
      { 
        const Moments<T>& weights = m_weights_node->GetMoments();
        const Moments<T>& parent2 = m_parent2_nodes.template get<n>()->GetMoments();
        const Moments<T>& child = m_child_node->GetMoments();
        return Model::CalcNP2Parent1(parent2,child)  * weights[n];
      }

      inline
      NaturalParameters<T>
      GetNaturalNot(typename boost::add_pointer<typename boost::mpl::at_c<typename parent2_t::type,n>::type>::type v) const
      { 
        const Moments<T>& weights = m_weights_node->GetMoments();
        const Moments<T>& parent1 = m_parent1_nodes.template get<n>()->GetMoments();
        const Moments<T>& child = m_child_node->GetMoments();
        return Model::CalcNP2Parent2(parent1,child) * weights[n];
      }

#undef n
