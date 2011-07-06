

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
#ifndef MACRO_DEFAULTS_HPP
#define MACRO_DEFAULTS_HPP




//The maximum number of components is given by FUSION_MAX_VECTOR_SIZE (default is 20)
//If we want to increase the number of componets beyond this, then need to increase this size.
#ifdef ENSEMBLE_LEARNING_MAX_COMPONENTS
#  define FUSION_MAX_VECTOR_SIZE ENSEMBLE_LEARNING_MAX_COMPONENTS
#endif
#ifndef ENSEMBLE_LEARNING_COMPONENTS
#  define ENSEMBLE_LEARNING_COMPONENTS 5
#endif
#ifndef ENSEMBLE_LEARNING_PLACEHOLDERS
#  define ENSEMBLE_LEARNING_PLACEHOLDERS 5
#endif


#endif  // guard for MACRO_DEFAULTS_HPP
