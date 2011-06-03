
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


#ifndef ENSEMBLE_LEARNING_HPP
#define ENSEMBLE_LEARNING_HPP

//The following header files contain classes that will be used by the user directly.
#include "EnsembleLearning/Output.hpp"
#include "EnsembleLearning/Builder.hpp"
#include "EnsembleLearning/calculation_tree/Factory.hpp"

/** @defgroup UserInterface The user interface of the Ensemble Learning Library.
 */

#endif //ENSEMBLE_LEARNING_HPP guard.



/**
 * @mainpage Ensemble Learning Documentation
 * @author Tom Shorrock (t.h.shorrock@gmail.com)
 * @date  June, 2011
 * @section purpose_sec Purpose
 * 
 * Provide a cross platform c++ library for Inference by Variational Ensemble Learning.
 *
 * @section example_sec Example:
 * @code#include "EnsembleLearning.hpp"
 * #include <vector>
 * #include <iostream>
 * using namespace ICR::EnsembleLearning;
 * 
 * int
 * main  (int ac, char **av)
 * {
 *   //Create the data
 *   rng* random = Random::Instance(); //a random number generator.
 *   std::vector<double> data(10);
 *   for(size_t i = 0; i<10; ++i)
 *   {
 *     data[i] = random->Gaussian(2, 3); //mean = 3, standard deviation = 2
 *   }
 *
 *   //Build the model
 *   Builder<double> build;  
 *  
 *   typedef Builder<double>::Variable Variable;  //A convenient alias
 *
 *   //Nodes to hold the learnt mean and precision
 *   Variable mean = build.gaussian(0.01,0.01);
 *   Variable prec = build.gamma(0.01,0.1);
 *
 *   //Model The data as Gaussian distributed, 
 *   //  each data point is modelled indepenantly
 *   for(size_t i; i<data.size(); ++i)
 *   {
 *    build.join(mean,prec,data[i];
 *   }
 * 
 *   //The model is complete, now need to do the inference.
 *   //  Iterate until convergance to within 0.1% or 100 iterations
 *   build.run(0.1,100);
 *
 *   //output the inferred mean and precision
 *   std::cout<<"mean      = "<<Mean(mean)<<" +- "<<StandardDeviation(mean)<<std::endl;
 *   std::cout<<"precision = "<<Mean(prec)<<" +- "<<StandardDeviation(prec)<<std::endl;
 *
 * }
 * @endcode
 *
 * \section feature_sec Feature
 * A swig interface file is provided for this library.
 * This means that a wrapper to the c++ library is easily created for any of the following languages:
 * - AllegroCL
 * - C# - Mono
 * - C# - MS .NET
 * - CFFI
 * - CHICKEN
 * - CLISP
 * - Go language
 * - Guile
 * - Java
 * - Lua
 * - MzScheme
 * - Ocaml
 * - Octave (and Matlab via Octave)
 * - Perl
 * - PHP
 * - Python
 * - R
 * - Ruby
 * - Tcl/Tk
 *
 * The following python script repeates the previous example.
 * @code
 * @endcode
 *
 *
 * \section dependancies_sec Dependancies.
 *  This library uses the Boost c++ libraries and the GNU Scientific Library.
 * 
 * \section git_repo_sec Location of the Repository
 * 
 * The library source code is stored in a git repository.
 * 
 * The Git repository is located at 
 * @verbatim  git@github.com:thshorrock/Ensemble-Learning.git  @endverbatim
 * To clone the repository type in a terminal ( or in the git bash terminal if you are using windows):
 * @code
 * git clone  git@github.com:thshorrock/Ensemble-Learning.git
 * @endcode
 * This will install a folder called "EnsembleLearning" at the location in which you typed the command.
 * 
 * 
 * 
 * \section installation_sec Installation
 * 
 * This library is most conveniently installed with cmake, although Bjam may also be used.
 * 
 *  - cmake is a meta-build tool.  It builds the local build tool for your system (for example linux make files, or Microsoft Visual Studio project files).
 *  - Bjam is a cross platform installation tool.
 *
 * @section licence_sec Licence
 * 
 * The library is licensed under the GPL, available for download: http://www.gnu.org/licenses/.
 * 
 * Since the library uses the Boost c++ libraries, the Boost lisence also applies.
 * 
 * @verbatim
 * Boost Software License - Version 1.0 - August 17th, 2003
 * 
 * Permission is hereby granted, free of charge, to any person or organization
 * obtaining a copy of the software and accompanying documentation covered by
 * this license (the "Software") to use, reproduce, display, distribute,
 * execute, and transmit the Software, and to prepare derivative works of the
 * Software, and to permit third-parties to whom the Software is furnished to
 * do so, all subject to the following:
 * 
 * The copyright notices in the Software and this entire statement, including
 * the above license grant, this restriction and the following disclaimer,
 * must be included in all copies of the Software, in whole or in part, and
 * all derivative works of the Software, unless such copies or derivative
 * works are solely in the form of machine-executable object code generated by
 * a source language processor.
 * 
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 * @endverbatim
 */

/** @file EnsembleLearning.hpp
 *  @brief The lecroy header files are specified
 *
 */
