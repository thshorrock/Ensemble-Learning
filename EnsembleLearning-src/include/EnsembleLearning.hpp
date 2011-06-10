
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
#include "EnsembleLearning/Input.hpp"
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
 *
 * @section example_sec Example:
 *
 *  See the list of examples.  The source code is found in the example folder.
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
 *
 *
 * @example InferData1.cpp
 *  In this example 20 data points are drawn at random from a Gaussian distribution
 *  with a mean equal to 3 and a precision (the inverse of the variance) equal to 10.
 *
 *  The distribution from which the data points were sampled is inferred from the data.
 *
 *  The generated data points are as plotted.
 * @image html data_points.png
 * @image latex  data_points.eps  "data points"
 *
 *  The data is modelled (correctly) as being Gaussian distrubuted. 
 *  The mean and precision of the Gaussian are inferred from the data.
 *  The prior knowledge assumed on the mean and precision is fairly non-commital,
 *  as is seen from the plots.
 * @image html priors.png
 * @image latex  priors.eps  "Prior distributions for the mean and precision"
 *  
 *  Nevertheless, the inferred mean and precision is dominated by the 20 data points.
 *   - The inferred mean is  \f$ 3.012 \pm 0.078 \f$
 *   - The inferred precision is \f$ 8.14 \pm 2.45 \f$
 *
 * @example InferMixture1.cpp
 *  In this example 500 data points are drawn at random from three different
 * Gaussian distributions, with
 *   - mean equal to 3 and precision 10,
 *   - mean equal to 6 and precision 6,
 *   - mean equal to 3 and precision 2
 *
 *  The distribution from which the data points were sampled is inferred from the data.
 *
 *  The generated data points are as plotted.
 * @image html data_points.png
 * @image latex  data_points.eps  "data points"
 *
 *  The data is modelled  as being taken from a mixture of 5 Gaussian distrubutions
 *  The mean, precision and weight of each distribution  are inferred from the data.
 *  The prior knowledge assumed on the mean and precision is fairly non-commital,
 *  as is seen from the plots.
 * @image html priors.png
 * @image latex  priors.eps  "Prior distributions for the mean and precision"
 *  
 *  The inferred weights for  500  data points are drawn in the table, the cost is  -1.40356.
 *  
 *  <table>
 *  <tr>
 *  <th></th> <th>weight</th> <th>mean</th><th>precision</th>
 *  </tr>
 *  <tr><th>Component 1</th> <td>45.8% +- 2.2%</td> <td>3.00 +- 0.02</td> <td>9.2 +- 0.9</td></tr>
 *  <tr><th>Component 2</th> <td>29.1% +- 2.0%</td> <td>6.00 +- 0.04</td> <td>5.1 +- 0.6</td></tr>
 *  <tr><th>Component 3</th> <td>21.0% +- 1.8%</td> <td>2.90 +- 0.07</td> <td>2.0 +- 0.3</td></tr>
 *  <tr><th>Component 4</th> <td>2.9%  +- 0.8%</td> <td>3.79 +- 0.02</td> <td>170 +- 61</td></tr>
 *  <tr><th>Component 5</th> <td>1.1%  +- 0.5%</td> <td>4.55 +- 0.06</td> <td>66  +- 35</td></tr>
 *  </table>
 *
 *
 */
/** @file EnsembleLearning.hpp
 *  @brief The lecroy header files are specified
 *
 */
