//
// File: OptimizationTools.h
// Created by: Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
// Created on: Sun Dec 14 09:43:32 2003
//

/*
Copyright ou © ou Copr. CNRS, (16 Novembre 2004) 

Julien.Dutheil@univ-montp2.fr

Ce logiciel est un programme informatique servant à fournir des classes
pour l'analyse de données phylogénétiques.

Ce logiciel est régi par la licence CeCILL soumise au droit français et
respectant les principes de diffusion des logiciels libres. Vous pouvez
utiliser, modifier et/ou redistribuer ce programme sous les conditions
de la licence CeCILL telle que diffusée par le CEA, le CNRS et l'INRIA 
sur le site "http://www.cecill.info".

En contrepartie de l'accessibilité au code source et des droits de copie,
de modification et de redistribution accordés par cette licence, il n'est
offert aux utilisateurs qu'une garantie limitée.  Pour les mêmes raisons,
seule une responsabilité restreinte pèse sur l'auteur du programme,  le
titulaire des droits patrimoniaux et les concédants successifs.

A cet égard  l'attention de l'utilisateur est attirée sur les risques
associés au chargement,  à l'utilisation,  à la modification et/ou au
développement et à la reproduction du logiciel par l'utilisateur étant 
donné sa spécificité de logiciel libre, qui peut le rendre complexe à 
manipuler et qui le réserve donc à des développeurs et des professionnels
avertis possédant  des  connaissances  informatiques approfondies.  Les
utilisateurs sont donc invités à charger  et  tester  l'adéquation  du
logiciel à leurs besoins dans des conditions permettant d'assurer la
sécurité de leurs systèmes et ou de leurs données et, plus généralement, 
à l'utiliser et l'exploiter dans les mêmes conditions de sécurité. 

Le fait que vous puissiez accéder à cet en-tête signifie que vous avez 
pris connaissance de la licence CeCILL, et que vous en avez accepté les
termes.
*/

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

Julien.Dutheil@univ-montp2.fr

This software is a computer program whose purpose is to provide classes
for phylogenetic data analysis.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef _OPTIMIZATIONTOOLS_H_
#define _OPTIMIZATIONTOOLS_H_

#include "TreeLikelihood.h"
#include "HomogeneousTreeLikelihood.h"
#include "AbstractHomogeneousTreeLikelihood.h"

// From the STL:
#include <iostream>
using namespace std;

/**
 * @brief Optimization methods for phylogenetic inference.
 *
 * NB: For now, only numerical parameter optimization is performed.
 *
 * This class provides optimization methods for phylogenetics.
 * They are combinations of classical optimization methods that can be found
 * in the NumCalc library.
 *
 * Parameters of the optimization methods are set to work with TreeLikelihood
 * object. Some non trivial parameters are left to the user choice (tolerance, maximum
 * number of function evaluation, output streams).
 */
class OptimizationTools
{
	public:
		OptimizationTools();
		virtual ~OptimizationTools();
	
	public:
		
		/**
		 * @brief Optimize a TreeLikelihood object with the Downhill Simplex method.
		 *
		 * This method optimize all parameters with the Downhill Simple method.
		 * The constraint mode is set to 'AUTO', meaning that all parameters are
		 * replaced by parameters of class AutoParameter.
		 * Indeed, this method does not work on likelihood objects with values 
		 * out of constraint: some null value for branch length for instance lead
		 * to NaN values for likelihood,  which 'attract' the algorithm.
		 *
		 * A condition over parameters is used as a stop condition for the algorithm.
		 *
		 * @param tl             A pointer toward the TreeLikelihood object to optimize.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @param verbose        The verbose level.
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static int optimizeWithDownhillSimplexMethod(
			TreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			unsigned int verbose = 1
			)	throw (Exception);

		/**
		 * @brief Optimize a TreeLikelihood object with Powell's multidimensions method.
		 *
		 * This method optimize all parameters with Powell's multidimensions method.
		 * The constraint mode is set to 'IGNORE', meaning that all constraints on
		 * parameters are removed.
		 * This is because Powell's algorithm may lead to transitionnal value out of
		 * parameter constraints. AutoParameter objects hence may lead to values that
		 * does not match the local - i.e. one-dimensional - optimum.
		 *
		 * A condition over parameters is used as a stop condition for the algorithm.
		 *
		 * @param tl             A pointer toward the TreeLikelihood object to optimize.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @param verbose        The verbose level.
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static int optimizeWithPowellMethod(
			TreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			unsigned int verbose = 1
			)	throw (Exception);
		
		/**
		 * @brief Optimize a TreeLikelihood object with Newton's method.
		 *
		 * This method optimize all parameters with Galtier's modified Newton's method.
		 * The constraint mode is set to 'AUTO', meaning that all constraints on
		 * parameters are taken into acount (see the AutoParameter class).
		 *
		 * A condition over function values is used as a stop condition for the algorithm.
		 *
		 * @param tl             A pointer toward the TreeLikelihood object to optimize.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @param verbose        The verbose level.
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static int optimizeWithNewtonMethod(
			TreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			unsigned int verbose = 1
			)	throw (Exception);
	
		/**
		 * @brief Optimize a TreeLikelihood object with Newton's method for branch length
		 * and Brent's one dimensional method for other parameters.
		 *
		 * A condition over function values is used as a stop condition for the algorithm.
		 *
		 * @param tl             A pointer toward the TreeLikelihood object to optimize.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @param verbose        The verbose level.
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static int optimizeWithNewtonBrentMethod(
			AbstractHomogeneousTreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			unsigned int verbose = 1)
			throw (Exception);
	
		static int optimizeWithDownhillSimplexAndPowellMethod(
			TreeLikelihood * tl,
			double ratio,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			unsigned int verbose = 1
			)	throw (Exception);
		
	private:
		
		class ScaleFunction: public Function {
				
			protected:
				TreeLikelihood * _tl;
				mutable ParameterList _brLen, _lambda;
				
			public:
				ScaleFunction(TreeLikelihood * tl);
				~ScaleFunction();
				
			public:
				void setParameters(const ParameterList & lambda) throw (ParameterNotFoundException, ConstraintException);
				double getValue() const throw (ParameterException);
				ParameterList getParameters() const throw (Exception) { return _lambda; }
				double getParameter(const string & name) const throw (ParameterNotFoundException) { return _lambda.getParameter(name) -> getValue(); };
				void setAllParametersValues(const ParameterList & params) 
					throw (ParameterNotFoundException, ConstraintException) {}
				void setParameterValue(const string & name, double value) 
					throw (ParameterNotFoundException, ConstraintException) {}
				void setParametersValues(const ParameterList & params)
					throw (ParameterNotFoundException, ConstraintException) {}
				void matchParametersValues(const ParameterList & params)
					throw (ConstraintException) {};
		};
	
	public:

		/**
		 * @brief Optimize the scale of a TreeLikelihood.
		 *
		 * This method only works on branch lengths parameters.
		 * It multiply all branch length by a factor 'x' which is optimized
		 * using Brent's algorithm in one dimension.
		 * This method may be usefull for scaling a tree whose branch lengths
		 * come from the Neighbor-Joining algorithm for instance.
		 *
		 * Practically, and contrarily to what one may expect, this does not
		 * speed up the optimization!
		 *
		 * A condition over parameters is used as a stop condition for the algorithm.
		 *
		 * @param tl             A pointer toward the TreeLikelihood object to optimize.
		 * @param tolerance      The tolerance to use in the algorithm.
		 * @param tlEvalMax      The maximum number of function evaluations.
		 * @param messageHandler The massage handler.
		 * @param profiler       The profiler.
		 * @param verbose        The verbose level.
		 * @throw Exception any exception thrown by the Optimizer.
		 */
		static int optimizeTreeScale(
			TreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			unsigned int verbose = 1
			)	throw (Exception);
	
		static int optimizeWithDownhillSimplexMethodAlphaSeparately(
			TreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			ostream * profilerAlpha  = &cout,
			unsigned int verbose = 1
			)	throw (Exception);	

		static int optimizeWithPowellMethodAlphaSeparately(
			TreeLikelihood * tl,
			double tolerance = 0.000001,
			int tlEvalMax = 1000000,
			ostream * messageHandler = &cout,
			ostream * profiler       = &cout,
			ostream * profilerAlpha  = &cout,
			unsigned int verbose = 1
			)	throw (Exception);	
};


#endif	//_OPTIMIZATIONTOOLS_H_
