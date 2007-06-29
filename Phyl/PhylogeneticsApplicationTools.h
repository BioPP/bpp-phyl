//
// File: PhylogeneticsApplicationTools.h
// Created by: Julien Dutheil
// Created on: Fri Oct 21 16:49 2005
// from old file ApplicationTools.h created on Sun Dec 14 09:36:26 2003
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#ifndef _PHYLOGENETICSAPPLICATIONTOOLS_H_
#define _PHYLOGENETICSAPPLICATIONTOOLS_H_

#include "Tree.h"
#include "SubstitutionModel.h"
#include "MarkovModulatedSubstitutionModel.h"
#include "HomogeneousTreeLikelihood.h"
#include "ClockTreeLikelihood.h"

// From Utils:
#include <Utils/TextTools.h>
#include <Utils/StringTokenizer.h>

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>

// From SeqLib:
#include <Seq/SiteContainer.h>
#include <Seq/VectorSiteContainer.h>

// From the STL:
#include <string>
#include <map>
#include <iostream>
using namespace std;

/**
 * @brief This class provides some common tools for applications.
 *
 * The functions parse some option file, create corresponding objects and send
 * a pointer toward it.
 * 
 * The option files are supposed to follow this simple format:
 * @code
 * parameterName = parameterContent
 * @endcode
 * with one parameter per line.
 *
 * @see ApplicationTools
 */
class PhylogeneticsApplicationTools
{
	
	public:
		PhylogeneticsApplicationTools();
		virtual ~PhylogeneticsApplicationTools();
	

		/**
		 * @brief Build a Tree object according to options.
		 *
		 * Only newick format is supported for now.
		 * Options used are:
		 * - tree.file = file_path, the path of the file to parse.
		 *
		 * @param params  The attribute map where options may be found.
		 * @param suffix  A suffix to be applied to each attribute name.
		 * @param suffixIsOptional Tell if the suffix is absolutely required.
		 * @param verbose Print some info to the 'message' output stream.
		 * @return A new Tree object according to options specified.
     * @throw Exception if an error occured.
		 */
		static TreeTemplate<Node> * getTree(
			map<string, string> & params,
			const string & suffix = "",
			bool suffixIsOptional = true,
			bool verbose = true) throw (Exception);
	
		/**
		 * @brief This function prints the options available for tree reading.
		 */
		static void printInputTreeHelp();

    /**
		 * @brief Build a SubstitutionModel object according to options.
		 *
		 * Options used are:
		 * - model = [JCnuc|K80|T92|HKY85|F84|TN93|GTR|JCprot|DSO78|JTT92|empirical], the substitution model to use.
		 * Options depending on the model specified:
		 * - If K80, T92 or HKY85 or F84 is to be used:
		 *   + kappa The transition/transversion ratio.
		 * - If T92 format is to be used:
		 *   + theta The GC ratio, or
		 * - If HKY85, F84 or TN93 is to be used:
		 *   + piA, piT, piC and piG: equilibrum frequencies.
		 * - If TN93 is to be used:
		 *   + kappa1, kappa2 The transition/transversion ratios.
		 * - If GTR is to be used:
		 *   + a, b, c, d, e rate parameters.
		 * - If GTR, TN93, HKY85, F84, T92, JTT92, DSO78 or empirical is to be used:
		 *   + model.use_observed_freq Tell if we must use the observed frequences. 
		 * - If empirical is to be used;
		 *   + model_empirical.file Give the path toward the data file to use, in PAML format.
		 *
		 * @param alphabet The alpabet to use in the model.
		 * @param data     A pointer toward the SiteContainer for which the substitution model is designed.
		 * 								 The alphabet associated to the data must be of the same type as the one specified for the model.
		 *                 May be equal to NULL, but in this cas use_observed_freq option will be unavailable.
		 * @param params   The attribute map where options may be found.
		 * @param suffix   A suffix to be applied to each attribute name.
		 * @param suffixIsOptional Tell if the suffix is absolutely required.
		 * @param verbose Print some info to the 'message' output stream.
		 * @return A new SubstitutionModel object according to options specified.
     * @throw Exception if an error occured.
		 */
		static SubstitutionModel * getSubstitutionModel(
			const Alphabet * alphabet,
			const SiteContainer * data, 
			map<string, string> & params,
			const string & suffix = "",
			bool suffixIsOptional = true,
			bool verbose = true) throw (Exception);
	
		/**
		 * @brief Build a DiscreteDistribution object according to options.
		 *
		 * Options used are:
		 * - rate_distribution = [constant|gamma], the distribution to use.
		 * Options depending on other options:
		 * - If 'gamma' is selected:
		 *   - rate_distribution_gamma.alpha = the shape of the distribution
		 *   - rate_distribution.classes_number = the number of categories to be used.
		 *
		 * @param params  The attribute map where options may be found.
		 * @param suffix  A suffix to be applied to each attribute name.
		 * @param suffixIsOptional Tell if the suffix is absolutely required.
		 * @param verbose Print some info to the 'message' output stream.
		 * @return A new DiscreteDistribution object according to options specified.
     * @throw Exception if an error occured.
		 */
		static DiscreteDistribution * getRateDistribution(
			map<string, string> & params,
			string suffix = "",
			bool suffixIsOptional = true,
			bool verbose = true) throw (Exception);
	
    /**
     * @brief Build a MarkovModulatedSubstitutionModel object, according to options.
     *
     * - covarion = [none|G2001|TS98]
     * - for the G2001 model, parameter of the rate distribution are the same as for getRateDistribution() function.
     *   The distribution must not be uniform.
     * - covarion_G2001.nu = the nu paremeter.
     * - covarion_TS98.s1, covarion_TS98.s2 = T&S rate change parameters.
     *   
     * @param model   The substitution model to use.
     * @param rDist   The rate distribution to use. Depending on the model, this may not be used.
  	 * @param params  The attribute map where options may be found.
		 * @param suffix  A suffix to be applied to each attribute name.
		 * @param suffixIsOptional Tell if the suffix is absolutely required.
		 * @param verbose Print some info to the 'message' output stream.
		 * @return A new MarkovModulatedSubstitutionModel object according to options specified.
     * @throw Exception if an error occured.
     */
    static MarkovModulatedSubstitutionModel * getCovarionProcess(
      SubstitutionModel * model,
      DiscreteDistribution * rDist,
			map<string, string> & params,
			string suffix = "",
			bool suffixIsOptional = true,
			bool verbose = true) throw (Exception); 
    
		/**
		 * @brief Optimize parameters according to options.
		 *
		 * Options used are:
     * - optimization = Tell if optimization must be performed.
		 * - optimization.message_handler = [std, file_path]
		 *   A path to a specific path (existing will be overwritten) or std for use
		 *   of the standard output.
		 * - optimization.profiler = [std, file_path], idem for the profiling (history
		 *   of all functions evaluations).
		 * - optimization.max_number_f_eval = The maximum number of function evaluation.
		 * - optimization.tolerance = The tolerance parameter (when to stop the optimization).
		 * - optimization.scale_first = Tell if we must scale the tree first.
		 * - optimization.ignore_parameter = A coma-separated list of parameter
		 *   names to ignore in the optimizing process.
     * - optimization.method = [DB|fullD] Algorithm to use: Derivatives+Brent or full derivatives, with numerical derivatives when required.
     * - optimization.method.derivatives = [gradient|newton] Use Conjugate Grandient or Newton-Rhaphson algorithm.
     * - optimization.final = [none|simplex|powell] Perform a downhill simplex or a Powell multidimensions optimization
     * - optimization.topology = Tell if we must optimize tree topology. Toplogy estimation uses the DB algorithm with Newton-Raphson during estimation.
     *   The previous options will be used only for final estimation of numerical parameters.
		 * Options depending on other options:
		 * - If optimization.scale_first is set to true:
		 *   - optimization.scale_first.tolerance = The tolerance of the scaling alogrithm.
		 *   - optimization.scale_first.max_number_f_eval = the maximum number of function evaluations
		 *     for the scaling algorithm.
     * - optimization.method_DB.nstep = number of progressive steps to use in DB algorithm.
     * - optimization.topology.algorithm = [nni] algorithm to use (for now, only Nearest Neighbor Interchanges (NNI) are implemented). 
     * - optimization.topology.nstep = Estimate numerical parameters every 'n' NNI rounds.
     * - optimization.topology.numfirst = Tell if numerical parameters must be estimated prior to topology search.
     * - optimization.topology.tolerance.before = Numerical parameters estimation prior to topology search.
     * - optimization.topology.tolerance.during = Numerical parameters estimation during topology search.
		 *
		 * @param tl      The TreeLikelihood function to optimize.
		 * @param params  The attribute map where options may be found.
		 * @param suffix  A suffix to be applied to each attribute name.
		 * @param suffixIsOptional Tell if the suffix is absolutely required.
		 * @param verbose Print some info to the 'message' output stream.
		 * @throw Exception Any exception that may happen during the optimization process.
     * @return A pointer toward the final likelihood object.
     * This pointer may be the same as passed in argument (tl), but in some cases the algorithm
     * clone this object. We may change this bahavior in the future...
     * You hence should write something like
     * @code
     * tl = PhylogeneticsApplicationTools::optimizeParameters(tl, ...);
     * @endcode
		 */
		static TreeLikelihood * optimizeParameters(
			TreeLikelihood * tl,
			map<string, string> & params,
			const string & suffix = "",
			bool suffixIsOptional = true,
			bool verbose = true)
			throw (Exception);
		
		/**
		 * @brief Optimize parameters according to options, with a molecular clock.
		 *
		 * Options used are:
     * - optimization = Tell if optimization must be performed.
		 * - optimization.message_handler = [std, file_path]
		 *   A path to a specific path (existing will be overwritten) or std for use
		 *   of the standard output.
		 * - optimization.profiler = [std, file_path], idem for the profiling (history
		 *   of all functions evaluations).
		 * - optimization.max_number_f_eval = The maximum number of function evaluation.
		 * - optimization.tolerance = The tolerance parameter (when to stop the optimization).
		 * - optimization.ignore_parameter = A coma-separated list of parameter
		 *   names to ignore in the optimizing process.
     * - optimization.method = [DB|fullD] Algorithm to use: Derivatives+Brent or full derivatives, with numerical derivatives when required.
     * - optimization.method.derivatives = [gradient|newton] Use Conjugate Grandient or Newton-Rhaphson algorithm.
     * - optimization.final = [none|simplex|powell] Perform a downhill simplex or a Powell multidimensions optimization
     * - optimization.method_DB.nstep = number of progressive steps to use in DB algorithm.
		 *
		 * @param tl      The ClockTreeLikelihood function to optimize.
		 * @param params  The attribute map where options may be found.
		 * @param suffix  A suffix to be applied to each attribute name.
		 * @param suffixIsOptional Tell if the suffix is absolutely required.
		 * @param verbose Print some info to the 'message' output stream.
		 * @throw Exception Any exception that may happen during the optimization process.
		 */
		static void optimizeParameters(
			ClockTreeLikelihood * tl,
			map<string, string> & params,
			const string & suffix = "",
			bool suffixIsOptional = true,
			bool verbose = true)
			throw (Exception);

    /**
		 * @brief This function prints the options available for substitution models.
		 */
		static void printSubstitutionModelHelp();
		
    /**
		 * @brief This function prints the options available for covarion models.
		 */
		static void printCovarionModelHelp();
		
		/**
		 * @brief This function prints the options available for rate distributions.
		 */
		static void printRateDistributionHelp();
		
		/**
		 * @brief This function prints the options available for optimization.
     *
     * @param topo Display options for topology estimation
     * @param clock Restrict to options available for molecular clock estimation
		 */
		static void printOptimizationHelp(bool topo = true, bool clock = false);
		
		/**
		 * @brief Write a tree according to options.
		 *
		 * Options used are:
		 * - output.tree.file = file_path, the file where to put the tree.
		 *
		 * NB: only Newick format is supported for now. 
		 *
		 * @param tree    The tree to write.
		 * @param params  The attribute map where options may be found.
		 * @param suffix  A suffix to be applied to each attribute name.
		 * @param verbose Print some info to the 'message' output stream.
     * @throw Exception if an error occured.
		 */
		static void writeTree(
			const TreeTemplate<Node> & tree,
			map<string, string> & params,
			const string & suffix = "",
			bool verbose = true) throw (Exception);
		
		/**
		 * @brief This function prints the options available for tree writing.
		 */
		static void printOutputTreeHelp();
		
};

#endif	//_PHYLOGENETICSAPPLICATIONTOOLS_H_

