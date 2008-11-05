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
#include "SubstitutionModelSet.h"
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

namespace bpp
{

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
     * @brief Build a SubstitutionModel object with default parameter values according to options.
     *
     * Options used are:
     * - %{prefix}name = [JCnuc|K80|T92|HKY85|F84|TN93|GTR|L95|JCprot[+F]|DSO78[+F]|JTT92[+F]|empirical[+F]], the substitution model to use.
     *
     * This function is mainly for internal usage, you're probably looking for the getSubstitutionModel or getSubstitutionModelSet function.
     *
     * @param alphabet The alpabet to use in the model.
     * @param params   The attribute map where options may be found.
     * @param prefix   A prefix to be applied to each attribute name (model number).
     * @param suffix   A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose Print some info to the 'message' output stream.
     * @return A new SubstitutionModel object according to options specified.
     * @throw Exception if an error occured.
     */
    static SubstitutionModel * getSubstitutionModelDefaultInstance(
      const Alphabet * alphabet,
      map<string, string> & params,
      const string & prefix,
      const string & suffix,
      bool suffixIsOptional,
      bool verbose) throw (Exception);

    /**
     * @brief Set parameter initial values of a given model according to options.
     *
     * Parameters actually depends on the model passed as argument.
     * See getSubstitutionModel for more information.
     *
     * This function is mainly for internal usage, you're probably looking for the getSubstitutionModel or getSubstitutionModelSet function.
     *
     * @param model  The model to set.
     * @param data   A pointer toward the SiteContainer for which the substitution model is designed.
     *               The alphabet associated to the data must be of the same type as the one specified for the model.
     *               May be equal to NULL, but in this case use_observed_freq option will be unavailable.
     * @param params The attribute map where options may be found.
     * @param prefix A prefix to be applied to each attribute name.
     * @param suffix A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose Print some info to the 'message' output stream.
     * @throw Exception if an error occured.
     */
    static void setSubstitutionModelParametersInitialValues(
      SubstitutionModel * model,
      const SiteContainer * data,
      map<string, string> & params,
      const string & prefix,
      const string & suffix,
      bool suffixIsOptional,
      bool verbose) throw (Exception);

    /**
     * @brief Build a SubstitutionModel object according to options.
     *
     * Options used are:
     * - model.name = [JCnuc|K80|T92|HKY85|F84|TN93|GTR|JCprot|DSO78|JTT92|empirical], the substitution model to use,
     *   +G2001 or +TS98 for enabling covarion models.
     *
     * Options depending on the model specified:
     * - If K80, T92 or HKY85 or F84 is to be used:
     *   + model.kappa The transition/transversion ratio.
     * - If T92 or L95 format is to be used:
     *   + model.theta The GC ratio, or
     * - If HKY85, F84 or TN93 is to be used:
     *   + model.piA, piT, piC and piG: equilibrum frequencies.
     * - If TN93 is to be used:
     *   + model.kappa1, kappa2 The transition/transversion ratios.
     * - If GTR is to be used:
     *   + model.a, b, c, d, e rate parameters.
     * - If L95 is to be used:
     *   + model.beta, gamma, delta rate parameters.
     * - If GTR, TN93, HKY85, F84, T92, JTT92, DSO78, L95 or empirical is to be used:
     *   + model.use_observed_freq Tell if we must use the observed frequences. 
     * - If empirical is to be used;
     *   + model_empirical.file Give the path toward the data file to use, in PAML format.
     * - For the T92, HKY85, F84, TN93, GTR, + all prot models, the model.use_observed_freq option lets you
     *   use the actual frequencies estimated from a given data set.
     *
     * Covarions options:
     * - for the G2001 model, parameter of the rate distribution are the same as for getRateDistribution() function.
     *   The distribution must not be uniform.
     * - model.nu = the nu paremeter.
     * - model.rate_distribution, etc (see getRateDistribution) for detailled options:
     *   The rate distribution associated to this covarion model.
     * - model.s1, model.s2 = T&S rate change parameters.
     *
     * @param alphabet The alphabet to use in the model.
     * @param data     A pointer toward the SiteContainer for which the substitution model is designed.
     *                 The alphabet associated to the data must be of the same type as the one specified for the model.
     *                 May be equal to NULL, but in this case use_observed_freq option will be unavailable.
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
     * @brief Set parameter initial values of a given model in a set according to options.
     *
     * Parameters actually depends on the model passed as argument.
     * See getSubstitutionModelSet for more information.
     *
     * This function is mainly for internal usage, you're probably looking for the getSubstitutionModel or getSubstitutionModelSet function.
     *
     * @param model The model to set.
     * @param data     A pointer toward the SiteContainer for which the substitution model is designed.
     *                 The alphabet associated to the data must be of the same type as the one specified for the model.
     *                 May be equal to NULL, but in this case use_observed_freq option will be unavailable.
     * @param existingParams (in/out) A map with already existing value that have been found in previous calls, and may be recalled here.
     *                       New parameters found here will be added.
     * @param specificParams (out) Parameters specific to this model will be recorded here.
     * @param sharedParams (out) remote parameters will be recorded here.
     * @param params The attribute map where options may be found.
     * @param prefix A prefix to be applied to each attribute name (model number in the set).
     * @param suffix A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose Print some info to the 'message' output stream.
     * @throw Exception if an error occured.
     */
    static void setSubstitutionModelParametersInitialValues(
      SubstitutionModel * model,
      const SiteContainer * data,
      map<string, double> & existingParams,
      vector<string> & specificParams,
      vector<string> & sharedParams,
      map<string, string> & params,
      const string & prefix,
      const string & suffix,
      bool suffixIsOptional,
      bool verbose) throw (Exception);

    /**
     * @brief Get A FrequenciesSet object according to options.
     *
     * Available options are:
     * nonhomogeneous.root_freq one of the following:
     * - balanced: All frequencies set to 1/n, where n is the size of the alphabet.
     * - observed: Use observed counts
     * - init: Specify each frequency using parameters with names ancA, ancC, ... ancA, ancR, ancN...
     * - balancedGC: use the GC content parametrization
     * - observedGC: use the observed GC content
     * - initGC: use ancTheta to set the ancestral GC content.
     *
     * @param alphabet The alpabet to use.
     * @param data      A pointer toward the SiteContainer for which the substitution model is designed.
     *                  The alphabet associated to the data must be of the same type as the one specified for the model.
     *                  May be equal to NULL, but in this cas use_observed_freq option will be unavailable.
     * @param params    The attribute map where options may be found.
     * @param rateFreqs A vector of rate catégories frequencies in case of a Markov Modulated Markov Model.
     *                  Ignored if a vector with size 0 is passed.
     * @param suffix    A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose   Print some info to the 'message' output stream.
     * @return A new FrequenciesSet object according to options specified.
     * @throw Exception if an error occured.
     */
    static FrequenciesSet * getFrequenciesSet(
      const Alphabet * alphabet,
      const SiteContainer * data, 
      map<string, string> & params,
      const vector<double> & rateFreqs,
      const string & suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true) throw (Exception);

    /**
     * @brief Build a SubstitutionModelSet object according to options.
     *
     * This model set is meant to be used with non-homogeneous substitution models of sequence evolution.
     *
     * Recognized options are:
     * - number_of_models: the number of distinct SubstitutionModel to use.
     *
     * Then, for each of the models, the following information must be provided:
     * - model1='model name'
     * - model1.'parameters'='value'
     * Model names and parameters follow the same syntaxe as for the getSubstitutionModel method.
     * - model1.nodes='list of nodes id, separated by comas'.
     * And then
     * - model2=...
     * etc.
     *
     * All models must be fully specified, and at the end of the description, all nodes must be attributed to a model,
     * otherwise an exception is thrown.
     * 
     * Finally, this is also allowed for models to share one or several parameters.
     * for instance:
     * @code
     * model1=T92
     * model1.kappa=2.0
     * model1.theta=0.5
     * model2=T92
     * model2.kappa=model1.kappa
     * model2.theta=0.5
     * @endcode
     * In this case model1 and model2 with have their own and independent theta parameter, but only one kappa parameter will be used for both models.
     * Note that
     * @code
     * model1=T92
     * model1.kappa=2.0
     * model1.theta=0.5
     * model1.nodes=1,2,3
     * model2=T92
     * model2.kappa=model1.kappa
     * model2.theta=model1.theta
     * model2.nodes=4,5,6
     * @endcode
     * is equivalent to
     * @code
     * model1=T92
     * model1.kappa=2.0
     * model1.theta=0.5
     * model1.nodes=1,2,3,4,5,6
     * @endcode
     * but will require more memory and use more CPU, since some calculations will be performed twice.
     *
     * @param alphabet The alpabet to use in all models.
     * @param data     A pointer toward the SiteContainer for which the substitution model is designed.
     *                  The alphabet associated to the data must be of the same type as the one specified for the model.
     *                 May be equal to NULL, but in this cas use_observed_freq option will be unavailable.
     * @param params   The attribute map where options may be found.
     * @param suffix   A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose Print some info to the 'message' output stream.
     * @return A new SubstitutionModelSet object according to options specified.
     * @throw Exception if an error occured.
     */
    static SubstitutionModelSet * getSubstitutionModelSet(
      const Alphabet * alphabet,
      const SiteContainer * data, 
      map<string, string> & params,
      const string & suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true) throw (Exception);


    /**
     * @brief Build a DiscreteDistribution object according to options, with default parameter values.
     *
     * Options used are:
     * - rate_distribution = [constant|gamma], the distribution to use.
     *   Add '+invariant' for invariants.
     * Options depending on other options:
     * - If 'gamma' is selected:
     *   - rate_distribution.classes_number = the number of categories to be used.
     * - If'+invariant' is selected:
     *   - rate_distribution.p = the proportion of invariants.
     *
     * This function is mainly for internal usage, you're probably looking for the getRateDistribution function.
     *
     * @param params  The attribute map where options may be found.
     * @param constDistAllowed Tell if the constant distribution is allowed as a choice.
     * @param prefix  A prefix to be applied to each attribute name.
     * @param suffix  A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose Print some info to the 'message' output stream.
     * @return A new DiscreteDistribution object according to options specified.
     * @throw Exception if an error occured.
     */
    static DiscreteDistribution * getRateDistributionDefaultInstance(
      map<string, string> & params,
      bool constDistAllowed = true,
      const string & prefix = "",
      const string & suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true) throw (Exception);

    /**
     * @brief Set parameter initial values of a given rate distribution according to options.
     *
     * Parameters actually depends on the rate distirbution passed as argument.
     * See getRateDistribution for more information.
     *
     * This function is mainly for internal usage, you're probably looking for the getRateDistribution function.
     *
     * @param rDist The distribution to set.
     * @param params  The attribute map where options may be found.
     * @param prefix  A prefix to be applied to each attribute name.
     * @param suffix  A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose Print some info to the 'message' output stream.
     * @throw Exception if an error occured.
     */
    static void setRateDistributionParametersInitialValues(
      DiscreteDistribution * rDist,
      map<string, string> & params,
      const string & prefix = "",
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
      const string & suffix = "",
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
     * - optimization.topology.algorithm_nni.method = [fast,better,phyml]
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
      DiscreteRatesAcrossSitesClockTreeLikelihood * tl,
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

} //end of namespace bpp.

#endif  //_PHYLOGENETICSAPPLICATIONTOOLS_H_

