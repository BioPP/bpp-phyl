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
     * @brief Build a SubstitutionModel object with default parameter values according to a keyval description.
     *
     * Check the Bio++ Program Suite documentation for a description of the syntax.
     * This function will resolve parameter aliasing, but will not assign initial values.
     * It is mainly for internal usage, you're probably looking for the getSubstitutionModel or getSubstitutionModelSet function.
     *
     * @param alphabet         The alpabet to use in the model.
     * @param modelDescription A string describing the model in the keyval syntax.
     * @param allowCovarions   Tell is a covarion model can be returned.
     * @param allowGaps        Tell is a gap model can be returned.
     * @param unparsedParameterValues [out] a map that will contain all the model parameters
     *                                names and their corresponding unparsed value, if they were found.
     * @param verbose Print some info to the 'message' output stream.
     * @return A new SubstitutionModel object according to options specified.
     * @throw Exception if an error occured.
     */
    static SubstitutionModel * getSubstitutionModelDefaultInstance(
      const Alphabet* alphabet,
      const string& modelDescription,
      map<string, string>& unparsedParameterValues,
      bool allowCovarions,
      bool allowGaps,
      bool verbose) throw (Exception);


    /**
     * @brief Set parameter initial values of a given model according to options.
     *
     * Parameters actually depends on the model passed as argument.
     * See getSubstitutionModel for more information.
     *
     * This function is mainly for internal usage, you're probably looking for the getSubstitutionModel or getSubstitutionModelSet function.
     *
     * @param model                   The model to set.
     * @param unparsedParameterValues A map that contains all the model parameters
     *                                names and their corresponding unparsed value, if they were found.
     * @param data   A pointer toward the SiteContainer for which the substitution model is designed.
     *               The alphabet associated to the data must be of the same type as the one specified for the model.
     *               May be equal to NULL, but in this case use_observed_freq option will be unavailable.
     * @param verbose Print some info to the 'message' output stream.
     * @throw Exception if an error occured.
     */
    static void setSubstitutionModelParametersInitialValues(
      SubstitutionModel* model,
      map<string, string>& unparsedParameterValues,
      const SiteContainer* data,
      bool verbose) throw (Exception);

    /**
     * @brief Build a SubstitutionModel object according to options.
     *
     * Creates a new substitution model object according to model description syntax
     * (see the Bio++ Progam Suite manual for a detailed description of this syntax). The
     * function also parses the parameter values and set them accordingly.
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
      const Alphabet* alphabet,
      const SiteContainer* data, 
      map<string, string>& params,
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
     * @param model                   The model to set.
     * @param unparsedParameterValues A map that contains all the model parameters
     *                                names and their corresponding unparsed value, if they were found.
     * @param modelPrefix The prefix to use to record prameters from this model.
     * @param data   A pointer toward the SiteContainer for which the substitution model is designed.
     *               The alphabet associated to the data must be of the same type as the one specified for the model.
     *               May be equal to NULL, but in this case use_observed_freq option will be unavailable.
     * @param existingParams (in/out) A map with already existing value that have been found in previous calls, and may be recalled here.
     *                       New parameters found here will be added.
     * @param specificParams (out) Parameters specific to this model will be recorded here.
     * @param sharedParams (out) remote parameters will be recorded here.
     * @param verbose Print some info to the 'message' output stream.
     * @throw Exception if an error occured.
     */
    static void setSubstitutionModelParametersInitialValues(
      SubstitutionModel* model,
      map<string, string>& unparsedParameterValues,
      const string& modelPrefix,
      const SiteContainer* data,
      map<string, double>& existingParams,
      vector<string>& specificParams,
      vector<string>& sharedParams,
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
     * @brief Build a rate distribution as a DiscreteDistribution object with default parameter values according to a keyval description.
     *
     * Any kind of discrete distribution with a mean of 1 can be used as a rate distribution.
     * Check the Bio++ Program Suite documentation for a description of the syntax.
     * This function will resolve parameter aliasing, but will not assign initial values.
     * It is mainly for internal usage, you're probably looking for the getRateDistribution function.
     *
     * @param distDescription         A string describing the model in the keyval syntax.
     * @param constDistAllowed        Tell if the constant distribution is allowed as a choice.
     * @param unparsedParameterValues [out] a map that will contain all the distribution parameters
     *                                names and their corresponding unparsed value, if they were found.
     * @param verbose                 Print some info to the 'message' output stream.
     * @return A new DiscreteDistribution object according to options specified.
     * @throw Exception if an error occured.
     */
    static DiscreteDistribution* getRateDistributionDefaultInstance(
      const string& distDescription,
      map<string, string>& unparsedParameterValues,
      bool constDistAllowed = true,
      bool verbose = true) throw (Exception);

    /**
     * @brief Set parameter initial values of a given rate distribution according to options.
     *
     * Parameters actually depends on the rate distribution passed as argument.
     * See getRateDistribution for more information.
     *
     * This function is mainly for internal usage, you're probably looking for the getRateDistribution function.
     *
     * @param rDist                   The distribution to set.
     * @param unparsedParameterValues a map that contains all the distribution parameters
     *                                names and their corresponding unparsed value, if they were found.
     * @param verbose                 Print some info to the 'message' output stream.
      * @throw Exception if an error occured.
     */
    static void setRateDistributionParametersInitialValues(
      DiscreteDistribution * rDist,
      map<string, string> & unparsedParameterValues,
      bool verbose = true) throw (Exception);

    /**
     * @brief Build a DiscreteDistribution object according to options.
     *
     * Creates a new rate distribution object according to distribution description syntax
     * (see the Bio++ Progam Suite manual for a detailed description of this syntax). The
     * function also parses the parameter values and set them accordingly.
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
     * @param tl               The TreeLikelihood function to optimize.
     * @param parameters       The initial list of parameters to optimize.
     *                         Use tl->getIndependentParameters() in order to estimate all parameters.
     * @param params           The attribute map where options may be found.
     * @param suffix           A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose          Print some info to the 'message' output stream.
     * @throw Exception        Any exception that may happen during the optimization process.
     * @return A pointer toward the final likelihood object.
     * This pointer may be the same as passed in argument (tl), but in some cases the algorithm
     * clone this object. We may change this bahavior in the future...
     * You hence should write something like
     * @code
     * tl = PhylogeneticsApplicationTools::optimizeParameters(tl, ...);
     * @endcode
     */
    static TreeLikelihood* optimizeParameters(
      TreeLikelihood* tl,
      const ParameterList& parameters,
      map<string, string>& params,
      const string& suffix = "",
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
     * @param tl               The ClockTreeLikelihood function to optimize.
     * @param parameters       The initial list of parameters to optimize.
     *                         Use tl->getIndependentParameters() in order to estimate all parameters.
     * @param params           The attribute map where options may be found.
     * @param suffix           A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose          Print some info to the 'message' output stream.
     * @throw Exception        Any exception that may happen during the optimization process.
     */
    static void optimizeParameters(
      DiscreteRatesAcrossSitesClockTreeLikelihood * tl,
      const ParameterList& parameters,
      map<string, string> & params,
      const string & suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true)
      throw (Exception);
  
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
   

    /**
     * @brief Output a SubstitutionModel description to a file.
     *
     * @param model The model to serialize.
     * @param out   The stream where to print.
     */
    static void printParameters(const SubstitutionModel* model, ostream& out);



    /**
     * @brief Output a SubstitutionModelSet description to a file.
     *
     * @param model The model set to serialize.
     * @param out   The stream where to print.
     */
    static void printParameters(const SubstitutionModelSet* modelSet, ostream& out);



    /**
     * @brief Output a DiscreteDistribution description to a file.
     *
     * @param model The model set to serialize.
     * @param out   The stream where to print.
     */
    static void printParameters(const DiscreteDistribution* rDist, ostream& out);

};

} //end of namespace bpp.

#endif  //_PHYLOGENETICSAPPLICATIONTOOLS_H_

