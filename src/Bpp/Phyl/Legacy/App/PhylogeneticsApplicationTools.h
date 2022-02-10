//
// File: PhylogeneticsApplicationTools.h
// Authors:
//   Julien Dutheil
// Created: 2005-10-21 16:49:00
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#ifndef BPP_PHYL_LEGACY_APP_PHYLOGENETICSAPPLICATIONTOOLS_H
#define BPP_PHYL_LEGACY_APP_PHYLOGENETICSAPPLICATIONTOOLS_H

#include <Bpp/Io/OutputStream.h>
#include <Bpp/Numeric/Function/Optimizer.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/MultipleDiscreteDistribution.h>
#include <Bpp/Text/TextTools.h>

#include "../../Model/MarkovModulatedSubstitutionModel.h"
#include "../../Model/SubstitutionModel.h"
#include "../../Tree/Tree.h"
#include "../Likelihood/ClockTreeLikelihood.h"
#include "../Likelihood/HomogeneousTreeLikelihood.h"
#include "../Model/MixedSubstitutionModelSet.h"
#include "../Model/SubstitutionModelSet.h"

// From SeqLib:
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/VectorProbabilisticSiteContainer.h>
#include <Bpp/Seq/AlphabetIndex/AlphabetIndex2.h>
#include <Bpp/Seq/Io/BppOAlphabetIndex2Format.h>

// From the STL:
#include <string>
#include <map>
#include <memory>
#include <vector>


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

class PhylogeneticsApplicationToolsOld
{
public:
  PhylogeneticsApplicationToolsOld();
  virtual ~PhylogeneticsApplicationToolsOld();

  /**
   * @brief Build a list ofTree objects according to options.
   *
   * See the Bio++ Program Suite manual for a description of available options.
   *
   * @param params           The attribute map where options may be
   * found.
   * @param mSeq             A map of pointers of AlignedValuesContainers,
   * necessary in case of random trees
   * @param unparsedParams   A map of parameters (BrLen) that will
   * be aliased after.
   * @param prefix           A prefix to be applied to each attribute name.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   * @return A new vector of Tree objects according to the specified options.
   */

  static std::map<size_t, Tree*> getTrees(
    const std::map<std::string, std::string>& params,
    const std::map<size_t, AlignedValuesContainer*>& mSeq,
    std::map<std::string, std::string>& unparsedParams,
    const std::string& prefix = "input.",
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

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
   * @param modelNumber The number of this model in the SubstitutionModelSet.
   * @param data A pointer toward an AlignedValuesContainer used for
   *             the initialization of substitution model when this
   *             data is needed (typically use_observed_freq option).
   *             The alphabet associated to the data must be of the
   *             same type as the one specified for the model. May be
   *             equal to NULL, but in this case will be unavailable.
   * @param sharedParams (out) remote parameters will be recorded here.
   * @param verbose Print some info to the 'message' output stream.
   */

  static void setSubstitutionModelParametersInitialValuesWithAliases(
    BranchModel& model,
    std::map<std::string, std::string>& unparsedParameterValues,
    size_t modelNumber,
    const AlignedValuesContainer* data,
    std::map<std::string, std::string>& sharedParams,
    bool verbose);

  /**
   * @brief Gets a SubstitutionModelSet object according to options.
   *
   * See setSubstitutionModelSet and setMixedSubstitutionModelSet
   * methods.
   */

  static SubstitutionModelSet* getSubstitutionModelSet(
    const Alphabet* alphabet,
    const GeneticCode* gcode,
    const AlignedValuesContainer* data,
    const std::map<std::string, std::string>& params,
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);


  /**
   * @brief Sets a SubstitutionModelSet object according to options.
   *
   * This model set is meant to be used with non-homogeneous substitution models of sequence evolution.
   *
   * Recognized options are:
   * - number_of_models: the number of distinct SubstitutionModel to use.
   *
   * Then, for each of the models, the following information must be provided:
   * - model1='model name(parameters'='value',...)
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
   * model1=T92(kappa=2.0, theta=0.5)
   * model2=T92(kappa=model1.kappa, theta=0.5)
   * @endcode
   * In this case model1 and model2 with have their own and independent theta parameter, but only one kappa parameter will be used for both models.
   * Note that
   * @code
   * model1=T92(kappa=2.0, theta=0.5)
   * model1.nodes=1,2,3
   * model2=T92(kappa= model1.kappa, theta=model1.theta)
   * model2.nodes=4,5,6
   * @endcode
   * is equivalent to
   * @code
   * model1=T92(kappa=2.0, theta=0.5)
   * model1.nodes=1,2,3,4,5,6
   * @endcode
   * but will require more memory and use more CPU, since some calculations will be performed twice.
   *
   * @param modelSet         The modified SubstitutionModelSet object according to options specified.
   * @param alphabet         The alpabet to use in all models.
   * @param gcode            The genetic code to use (only for codon models, otherwise can be set to 0).
   *                         If set to NULL and a codon model is requested, an Exception will be thrown.
   * @param data             A pointer toward the AlignedValuesContainer for which the substitution model is designed.
   *                         The alphabet associated to the data must be of the same type as the one specified for the model.
   *                         May be equal to NULL, but in this cas use_observed_freq option will be unavailable.
   * @param params           The attribute map where options may be found.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   */

  static void setSubstitutionModelSet(
    SubstitutionModelSet& modelSet,
    const Alphabet* alphabet,
    const GeneticCode* gcode,
    const AlignedValuesContainer* data,
    const std::map<std::string, std::string>& params,
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

  /**
   * @brief Complete a MixedSubstitutionModelSet object according to
   * options, given this model has already been filled through
   * setSubstitutionModelSet method.
   *
   * In addition, this method builds the allowed combinations of
   * submodels of the different mixed models.
   *
   * If none combination is given, then all possible submodels
   * combinations will be considered.
   *
   * The submodels dependencies are given a sets of combinations of
   * the mixed variables of the mixed models. For instance, if we
   * have:
   *
   * @code
   * model1=MixedModel(model=T92(kappa=Gamma(n=4), theta=0.5))
   * model2=MixedModel(model=T92(kappa=Gaussian(n=5), theta=Beta(n=3)))
   * @endcode
   *
   * In this case model1 is a mixture of 4 T92 submodels and model2
   * a mixture of 15 T92 submodels. These submodels are denoted with
   * the parameter name and the class number. For example, the
   * submodels of model1 are denoted model1[kappa_1], ...,
   * model1[kappa_4], and the submodels of model2 are denoted
   * model2[kappa_1,theta_1], ..., model2[kappa_5, theta_3].
   * Additionnaly, for instance, model2[kappa_2] denotes all the
   * submodels whose description has kappa_2.
   *
   * By default, when switching from model1 to model2, a site is
   * allowed to switch between any submodel of model1 and any
   * submodel of model2. If the only allowed combination is that a
   * site follows submodels model1(kappa_1) and
   * model2(kappa_3,theta_2), it is denoted:
   *
   * @code
   * site.allowedpaths= model1[kappa_1] & model2[kappa_3,theta_2]
   * @endcode
   *
   * With additional co()mbination saying that a site can follow
   * submodels model1[kappa_2] and any submodel of model2[kappa_3]
   * is denoted:
   *
   * @code
   * site.allowedpaths= model1[kappa_1] & model2[kappa_3,theta_2] |
   *                    model1[kappa_2] & model2[kappa_3]
   * @endcode
   *
   * See MixedSubstitutionModelSet description for further
   * information.
   *
   * @param mixedModelSet    The modified MixedSubstitutionModelSet object according to options specified.
   * @param alphabet         The alpabet to use in all models.
   * @param data             A pointer toward the AlignedValuesContainer for which the substitution model is designed.
   *                         The alphabet associated to the data must be of the same type as the one specified for the model.
   *                         May be equal to NULL, but in this cas use_observed_freq option will be unavailable.
   * @param params           The attribute map where options may be found.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   */

  static void completeMixedSubstitutionModelSet(
    MixedSubstitutionModelSet& mixedModelSet,
    const Alphabet* alphabet,
    const AlignedValuesContainer* data,
    const std::map<std::string, std::string>& params,
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

  /**
   * @brief Optimize parameters according to options.
   *
   * @param tl               The TreeLikelihood function to optimize.
   * @param parameters       The initial list of parameters to optimize.
   *                         Use tl->getIndependentParameters() in order to estimate all parameters.
   * @param params           The attribute map where options may be found.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   * @return A pointer toward the final likelihood object.
   * This pointer may be the same as passed in argument (tl), but in some cases the algorithm
   * clone this object. We may change this bahavior in the future...
   * You hence should write something like
   * @code
   * tl = PhylogeneticsApplicationToolsOld::optimizeParameters(tl, ...);
   * @endcode
   */

  static TreeLikelihood* optimizeParameters(
    TreeLikelihood* tl,
    const ParameterList& parameters,
    const std::map<std::string, std::string>& params,
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

  /**
   * @brief Optimize parameters according to options, with a molecular clock.
   *
   * @param tl               The ClockTreeLikelihood function to optimize.
   * @param parameters       The initial list of parameters to optimize.
   *                         Use tl->getIndependentParameters() in order to estimate all parameters.
   * @param params           The attribute map where options may be found.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   */

  static void optimizeParameters(
    DiscreteRatesAcrossSitesClockTreeLikelihood* tl,
    const ParameterList& parameters,
    const std::map<std::string, std::string>& params,
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

  /**
   * @brief Write a tree according to options.
   *
   * See the Bio++ Program Suite manual for a descriptio of all available options.
   *
   * @param trees            The trees to write.
   * @param params           The attribute map where options may be found.
   * @param prefix           A prefix to be applied to each attribute name.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param checkOnly        If this parameter is set to true, then all
   *                         options are checked and error messages
   *                         sent, but no file is written.
   * @param warn Set the warning level (0: always display warnings,
   * >0 display warnings on demand).
   */

  static void writeTrees(
    const std::vector<const Tree*>& trees,
    const std::map<std::string, std::string>& params,
    const std::string& prefix = "output.",
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    bool checkOnly = false,
    int warn = 1);

  static void writeTrees(
    const std::vector<const TreeTemplate<Node>* >& trees,
    const std::map<std::string, std::string>& params,
    const std::string& prefix = "output.",
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    bool checkOnly = false,
    int warn = 1);

  /**
   * @brief Output a SubstitutionModelSet description to a file.
   *
   * @param modelSet The model set to serialize.
   * @param out      The stream where to print.
   * @param warn  Set the warning level (0: always display warnings, >0 display warnings on demand).
   * @param withAlias outputs the alias names of the aliased
   *                      Parameters instead of the values (default
   *                      : true).
   */
  static void printParameters(const SubstitutionModelSet* modelSet, OutputStream& out, int warn = 1, bool withAlias = true);

};
} // end of namespace bpp.
#endif // BPP_PHYL_LEGACY_APP_PHYLOGENETICSAPPLICATIONTOOLS_H
