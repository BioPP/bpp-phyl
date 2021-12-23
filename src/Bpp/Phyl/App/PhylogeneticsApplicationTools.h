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

#ifndef BPP_PHYL_APP_PHYLOGENETICSAPPLICATIONTOOLS_H
#define BPP_PHYL_APP_PHYLOGENETICSAPPLICATIONTOOLS_H

#include <Bpp/Io/OutputStream.h>
#include <Bpp/Numeric/Function/Optimizer.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/MultipleDiscreteDistribution.h>
#include <Bpp/Phyl/Likelihood/DataFlow/DataFlow.h>
#include <Bpp/Text/TextTools.h>

#include "../Likelihood/ModelScenario.h"
#include "../Likelihood/PhyloLikelihoods/PhyloLikelihoodContainer.h"
#include "../Likelihood/PhyloLikelihoods/SetOfAlignedPhyloLikelihood.h"
#include "../Likelihood/PhyloLikelihoods/SingleDataPhyloLikelihood.h"
#include "../Likelihood/SequenceEvolution.h"
#include "../Likelihood/AutonomousSubstitutionProcess.h"
#include "../Likelihood/SubstitutionProcessCollection.h"
#include "../Likelihood/SubstitutionProcessCollectionMember.h"
#include "../Mapping/SubstitutionCount.h"
#include "../Model/MarkovModulatedSubstitutionModel.h"
#include "../Model/SubstitutionModel.h"
#include "../Tree/PhyloTree.h"

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

class PhylogeneticsApplicationTools
{
public:
  PhylogeneticsApplicationTools();
  virtual ~PhylogeneticsApplicationTools();


  /**
   * @brief Build a Tree object according to options.
   *
   * See the Bio++ Program Suite manual for a description of available options.
   *
   * @param params           The attribute map where options may be found.
   * @param prefix           A prefix to be applied to each attribute name.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   * @return A new Tree object according to the specified options.
   */
  static Tree* getTree(
    const std::map<std::string, std::string>& params,
    const std::string& prefix = "input.",
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

  /**
   * @brief Build a list ofTree objects according to options.
   *
   * See the Bio++ Program Suite manual for a description of available options.
   *
   * @param params           The attribute map where options may be found.
   * @param prefix           A prefix to be applied to each attribute name.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   * @return A new vector of Tree objects according to the specified options.
   */
  static std::vector<Tree*> getTrees(
    const std::map<std::string, std::string>& params,
    const std::string& prefix = "input.",
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

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

  static std::map<size_t, std::shared_ptr<PhyloTree> > getPhyloTrees(
    const std::map<std::string, std::string>& params,
    const std::map<size_t, AlignedValuesContainer*>& mSeq,
    std::map<std::string, std::string>& unparsedParams,
    const std::string& prefix = "input.",
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

  /**
   * @brief Build a SubstitutionModel object according to options.
   *
   * Creates a new substitution model object according to model description syntax
   * (see the Bio++ Progam Suite manual for a detailed description of this syntax). The
   * function also parses the parameter values and set them accordingly.
   *
   * @param alphabet         The alphabet to use in the model.
   *
   * @param gCode The genetic code to use (only for codon models,
   *              otherwise can be set to 0). If set to NULL and a
   *              codon model is requested, an Exception will be
   *              thrown.
   *
   * @param data A pointer toward an AlignedValuesContainer used for
   *             the initialization of substitution model when this
   *             data is needed (typically use_observed_freq option).
   *             The alphabet associated to the data must be of the
   *             same type as the one specified for the model. May be
   *             equal to NULL, but in this case will be unavailable.
   *
   * @param params The attribute map where options may be found.
   *
   * @param unparsedparams   the map of the aliases between
   * parameters names.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   * @return A new SubstitutionModel object according to options specified.
   */

  static SubstitutionModel* getSubstitutionModel(
    const Alphabet* alphabet,
    const GeneticCode* gCode,
    const AlignedValuesContainer* data,
    const std::map<std::string, std::string>& params,
    std::map<std::string, std::string>& unparsedparams,
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);


  static BranchModel* getBranchModel(
    const Alphabet* alphabet,
    const GeneticCode* gCode,
    const AlignedValuesContainer* data,
    const std::map<std::string, std::string>& params,
    std::map<std::string, std::string>& unparsedparams,
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

  /*
   * @brief The same as before, but with several models.
   *
   */

  static std::map<size_t, std::shared_ptr<BranchModel> > getBranchModels(
    const Alphabet* alphabet,
    const GeneticCode* gCode,
    const std::map<size_t, AlignedValuesContainer*>& mData,
    const std::map<std::string, std::string>& params,
    std::map<std::string, std::string>& unparsedparams,
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
   * @brief Get A FrequencySet object for root frequencies (NH models) according to options.
   *
   * @param alphabet         The alpabet to use.
   * @param gCode            The genetic code to use (only for codon alphabets, otherwise can be set to 0).
   *                         If set to NULL and a codon frequencies set is requested, an Exception will be thrown.
   * @param data   A pointer toward the AlignedValuesContainer for which the root frequencies are designed.
   *               The alphabet associated to the data must be of the same type as the one specified for the model.
   *               May be equal to NULL, but in this case use_observed_freq option will be unavailable.
   * @param params           The attribute map where options may be
   * found.
   * @param sharedparams     The attribute map of aliases (out).
   * @param rateFreqs        A vector of rate catÃÂ©gories frequencies in case of a Markov Modulated Markov Model.
   *                         Ignored if a vector with size 0 is passed.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   * @return A new FrequencySet object according to options specified.
   */

  static std::shared_ptr<FrequencySet> getRootFrequencySet(
    const Alphabet* alphabet,
    const GeneticCode* gCode,
    const AlignedValuesContainer* data,
    const std::map<std::string, std::string>& params,
    std::map<std::string, std::string>& sharedparams,
    const std::vector<double>& rateFreqs,
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

  /*
   * @brief The same, but with several FrequencySet.
   *
   */

  static std::map<size_t, std::shared_ptr<FrequencySet> > getRootFrequencySets(
    const Alphabet* alphabet,
    const GeneticCode* gCode,
    const std::map<size_t, AlignedValuesContainer*>& mData,
    const std::map<std::string, std::string>& params,
    std::map<std::string, std::string>& sharedparams,
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

  /**
   * @brief Get A FrequencySet object according to options.
   *
   * @param alphabet         The alpabet to use.
   * @param gCode            The genetic code to use (only for codon alphabets, otherwise can be set to 0).
   *                         If set to NULL and a codon frequencies set is requested, an Exception will be thrown.
   * @param freqDescription  A string in the keyval syntaxe describing the frequency set to use.:if expand("%") == ""|browse confirm w|else|confirm w|endif
   *
   * @param data A pointer toward an AlignedValuesContainer used for
   *             the initialization of the frequency set when this
   *             data is needed (typically use_observed_freq option).
   *             The alphabet associated to the data must be of the
   *             same type as the one specified for the frequency set. May be
   *             equal to NULL, but in this case will be unavailable.
   * @param sharedParams     (out) remote parameters will be recorded here.
   * @param rateFreqs        A vector of rate catÃÂ©gories frequencies in case of a Markov Modulated Markov Model.
   *                         Ignored if a vector with size 0 is passed.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   * @return A new FrequencySet object according to options specified.
   */

  static std::shared_ptr<FrequencySet> getFrequencySet(
    const Alphabet* alphabet,
    const GeneticCode* gCode,
    const std::string& freqDescription,
    const AlignedValuesContainer* data,
    std::map<std::string, std::string>& sharedParams,
    const std::vector<double>& rateFreqs,
    bool verbose = true,
    int warn = 1);

  /**
   * @brief Get A FrequencySet object according to options.
   *
   * @param alphabet         The alpabet to use.
   * @param gCode            The genetic code to use (only for codon alphabets, otherwise can be set to 0).
   *                         If set to NULL and a codon frequencies set is requested, an Exception will be thrown.
   * @param freqDescription  A string in the keyval syntaxe describing the frequency set to use.:if expand("%") == ""|browse confirm w|else|confirm w|endif
   *
   * @param data A pointer toward an AlignedValuesContainer used for
   *             the initialization of the frequency set when this
   *             data is needed (typically use_observed_freq option).
   *             The alphabet associated to the data must be of the
   *             same type as the one specified for the frequency set. May be
   *             equal to NULL, but in this case will be unavailable.
   * @param data             A pointer toward the AlignedValuesContainer for which the substitution model is designed.
   *                         The alphabet associated to the data must be of the same type as the one specified for the model.
   *                         May be equal to NULL, but in this cas use_observed_freq option will be unavailable.
   * @param rateFreqs        A vector of rate catÃÂ©gories frequencies in case of a Markov Modulated Markov Model.
   *                         Ignored if a vector with size 0 is passed.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   * @return A new FrequencySet object according to options specified.
   */
  static std::shared_ptr<FrequencySet> getFrequencySet(
    const Alphabet* alphabet,
    const GeneticCode* gCode,
    const std::string& freqDescription,
    const AlignedValuesContainer* data,
    const std::vector<double>& rateFreqs,
    bool verbose = true,
    int warn = 1)
  {
    std::map<std::string, std::string> sharedParams;
    return getFrequencySet(alphabet, gCode, freqDescription, data, sharedParams, rateFreqs, verbose, warn);
  }

  /**
   * @brief Build map of ModelPaths from a map of BranchModel.
   *
   * See the Bio++ Program Suite manual for a description of
   * available options.
   *
   * @param params           The attribute map where options may be found.
   * @param mModel           A map of shared pointers of BranchModels.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   *
   * Warning: the model FIRST described in a ModelPath will be the
   * leading model (ie which "decides" of the submodel
   * probabilities).
   *
   * @return A map of ModelPath shared_pointers objects according to
   * the specified options.
   */

  static map<size_t, std::shared_ptr<ModelPath> > getModelPaths(
    const std::map<std::string, std::string>& params,
    const map<size_t, std::shared_ptr<BranchModel> >& mModel,
    const string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

  /**
   * @brief Build map of ModelScenarios from a map of ModelPaths.
   *
   * See the Bio++ Program Suite manual for a description of
   * available options.
   *
   * @param params           The attribute map where options may be found.
   * @param mModelPath       A map of shared pointers of ModelPaths.
   * @param mModel           A map of shared pointers of BranchModels.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   *
   * @return A map of ModelScenarios shared_pointers objects
   * according to the specified options.
   */

  static map<size_t, std::shared_ptr<ModelScenario> > getModelScenarios(
    const std::map<std::string, std::string>& params,
    const map<size_t, std::shared_ptr<ModelPath> >& mModelPath,
    const map<size_t, std::shared_ptr<BranchModel> >& mModel,
    const string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

  /**
   * @brief Get a Register instance.
   *
   * @param regTypeDesc The description of the register.
   * @param stateMap The stateMap to use.
   * @param genCode when codon Alphabet, the genetic Code (otherwise,
   *                default : 0)
   * @param weights [out] AlphabetIndex2 pointer if "weights" argument
   *        is provided, null otherwise
   * @param distances [out] AlphabetIndex2 pointer if "distances"
   *        argument is provided, null otherwise
   * @param verbose if outputs  reading
   * @return A SubstitutionRegister object.
   */

  static SubstitutionRegister* getSubstitutionRegister(
    const std::string& regTypeDesc,
    const StateMap& stateMap,
    const GeneticCode* genCode,
    AlphabetIndex2*& weights,
    AlphabetIndex2*& distances,
    bool verbose = true);

  /**
   * @brief Sets a SubstitutionProcess object according to options.
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
   * All models must be fully specified, and at the end of the
   * description, all nodes must be attributed to a model, otherwise
   * an exception is thrown.
   *
   * Finally, this is also allowed for models to share one or
   * several parameters. for instance:
   * @code
   * model1=T92(kappa=2.0, theta=0.5)
   * model2=T92(kappa=model1.kappa, theta=0.5)
   * @endcode
   *
   * In this case model1 and model2 with have their own and
   * independent theta parameter, but only one kappa parameter will
   * be used for both models.
   *
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
   *
   * but will require more memory and use more CPU, since some
   * calculations will be performed twice.
   *
   * @param alphabet The alpabet to use in all models.
   * @param gCode            The genetic code to use (only for codon alphabets, otherwise can be set to 0).
   *                         If set to NULL and a codon frequencies set is requested, an Exception will be thrown.
   * @param data A pointer toward an AlignedValuesContainer used for
   *             the initialization of process set when this
   *             data is needed (typically use_observed_freq option).
   *             The alphabet associated to the data must be of the
   *             same type as the one specified for the process set. May be
   *             equal to NULL, but in this case will be unavailable.
   *
   * @param vTree A vector of pointers of Trees, used to set processes.
   * @param params   The attribute map where options may be found.
   * @param suffix   A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   */

  static AutonomousSubstitutionProcess* getSubstitutionProcess(
    const Alphabet* alphabet,
    const GeneticCode* gCode,
    const AlignedValuesContainer* pData,
    const vector<PhyloTree*>& vTree,
    const std::map<std::string, std::string>& params,
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

  static SubstitutionProcessCollection* getSubstitutionProcessCollection(
    const Alphabet* alphabet,
    const GeneticCode* gCode,
    const map<size_t, std::shared_ptr<PhyloTree> >& mTree,
    const map<size_t, std::shared_ptr<BranchModel> >& mMod,
    const map<size_t, std::shared_ptr<FrequencySet> >& mRootFreq,
    const map<size_t, std::shared_ptr<DiscreteDistribution> >& mDist,
    const map<size_t, std::shared_ptr<ModelScenario> >& mScen,
    const std::map<std::string, std::string>& params,
    map<string, string>& unparsedparams,
    const string& suffix = "",
    bool suffixIsOptional  = true,
    bool verbose = true,
    int warn = 1);

  static bool addSubstitutionProcessCollectionMember(
    SubstitutionProcessCollection* SubProColl,
    size_t procNum,
    const std::map<std::string, std::string>& params,
    bool verbose = true,
    int warn = 1);


  /**
   * @brief Build a map of Sequence Evolution, ie ways how sequence
   * evolve, which may use several processes.
   *
   */

  static std::map<size_t, SequenceEvolution*> getSequenceEvolutions(
    SubstitutionProcessCollection& SPC,
    const std::map<std::string, std::string>& params,
    map<string, string>& unparsedParams,
    const string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);


  /**
   * @brief Build a Phylogeny container from parameters map.
   *
   * Default data compression is simple.
   *
   * The PhyloLikelihood with number 0 is a specific
   * PhyloLikelihood, with name "result" in the BppO file.
   *
   */

  static std::shared_ptr<PhyloLikelihoodContainer> getPhyloLikelihoodContainer(
    Context& context,
    SubstitutionProcessCollection& SPC,
    std::map<size_t, SequenceEvolution*>& mSeqEvol,
    const std::map<size_t, AlignedValuesContainer*>& mData,
    const std::map<std::string, std::string>& params,
    const string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

  /**
   * @brief Build a multi-dimension distribution as a
   * MultipleDiscreteDistribution object with default parameter
   * values according to a keyval description.
   *
   * Check the Bio++ Program Suite documentation for a description of the syntax.
   * It is mainly for internal usage, you're probably looking for the getRateDistribution function.
   *
   * @param distDescription         A string describing the model in the keyval syntax.
   * @param unparsedParameterValues [out] a map that will contain all the distribution parameters
   *                                names and their corresponding unparsed value, if they were found.
   * @param verbose                 Print some info to the 'message' output stream.
   * @return A new MultipleDiscreteDistribution object according to options specified.
   */
  static MultipleDiscreteDistribution* getMultipleDistributionDefaultInstance(
    const std::string& distDescription,
    std::map<std::string, std::string>& unparsedParameterValues,
    bool verbose = true);

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
   */

  static DiscreteDistribution* getRateDistribution(
    const std::map<std::string, std::string>& params,
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true);

  /**
   * brief Same as before, but using several distributions.
   *
   */

  static std::map<size_t, std::shared_ptr<DiscreteDistribution> > getRateDistributions(
    const std::map<std::string, std::string>& params,
    const string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true);

  /**
   * @brief Optimize parameters according to options.
   *
   * @param lik              The PhyloLikelihood function to optimize.
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
   * tl = PhylogeneticsApplicationTools::optimizeParameters(tl, ...);
   * @endcode
   */

  static PhyloLikelihood* optimizeParameters(
    PhyloLikelihood* lik,
    const ParameterList& parameters,
    const std::map<std::string, std::string>& params,
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    int warn = 1);

  /**
   * @brief Check if parameter values are close to their definition boundary.
   *
   * This allows the detection of potential optimization issues.
   * A warning message will be output for each problematic parameter.
   *
   * @param pl A list of parameters. Parameters without constraint will be ignored.
   */
  static void checkEstimatedParameters(const ParameterList& pl);

  /**
   * @brief Get a SubstitutionCount instance.
   *
   * @param alphabet The alphabet to use.
   * @param model The model to use.
   * @param params The attribute map where options may be found.
   * @param suffix Optional suffix for command name.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn Set the warning level (0: always display warnings, >0 display warnings on demand).
   * @return A SubstitutionCount object.
   */
  static SubstitutionCount* getSubstitutionCount(
    const Alphabet* alphabet,
    const SubstitutionModel* model,
    const std::map<std::string, std::string>& params,
    string suffix = "",
    bool verbose = true,
    int warn = 1);

  /**
   * @brief Write a tree according to options.
   *
   * See the Bio++ Program Suite manual for a descriptio of all available options.
   *
   * @param tree    The tree to write.
   * @param params  The attribute map where options may be found.
   * @param prefix  A prefix to be applied to each attribute name.
   * @param suffix  A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose Print some info to the 'message' output stream.
   * @param checkOnly If this parameter is set to true, then all options are
   * checked and error messages sent, but no file is written.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   */

  static void writeTree(
    const TreeTemplate<Node>& tree,
    const std::map<std::string, std::string>& params,
    const std::string& prefix = "output.",
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    bool checkOnly = false,
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

  static void writePhyloTrees(
    const std::vector<const PhyloTree*>& trees,
    const std::map<std::string, std::string>& params,
    const std::string& prefix = "output.",
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    bool checkOnly = false,
    int warn = 1);

  /**
   * @brief Write a tree according to options.
   *
   * See the Bio++ Program Suite manual for a descriptio of all available options.
   *
   * @param spc              The SubstitutionProcessCollection of all objects
   * @param params           The attribute map where options may be found.
   * @param prefix           A prefix to be applied to each attribute name.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param checkOnly        If this parameter is set to true, then all
   *                         options are checked and error messages
   *                         sent, but no file is written.
   * @param withIds          If this parameter is set to true, ids
   *                         are added to the names of the nodes in output.
   * @param warn Set the warning level (0: always display warnings,
   * >0 display warnings on demand).
   */

  static void writePhyloTrees(
    const SubstitutionProcessCollection& spc,
    const std::map<std::string, std::string>& params,
    const std::string& prefix = "output.",
    const std::string& suffix = "",
    bool suffixIsOptional = true,
    bool verbose = true,
    bool checkOnly = false,
    bool withIds = false,
    int warn = 1);


  /**
   * @brief Output a SubstitutionModel description to a file.
   *
   * @param model The model to serialize.
   * @param out   The stream where to print.
   * @param warn  Set the warning level (0: always display warnings,
   *                      >0 display warnings on demand).
   */

  static void printParameters(const BranchModel* model, OutputStream& out, int warn = 1);

  /**
   * @brief Output a SubstitutionProcess description to a file.
   *
   * @param process The process to serialize.
   * @param out      The stream where to print.
   * @param warn  Set the warning level (0: always display warnings, >0 display warnings on demand).
   */

  static void printParameters(const SubstitutionProcess* process, OutputStream& out, int warn = 1);

  /**
   * @brief Output a SubstitutionProcessCollection description to a file.
   *
   * @param collection The process collection to serialize.
   * @param out      The stream where to print.
   * @param warn  Set the warning level (0: always display warnings, >0 display warnings on demand).
   * @param withAlias outputs the alias names of the aliased
   *                      Parameters instead of the values (default
   *                      : true).
   */

  static void printParameters(const SubstitutionProcessCollection* collection, OutputStream& out, int warn = 1, bool withAlias = true);

  /**
   * @brief Output the description of the  PhyloLikelihoods IN USE a
   * a PhyloLikelihoodContainer to a file.
   *
   * Starting from PhyloLikelihood with number 0 (with name
   * "result"), PhyloLikelihoods involved in the computation are
   * successively serialized.
   *
   * @param phylocont The PhyloLikelihoodContainer to serialize.
   * @param out       The stream where to print.
   * @param warn  Set the warning level (0: always display warnings, >0 display warnings on demand).
   */

  static void printParameters(const PhyloLikelihoodContainer& phylocont, OutputStream& out, int warn = 1);


  static void printParameters(const SingleDataPhyloLikelihood& phylolike, OutputStream& out, size_t nPhylo = 1, int warn = 1);

  /**
   * @brief Output a Sequence Evolution description to a file.
   *
   * @param evol The Sequence Evolution to serialize.
   * @param out       The stream where to print.
   * @param nEvol    The number of the evolution
   * @param warn  Set the warning level (0: always display warnings, >0 display warnings on demand).
   */

  static void printParameters(const SequenceEvolution* evol, OutputStream& out, size_t nEvol = 1, int warn = 1);

  /**
   * @brief Output a DiscreteDistribution description to a file.
   *
   * @param rDist The rate distribution to serialize.
   * @param out   The stream where to print.
   * @param withAlias outputs the alias names of the aliased
   *                      Parameters instead of the values (default
   *                      : true).
   */
  static void printParameters(const DiscreteDistribution* rDist, OutputStream& out, bool withAlias = true);

  /**
   * @brief Output information on the computation to a file.
   *
   * @param phylocont The phylolikelihood to serialize.
   * @param infosFile   The name of the file where to print.
   * @param warn  Set the warning level (0: always display warnings, >0 display warnings on demand).
   */

  static void printAnalysisInformation(const PhyloLikelihoodContainer& phylocont, const std::string& infosFile, int warn = 1);

  static void printAnalysisInformation(const SingleDataPhyloLikelihood& phylolike, const std::string& infosFile, int warn = 1);

  static void printAnalysisInformation(const SetOfAlignedPhyloLikelihood& sOAP, const std::string& infosFile, int warn = 1);
};
} // end of namespace bpp.
#endif // BPP_PHYL_APP_PHYLOGENETICSAPPLICATIONTOOLS_H
