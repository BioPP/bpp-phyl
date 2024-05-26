// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_APP_PHYLOGENETICSAPPLICATIONTOOLS_H
#define BPP_PHYL_APP_PHYLOGENETICSAPPLICATIONTOOLS_H

#include "../Likelihood/ModelScenario.h"
#include "../Likelihood/PhyloLikelihoods/PhyloLikelihoodContainer.h"
#include "../Likelihood/PhyloLikelihoods/AlignedPhyloLikelihoodSet.h"
#include "../Likelihood/PhyloLikelihoods/SingleDataPhyloLikelihood.h"
#include "../Likelihood/SequenceEvolution.h"
#include "../Likelihood/AutonomousSubstitutionProcess.h"
#include "../Likelihood/SubstitutionProcessCollection.h"
#include "../Likelihood/SubstitutionProcessCollectionMember.h"
#include "../Mapping/SubstitutionCount.h"
#include "../Model/MarkovModulatedSubstitutionModel.h"
#include "../Model/SubstitutionModel.h"
#include "../Tree/PhyloTree.h"
#include "../Likelihood/DataFlow/DataFlow.h"

// From bpp-core:
#include <Bpp/Io/OutputStream.h>
#include <Bpp/Numeric/Function/Optimizer.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/MultipleDiscreteDistribution.h>
#include <Bpp/Text/TextTools.h>

// From bpp-seq:
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
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
   *
   * @return A new Tree object according to the specified options.
   */
  static std::unique_ptr<Tree> getTree(
      const std::map<std::string, std::string>& params,
      const std::string& prefix = "input.",
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1);

  /**
   * @brief Build a list of Tree objects according to options.
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

  static std::vector<std::unique_ptr<Tree>> getTrees(
      const std::map<std::string, std::string>& params,
      const std::string& prefix = "input.",
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1);

  /**
   * @brief Build a map of <number,PhyloTree> according to options.
   *
   * See the Bio++ Program Suite manual for a description of available options.
   *
   * @param params           The attribute map where options may be found.
   * @param mSeq             A map of pointers of AlignmentDataInterface<std::string>s,
   * necessary in case of random trees.
   * @param unparsedParams A map of parameters (BrLen) that will be
   * aliased after.
   * @param prefix           A prefix to be applied to each attribute name.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   *
   * @return A new vector of Tree objects according to the specified options.
   */

  static std::map<size_t, std::shared_ptr<PhyloTree>> getPhyloTrees(
      const std::map<std::string, std::string>& params,
      const std::map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mSeq,
      std::map<std::string, std::string>& unparsedParams,
      const std::string& prefix = "input.",
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1);

  /**
   * @brief Build a BranchModel (most general class of branch models)
   * object according to options.
   *
   * Creates a new branch model object according to model description
   * syntax (see the Bio++ Program Suite manual for a detailed
   * description of this syntax). The function also parses the
   * parameter values and set them accordingly.
   *
   * @param alphabet   The alphabet to use in the model.
   * @param gCode The genetic code to use (only for codon models,
   *              otherwise can be set to 0). If set to NULL and a
   *              codon model is requested, an Exception will be
   *              thrown.
   * @param data A pointer toward an AlignmentDataInterface<std::string> used for
   *             the initialization of substitution model when this
   *             data is needed (typically use_observed_freq option).
   *             The alphabet associated to the data must be of the
   *             same type as the one specified for the model. May be
   *             equal to NULL, but in this case will be unavailable.
   * @param params The attribute map where options may be found.
   * @param unparsedparams   the map of the aliases between
   * parameters names.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   *
   * @return A new BranchModel object.
   */
  static std::unique_ptr<BranchModelInterface> getBranchModel(
      std::shared_ptr<const Alphabet> alphabet,
      std::shared_ptr<const GeneticCode> gCode,
      std::shared_ptr<const AlignmentDataInterface> data,
      const std::map<std::string, std::string>& params,
      std::map<std::string, std::string>& unparsedparams,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1);

  /**
   * @brief Build a SubstitutionModel object according to options (ie
   * a BranchModel with a generator).
   *
   * Creates a new substitution model object according to model description syntax
   * (see the Bio++ Program Suite manual for a detailed description of this syntax). The
   * function also parses the parameter values and set them accordingly.
   *
   * @param alphabet   The alphabet to use in the model.
   * @param gCode The genetic code to use (only for codon models,
   *              otherwise can be set to 0). If set to NULL and a
   *              codon model is requested, an Exception will be
   *              thrown.
   * @param data A pointer toward an AlignmentDataInterface<std::string> used for
   *             the initialization of substitution model when this
   *             data is needed (typically use_observed_freq option).
   *             The alphabet associated to the data must be of the
   *             same type as the one specified for the model. May be
   *             equal to NULL, but in this case will be unavailable.
   * @param params The attribute map where options may be found.
   * @param unparsedparams   the map of the aliases between
   * parameters names.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   *
   * @return A new SubstitutionModel object.
   */
  static std::unique_ptr<SubstitutionModelInterface> getSubstitutionModel(
      std::shared_ptr<const Alphabet> alphabet,
      std::shared_ptr<const GeneticCode> gCode,
      std::shared_ptr<const AlignmentDataInterface> data,
      const std::map<std::string, std::string>& params,
      std::map<std::string, std::string>& unparsedparams,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1);


  /**
   * @brief The same as before, but returns a map <number, branch model>.
   */
  static std::map<size_t, std::unique_ptr<BranchModelInterface>> getBranchModels(
      std::shared_ptr<const Alphabet> alphabet,
      std::shared_ptr<const GeneticCode> gCode,
      const std::map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mData,
      const std::map<std::string, std::string>& params,
      std::map<std::string, std::string>& unparsedparams,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1);

  /**
   * @brief Get A FrequencySet object for root frequencies (NH models) according to options.
   *
   * @param alphabet         The alphabet to use.
   * @param gCode            The genetic code to use (only for codon alphabets, otherwise can be set to 0).
   *                         If set to NULL and a codon frequencies set is requested, an Exception will be thrown.
   * @param mData A map <size_t, AlignmentDataInterface> used for
   *             the initialization of the frequency set when this
   *             data is needed (typically use_observed_freq option).
   *             The alphabet associated to the data must be of the
   *             same type as the one specified for the frequency set. May be
   *             equal to NULL, but in this case will be unavailable.
   * @param nData the number of the AlignmentDataInterface used in the mData map.
   *             0 if this number has not been defined.
   * @param params           The attribute map where options may be found.
   * @param sharedparams     The attribute map of aliases (out).
   * @param rateFreqs        A vector of rate categories frequencies in case of a Markov Modulated Markov Model.
   *                         Ignored if a vector with size 0 is passed.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   *
   * @return A new FrequencySet object according to options specified.
   */
  static std::unique_ptr<FrequencySetInterface> getRootFrequencySet(
      std::shared_ptr<const Alphabet> alphabet,
      std::shared_ptr<const GeneticCode> gCode,
      const std::map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mData,
      size_t nData, 
      const std::map<std::string, std::string>& params,
      std::map<std::string, std::string>& sharedparams,
      const std::vector<double>& rateFreqs,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1);

  /**
   * @brief The same, but returns a map <number, shared_ptr<FrequencySetInterface>>.
   */
  static std::map<size_t, std::unique_ptr<FrequencySetInterface>> getRootFrequencySets(
      std::shared_ptr<const Alphabet> alphabet,
      std::shared_ptr<const GeneticCode> gCode,
      const std::map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mData,
      const std::map<std::string, std::string>& params,
      std::map<std::string, std::string>& sharedparams,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1);

  /**
   * @brief Get A FrequencySet object according to options.
   *
   * @param alphabet         The alphabet to use.
   * @param gCode            The genetic code to use (only for codon alphabets, otherwise can be set to 0).
   *                         If set to NULL and a codon frequencies set is requested, an Exception will be thrown.
   * @param freqDescription  A string in the keyval syntax describing the frequency set to use.
   * @param data A pointer toward an AlignmentDataInterface used for
   *             the initialization of the frequency set when this
   *             data is needed (typically use_observed_freq option).
   *             The alphabet associated to the data must be of the
   *             same type as the one specified for the frequency set. May be
   *             equal to NULL, but in this case will be unavailable.
   * @param sharedParams     (out) remote parameters will be recorded here.
   * @param rateFreqs        A vector of rate categories frequencies in case of a Markov Modulated Markov Model.
   *                         Ignored if a vector with size 0 is passed.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   *
   * @return A new FrequencySet object according to options specified.
   */
  
  static std::unique_ptr<FrequencySetInterface> getFrequencySet(
      std::shared_ptr<const Alphabet> alphabet,
      std::shared_ptr<const GeneticCode> gCode,
      const std::string& freqDescription,
      std::shared_ptr<const AlignmentDataInterface> data,
      std::map<std::string, std::string>& sharedParams,
      const std::vector<double>& rateFreqs,
      bool verbose = true,
      int warn = 1)
  {
    std::map<size_t, std::shared_ptr<const AlignmentDataInterface>> mData;
    mData[1]=data;

    return getFrequencySet(alphabet, gCode, freqDescription, mData, 1, sharedParams, rateFreqs, verbose, warn);
  }

  /**
   * @brief Get A FrequencySet object according to options.
   *
   * @param alphabet         The alphabet to use.
   * @param gCode            The genetic code to use (only for codon alphabets, otherwise can be set to 0).
   *                         If set to NULL and a codon frequencies set is requested, an Exception will be thrown.
   * @param freqDescription  A string in the keyval syntax describing the frequency set to use.
   * @param mData A map <size_t, AlignmentDataInterface> used for
   *             the initialization of the frequency set when this
   *             data is needed (typically use_observed_freq option).
   *             The alphabet associated to the data must be of the
   *             same type as the one specified for the frequency set. May be
   *             equal to NULL, but in this case will be unavailable.
   * @param nData the number of the AlignmentDataInterface used in the mData map.
   *             0 if this number has not been defined.
   * @param sharedParams     (out) remote parameters will be recorded here.
   * @param rateFreqs        A vector of rate categories frequencies in case of a Markov Modulated Markov Model.
   *                         Ignored if a vector with size 0 is passed.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   *
   * @return A new FrequencySet object according to options specified.
   */
  static std::unique_ptr<FrequencySetInterface> getFrequencySet(
    std::shared_ptr<const Alphabet> alphabet,
    std::shared_ptr<const GeneticCode> gCode,
    const std::string& freqDescription,
    const std::map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mData,
    size_t nData,
    std::map<std::string, std::string>& sharedParams,
    const std::vector<double>& rateFreqs,
    bool verbose = true,
    int warn = 1);

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
   * leading model (ie on which the submodel probabilities are
   * chosen).
   *
   * @return A map of ModelPath shared_pointers objects according to
   * the specified options.
   */
  static map<size_t, std::unique_ptr<ModelPath>> getModelPaths(
      const std::map<std::string, std::string>& params,
      const map<size_t, std::shared_ptr<BranchModelInterface>>& mModel,
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
  static map<size_t, std::unique_ptr<ModelScenario>> getModelScenarios(
      const std::map<std::string, std::string>& params,
      const map<size_t, std::shared_ptr<ModelPath>>& mModelPath,
      const map<size_t, std::shared_ptr<BranchModelInterface>>& mModel,
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
  static std::unique_ptr<SubstitutionRegisterInterface> getSubstitutionRegister(
      const std::string& regTypeDesc,
      std::shared_ptr<const StateMapInterface> stateMap,
      std::shared_ptr<const GeneticCode> genCode,
      std::shared_ptr<AlphabetIndex2>& weights,
      std::shared_ptr<AlphabetIndex2>& distances,
      bool verbose = true);

  /**
   * @brief Sets a SubstitutionProcess object according to options.
   *
   * See the Bio++ Program Suite manual for a description of
   * available options.
   *
   * @param alphabet The alphabet to use in all models.
   * @param gCode The genetic code to use (only for codon alphabets,
   *              otherwise can be set to 0). If set to NULL and a
   *              codon frequencies set is requested, an Exception
   *              will be thrown.
   * @param data an AlignmentDataInterface<std::string> used for
   *             the initialization of process set when this
   *             data is needed (typically use_observed_freq option).
   *             The alphabet associated to the data must be of the
   *             same type as the one specified for the process set. May be
   *             equal to NULL, but in this case will be unavailable.
   * @param vTree A vector of pointers of Trees, used to set processes.
   * @param params   The attribute map where options may be found.
   * @param suffix   A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   *
   * @return A pointer to a AutonomousSubstitutionProcess.
   */
  static std::unique_ptr<AutonomousSubstitutionProcessInterface> getSubstitutionProcess(
      std::shared_ptr<const Alphabet> alphabet,
      std::shared_ptr<const GeneticCode> gCode,
      std::shared_ptr<const AlignmentDataInterface> data,
      const vector<std::shared_ptr<PhyloTree>>& vTree,
      const std::map<std::string, std::string>& params,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1);

  /**
   * @brief Builds a  SubstitutionProcessCollection from many maps of relevant objects.
   *
   * @param alphabet The alphabet to use in all models.
   * @param gCode The genetic code to use (only for codon alphabets,
   *              otherwise can be set to 0). If set to NULL and a
   *              codon frequencies set is requested, an Exception
   *              will be thrown.
   * @param mTree Map of the PhyloTrees
   * @param mMod Map of Branch Models
   * @param mRootFreq Map of FrequencySet
   * @param mDist Map of Rate Distributions
   * @param mScen Map of Scenarios
   * @param params   The attribute map where options may be found.
   * @param unparsedparams   map of the attributes that will not be parsed there
   * @param suffix   A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   *
   * @return A shared pointer to a SubstitutionProcessCollection.
   */
  static std::unique_ptr<SubstitutionProcessCollection> getSubstitutionProcessCollection(
      std::shared_ptr<const Alphabet> alphabet,
      std::shared_ptr<const GeneticCode> gCode,
      const map<size_t, std::shared_ptr<PhyloTree>>& mTree,
      const map<size_t, std::shared_ptr<BranchModelInterface>>& mMod,
      const map<size_t, std::shared_ptr<FrequencySetInterface>>& mRootFreq,
      const map<size_t, std::shared_ptr<DiscreteDistributionInterface>>& mDist,
      const map<size_t, std::shared_ptr<ModelScenario>>& mScen,
      const std::map<std::string, std::string>& params,
      map<string, string>& unparsedparams,
      const string& suffix = "",
      bool suffixIsOptional  = true,
      bool verbose = true,
      int warn = 1);


  /**
   * @brief Adds a SubstitutionProcessCollectionMember to a Collection, under a given number.
   *
   * @param SubProColl SubstitutionProcessCollection where the
   * SubstitutionProcessCollectionMember will be added.
   * @param procNum the number of the process
   * @param params   The attribute map where options may be found.
   * @param verbose Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   *
   * @return if the addition has been successful.
   */
  static bool addSubstitutionProcessCollectionMember(
      SubstitutionProcessCollection& SubProColl,
      size_t procNum,
      const std::map<std::string, std::string>& params,
      bool verbose = true,
      int warn = 1);


  /**
   * @brief Build a map of Sequence Evolution, ie ways how sequence
   * evolve, which may use several processes.
   */
  static std::map<size_t, std::unique_ptr<SequenceEvolution>> getSequenceEvolutions(
      std::shared_ptr<SubstitutionProcessCollection> SPC,
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
   */
  static std::shared_ptr<PhyloLikelihoodContainer> getPhyloLikelihoodContainer(
      Context& context,
      std::shared_ptr<SubstitutionProcessCollection> SPC,
      std::map<size_t, std::shared_ptr<SequenceEvolution>>& mSeqEvol,
      const std::map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mData,
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
   * (see the Bio++ Program Suite manual for a detailed description of this syntax). The
   * function also parses the parameter values and set them accordingly.
   *
   * @param params  The attribute map where options may be found.
   * @param suffix  A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose Print some info to the 'message' output stream.
   * @return A new DiscreteDistribution object according to options specified.
   */
  static std::unique_ptr<DiscreteDistributionInterface> getRateDistribution(
      const std::map<std::string, std::string>& params,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true);

  /**
   * brief Same as before, but using several distributions.
   *
   */

  static std::map<size_t, std::shared_ptr<DiscreteDistributionInterface>> getRateDistributions(
      const std::map<std::string, std::string>& params,
      const string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true);

  /**
   * @brief Optimize parameters according to options.
   *
   * @param lik              The PhyloLikelihood function to optimize.
   * @param params           The attribute map where options may be found.
   * @param suffix           A suffix to be applied to each attribute name.
   * @param suffixIsOptional Tell if the suffix is absolutely required.
   * @param verbose          Print some info to the 'message' output stream.
   * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
   *
   * @return A pointer toward the final likelihood object.
   *
   * This pointer may be the same as passed in argument (tl), but in some cases the algorithm
   * clone this object. We may change this behavior in the future...
   * You hence should write something like
   * @code
   * tl = PhylogeneticsApplicationTools::optimizeParameters(tl, ...);
   * @endcode
   */
  static std::shared_ptr<PhyloLikelihoodInterface> optimizeParameters(
      std::shared_ptr<PhyloLikelihoodInterface> lik,
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
  static std::unique_ptr<SubstitutionCountInterface> getSubstitutionCount(
      std::shared_ptr<const Alphabet> alphabet,
      std::shared_ptr<const SubstitutionModelInterface> model,
      const std::map<std::string, std::string>& params,
      string suffix = "",
      bool verbose = true,
      int warn = 1);

  /**
   * @brief Output methods:
   *
   * @{
   *
   */

  /**
   * @brief Write a tree according to options.
   *
   * See the Bio++ Program Suite manual for a description of all available options.
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
   * See the Bio++ Program Suite manual for a description of all available options.
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
   * See the Bio++ Program Suite manual for a description of all available options.
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
  static void printParameters(const BranchModelInterface& model, OutputStream& out, int warn = 1);

  /**
   * @brief Output a SubstitutionProcess description to a file.
   *
   * @param process The process to serialize.
   * @param out      The stream where to print.
   * @param warn  Set the warning level (0: always display warnings, >0 display warnings on demand).
   */
  static void printParameters(const SubstitutionProcessInterface& process, OutputStream& out, int warn = 1);

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
  static void printParameters(const SubstitutionProcessCollection& collection, OutputStream& out, int warn = 1, bool withAlias = true);

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


  static void printParameters(const SingleDataPhyloLikelihoodInterface& phylolike, OutputStream& out, size_t nPhylo = 1, int warn = 1);

  /**
   * @brief Output a Sequence Evolution description to a file.
   *
   * @param evol The Sequence Evolution to serialize.
   * @param out       The stream where to print.
   * @param nEvol    The number of the evolution
   * @param warn  Set the warning level (0: always display warnings, >0 display warnings on demand).
   */

  static void printParameters(const SequenceEvolution& evol, OutputStream& out, size_t nEvol = 1, int warn = 1);

  /**
   * @brief Output a DiscreteDistribution description to a file.
   *
   * @param rDist The rate distribution to serialize.
   * @param out   The stream where to print.
   * @param withAlias outputs the alias names of the aliased
   *                      Parameters instead of the values (default
   *                      : true).
   */
  static void printParameters(const DiscreteDistributionInterface& rDist, OutputStream& out, bool withAlias = true);

  /**
   * @brief Output information on the computation to a file.
   *
   * @param phylocont The phylolikelihood to serialize.
   * @param infosFile   The name of the file where to print.
   * @param warn  Set the warning level (0: always display warnings, >0 display warnings on demand).
   */

  static void printAnalysisInformation(const PhyloLikelihoodContainer& phylocont, const std::string& infosFile, int warn = 1);

  static void printAnalysisInformation(const SingleDataPhyloLikelihoodInterface& phylolike, const std::string& infosFile, int warn = 1);

  static void printAnalysisInformation(const AlignedPhyloLikelihoodSetInterface& sOAP, const std::string& infosFile, int warn = 1);

  /**
   * @}
   */
  template<class T>
  static std::map<size_t, std::shared_ptr<T>> uniqueToSharedMap(std::map<size_t, std::unique_ptr<T>>& mu)
  {
    std::map<size_t, std::shared_ptr<T>> ms;
    for (auto& it : mu)
    {
      std::shared_ptr<T> pt = std::move(it.second);
      ms.emplace(it.first, pt);
    }
    return ms;
  }
};
} // end of namespace bpp.
#endif // BPP_PHYL_APP_PHYLOGENETICSAPPLICATIONTOOLS_H
