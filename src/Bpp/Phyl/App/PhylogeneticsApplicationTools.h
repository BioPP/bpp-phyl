//
// File: PhylogeneticsApplicationTools.h
// Created by: Julien Dutheil
// Created on: Fri Oct 21 16:49 2005
// from old file ApplicationTools.h created on Sun Dec 14 09:36:26 2003
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "../Tree/Tree.h"
#include "../Model/SubstitutionModel.h"
#include "../Model/SubstitutionModelSet.h"
#include "../Model/MixedSubstitutionModelSet.h"
#include "../Model/MarkovModulatedSubstitutionModel.h"
#include "../Likelihood/HomogeneousTreeLikelihood.h"
#include "../Likelihood/ClockTreeLikelihood.h"
#include "../Mapping/SubstitutionCount.h"
#include <Bpp/Text/TextTools.h>
#include <Bpp/Io/OutputStream.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/MultipleDiscreteDistribution.h>
#include <Bpp/Numeric/Function/Optimizer.h>

#include "../NewLikelihood/PhyloLikelihoods/SingleDataPhyloLikelihood.h"
#include "../NewLikelihood/PhyloLikelihoods/PhyloLikelihoodContainer.h"
#include "../NewLikelihood/SubstitutionProcessCollection.h"
#include "../NewLikelihood/SubstitutionProcessCollectionMember.h"
#include "../NewLikelihood/SubstitutionProcess.h"
#include "../NewLikelihood/SequenceEvolution.h"

// From SeqLib:
#include <Bpp/Seq/Container/SiteContainer.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

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
     * @throw Exception if an error occured.
     */
    static Tree* getTree(
      std::map<std::string, std::string>& params,
      const std::string& prefix = "input.",
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1) throw (Exception);
 
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
     * @throw Exception if an error occured.
     */
    static std::vector<Tree*> getTrees(
      std::map<std::string, std::string>& params,
      const std::string& prefix = "input.",
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1) throw (Exception);

    /**
     * @brief Build a list ofTree objects according to options.
     *
     * See the Bio++ Program Suite manual for a description of available options.
     *
     * @param params           The attribute map where options may be
     * found.
     * @param mSeq             A map of pointers of SiteContainers,
     * necessary in case of random trees
     * @param unparsedParams   A map of parameters (BrLen) that will
     * be aliased after.
     * @param prefix           A prefix to be applied to each attribute name.
     * @param suffix           A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose          Print some info to the 'message' output stream.
     * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
     * @return A new vector of Tree objects according to the specified options.
     * @throw Exception if an error occured.
     */
    
    static std::map<size_t, Tree*> getTrees(
      std::map<std::string, std::string>& params,
      const std::map<size_t, SiteContainer*>& mSeq,
      std::map<std::string, std::string>& unparsedParams,
      const std::string& prefix = "input.",
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1) throw (Exception);
    

    /**
     * @brief Build a SubstitutionModel object according to options.
     *
     * Creates a new substitution model object according to model description syntax
     * (see the Bio++ Progam Suite manual for a detailed description of this syntax). The
     * function also parses the parameter values and set them accordingly.
     *
     * @param alphabet         The alphabet to use in the model.
     * @param gCode            The genetic code to use (only for codon models, otherwise can be set to 0).
     *                         If set to NULL and a codon model is requested, an Exception will be thrown.
     * @param data             A pointer toward the SiteContainer for which the substitution model is designed.
     *                         The alphabet associated to the data must be of the same type as the one specified for the model.
     *                         May be equal to NULL, but in this case use_observed_freq option will be unavailable.
     * @param params           The attribute map where options may be
     * found.
     * @param unparsedparams   the map of the aliases between
     * parameters names.
     * @param suffix           A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose          Print some info to the 'message' output stream.
     * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
     * @return A new SubstitutionModel object according to options specified.
     * @throw Exception if an error occured.
     */
    
    static SubstitutionModel* getSubstitutionModel(
      const Alphabet* alphabet,
      const GeneticCode* gCode,
      const SiteContainer* data, 
      std::map<std::string, std::string>& params,
      std::map<std::string, std::string>& unparsedparams,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1) throw (Exception);
  
  

    /*
     * @brief The same as before, but with several models.
     *
     */
    
    static std::map<size_t, SubstitutionModel*> getSubstitutionModels(
      const Alphabet* alphabet,
      const GeneticCode* gCode,
      const std::map<size_t, SiteContainer*>& mData, 
      std::map<std::string, std::string>& params,
      std::map<std::string, std::string>& unparsedparams,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1) throw (Exception);


    /**
     * @brief Build a map of SubstitutionProcess objects according to options.
     *
     * Creates a new substitution model object according to model description syntax
     * (see the Bio++ Progam Suite manual for a detailed description of this syntax). The
     * function also parses the parameter values and set them accordingly.
     *
     * @param alphabet         The alphabet to use in the model.
     * @param gCode            The genetic code to use (only for codon models, otherwise can be set to 0).
     *                         If set to NULL and a codon model is requested, an Exception will be thrown.
     * @param mData            A map of pointers toward the
     * SiteContainers for which the substitution process wil be designed.
     *                         The alphabet associated to the data
     * must be of the same type as the one specified for the process. 
     * @param  params           The attribute map where options may be
     * found.
     $ @param unparsedparams   Unparsed params waiting for aliasing.
     * @param suffix           A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose          Print some info to the 'message' output stream.
     * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
     * @return A new SubstitutionModel object according to options specified.
     * @throw Exception if an error occured.
     */
    
    static std::map<size_t, SubstitutionProcess*> getSubstitutionProcesses(
      const Alphabet* alphabet,
      const GeneticCode* gCode,
      const std::map<size_t, SiteContainer*>& mData, 
      std::map<std::string, std::string>& params,
      std::map<std::string, std::string>& unparsedparams,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1) throw (Exception);


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
     * @param data   A pointer toward the SiteContainer for which the substitution model is designed.
     *               The alphabet associated to the data must be of the same type as the one specified for the model.
     *               May be equal to NULL, but in this case use_observed_freq option will be unavailable.
     * @param sharedParams (out) remote parameters will be recorded here.
     * @param verbose Print some info to the 'message' output stream.
     * @throw Exception if an error occured.
     */

    static void setSubstitutionModelParametersInitialValuesWithAliases(
      SubstitutionModel& model,
      std::map<std::string, std::string>& unparsedParameterValues,
      size_t modelNumber,
      const SiteContainer* data,
      std::map<std::string, std::string>& sharedParams,
      bool verbose) throw (Exception);

    /**
     * @brief Get A FrequenciesSet object for root frequencies (NH models) according to options.
     *
     * @param alphabet         The alpabet to use.
     * @param gCode            The genetic code to use (only for codon alphabets, otherwise can be set to 0).
     *                         If set to NULL and a codon frequencies set is requested, an Exception will be thrown.
     * @param data             A pointer toward the SiteContainer for which the substitution model is designed.
     *                         The alphabet associated to the data must be of the same type as the one specified for the model.
     *                         May be equal to NULL, but in this cas use_observed_freq option will be unavailable.
     * @param params           The attribute map where options may be
     * found.
     * @param sharedparams     The attribute map of aliases (out).
     * @param rateFreqs        A vector of rate catégories frequencies in case of a Markov Modulated Markov Model.
     *                         Ignored if a vector with size 0 is passed.
     * @param suffix           A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose          Print some info to the 'message' output stream.
     * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
     * @return A new FrequenciesSet object according to options specified.
     * @throw Exception if an error occured.
     */
    
    static FrequenciesSet* getRootFrequenciesSet(
      const Alphabet* alphabet,
      const GeneticCode* gCode,
      const SiteContainer* data, 
      std::map<std::string, std::string>& params,
      std::map<std::string, std::string>& sharedparams,
      const std::vector<double>& rateFreqs,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1) throw (Exception);

    /*
     * @brief The same, but with several FrequenciesSet.
     *
     */

    static std::map<size_t, FrequenciesSet*> getRootFrequenciesSets(
      const Alphabet* alphabet,
      const GeneticCode* gCode,
      const std::map<size_t, SiteContainer*>& mData, 
      std::map<std::string, std::string>& params,
      std::map<std::string, std::string>& sharedparams,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1) throw (Exception);

    /**
     * @brief Get A FrequenciesSet object according to options.
     *
     * @param alphabet         The alpabet to use.
     * @param gCode            The genetic code to use (only for codon alphabets, otherwise can be set to 0).
     *                         If set to NULL and a codon frequencies set is requested, an Exception will be thrown.
     * @param freqDescription  A string in the keyval syntaxe describing the frequency set to use.:if expand("%") == ""|browse confirm w|else|confirm w|endif
     * 
     * @param data             A pointer toward the SiteContainer for which the substitution model is designed.
     *                         The alphabet associated to the data must be of the same type as the one specified for the model.
     *                         May be equal to NULL, but in this cas use_observed_freq option will be unavailable.
     * @param sharedParams     (out) remote parameters will be recorded here.
     * @param rateFreqs        A vector of rate catégories frequencies in case of a Markov Modulated Markov Model.
     *                         Ignored if a vector with size 0 is passed.
     * @param verbose          Print some info to the 'message' output stream.
     * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
     * @return A new FrequenciesSet object according to options specified.
     * @throw Exception if an error occured.
     */

    static FrequenciesSet* getFrequenciesSet(
        const Alphabet* alphabet,
        const GeneticCode* gCode,
        const std::string& freqDescription,
        const SiteContainer* data, 
        std::map<std::string, std::string>& sharedParams,
        const std::vector<double>& rateFreqs,
        bool verbose = true,
        int warn = 1)
      throw (Exception);

    /**
     * @brief Get A FrequenciesSet object according to options.
     *
     * @param alphabet         The alpabet to use.
     * @param gCode            The genetic code to use (only for codon alphabets, otherwise can be set to 0).
     *                         If set to NULL and a codon frequencies set is requested, an Exception will be thrown.
     * @param freqDescription  A string in the keyval syntaxe describing the frequency set to use.:if expand("%") == ""|browse confirm w|else|confirm w|endif
     * 
     * @param data             A pointer toward the SiteContainer for which the substitution model is designed.
     *                         The alphabet associated to the data must be of the same type as the one specified for the model.
     *                         May be equal to NULL, but in this cas use_observed_freq option will be unavailable.
     * @param rateFreqs        A vector of rate catégories frequencies in case of a Markov Modulated Markov Model.
     *                         Ignored if a vector with size 0 is passed.
     * @param verbose          Print some info to the 'message' output stream.
     * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
     * @return A new FrequenciesSet object according to options specified.
     * @throw Exception if an error occured.
     */
    static FrequenciesSet* getFrequenciesSet(
        const Alphabet* alphabet,
        const GeneticCode* gCode,
        const std::string& freqDescription,
        const SiteContainer* data, 
        const std::vector<double>& rateFreqs,
        bool verbose = true,
        int warn = 1)
      throw (Exception)
    {
      std::map<std::string, std::string> sharedParams;
      return getFrequenciesSet(alphabet, gCode, freqDescription, data, sharedParams, rateFreqs, verbose, warn);
    }

    /**
     * @brief Gets a SubstitutionModelSet object according to options.
     *
     * See setSubstitutionModelSet and setMixedSubstitutionModelSet
     * methods.
     */
    
    static SubstitutionModelSet* getSubstitutionModelSet(
      const Alphabet* alphabet,
      const GeneticCode* gcode,
      const SiteContainer* data, 
      std::map<std::string, std::string>& params,
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
     * @param data             A pointer toward the SiteContainer for which the substitution model is designed.
     *                         The alphabet associated to the data must be of the same type as the one specified for the model.
     *                         May be equal to NULL, but in this cas use_observed_freq option will be unavailable.
     * @param params           The attribute map where options may be found.
     * @param suffix           A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose          Print some info to the 'message' output stream.
     * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
     * @throw Exception if an error occured.
     */
    static void setSubstitutionModelSet(
      SubstitutionModelSet& modelSet,
      const Alphabet* alphabet,
      const GeneticCode* gcode,
      const SiteContainer* data, 
      std::map<std::string, std::string>& params,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1);
    
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
     * @param pData A pointer toward the SiteContainer for which the
     *             substitution process is designed.
     *
     *             The alphabet associated to the data must be of the
     *             same type as the one specified for the model. May
     *             be equal to NULL, but in this cas use_observed_freq
     *             option will be unavailable.
     *
     * @param vTree A vector of pointers of Trees, used to set processes.
     * @param params   The attribute map where options may be found.
     * @param suffix   A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose Print some info to the 'message' output stream.
     * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
     * @throw Exception if an error occured.
     */
    
    static SubstitutionProcess* getSubstitutionProcess(
      const Alphabet* alphabet,
      const GeneticCode* gCode,
      const SiteContainer* pData, 
      const vector<Tree*>& vTree, 
      std::map<std::string, std::string>& params,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1);
    
    static SubstitutionProcessCollection* getSubstitutionProcessCollection(
      const Alphabet* alphabet,
      const GeneticCode* gCode,
      const map<size_t, Tree*>& mTree,
      const map<size_t, SubstitutionModel*>& mMod,
      const map<size_t, FrequenciesSet*>& mRootFreq,
      const map<size_t, DiscreteDistribution*>& mDist,
      map<string, string>& params,
      map<string, string>& unparsedparams,
      const string& suffix = "",
      bool suffixIsOptional  = true,
      bool verbose = true,
      int warn = 1);


    static void addSubstitutionProcessCollectionMember(
      SubstitutionProcessCollection* SubProColl, 
      size_t procNum,
      map<string, string>& params,
      bool verbose = true,
      int warn = 1);


    static std::map<size_t, SequenceEvolution*> getSequenceEvolutions(
      SubstitutionProcessCollection& SPC,
      map<string, string>& params,
      map<string, string>& unparsedParams,
      const string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1) throw (Exception);
    

    /**
     * @brief Build a Phylogeny container from parameters map.
     *
     * Default data compression is simple.
     *
     * The PhyloLikelihood with number 0 is a specific
     * PhyloLikelihood, with name "result" in the BppO file.
     *
     */
    
    static PhyloLikelihoodContainer* getPhyloLikelihoodContainer(
      SubstitutionProcessCollection& SPC,
      std::map<size_t, SequenceEvolution*>& mSeqEvol,
      const std::map<size_t, SiteContainer*>& mData,
      map<string, string>& params,
      const string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1) throw (Exception);
    
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
     * @param data             A pointer toward the SiteContainer for which the substitution model is designed.
     *                         The alphabet associated to the data must be of the same type as the one specified for the model.
     *                         May be equal to NULL, but in this cas use_observed_freq option will be unavailable.
     * @param params           The attribute map where options may be found.
     * @param suffix           A suffix to be applied to each attribute name.
     * @param suffixIsOptional Tell if the suffix is absolutely required.
     * @param verbose          Print some info to the 'message' output stream.
     * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
     * @throw Exception if an error occured.
     */
    static void completeMixedSubstitutionModelSet(
      MixedSubstitutionModelSet& mixedModelSet,
      const Alphabet* alphabet,
      const SiteContainer* data, 
      std::map<std::string, std::string>& params,
      const std::string& suffix = "",
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
     * @throw Exception if an error occured.
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
     * @throw Exception if an error occured.
     */
    static DiscreteDistribution* getRateDistribution(
      std::map<std::string, std::string>& params,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true)
      throw (Exception);

    /**
     * brief Same as before, but using several distributions.
     *
     */
    
    static std::map<size_t, DiscreteDistribution*> getRateDistributions(
      map<string, string>& params,
      const string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true)
      throw (Exception);

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
      std::map<std::string, std::string>& params,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1)
      throw (Exception);
    
    static PhyloLikelihood* optimizeParameters(
      PhyloLikelihood* lik,
      const ParameterList& parameters,
      std::map<std::string, std::string>& params,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1)
      throw (Exception);
    
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
     * @throw Exception        Any exception that may happen during the optimization process.
     */
    static void optimizeParameters(
      DiscreteRatesAcrossSitesClockTreeLikelihood* tl,
      const ParameterList& parameters,
      std::map<std::string, std::string>& params,
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      int warn = 1)
      throw (Exception);
    
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
      map<string, string>& params,
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
     * @throw Exception if an error occured.
     */
    static void writeTree(
      const TreeTemplate<Node>& tree,
      std::map<std::string, std::string>& params,
      const std::string& prefix = "output.",
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      bool checkOnly = false,
      int warn = 1) throw (Exception);
    
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
     * @param checkOnly        If this parameter is set to true, then all options are
     *                         checked and error messages sent, but no file is written.
     * @param warn             Set the warning level (0: always display warnings, >0 display warnings on demand).
     * @throw Exception if an error occured.
     */
    static void writeTrees(
      const std::vector<const Tree*>& trees,
      std::map<std::string, std::string>& params,
      const std::string& prefix = "output.",
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      bool checkOnly = false,
      int warn = 1) throw (Exception);

    static void writeTrees(
      const std::vector<const TreeTemplate<Node>* >& trees,
      std::map<std::string, std::string>& params,
      const std::string& prefix = "output.",
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      bool checkOnly = false,
      int warn = 1) throw (Exception);

    static void writeTrees(
      const SubstitutionProcessCollection& spc,
      std::map<std::string, std::string>& params,
      const std::string& prefix = "output.",
      const std::string& suffix = "",
      bool suffixIsOptional = true,
      bool verbose = true,
      bool checkOnly = false,
      int warn = 1) throw (Exception);

    
    /**
     * @brief Output a SubstitutionModel description to a file.
     *
     * @param model The model to serialize.
     * @param out   The stream where to print.
     * @param warn  Set the warning level (0: always display warnings, >0 display warnings on demand).
     */
    static void printParameters(const SubstitutionModel* model, OutputStream& out,int warn = 1);

    /**
     * @brief Output a SubstitutionModelSet description to a file.
     *
     * @param modelSet The model set to serialize.
     * @param out      The stream where to print.
     * @param warn  Set the warning level (0: always display warnings, >0 display warnings on demand).
     */
    static void printParameters(const SubstitutionModelSet* modelSet, OutputStream& out, int warn = 1);

    /**
     * @brief Output a SubstitutionProcess description to a file.
     *
     * @param process The process to serialize.
     * @param out      The stream where to print.
     * @param warn  Set the warning level (0: always display warnings, >0 display warnings on demand).
     */
    
    static void printParameters(const SubstitutionProcess* process, OutputStream& out, int warn = 1);

    /**
     * @brief Output information on the computation to a file.
     *
     * @param phylolike The phylolikelihood to serialize.
     * @param out      The stream where to print.
     * @param warn  Set the warning level (0: always display warnings, >0 display warnings on demand).
     */
    
    static void printAnalysisInformation(const PhyloLikelihood* phylolike, OutputStream& out, int warn = 1);

    static void printAnalysisInformation(const SingleDataPhyloLikelihood* phylolike, OutputStream& out, int warn = 1);

    /**
     * @brief Output a SubstitutionProcessCollection description to a file.
     *
     * @param collection The process collection to serialize.
     * @param out      The stream where to print.
     * @param warn  Set the warning level (0: always display warnings, >0 display warnings on demand).
     */

    static void printParameters(const SubstitutionProcessCollection* collection, OutputStream& out, int warn = 1);

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
     */
    static void printParameters(const DiscreteDistribution* rDist, OutputStream& out);

  };

} //end of namespace bpp.

#endif  //_PHYLOGENETICSAPPLICATIONTOOLS_H_

