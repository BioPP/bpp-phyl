//
// File: PhylogeneticsApplicationTools.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 21 16:49 2005
// from old file ApplicationTools.cpp created on Sun Dec 14 09:36:26 2003
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

#include "PhylogeneticsApplicationTools.h"
#include "../Model/SubstitutionModel.h"
#include "../Model/MixedTransitionModel.h"
#include "../Model/WrappedModel.h"
#include "../Model/Protein/Coala.h"
#include "../Model/FrequenciesSet/MvaFrequenciesSet.h"
#include "../Likelihood/TreeLikelihood.h"
#include "../Mapping/LaplaceSubstitutionCount.h"
#include "../Mapping/UniformizationSubstitutionCount.h"
#include "../Mapping/DecompositionSubstitutionCount.h"
#include "../Mapping/NaiveSubstitutionCount.h"
#include "../Mapping/OneJumpSubstitutionCount.h"
#include "../OptimizationTools.h"
#include "../Tree.h"
#include "../Io/BppOTreeReaderFormat.h"
#include "../Io/BppOMultiTreeReaderFormat.h"
#include "../Io/BppOTreeWriterFormat.h"
#include "../Io/BppOMultiTreeWriterFormat.h"
#include "../Io/BppOSubstitutionModelFormat.h"
#include "../Io/BppOTransitionModelFormat.h"
#include "../Io/BppOFrequenciesSetFormat.h"
#include "../Io/BppORateDistributionFormat.h"

// From bpp-core
#include <Bpp/Io/BppODiscreteDistributionFormat.h>
#include <Bpp/Io/BppOParametrizableFormat.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/Prob/DirichletDiscreteDistribution.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

using namespace bpp;

// From the STL:
#include <fstream>
#include <memory>
#include <set>
#include <vector>

using namespace std;


/*************************************************************/
/*****************  TREES ************************************/
/*************************************************************/


/******************************************************************************/

Tree* PhylogeneticsApplicationTools::getTree(
  map<string, string>& params,
  const string& prefix,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  string format = ApplicationTools::getStringParameter(prefix + "tree.format", params, "Newick", suffix, suffixIsOptional, warn);
  string treeFilePath = ApplicationTools::getAFilePath(prefix + "tree.file", params, true, true, suffix, suffixIsOptional, "none", warn);

  BppOTreeReaderFormat bppoReader(warn);
  unique_ptr<ITree> iTree(bppoReader.read(format));
  if (verbose)
  {
    ApplicationTools::displayResult("Input tree file " + suffix, treeFilePath);
    ApplicationTools::displayResult("Input tree format " + suffix, iTree->getFormatName());
  }
  Tree* tree = iTree->readTree(treeFilePath);
  return tree;
}

/******************************************************************************/

vector<Tree*> PhylogeneticsApplicationTools::getTrees(
  map<string, string>& params,
  const string& prefix,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  string format = ApplicationTools::getStringParameter(prefix + "trees.format", params, "Newick", suffix, suffixIsOptional, warn);
  string treeFilePath = ApplicationTools::getAFilePath(prefix + "trees.file", params, true, true, suffix, suffixIsOptional, "none", warn);

  BppOMultiTreeReaderFormat bppoReader(warn);
  unique_ptr<IMultiTree> iTrees(bppoReader.read(format));
  if (verbose)
  {
    ApplicationTools::displayResult("Input trees file " + suffix, treeFilePath);
    ApplicationTools::displayResult("Input trees format " + suffix, iTrees->getFormatName());
  }
  vector<Tree*> trees;
  iTrees->readTrees(treeFilePath, trees);

  if (verbose)
  {
    ApplicationTools::displayResult("Number of trees in file", trees.size());
  }
  return trees;
}


/*************************************************************/
/******* MODELS **********************************************/
/*************************************************************/

/******************************************************************************/

SubstitutionModel* PhylogeneticsApplicationTools::getSubstitutionModel(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const SiteContainer* data,
  std::map<std::string, std::string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  BppOSubstitutionModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, verbose, warn + 1);
  string modelDescription;
  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(alphabet);
  if (ca)
  {
    modelDescription = ApplicationTools::getStringParameter("model", params, "CodonRate(model=JC69)", suffix, suffixIsOptional, warn);
    if (!gCode)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModel(): a GeneticCode instance is required for instanciating a codon model.");
    bIO.setGeneticCode(gCode);
  }
  else if (AlphabetTools::isWordAlphabet(alphabet))
    modelDescription = ApplicationTools::getStringParameter("model", params, "Word(model=JC69)", suffix, suffixIsOptional, warn);
  else
    modelDescription = ApplicationTools::getStringParameter("model", params, "JC69", suffix, suffixIsOptional, warn);

  SubstitutionModel* model = bIO.readSubstitutionModel(alphabet, modelDescription, data, true);
  return model;
}

TransitionModel* PhylogeneticsApplicationTools::getTransitionModel(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const SiteContainer* data,
  std::map<std::string, std::string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  BppOTransitionModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, verbose, warn + 1);
  string modelDescription;
  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(alphabet);
  if (ca)
  {
    modelDescription = ApplicationTools::getStringParameter("model", params, "CodonRate(model=JC69)", suffix, suffixIsOptional, warn);
    if (!gCode)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModel(): a GeneticCode instance is required for instanciating a codon model.");
    bIO.setGeneticCode(gCode);
  }
  else if (AlphabetTools::isWordAlphabet(alphabet))
    modelDescription = ApplicationTools::getStringParameter("model", params, "Word(model=JC69)", suffix, suffixIsOptional, warn);
  else
    modelDescription = ApplicationTools::getStringParameter("model", params, "JC69", suffix, suffixIsOptional, warn);

  TransitionModel* model = bIO.readTransitionModel(alphabet, modelDescription, data, true);
  return model;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::setSubstitutionModelParametersInitialValuesWithAliases(
  TransitionModel& model,
  std::map<std::string, std::string>& unparsedParameterValues,
  size_t modelNumber,
  const SiteContainer* data,
  std::map<std::string, std::string>& sharedParams,
  bool verbose)
{
  string initFreqs = ApplicationTools::getStringParameter(model.getNamespace() + "initFreqs", unparsedParameterValues, "", "", true, 2);

  if (verbose)
    ApplicationTools::displayResult("Frequencies Initialization for model", (initFreqs == "") ? "None" : initFreqs);

  if (initFreqs != "")
  {
    if (initFreqs == "observed")
    {
      if (!data)
        throw Exception("Missing data for observed frequencies");
      unsigned int psi = ApplicationTools::getParameter<unsigned int>(model.getNamespace() + "initFreqs.observedPseudoCount", unparsedParameterValues, 0);
      model.setFreqFromData(*data, psi);
    }
    else if (initFreqs.substr(0, 6) == "values")
    {
      // Initialization using the "values" argument
      map<int, double> frequencies;

      string rf = initFreqs.substr(6);
      StringTokenizer strtok(rf.substr(1, rf.length() - 2), ",");
      int i = 0;
      while (strtok.hasMoreToken())
        frequencies[i++] = TextTools::toDouble(strtok.nextToken());
      model.setFreq(frequencies);
    }
    else
      throw Exception("Unknown initFreqs argument");
  }


  ParameterList pl = model.getIndependentParameters();
  for (size_t i = 0; i < pl.size(); ++i)
  {
    AutoParameter ap(pl[i]);
    ap.setMessageHandler(ApplicationTools::warning.get());
    pl.setParameter(i, ap);
  }

  for (size_t i = 0; i < pl.size(); ++i)
  {
    const string pName = pl[i].getName();
    size_t posp = model.getParameterNameWithoutNamespace(pName).rfind(".");
    string value;
    bool test1 = (initFreqs == "");
    bool test2 = (model.getParameterNameWithoutNamespace(pName).substr(posp + 1, 5) != "theta");
    bool test3 = (unparsedParameterValues.find(pName) != unparsedParameterValues.end());

    if (test1 || test2 || test3)
    {
      if (!test1 && !test2 && test3)
        ApplicationTools::displayWarning("Warning, initFreqs argument is set and a value is set for parameter " + pName);

      value = ApplicationTools::getStringParameter(pName, unparsedParameterValues, TextTools::toString(pl[i].getValue()));

      try
      {
        pl[i].setValue(TextTools::toDouble(value));
        if (verbose)
          ApplicationTools::displayResult("Parameter found", pName + +"_" + TextTools::toString(modelNumber) + "=" + TextTools::toString(pl[i].getValue()));
      }
      catch (Exception& e)
      {
        sharedParams[pl[i].getName() + "_" + TextTools::toString(modelNumber)] = value;
      }
    }
  }

  model.matchParametersValues(pl);
}


/******************************************************/
/**** FREQUENCIES SET *********************************/
/******************************************************/

/******************************************************************************/

FrequenciesSet* PhylogeneticsApplicationTools::getRootFrequenciesSet(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const SiteContainer* data,
  std::map<std::string, std::string>& params,
  std::map<std::string, std::string>& sharedparams,
  const std::vector<double>& rateFreqs,
  const std::string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  string freqDescription = ApplicationTools::getStringParameter("nonhomogeneous.root_freq", params, "Full(init=observed)", suffix, suffixIsOptional, warn);
  if (freqDescription == "None")
  {
    return 0;
  }
  else
  {
    map<string, string> unparams;

    FrequenciesSet* freq = getFrequenciesSet(alphabet, gCode, freqDescription, data, unparams, rateFreqs, verbose, warn + 1);
    freq->setNamespace("root." + freq->getNamespace());

    for (const auto& it:unparams)
      sharedparams["root." + it.first] = it.second;

    if (verbose)
      ApplicationTools::displayResult("Root frequencies ", freq->getName());
    return freq;
  }
}

/******************************************************************************/

FrequenciesSet* PhylogeneticsApplicationTools::getFrequenciesSet(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const std::string& freqDescription,
  const SiteContainer* data,
  std::map<std::string, std::string>& sharedparams,
  const std::vector<double>& rateFreqs,
  bool verbose,
  int warn)
{
  map<string, string> unparsedParameterValues;
  BppOFrequenciesSetFormat bIO(BppOFrequenciesSetFormat::ALL, verbose, warn);
  if (AlphabetTools::isCodonAlphabet(alphabet))
  {
    if (!gCode)
      throw Exception("PhylogeneticsApplicationTools::getFrequenciesSet(): a GeneticCode instance is required for instanciating a codon frequencies set.");
    bIO.setGeneticCode(gCode);
  }
  unique_ptr<FrequenciesSet> pFS(bIO.read(alphabet, freqDescription, data, true));

  std::map<std::string, std::string> unparsedparam = bIO.getUnparsedArguments();

  sharedparams.insert(unparsedparam.begin(), unparsedparam.end());

  // /////// To be changed for input normalization
  if (rateFreqs.size() > 0)
  {
    pFS.reset(new MarkovModulatedFrequenciesSet(pFS.release(), rateFreqs));
  }

  return pFS.release();
}

/******************************************************/
/**** SUBSTITUTION MODEL SET **************************/
/******************************************************/

/******************************************************************************/

SubstitutionModelSet* PhylogeneticsApplicationTools::getSubstitutionModelSet(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const SiteContainer* data,
  std::map<std::string, std::string>& params,
  const std::string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  if (!ApplicationTools::parameterExists("nonhomogeneous.number_of_models", params))
    throw Exception("A value is needed for this parameter: nonhomogeneous.number_of_models .");
  size_t nbModels = ApplicationTools::getParameter<size_t>("nonhomogeneous.number_of_models", params, 1, suffix, suffixIsOptional, warn);
  if (nbModels == 0)
    throw Exception("The number of models can't be 0 !");

  // bool nomix = true;
  // for (size_t i = 0; nomix &(i < nbModels); i++)
  // {
  //   string prefix = "model" + TextTools::toString(i + 1);
  //   string modelDesc;
  //   modelDesc = ApplicationTools::getStringParameter(prefix, params, "", suffix, suffixIsOptional, warn);

  //   if (modelDesc.find("Mixed") != string::npos)
  //     nomix = false;
  // }

  SubstitutionModelSet* modelSet, * modelSet1 = 0;
  modelSet1 = new SubstitutionModelSet(alphabet);
  setSubstitutionModelSet(*modelSet1, alphabet, gCode, data, params, suffix, suffixIsOptional, verbose, warn);

  if (modelSet1->hasMixedTransitionModel())
  {
    modelSet = new MixedSubstitutionModelSet(*modelSet1);
    completeMixedSubstitutionModelSet(*dynamic_cast<MixedSubstitutionModelSet*>(modelSet), alphabet, data, params, suffix, suffixIsOptional, verbose);
  }
  else
    modelSet = modelSet1;

  return modelSet;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::setSubstitutionModelSet(
  SubstitutionModelSet& modelSet,
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const SiteContainer* data,
  map<string, string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  modelSet.clear();
  if (!ApplicationTools::parameterExists("nonhomogeneous.number_of_models", params))
    throw Exception("You must specify this parameter: 'nonhomogeneous.number_of_models'.");
  size_t nbModels = ApplicationTools::getParameter<size_t>("nonhomogeneous.number_of_models", params, 1, suffix, suffixIsOptional, warn);
  if (nbModels == 0)
    throw Exception("The number of models can't be 0 !");

  if (verbose)
    ApplicationTools::displayResult("Number of distinct models", TextTools::toString(nbModels));

  BppOTransitionModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, false, warn);

  // ///////////////////////////////////////////
  // Build a new model set object:

  vector<double> rateFreqs;
  string tmpDesc;
  if (AlphabetTools::isCodonAlphabet(alphabet))
  {
    if (!gCode)
      throw Exception("PhylogeneticsApplicationTools::setSubstitutionModelSet(): a GeneticCode instance is required for instanciating a codon model.");
    bIO.setGeneticCode(gCode);
    tmpDesc = ApplicationTools::getStringParameter("model1", params, "CodonRate(model=JC69)", suffix, suffixIsOptional, warn);
  }
  else if (AlphabetTools::isWordAlphabet(alphabet))
    tmpDesc = ApplicationTools::getStringParameter("model1", params, "Word(model=JC69)", suffix, suffixIsOptional, warn);
  else
    tmpDesc = ApplicationTools::getStringParameter("model1", params, "JC69", suffix, suffixIsOptional, warn);

  unique_ptr<TransitionModel> tmp(bIO.readTransitionModel(alphabet, tmpDesc, data, false));


  if (tmp->getNumberOfStates() != alphabet->getSize())
  {
    // Markov-Modulated Markov Model...
    size_t n = static_cast<size_t>(tmp->getNumberOfStates() / alphabet->getSize());
    rateFreqs = vector<double>(n, 1. / static_cast<double>(n)); // Equal rates assumed for now, may be changed later (actually, in the most general case,
  }

  // ////////////////////////////////////
  // Deal with root frequencies

  map<string, string> unparsedParameters;

  bool stationarity = ApplicationTools::getBooleanParameter("nonhomogeneous.stationarity", params, false, "", true, warn);
  FrequenciesSet* rootFrequencies = 0;
  if (!stationarity)
  {
    rootFrequencies = getRootFrequenciesSet(alphabet, gCode, data, params, unparsedParameters, rateFreqs, suffix, suffixIsOptional, verbose);
    stationarity = !rootFrequencies;
    string freqDescription = ApplicationTools::getStringParameter("nonhomogeneous.root_freq", params, "", suffix, suffixIsOptional, warn);
    if (freqDescription.substr(0, 10) == "MVAprotein")
    {
      if (dynamic_cast<Coala*>(tmp.get()))
        dynamic_cast<MvaFrequenciesSet*>(rootFrequencies)->initSet(dynamic_cast<CoalaCore*>(tmp.get()));
      else
        throw Exception("The MVAprotein frequencies set at the root can only be used if a Coala model is used on branches.");
    }
  }

  ApplicationTools::displayBooleanResult("Stationarity assumed", stationarity);

  if (!stationarity)
    modelSet.setRootFrequencies(rootFrequencies);


  // //////////////////////////////////////
  // Now parse all models:

  bIO.setVerbose(true);

  for (size_t i = 0; i < nbModels; i++)
  {
    string prefix = "model" + TextTools::toString(i + 1);
    string modelDesc;
    if (AlphabetTools::isCodonAlphabet(alphabet))
      modelDesc = ApplicationTools::getStringParameter(prefix, params, "CodonRate(model=JC69)", suffix, suffixIsOptional, warn);
    else if (AlphabetTools::isWordAlphabet(alphabet))
      modelDesc = ApplicationTools::getStringParameter(prefix, params, "Word(model=JC69)", suffix, suffixIsOptional, warn);
    else
      modelDesc = ApplicationTools::getStringParameter(prefix, params, "JC69", suffix, suffixIsOptional, warn);

    unique_ptr<TransitionModel> model(bIO.readTransitionModel(alphabet, modelDesc, data, false));

    map<string, string> unparsedModelParameters = bIO.getUnparsedArguments();
    map<string, string> sharedParameters;

    setSubstitutionModelParametersInitialValuesWithAliases(
      *model,
      unparsedModelParameters, i + 1, data,
      sharedParameters,
      verbose);

    unparsedParameters.insert(sharedParameters.begin(), sharedParameters.end());

    vector<int> nodesId = ApplicationTools::getVectorParameter<int>(prefix + ".nodes_id", params, ',', ':', "", suffix, suffixIsOptional, warn);
    
    if (verbose)
      ApplicationTools::displayResult("Model" + TextTools::toString(i + 1) + " is associated to", TextTools::toString(nodesId.size()) + " node(s).");

    modelSet.addModel(model.release(), nodesId);
  }

  // Finally check parameter aliasing:
  string aliasDesc = ApplicationTools::getStringParameter("nonhomogeneous.alias", params, "", suffix, suffixIsOptional, warn);
  StringTokenizer st(aliasDesc, ",");
  while (st.hasMoreToken())
  {
    string alias = st.nextToken();
    string::size_type index = alias.find("->");
    if (index == string::npos)
      throw Exception("PhylogeneticsApplicationTools::setSubstitutionModelSet. Bad alias syntax, should contain `->' symbol: " + alias);
    unparsedParameters[alias.substr(0, index)] = alias.substr(index + 2);
  }

  // alias unparsedParameters

  modelSet.aliasParameters(unparsedParameters, verbose);
}

/******************************************************************************/
void PhylogeneticsApplicationTools::completeMixedSubstitutionModelSet(
  MixedSubstitutionModelSet& mixedModelSet,
  const Alphabet* alphabet,
  const SiteContainer* data,
  map<string, string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  // /////////////////////////////////////////
  // Looks for the allowed paths

  size_t numd;
  if (!ApplicationTools::parameterExists("site.number_of_paths", params))
    numd = 0;
  else
    numd = ApplicationTools::getParameter<size_t>("site.number_of_paths", params, 1, suffix, suffixIsOptional, warn);

  if (verbose)
    ApplicationTools::displayResult("Number of distinct paths", TextTools::toString(numd));

  vector<string> vdesc;
  size_t numi=0;
  while (numi<numd)
  {
    string desc = ApplicationTools::getStringParameter("site.path" + TextTools::toString(numi+1), params, "",  suffix, suffixIsOptional, warn);
    if (desc.size() == 0)
      break;
    else
      vdesc.push_back(desc);
    numi++;
  }

  if (vdesc.size() == 0)
  {
    mixedModelSet.complete();
    mixedModelSet.computeHyperNodesProbabilities();
    return;
  }

  for (auto& desc:vdesc)
  {
    mixedModelSet.addEmptyHyperNode();
    StringTokenizer st(desc, "&");
    while (st.hasMoreToken())
    {
      string submodel = st.nextToken();
      Vint submodelNb;
      string::size_type indexo = submodel.find("[");
      string::size_type indexf = submodel.find("]");
      if ((indexo == string::npos) | (indexf == string::npos))
        throw Exception("PhylogeneticsApplicationTools::completeMixedSubstitutionModelSet. Bad path syntax, should contain `[]' symbols: " + submodel);
      size_t num = TextTools::to<size_t>(submodel.substr(5, indexo - 5));
      const MixedTransitionModel* pSM = dynamic_cast<const MixedTransitionModel*>(mixedModelSet.getModel(num - 1));
      if (!pSM)
        throw BadIntegerException("PhylogeneticsApplicationTools::completeMixedSubstitutionModelSet: Wrong model for number", static_cast<int>(num - 1));

      string lp2 = submodel.substr(indexo + 1, indexf - indexo - 1);      
      StringTokenizer stp2(lp2, ",");
      while (stp2.hasMoreToken())
      {
        string p2=stp2.nextToken();
        
        try  {
          int n2=TextTools::toInt(p2);
          if (n2<=0 || n2>(int)(pSM->getNumberOfModels()))
            throw BadIntegerException("PhylogeneticsApplicationTools::completeMixedSubstitutionModelSet: Wrong model for number", static_cast<int>(n2));
          submodelNb.push_back(n2-1);
        }
        catch (Exception& e)
        {
          Vint submodnb = pSM->getSubmodelNumbers(p2);
          if (submodelNb.size()==0)
            submodelNb=submodnb;
          else
            submodelNb=VectorTools::vectorIntersection(submodelNb,submodnb);
        }
      }
      
      mixedModelSet.addToHyperNode(num - 1, submodelNb);
    }
    
    if (!mixedModelSet.getHyperNode(mixedModelSet.getNumberOfHyperNodes() - 1).isComplete())
      throw Exception("A path should own at least a submodel of each mixed model: " + desc);

    if (verbose)
      ApplicationTools::displayResult("Site Path", desc);
  }

  // / Checks if the paths are separate
  if (!mixedModelSet.hasExclusivePaths())
    throw Exception("All paths must be disjoint.");

  // / Puts all the remaining models in a new path
  string st;
  st = (mixedModelSet.complete()) ? "Yes" : "No";

  if (verbose)
    ApplicationTools::displayResult("Site Path Completion", st);

  mixedModelSet.computeHyperNodesProbabilities();

  if (!mixedModelSet.getHyperNode(mixedModelSet.getNumberOfHyperNodes() - 1).isComplete())
    throw Exception("The remaining submodels can not create a complete path.");
}


/******************************************************/
/*** DISTRIBUTIONS ********************************/
/******************************************************/


/******************************************************************************/

MultipleDiscreteDistribution* PhylogeneticsApplicationTools::getMultipleDistributionDefaultInstance(
  const std::string& distDescription,
  std::map<std::string, std::string>& unparsedParameterValues,
  bool verbose)
{
  string distName;
  MultipleDiscreteDistribution* pMDD  = 0;
  map<string, string> args;
  KeyvalTools::parseProcedure(distDescription, distName, args);

  if (distName == "Dirichlet")
  {
    if (args.find("classes") == args.end())
      throw Exception("Missing argument 'classes' (vector of number of classes) in " + distName
                      + " distribution");
    if (args.find("alphas") == args.end())
      throw Exception("Missing argument 'alphas' (vector of Dirichlet shape parameters) in Dirichlet distribution");
    vector<double> alphas;
    vector<size_t> classes;

    string rf = args["alphas"];
    StringTokenizer strtok(rf.substr(1, rf.length() - 2), ",");
    while (strtok.hasMoreToken())
      alphas.push_back(TextTools::toDouble(strtok.nextToken()));

    rf = args["classes"];
    StringTokenizer strtok2(rf.substr(1, rf.length() - 2), ",");
    while (strtok2.hasMoreToken())
      classes.push_back(TextTools::to<size_t>(strtok2.nextToken()));

    pMDD = new DirichletDiscreteDistribution(classes, alphas);
    vector<string> v = pMDD->getParameters().getParameterNames();

    for (size_t i = 0; i < v.size(); i++)
    {
      unparsedParameterValues[v[i]] = TextTools::toString(pMDD->getParameterValue(pMDD->getParameterNameWithoutNamespace(v[i])));
    }
  }
  else
    throw Exception("Unknown multiple distribution name: " + distName);

  return pMDD;
}

/******************************************************************************/

DiscreteDistribution* PhylogeneticsApplicationTools::getRateDistribution(
  map<string, string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose)
{
  string distDescription = ApplicationTools::getStringParameter("rate_distribution", params, "Constant()", suffix, suffixIsOptional);

  string distName;
  map<string, string> args;
  KeyvalTools::parseProcedure(distDescription, distName, args);

  BppORateDistributionFormat bIO(true);
  unique_ptr<DiscreteDistribution> rDist(bIO.read(distDescription, true));

  if (verbose)
  {
    ApplicationTools::displayResult("Rate distribution", distName);
    ApplicationTools::displayResult("Number of classes", TextTools::toString(rDist->getNumberOfCategories()));
  }

  return rDist.release();
}


/*************************************************************/
/*****  OPTIMIZATORS *****************************************/
/*************************************************************/

/******************************************************************************/

TreeLikelihood* PhylogeneticsApplicationTools::optimizeParameters(
  TreeLikelihood* tl,
  const ParameterList& parameters,
  std::map<std::string, std::string>& params,
  const std::string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  string optimization = ApplicationTools::getStringParameter("optimization", params, "FullD(derivatives=Newton)", suffix, suffixIsOptional, warn);
  if (optimization == "None")
    return tl;
  string optName;
  map<string, string> optArgs;
  KeyvalTools::parseProcedure(optimization, optName, optArgs);

  unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2, suffix, suffixIsOptional, warn + 1);

  string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
  OutputStream* messageHandler =
    (mhPath == "none") ? 0 :
    (mhPath == "std") ? ApplicationTools::message.get() :
    new StlOutputStream(new ofstream(mhPath.c_str(), ios::out));
  if (verbose)
    ApplicationTools::displayResult("Message handler", mhPath);

  string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
  OutputStream* profiler =
    (prPath == "none") ? 0 :
    (prPath == "std") ? ApplicationTools::message.get() :
    new StlOutputStream(new ofstream(prPath.c_str(), ios::out));

  if (profiler)
    profiler->setPrecision(20);
  if (verbose)
    ApplicationTools::displayResult("Profiler", prPath);

  bool scaleFirst = ApplicationTools::getBooleanParameter("optimization.scale_first", params, false, suffix, suffixIsOptional, warn + 1);
  if (scaleFirst)
  {
    // We scale the tree before optimizing each branch length separately:
    if (verbose)
      ApplicationTools::displayMessage("Scaling the tree before optimizing each branch length separately.");
    double tolerance = ApplicationTools::getDoubleParameter("optimization.scale_first.tolerance", params, .0001, suffix, suffixIsOptional, warn + 1);
    if (verbose)
      ApplicationTools::displayResult("Scaling tolerance", TextTools::toString(tolerance));
    unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.scale_first.max_number_f_eval", params, 1000000, suffix, suffixIsOptional, warn + 1);
    if (verbose)
      ApplicationTools::displayResult("Scaling max # f eval", TextTools::toString(nbEvalMax));
    OptimizationTools::optimizeTreeScale(
      tl,
      tolerance,
      nbEvalMax,
      messageHandler,
      profiler);
    if (verbose)
      ApplicationTools::displayResult("New tree likelihood", -tl->getValue());
  }

  // Should I ignore some parameters?
  ParameterList parametersToEstimate = parameters;
  vector<string> parNames = parametersToEstimate.getParameterNames();

  if (params.find("optimization.ignore_parameter") != params.end())
    throw Exception("optimization.ignore_parameter is deprecated, use optimization.ignore_parameters instead!");
  string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameters", params, "", suffix, suffixIsOptional, warn + 1);
  StringTokenizer st(paramListDesc, ",");
  while (st.hasMoreToken())
  {
    try
    {
      string param = st.nextToken();
      if (param == "BrLen")
      {
        vector<string> vs = tl->getBranchLengthsParameters().getParameterNames();
        parametersToEstimate.deleteParameters(vs);
        if (verbose)
          ApplicationTools::displayResult("Parameter ignored", string("Branch lengths"));
      }
      else if (param == "Ancient")
      {
        NonHomogeneousTreeLikelihood* nhtl = dynamic_cast<NonHomogeneousTreeLikelihood*>(tl);
        if (!nhtl)
          ApplicationTools::displayWarning("The 'Ancient' parameters do not exist in homogeneous models, and will be ignored.");
        else
        {
          vector<string> vs = nhtl->getRootFrequenciesParameters().getParameterNames();
          parametersToEstimate.deleteParameters(vs);
        }
        if (verbose)
          ApplicationTools::displayResult("Parameter ignored", string("Root frequencies"));
      }
      else if (param == "Model")
      {
        vector<string> vs;
        vector<string> vs1 = tl->getSubstitutionModelParameters().getParameterNames();
        NonHomogeneousTreeLikelihood* nhtl = dynamic_cast<NonHomogeneousTreeLikelihood*>(tl);
        if (nhtl != NULL)
        {
          vector<string> vs2 = nhtl->getRootFrequenciesParameters().getParameterNames();
          VectorTools::diff(vs1, vs2, vs);
        }
        else
          vs = vs1;

        parametersToEstimate.deleteParameters(vs);
        if (verbose)
          ApplicationTools::displayResult("Parameter ignored", string("Model"));
      }
      else if (param.find("*") != string::npos)
      {
        vector<string> vs = ApplicationTools::matchingParameters(param, parNames);

        for (vector<string>::iterator it = vs.begin(); it != vs.end(); it++)
        {
          parametersToEstimate.deleteParameter(*it);
          if (verbose)
            ApplicationTools::displayResult("Parameter ignored", *it);
        }
      }
      else
      {
        parametersToEstimate.deleteParameter(param);
        if (verbose)
          ApplicationTools::displayResult("Parameter ignored", param);
      }
    }
    catch (ParameterNotFoundException& pnfe)
    {
      ApplicationTools::displayWarning("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
    }
  }

  // Should I constrain some parameters?
  vector<string> parToEstNames = parametersToEstimate.getParameterNames();

  if (params.find("optimization.constrain_parameter") != params.end())
    throw Exception("optimization.constrain_parameter is deprecated, use optimization.constrain_parameters instead!");
  paramListDesc = ApplicationTools::getStringParameter("optimization.constrain_parameters", params, "", suffix, suffixIsOptional, warn + 1);

  string constraint = "";
  string pc, param = "";

  StringTokenizer st2(paramListDesc, ",");
  while (st2.hasMoreToken())
  {
    try
    {
      pc = st2.nextToken();
      string::size_type index = pc.find("=");
      if (index == string::npos)
        throw Exception("PhylogeneticsApplicationTools::optimizeParamaters. Bad constrain syntax, should contain `=' symbol: " + pc);
      param = pc.substr(0, index);
      constraint = pc.substr(index + 1);
      IntervalConstraint ic(constraint);

      vector<string> parNames2;

      if (param == "BrLen")
        parNames2  = tl->getBranchLengthsParameters().getParameterNames();
      else if (param == "Ancient")
      {
        NonHomogeneousTreeLikelihood* nhtl = dynamic_cast<NonHomogeneousTreeLikelihood*>(tl);
        if (!nhtl)
          ApplicationTools::displayWarning("The 'Ancient' parameters do not exist in homogeneous models, and will be ignored.");
        else
        {
          parNames2 = nhtl->getRootFrequenciesParameters().getParameterNames();
          ApplicationTools::displayResult("Parameter ignored", string("Root frequencies"));
        }
      }
      else if (param == "Model")
      {
        vector<string> vs1 = tl->getSubstitutionModelParameters().getParameterNames();
        NonHomogeneousTreeLikelihood* nhtl = dynamic_cast<NonHomogeneousTreeLikelihood*>(tl);
        if (nhtl != NULL)
        {
          vector<string> vs2 = nhtl->getRootFrequenciesParameters().getParameterNames();
          VectorTools::diff(vs1, vs2, parNames2);
        }
        else
          parNames2 = vs1;
      }
      else if (param.find("*") != string::npos)
        parNames2 = ApplicationTools::matchingParameters(param, parToEstNames);
      else
        parNames2.push_back(param);


      for (size_t i = 0; i < parNames2.size(); i++)
      {
        Parameter& par = parametersToEstimate.getParameter(parNames2[i]);
        if (par.hasConstraint())
        {
          par.setConstraint(std::shared_ptr<Constraint>(ic & (*par.getConstraint())));
          if (par.getConstraint()->isEmpty())
            throw Exception("Empty interval for parameter " + parNames[i] + par.getConstraint()->getDescription());
        }
        else
          par.setConstraint(std::shared_ptr<Constraint>(ic.clone()));

        if (verbose)
          ApplicationTools::displayResult("Parameter constrained " + par.getName(), par.getConstraint()->getDescription());
      }
    }
    catch (ParameterNotFoundException& pnfe)
    {
      ApplicationTools::displayWarning("Parameter '" + pnfe.getParameter() + "' not found, and so can't be constrained!");
    }
    catch (ConstraintException& pnfe)
    {
      throw Exception("Parameter '" + param + "' does not fit the constraint " + constraint);
    }
  }


  // /////
  // / optimization options

  unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000, suffix, suffixIsOptional, warn + 1);
  if (verbose)
    ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));

  double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", params, .000001, suffix, suffixIsOptional, warn + 1);
  if (verbose)
    ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));

  // Backing up or restoring?
  unique_ptr<BackupListener> backupListener;
  string backupFile = ApplicationTools::getAFilePath("optimization.backup.file", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
  if (backupFile != "none")
  {
    ApplicationTools::displayResult("Parameters will be backup to", backupFile);
    backupListener.reset(new BackupListener(backupFile));
    if (FileTools::fileExists(backupFile))
    {
      ApplicationTools::displayMessage("A backup file was found! Try to restore parameters from previous run...");
      ifstream bck(backupFile.c_str(), ios::in);
      vector<string> lines = FileTools::putStreamIntoVectorOfStrings(bck);
      double fval = TextTools::toDouble(lines[0].substr(5));
      ParameterList pl = tl->getParameters();
      for (size_t l = 1; l < lines.size(); ++l)
      {
        if (!TextTools::isEmpty(lines[l]))
        {
          StringTokenizer stp(lines[l], "=");
          if (stp.numberOfRemainingTokens() != 2)
          {
            cerr << "Corrupted backup file!!!" << endl;
            cerr << "at line " << l << ": " << lines[l] << endl;
          }
          string pname  = stp.nextToken();
          string pvalue = stp.nextToken();
          try {
            size_t p = pl.whichParameterHasName(pname);
            pl.setParameter(p, AutoParameter(pl[p]));
            pl[p].setValue(TextTools::toDouble(pvalue));
          }
          catch(Exception& e)
          {
          }
        }
      }
      bck.close();
      tl->setParameters(pl);
      if (abs(tl->getValue() - fval) > 0.000001)
        ApplicationTools::displayWarning("Warning, incorrect likelihood value after restoring from backup file.");
      ApplicationTools::displayResult("Restoring log-likelihood", -fval);
    }
  }

  // There it goes...
  bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", params, false, suffix, suffixIsOptional, warn + 1);
  if (verbose)
    ApplicationTools::displayResult("Optimize topology", optimizeTopo ? "yes" : "no");
  string nniMethod = ApplicationTools::getStringParameter("optimization.topology.algorithm_nni.method", params, "phyml", suffix, suffixIsOptional, warn + 1);
  string nniAlgo;
  if (nniMethod == "fast")
  {
    nniAlgo = NNITopologySearch::FAST;
  }
  else if (nniMethod == "better")
  {
    nniAlgo = NNITopologySearch::BETTER;
  }
  else if (nniMethod == "phyml")
  {
    nniAlgo = NNITopologySearch::PHYML;
  }
  else
    throw Exception("Unknown NNI algorithm: '" + nniMethod + "'.");


  string order = ApplicationTools::getStringParameter("derivatives", optArgs, "Newton", "", true, warn + 1);
  string optMethodDeriv;
  if (order == "Gradient")
  {
    optMethodDeriv = OptimizationTools::OPTIMIZATION_GRADIENT;
  }
  else if (order == "Newton")
  {
    optMethodDeriv = OptimizationTools::OPTIMIZATION_NEWTON;
  }
  else if (order == "BFGS")
  {
    optMethodDeriv = OptimizationTools::OPTIMIZATION_BFGS;
  }
  else
    throw Exception("Unknown derivatives algorithm: '" + order + "'.");
  if (verbose)
    ApplicationTools::displayResult("Optimization method", optName);
  if (verbose)
    ApplicationTools::displayResult("Algorithm used for derivable parameters", order);

  // See if we should reparametrize:
  bool reparam = ApplicationTools::getBooleanParameter("optimization.reparametrization", params, false, suffix, suffixIsOptional, warn + 1);
  if (verbose)
    ApplicationTools::displayResult("Reparametrization", (reparam ? "yes" : "no"));

  // See if we should use a molecular clock constraint:
  string clock = ApplicationTools::getStringParameter("optimization.clock", params, "None", suffix, suffixIsOptional, warn + 1);
  if (clock != "None" && clock != "Global")
    throw Exception("Molecular clock option not recognized, should be one of 'Global' or 'None'.");
  bool useClock = (clock == "Global");
  if (useClock && optimizeTopo)
    throw Exception("PhylogeneticsApplicationTools::optimizeParameters. Cannot optimize topology with a molecular clock.");
  if (verbose)
    ApplicationTools::displayResult("Molecular clock", clock);

  unsigned int n = 0;
  if ((optName == "D-Brent") || (optName == "D-BFGS"))
  {
    // Uses Newton-Brent method or Newton-BFGS method
    string optMethodModel;
    if (optName == "D-Brent")
      optMethodModel = OptimizationTools::OPTIMIZATION_BRENT;
    else
      optMethodModel = OptimizationTools::OPTIMIZATION_BFGS;

    unsigned int nstep = ApplicationTools::getParameter<unsigned int>("nstep", optArgs, 1, "", true, warn + 1);

    if (optimizeTopo)
    {
      bool optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, warn + 1);
      unsigned int topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, suffix, suffixIsOptional, warn + 1);
      double tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional, warn + 1);
      double tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, suffixIsOptional, warn + 1);
      tl = OptimizationTools::optimizeTreeNNI(
        dynamic_cast<NNIHomogeneousTreeLikelihood*>(tl), parametersToEstimate,
        optNumFirst, tolBefore, tolDuring, nbEvalMax, topoNbStep, messageHandler, profiler,
        reparam, optVerbose, optMethodDeriv, nstep, nniAlgo);
    }

    if (verbose && nstep > 1)
      ApplicationTools::displayResult("# of precision steps", TextTools::toString(nstep));
    parametersToEstimate.matchParametersValues(tl->getParameters());
    n = OptimizationTools::optimizeNumericalParameters(
      dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(tl), parametersToEstimate,
      backupListener.get(), nstep, tolerance, nbEvalMax, messageHandler, profiler, reparam, optVerbose, optMethodDeriv, optMethodModel);
  }
  else if (optName == "FullD")
  {
    // Uses Newton-raphson algorithm with numerical derivatives when required.

    if (optimizeTopo)
    {
      bool optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, warn + 1);
      unsigned int topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, suffix, suffixIsOptional, warn + 1);
      double tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional, warn + 1);
      double tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, suffixIsOptional, warn + 1);
      tl = OptimizationTools::optimizeTreeNNI2(
        dynamic_cast<NNIHomogeneousTreeLikelihood*>(tl), parametersToEstimate,
        optNumFirst, tolBefore, tolDuring, nbEvalMax, topoNbStep, messageHandler, profiler,
        reparam, optVerbose, optMethodDeriv, nniAlgo);
    }

    parametersToEstimate.matchParametersValues(tl->getParameters());
    n = OptimizationTools::optimizeNumericalParameters2(
      dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(tl), parametersToEstimate,
      backupListener.get(), tolerance, nbEvalMax, messageHandler, profiler, reparam, useClock, optVerbose, optMethodDeriv);
  }
  else
    throw Exception("Unknown optimization method: " + optName);

  string finalMethod = ApplicationTools::getStringParameter("optimization.final", params, "none", suffix, suffixIsOptional, warn + 1);
  Optimizer* finalOptimizer  = 0;
  if (finalMethod == "none")
  {}
  else if (finalMethod == "simplex")
  {
    finalOptimizer = new DownhillSimplexMethod(tl);
  }
  else if (finalMethod == "powell")
  {
    finalOptimizer = new PowellMultiDimensions(tl);
  }
  else
    throw Exception("Unknown final optimization method: " + finalMethod);

  if (finalOptimizer)
  {
    parametersToEstimate.matchParametersValues(tl->getParameters());
    if (verbose)
      ApplicationTools::displayResult("Final optimization step", finalMethod);
    finalOptimizer->setProfiler(profiler);
    finalOptimizer->setMessageHandler(messageHandler);
    finalOptimizer->setMaximumNumberOfEvaluations(nbEvalMax);
    finalOptimizer->getStopCondition()->setTolerance(tolerance);
    finalOptimizer->setVerbose(verbose);
    finalOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    finalOptimizer->init(parametersToEstimate);
    finalOptimizer->optimize();
    n += finalOptimizer->getNumberOfEvaluations();
    delete finalOptimizer;
  }

  if (prPath != "none" && prPath != "std")
    delete profiler;
  if (mhPath != "none" && mhPath != "std")
    delete messageHandler;

  if (verbose)
    ApplicationTools::displayResult("Performed", TextTools::toString(n) + " function evaluations.");
  if (backupFile != "none")
  {
    string bf=backupFile+".def";
    rename(backupFile.c_str(),bf.c_str());
  }
  return tl;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::optimizeParameters(
  DiscreteRatesAcrossSitesClockTreeLikelihood* tl,
  const ParameterList& parameters,
  map<string, string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  string optimization = ApplicationTools::getStringParameter("optimization", params, "FullD(derivatives=Newton)", suffix, suffixIsOptional, warn);
  if (optimization == "None")
    return;
  string optName;
  map<string, string> optArgs;
  KeyvalTools::parseProcedure(optimization, optName, optArgs);

  unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2, suffix, suffixIsOptional, warn + 1);

  string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
  OutputStream* messageHandler =
    (mhPath == "none") ? 0 :
    (mhPath == "std") ? ApplicationTools::message.get() :
    new StlOutputStream(new ofstream(mhPath.c_str(), ios::out));
  if (verbose)
    ApplicationTools::displayResult("Message handler", mhPath);

  string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
  OutputStream* profiler =
    (prPath == "none") ? 0 :
    (prPath == "std") ? ApplicationTools::message.get() :
    new StlOutputStream(new ofstream(prPath.c_str(), ios::out));
  if (profiler)
    profiler->setPrecision(20);
  if (verbose)
    ApplicationTools::displayResult("Profiler", prPath);

  ParameterList parametersToEstimate = parameters;

  // Should I ignore some parameters?
  if (params.find("optimization.ignore_parameter") != params.end())
    throw Exception("optimization.ignore_parameter is deprecated, use optimization.ignore_parameters instead!");
  string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameters", params, "", suffix, suffixIsOptional, warn + 1);
  StringTokenizer st(paramListDesc, ",");
  while (st.hasMoreToken())
  {
    try
    {
      string param = st.nextToken();
      if (param == "BrLen")
      {
        vector<string> vs = tl->getBranchLengthsParameters().getParameterNames();
        parametersToEstimate.deleteParameters(vs);
        if (verbose)
          ApplicationTools::displayResult("Parameter ignored", string("Branch lengths"));
      }
      else if (param == "Ancient")
      {
        NonHomogeneousTreeLikelihood* nhtl = dynamic_cast<NonHomogeneousTreeLikelihood*>(tl);
        if (!nhtl)
          ApplicationTools::displayWarning("The 'Ancient' parameters do not exist in homogeneous models, and will be ignored.");
        else
        {
          vector<string> vs = nhtl->getRootFrequenciesParameters().getParameterNames();
          parametersToEstimate.deleteParameters(vs);
        }
        if (verbose)
          ApplicationTools::displayResult("Parameter ignored", string("Root frequencies"));
      }
      else
      {
        parametersToEstimate.deleteParameter(param);
        if (verbose)
          ApplicationTools::displayResult("Parameter ignored", param);
      }
    }
    catch (ParameterNotFoundException& pnfe)
    {
      ApplicationTools::displayError("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
    }
  }

  unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000, suffix, suffixIsOptional, warn + 1);
  if (verbose)
    ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));

  double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", params, .000001, suffix, suffixIsOptional, warn + 1);
  if (verbose)
    ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));

  string order  = ApplicationTools::getStringParameter("derivatives", optArgs, "Gradient", "", true, warn + 1);
  string optMethod, derMethod;
  if (order == "Gradient")
  {
    optMethod = OptimizationTools::OPTIMIZATION_GRADIENT;
  }
  else if (order == "Newton")
  {
    optMethod = OptimizationTools::OPTIMIZATION_NEWTON;
  }
  else
    throw Exception("Option '" + order + "' is not known for 'optimization.method.derivatives'.");
  if (verbose)
    ApplicationTools::displayResult("Optimization method", optName);
  if (verbose)
    ApplicationTools::displayResult("Algorithm used for derivable parameters", order);

  // Backing up or restoring?
  unique_ptr<BackupListener> backupListener;
  string backupFile = ApplicationTools::getAFilePath("optimization.backup.file", params, false, false, suffix, suffixIsOptional, "none", warn + 1);
  if (backupFile != "none")
  {
    ApplicationTools::displayResult("Parameters will be backup to", backupFile);
    backupListener.reset(new BackupListener(backupFile));
    if (FileTools::fileExists(backupFile))
    {
      ApplicationTools::displayMessage("A backup file was found! Try to restore parameters from previous run...");
      ifstream bck(backupFile.c_str(), ios::in);
      vector<string> lines = FileTools::putStreamIntoVectorOfStrings(bck);
      double fval = TextTools::toDouble(lines[0].substr(5));
      ParameterList pl = tl->getParameters();
      for (size_t l = 1; l < lines.size(); ++l)
      {
        if (!TextTools::isEmpty(lines[l]))
        {
          StringTokenizer stp(lines[l], "=");
          if (stp.numberOfRemainingTokens() != 2)
          {
            cerr << "Corrupted backup file!!!" << endl;
            cerr << "at line " << l << ": " << lines[l] << endl;
          }
          string pname  = stp.nextToken();
          string pvalue = stp.nextToken();
          size_t p = pl.whichParameterHasName(pname);
          pl.setParameter(p, AutoParameter(pl[p]));
          pl[p].setValue(TextTools::toDouble(pvalue));
        }
      }
      bck.close();
      tl->setParameters(pl);
      if (abs(tl->getValue() - fval) > 0.000001)
        throw Exception("Incorrect likelihood value after restoring, from backup file. Remove backup file and start from scratch :s");
      ApplicationTools::displayResult("Restoring log-likelihood", -fval);
    }
  }

  size_t n = 0;
  if (optName == "D-Brent")
  {
    // Uses Newton-Brent method:
    unsigned int nstep = ApplicationTools::getParameter<unsigned int>("nstep", optArgs, 1, "", true, warn + 1);
    if (verbose && nstep > 1)
      ApplicationTools::displayResult("# of precision steps", TextTools::toString(nstep));
    n = OptimizationTools::optimizeNumericalParametersWithGlobalClock(
      tl,
      parametersToEstimate,
      backupListener.get(),
      nstep,
      tolerance,
      nbEvalMax,
      messageHandler,
      profiler,
      optVerbose,
      optMethod);
  }
  else if (optName == "FullD")
  {
    // Uses Newton-raphson alogrithm with numerical derivatives when required.
    n = OptimizationTools::optimizeNumericalParametersWithGlobalClock2(
      tl,
      parametersToEstimate,
      backupListener.get(),
      tolerance,
      nbEvalMax,
      messageHandler,
      profiler,
      optVerbose,
      optMethod);
  }
  else
    throw Exception("Unknown optimization method: " + optName);

  string finalMethod = ApplicationTools::getStringParameter("optimization.final", params, "none", suffix, suffixIsOptional, warn + 1);
  Optimizer* finalOptimizer  = 0;
  if (finalMethod == "none")
  {}
  else if (finalMethod == "simplex")
  {
    finalOptimizer = new DownhillSimplexMethod(tl);
  }
  else if (finalMethod == "powell")
  {
    finalOptimizer = new PowellMultiDimensions(tl);
  }
  else
    throw Exception("Unknown final optimization method: " + finalMethod);

  if (finalOptimizer)
  {
    parametersToEstimate.matchParametersValues(tl->getParameters());
    ApplicationTools::displayResult("Final optimization step", finalMethod);
    finalOptimizer->setProfiler(profiler);
    finalOptimizer->setMessageHandler(messageHandler);
    finalOptimizer->setMaximumNumberOfEvaluations(nbEvalMax);
    finalOptimizer->getStopCondition()->setTolerance(tolerance);
    finalOptimizer->setVerbose(verbose);
    finalOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    finalOptimizer->init(parametersToEstimate);
    finalOptimizer->optimize();
    n += finalOptimizer->getNumberOfEvaluations();
    delete finalOptimizer;
  }

  if (verbose)
    ApplicationTools::displayResult("Performed", TextTools::toString(n) + " function evaluations.");
  if (backupFile != "none")
  {
    string bf=backupFile+".def";
    rename(backupFile.c_str(),bf.c_str());
  }
}

/******************************************************************************/

void PhylogeneticsApplicationTools::checkEstimatedParameters(const ParameterList& pl)
{
  for (size_t i = 0; i < pl.size(); ++i)
  {
    auto constraint = pl[i].getConstraint();
    if (constraint)
    {
      double value = pl[i].getValue();
      if (!constraint->isCorrect(value - 1e-6) || !constraint->isCorrect(value + 1e-6))
      {
        ApplicationTools::displayWarning("This parameter has a value close to the boundary: " + pl[i].getName() + "(" + TextTools::toString(value) + ").");
      }
    }
  }
}

/******************************************************************************/

void PhylogeneticsApplicationTools::writeTree(
  const TreeTemplate<Node>& tree,
  map<string, string>& params,
  const string& prefix,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  bool checkOnly,
  int warn)
{
  string format = ApplicationTools::getStringParameter(prefix + "tree.format", params, "Newick", suffix, suffixIsOptional, warn);
  string file = ApplicationTools::getAFilePath(prefix + "tree.file", params, false, false, suffix, suffixIsOptional, "none", warn);
 
  BppOTreeWriterFormat bppoWriter(warn);
  unique_ptr<OTree> oTree(bppoWriter.read(format));
  if (verbose)
  {
    ApplicationTools::displayResult("Output tree file " + suffix, file);
    ApplicationTools::displayResult("Output tree format " + suffix, oTree->getFormatName());
  }
  if (!checkOnly && file != "none")
    oTree->writeTree(tree, file, true);  
}

/******************************************************************************/

void PhylogeneticsApplicationTools::writeTrees(
  const vector<Tree*>& trees,
  map<string, string>& params,
  const string& prefix,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  bool checkOnly,
  int warn)
{
  string format = ApplicationTools::getStringParameter(prefix + "trees.format", params, "Newick", suffix, suffixIsOptional, warn);
  string file = ApplicationTools::getAFilePath(prefix + "trees.file", params, true, false, suffix, suffixIsOptional, "none", warn);

  BppOMultiTreeWriterFormat bppoWriter(warn);
  unique_ptr<OMultiTree> oTrees(bppoWriter.read(format));
  if (verbose)
  {
    ApplicationTools::displayResult("Output trees file " + suffix, file);
    ApplicationTools::displayResult("Output trees format " + suffix, oTrees->getFormatName());
  }
  if (!checkOnly && file != "none")
    oTrees->writeTrees(trees, file, true);
  
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printParameters(const TransitionModel* model, OutputStream& out, int warn, bool withAlias)
{
  out << "model=";
  map<string, string> globalAliases;
  vector<string> writtenNames;
  BppOSubstitutionModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, false, warn);
  bIO.write(*model, out, globalAliases, writtenNames);
  out.endLine();
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printParameters(const SubstitutionModelSet* modelSet, OutputStream& out, int warn, bool withAlias)
{
  (out << "nonhomogeneous=general").endLine();
  (out << "nonhomogeneous.number_of_models=" << modelSet->getNumberOfModels()).endLine();

  if (modelSet->isStationary())
    (out << "nonhomogeneous.stationarity = yes");
    
  // Get the parameter links:
  map< size_t, vector<string> > modelLinks; // for each model index, stores the list of global parameters.
  map< string, set<size_t> > parameterLinks; // for each parameter name, stores the list of model indices, wich should be sorted.
  vector<string> writtenNames;

  // Loop over all models:
  for (size_t i = 0; i < modelSet->getNumberOfModels(); i++)
  {
    const TransitionModel* model = modelSet->getModel(i);

    map<string, string> aliases;

    // First get the aliases for this model:

    if (withAlias)
    {
      ParameterList pl = model->getParameters();

      for (size_t np = 0; np < pl.size(); np++)
      {
        string nfrom = modelSet->getFrom(pl[np].getName() + "_" + TextTools::toString(i + 1));
        if (nfrom != "")
          aliases[pl[np].getName()] = nfrom;
      }
    }

    // Now print it:
    writtenNames.clear();
    out.endLine() << "model" << (i + 1) << "=";
    BppOSubstitutionModelFormat bIOsm(BppOSubstitutionModelFormat::ALL, true, true, true, false, warn);
    bIOsm.write(*model, out, aliases, writtenNames);

    out.endLine();
    vector<int> ids = modelSet->getNodesWithModel(i);
    out << "model" << (i + 1) << ".nodes_id=" << ids[0];
    for (size_t j = 1; j < ids.size(); ++j)
    {
      out << "," << ids[j];
    }
    out.endLine();
  }

  // First get the aliases for this frequencies set

  if (!modelSet->isStationary())
  {
    
    const FrequenciesSet* pFS = modelSet->getRootFrequenciesSet();

    ParameterList plf = pFS->getParameters();

    map<string, string> aliases;
    
    if (withAlias)
    {
      for (size_t np = 0; np < plf.size(); np++)
      {
        string nfrom = modelSet->getFrom(plf[np].getName());
        if (nfrom != "")
          aliases[plf[np].getName()] = nfrom;
      }
    }
    
    // Root frequencies:
    out.endLine();
    (out << "# Root frequencies:").endLine();
    out << "nonhomogeneous.root_freq=";
    
    BppOFrequenciesSetFormat bIO(BppOFrequenciesSetFormat::ALL, false, warn);
    bIO.write(pFS, out, aliases, writtenNames);
  }
  
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printParameters(const DiscreteDistribution* rDist, OutputStream& out, bool withAlias)
{
  out << "rate_distribution=";
  map<string, string> globalAliases;
  vector<string> writtenNames;
  const BppORateDistributionFormat* bIO = new BppORateDistributionFormat(true);

  bIO->write(*rDist, out, globalAliases, writtenNames);
  delete bIO;
  out.endLine();
}

/************************
* Substitution Mapping *
************************/

SubstitutionCount* PhylogeneticsApplicationTools::getSubstitutionCount(
  const Alphabet* alphabet,
  const SubstitutionModel* model,
  map<string, string>& params,
  string suffix,
  bool verbose,
  int warn)
{
  SubstitutionCount* substitutionCount = 0;
  string nijtOption;
  map<string, string> nijtParams;
  string nijtText = ApplicationTools::getStringParameter("nijt", params, "Uniformization", suffix, true, warn);
  KeyvalTools::parseProcedure(nijtText, nijtOption, nijtParams);

  if (nijtOption == "Laplace")
  {
    size_t trunc = ApplicationTools::getParameter<size_t>("trunc", nijtParams, 10, suffix, true, warn + 1);
    substitutionCount = new LaplaceSubstitutionCount(model, trunc);
  }
  else if (nijtOption == "Uniformization")
  {
    string weightOption = ApplicationTools::getStringParameter("weight", nijtParams, "None", "", true, warn + 1);
    AlphabetIndex2* weights = SequenceApplicationTools::getAlphabetIndex2(alphabet, weightOption, "Substitution weight scheme:");
    substitutionCount = new UniformizationSubstitutionCount(model, new TotalSubstitutionRegister(model), weights);
  }
  else if (nijtOption == "Decomposition")
  {
    string weightOption = ApplicationTools::getStringParameter("weight", nijtParams, "None", "", true, warn + 1);
    AlphabetIndex2* weights = SequenceApplicationTools::getAlphabetIndex2(alphabet, weightOption, "Substitution weight scheme:");
    const ReversibleSubstitutionModel* revModel = dynamic_cast<const ReversibleSubstitutionModel*>(model);
    if (revModel)
      substitutionCount = new DecompositionSubstitutionCount(revModel, new TotalSubstitutionRegister(model), weights);
    else
      throw Exception("Decomposition method can only be used with reversible substitution models.");
  }
  else if (nijtOption == "Naive")
  {
    string weightOption = ApplicationTools::getStringParameter("weight", nijtParams, "None", "", true, warn + 1);
    AlphabetIndex2* weights = SequenceApplicationTools::getAlphabetIndex2(alphabet, weightOption, "Substitution weight scheme:");
    substitutionCount = new NaiveSubstitutionCount(model, new TotalSubstitutionRegister(model), false, weights);
  }
  else if (nijtOption == "Label")
  {
    substitutionCount = new LabelSubstitutionCount(model);
  }
  else if (nijtOption == "ProbOneJump")
  {
    substitutionCount = new OneJumpSubstitutionCount(model);
  }
  else
  {
    ApplicationTools::displayError("Invalid option '" + nijtOption + ", in 'nijt' parameter.");
    exit(-1);
  }
  ApplicationTools::displayResult("Substitution count procedure", nijtOption);

  // Send results:
  return substitutionCount;
}

/****************************************************************************/


SubstitutionRegister* PhylogeneticsApplicationTools::getSubstitutionRegister(const std::string& regTypeDesc, const SubstitutionModel* model, bool verbose)
{
  string regType = "";
  map<string, string> regArgs;
  KeyvalTools::parseProcedure(regTypeDesc, regType, regArgs);
  
  SubstitutionRegister* reg = 0;

  if (regType=="Combination")
  {
    VectorOfSubstitionRegisters* vreg= new VectorOfSubstitionRegisters(model);

    size_t i = 0;
    while (++i)
    {
      string regDesc = ApplicationTools::getStringParameter("reg" + TextTools::toString(i), regArgs, "", "", false, 1);
      if (regDesc=="")
        break;
      
      SubstitutionRegister* sreg=getSubstitutionRegister(regDesc, model);

      vreg->addRegister(sreg);
    }
    
    reg=vreg;
  }
  else if (regType == "All")
  {
    reg = new ComprehensiveSubstitutionRegister(model, false);
  }
  else if (regType == "Total")
  {
    reg = new TotalSubstitutionRegister(model);
  }    
  else if (regType == "Selected"){  
    string subsList = ApplicationTools::getStringParameter("substitution.list", regArgs, "All", "", true, false);
    reg = new SelectedSubstitutionRegister(model, subsList);  
  }

  
  // Alphabet dependent registers

  
  else if (AlphabetTools::isNucleicAlphabet(model->getAlphabet()))
  {
    const NucleotideSubstitutionModel* nmodel=dynamic_cast<const NucleotideSubstitutionModel*>(model);
    if (!nmodel)
    {
      const WrappedSubstitutionModel* wmodel=dynamic_cast<const WrappedSubstitutionModel*>(model);
      if (wmodel)
      {
        nmodel=dynamic_cast<const NucleotideSubstitutionModel*>(&wmodel->getSubstitutionModel());
      }
    }

    if (!nmodel)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionRegister : model and alphabet do not fit " + model->getAlphabet()->getAlphabetType() + " vs " + model->getName());

    
    if (regType == "GC")
      reg = new GCSubstitutionRegister(nmodel, false);
    else if (regType == "TsTv")
      reg = new TsTvSubstitutionRegister(nmodel);
    else if (regType == "SW")
      reg = new SWSubstitutionRegister(nmodel);
    else
      throw Exception("Unsupported substitution categorization: " + regType + " for alphabet " + model->getAlphabet()->getAlphabetType());
  }
  
  else if (AlphabetTools::isCodonAlphabet(model->getAlphabet()))
  {
    const CodonSubstitutionModel* cmodel=dynamic_cast<const CodonSubstitutionModel*>(model);
    if (!cmodel)
    {
      const WrappedSubstitutionModel* wmodel=dynamic_cast<const WrappedSubstitutionModel*>(model);
      if (wmodel)
      {
        cmodel=dynamic_cast<const CodonSubstitutionModel*>(&wmodel->getSubstitutionModel());
      }
    }

    if (!cmodel)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionRegister : model and alphabet do not fit " + model->getAlphabet()->getAlphabetType() + " vs " + model->getName());
  
    if (regType == "IntraAA")
      reg = new AAInteriorSubstitutionRegister(cmodel);
    else if (regType == "InterAA")
      reg = new AAExteriorSubstitutionRegister(cmodel);
    else if (regType == "GC")
      reg = new GCSynonymousSubstitutionRegister(cmodel);
    else if (regType == "TsTv")
      reg = new TsTvSubstitutionRegister(cmodel);
    else if (regType == "SW")
      reg = new SWSubstitutionRegister(cmodel);
    else if (regType == "KrKc")
      reg = new KrKcSubstitutionRegister(cmodel);
    else if (regType == "DnDs")
      reg = new DnDsSubstitutionRegister(cmodel, false);
    else
      throw Exception("Unsupported substitution categorization: " + regType + " for alphabet " + model->getAlphabet()->getAlphabetType());
  }
  
  else if (AlphabetTools::isProteicAlphabet(model->getAlphabet()))
  {
    const ProteinSubstitutionModel* pmodel=dynamic_cast<const ProteinSubstitutionModel*>(model);
    if (!pmodel)
    {
      const WrappedSubstitutionModel* wmodel=dynamic_cast<const WrappedSubstitutionModel*>(model);
      if (wmodel)
      {
        pmodel=dynamic_cast<const ProteinSubstitutionModel*>(&wmodel->getSubstitutionModel());
      }
    }

    if (!pmodel)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionRegister : model and alphabet do not fit " + model->getAlphabet()->getAlphabetType() + " vs " + model->getName());
  
    if (regType == "KrKc")
      reg = new KrKcSubstitutionRegister(pmodel);
    else
      throw Exception("Unsupported substitution categorization: " + regType + " for alphabet " + model->getAlphabet()->getAlphabetType());
  }
  
  CategorySubstitutionRegister* csr=dynamic_cast<CategorySubstitutionRegister*>(reg);
  if (csr)
    csr->setStationarity(ApplicationTools::getBooleanParameter("stationarity", regArgs, true));

  if (verbose)
    ApplicationTools::displayResult("Substitution Register", regType);

  return reg;
}

