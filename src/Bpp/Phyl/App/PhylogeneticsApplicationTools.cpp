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
#include "../Model.all"
#include "../Likelihood.all"
#include "../OptimizationTools.h"
#include "../Tree.h"
#include "../Io/Newick.h"
#include "../Io/NexusIoTree.h"
#include "../Io/Nhx.h"
#include "../Io/BppOSubstitutionModelFormat.h"
#include "../Io/BppOFrequenciesSetFormat.h"

#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/Prob.all>
#include <Bpp/Numeric/Function.all>

// From SeqLib:
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
  bool verbose) throw (Exception)
{
  string format = ApplicationTools::getStringParameter(prefix + "tree.format", params, "Newick", suffix, suffixIsOptional, true);
  string treeFilePath = ApplicationTools::getAFilePath(prefix + "tree.file", params, true, true, suffix, suffixIsOptional);

  ITree* treeReader;
  if (format == "Newick")
    treeReader = new Newick(true);
  else if (format == "Nexus")
    treeReader = new NexusIOTree();
  else if (format == "NHX")
    treeReader = new Nhx();
  else
    throw Exception("Unknow format for tree reading: " + format);
  Tree* tree = treeReader->read(treeFilePath);
  delete treeReader;

  if (verbose)
    ApplicationTools::displayResult("Tree file", treeFilePath);
  return tree;
}

/******************************************************************************/

vector<Tree*> PhylogeneticsApplicationTools::getTrees(
  map<string, string>& params,
  const string& prefix,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  string format = ApplicationTools::getStringParameter(prefix + "trees.format", params, "Newick", suffix, suffixIsOptional, true);
  string treeFilePath = ApplicationTools::getAFilePath(prefix + "trees.file", params, true, true, suffix, suffixIsOptional);

  IMultiTree* treeReader;
  if (format == "Newick")
    treeReader = new Newick(true);
  else if (format == "Nexus")
    treeReader = new NexusIOTree();
  else if (format == "NHX")
    treeReader = new Nhx();
  else
    throw Exception("Unknow format for tree reading: " + format);
  vector<Tree*> trees;
  treeReader->read(treeFilePath, trees);
  delete treeReader;

  if (verbose)
  {
    ApplicationTools::displayResult("Tree file", treeFilePath);
    ApplicationTools::displayResult("Number of trees in file", trees.size());
  }
  return trees;
}


/*************************************************************/
/******* MODELS **********************************************/
/*************************************************************/

/******************************************************************************/

SubstitutionModel* PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance(
  const Alphabet* alphabet,
  const string& modelDescription,
  map<string, string>& unparsedParameterValues,
  bool allowCovarions,
  bool allowMixed,
  bool allowGaps,
  bool verbose) throw (Exception)
{
  BppOSubstitutionModelFormat* bIO = new BppOSubstitutionModelFormat();
  SubstitutionModel* pS = bIO->read(alphabet, modelDescription, unparsedParameterValues, allowCovarions, allowMixed, allowMixed, verbose);
  delete bIO;
  return pS;
  
}

/******************************************************************************/

SubstitutionModel* PhylogeneticsApplicationTools::getSubstitutionModel(
  const Alphabet* alphabet,
  const SiteContainer* data,
  std::map<std::string, std::string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  string modelDescription;
  if (AlphabetTools::isCodonAlphabet(alphabet))
    modelDescription = ApplicationTools::getStringParameter("model", params, "CodonRate(model=JC69)", suffix, suffixIsOptional, verbose);
  else if (AlphabetTools::isWordAlphabet(alphabet))
    modelDescription = ApplicationTools::getStringParameter("model", params, "Word(model=JC69)", suffix, suffixIsOptional, verbose);
  else
    modelDescription = ApplicationTools::getStringParameter("model", params, "JC69", suffix, suffixIsOptional, verbose);

  map<string, string> unparsedParameterValues;
  SubstitutionModel* model = getSubstitutionModelDefaultInstance(alphabet, modelDescription, unparsedParameterValues, true, true, true, verbose);
  setSubstitutionModelParametersInitialValues(model, unparsedParameterValues, data, verbose);
  return model;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::setSubstitutionModelParametersInitialValues(
  SubstitutionModel* model,
  std::map<std::string, std::string>& unparsedParameterValues,
  const SiteContainer* data,
  bool verbose) throw (Exception)
{
  string initFreqs = ApplicationTools::getStringParameter(model->getNamespace() + "initFreqs", unparsedParameterValues, "", "", true, false);
  if (verbose)
    ApplicationTools::displayResult("Frequencies initialization for model", (initFreqs == "") ? "None" : initFreqs);

  if (initFreqs != "")
  {
    if (initFreqs == "observed")
    {
      if (!data)
        throw Exception("Missing data for observed frequencies");
      unsigned int psi = ApplicationTools::getParameter<unsigned int>(model->getNamespace() + "initFreqs.observedPseudoCount", unparsedParameterValues, 0);
      model->setFreqFromData(*data, psi);
    }
    else if (initFreqs.substr(0, 6) == "values")
    {
      // Initialization using the "values" argument
      map<int, double> frequencies;

      string rf = initFreqs.substr(6);
      StringTokenizer strtok(rf.substr(1, rf.length() - 2), ",");
      unsigned int i = 0;
      while (strtok.hasMoreToken())
        frequencies[i++] = TextTools::toDouble(strtok.nextToken());
      model->setFreq(frequencies);
    }
    else
      throw Exception("Unknown initFreqs argument");
  }

  ParameterList pl = model->getIndependentParameters();
  for (unsigned int i = 0; i < pl.size(); i++)
  {
    AutoParameter ap(pl[i]);
    ap.setMessageHandler(ApplicationTools::warning);
    pl.setParameter(i, ap);
  }
  unsigned int posp;
  for (unsigned int i = 0; i < pl.size(); i++)
  {
    const string pName = pl[i].getName();
    posp = pName.rfind(".");
    bool test1 = (initFreqs == "");
    bool test2 = (model->getParameterNameWithoutNamespace(pName).size() < posp + 6) || (model->getParameterNameWithoutNamespace(pName).substr(posp + 1, 5) != "theta");
    bool test3 = (unparsedParameterValues.find(pName) != unparsedParameterValues.end());
    if (test1 || test2 || test3)
    {
      if (!test1 && !test2 && test3)
        ApplicationTools::displayWarning("Warning, initFreqs argument is set and a value is set for parameter " + pName);
      double value = ApplicationTools::getDoubleParameter(pName, unparsedParameterValues, pl[i].getValue());
      pl[i].setValue(value);
    }
    if (verbose)
      ApplicationTools::displayResult("Parameter found", pName + "=" + TextTools::toString(pl[i].getValue()));
  }

  model->matchParametersValues(pl);
}

/******************************************************************************/

void PhylogeneticsApplicationTools::setSubstitutionModelParametersInitialValues(
  SubstitutionModel* model,
  std::map<std::string, std::string>& unparsedParameterValues,
  const std::string& modelPrefix,
  const SiteContainer* data,
  std::map<std::string, double>& existingParams,
  std::vector<std::string>& specificParams,
  std::vector<std::string>& sharedParams,
  bool verbose) throw (Exception)
{
  string initFreqs = ApplicationTools::getStringParameter(model->getNamespace() + "initFreqs", unparsedParameterValues, "", "", true, false);

  if (verbose)
    ApplicationTools::displayResult("Frequencies Initialization for model", (initFreqs == "") ? "None" : initFreqs);

  if (initFreqs != "")
  {
    if (initFreqs == "observed")
    {
      if (!data)
        throw Exception("Missing data for observed frequencies");
      unsigned int psi = ApplicationTools::getParameter<unsigned int>(model->getNamespace() + "initFreqs.observedPseudoCount", unparsedParameterValues, 0);
      model->setFreqFromData(*data, psi);
    }
    else if (initFreqs.substr(0, 6) == "values")
    {
      // Initialization using the "values" argument
      map<int, double> frequencies;

      string rf = initFreqs.substr(6);
      StringTokenizer strtok(rf.substr(1, rf.length() - 2), ",");
      unsigned int i = 0;
      while (strtok.hasMoreToken())
        frequencies[i++] = TextTools::toDouble(strtok.nextToken());
      model->setFreq(frequencies);
    }
    else
      throw Exception("Unknown initFreqs argument");
  }

  ParameterList pl = model->getIndependentParameters();
  for (unsigned int i = 0; i < pl.size(); i++)
  {
    AutoParameter ap(pl[i]);
    ap.setMessageHandler(ApplicationTools::warning);
    pl.setParameter(i, ap);
  }

  for (unsigned int i = 0; i < pl.size(); i++)
  {
    const string pName = pl[i].getName();
    unsigned int posp = model->getParameterNameWithoutNamespace(pName).rfind(".");
    string value;
    bool test1 = (initFreqs == "");
    bool test2 = (model->getParameterNameWithoutNamespace(pName).substr(posp + 1, 5) != "theta");
    bool test3 = (unparsedParameterValues.find(pName) != unparsedParameterValues.end());
    if (test1 || test2 || test3)
    {
      if (!test1 && !test2 && test3)
        ApplicationTools::displayWarning("Warning, initFreqs argument is set and a value is set for parameter " + pName);
      value = ApplicationTools::getStringParameter(pName, unparsedParameterValues, TextTools::toString(pl[i].getValue()));
      if (value.size() > 5 && value.substr(0, 5) == "model")
      {
        if (existingParams.find(value) != existingParams.end())
        {
          pl[i].setValue(existingParams[value]);
          sharedParams.push_back(value);
        }
        else
          throw Exception("Error, unknown parameter " + value);
      }
      else
      {
        double value2 = TextTools::toDouble(value);
        existingParams[modelPrefix + pName] = value2;
        specificParams.push_back(pName);
        pl[i].setValue(value2);
      }
    }
    else
    {
      existingParams[modelPrefix + pName] = pl[i].getValue();
      specificParams.push_back(pName);
    }
    if (verbose)
      ApplicationTools::displayResult("Parameter found", modelPrefix + pName + "=" + TextTools::toString(pl[i].getValue()));
  }
  model->matchParametersValues(pl);
}


/******************************************************/
/**** FREQUENCIES SET *********************************/
/******************************************************/

/******************************************************************************/

FrequenciesSet* PhylogeneticsApplicationTools::getRootFrequenciesSet(
  const Alphabet* alphabet,
  const SiteContainer* data,
  std::map<std::string, std::string>& params,
  const std::vector<double>& rateFreqs,
  const std::string& suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  string freqDescription = ApplicationTools::getStringParameter("nonhomogeneous.root_freq", params, "Full(init=observed)", suffix, suffixIsOptional);
  if (freqDescription == "None")
  {
    return 0;
  }
  else
  {
    FrequenciesSet* freq = getFrequenciesSet(alphabet, freqDescription, data, rateFreqs, verbose);
    if (verbose)
      ApplicationTools::displayResult("Root frequencies ", freq->getName());
    return freq;
  }
}

/******************************************************************************/

FrequenciesSet* PhylogeneticsApplicationTools::getFrequenciesSet(
  const Alphabet* alphabet,
  const std::string& freqDescription,
  const SiteContainer* data,
  const std::vector<double>& rateFreqs,
  bool verbose) throw (Exception)
{
  map<string, string> unparsedParameterValues;
  FrequenciesSet* pFS = getFrequenciesSetDefaultInstance(alphabet, freqDescription, unparsedParameterValues);

  // Now we set the initial frequencies according to options:
  if (unparsedParameterValues.find("init") != unparsedParameterValues.end())
  {
    // Initialization using the "init" option
    string init = unparsedParameterValues["init"];
    if (init == "observed")
    {
      if (!data)
        throw Exception("Missing data for observed frequencies");
      unsigned int psc = 0;
      if (unparsedParameterValues.find("init.observedPseudoCount") != unparsedParameterValues.end())
        psc = TextTools::toInt(unparsedParameterValues["init.observedPseudoCount"]);

      map<int, double> freqs;
      SequenceContainerTools::getFrequencies(*data, freqs, psc);

      pFS->setFrequenciesFromMap(freqs);
    }
    else if (init.substr(0, 6) == "values")
    {
      // Initialization using the "values" argument
      vector<double> frequencies;
      string rf = init.substr(6);

      StringTokenizer strtok(rf.substr(1, rf.length() - 2), ",");
      while (strtok.hasMoreToken())
        frequencies.push_back(TextTools::toDouble(strtok.nextToken()));
      pFS->setFrequencies(frequencies);
    }
    else if (init == "balanced")
    {
      // Nothing to do here, this is the default instanciation.
    }
    else
      throw Exception("Unknown init argument");
  }
  else
  {
    // Explicit initialization of each parameter
    ParameterList pl = pFS->getParameters();

    for (unsigned int i = 0; i < pl.size(); i++)
    {
      AutoParameter ap(pl[i]);
      if (verbose)
        ap.setMessageHandler(ApplicationTools::warning);
      pl.setParameter(i, ap);
    }

    for (unsigned int i = 0; i < pl.size(); i++)
    {
      const string pName = pl[i].getName();
      double value = ApplicationTools::getDoubleParameter(pName, unparsedParameterValues, pl[i].getValue());

      pl[i].setValue(value);
      if (verbose)
        ApplicationTools::displayResult("Parameter found", pName + "=" + TextTools::toString(pl[i].getValue()));
    }

    pFS->matchParametersValues(pl);
  }

  // /////// To be changed for input normalization
  if (rateFreqs.size() > 0)
  {
    pFS = new MarkovModulatedFrequenciesSet(pFS, rateFreqs);
  }

  return pFS;
}

/******************************************************************************/


FrequenciesSet* PhylogeneticsApplicationTools::getFrequenciesSetDefaultInstance(
  const Alphabet* alphabet,
  const std::string& freqDescription,
  std::map<std::string, std::string>& unparsedParameterValues) throw (Exception)
{
  BppOFrequenciesSetFormat* bIO=new BppOFrequenciesSetFormat();
  FrequenciesSet* pFS= bIO->read(alphabet, freqDescription, unparsedParameterValues);
  delete bIO;
  return pFS;
}


/******************************************************/
/**** SUBSTITUTION MODEL SET **************************/
/******************************************************/

/******************************************************************************/

SubstitutionModelSet* PhylogeneticsApplicationTools::getSubstitutionModelSet(
  const Alphabet* alphabet,
  const SiteContainer* data,
  std::map<std::string, std::string>& params,
  const std::string& suffix,
  bool suffixIsOptional,
  bool verbose)
{
  if (!ApplicationTools::parameterExists("nonhomogeneous.number_of_models", params))
    throw Exception("You must specify this parameter: nonhomogeneous.number_of_models .");
  unsigned int nbModels = ApplicationTools::getParameter<unsigned int>("nonhomogeneous.number_of_models", params, 1, suffix, suffixIsOptional, false);
  if (nbModels == 0)
    throw Exception("The number of models can't be 0 !");

  bool nomix = true;
  for (unsigned int i = 0; nomix &(i < nbModels); i++)
  {
    string prefix = "model" + TextTools::toString(i + 1);
    string modelDesc;
    modelDesc = ApplicationTools::getStringParameter(prefix, params, "", suffix, suffixIsOptional, verbose);

    if (modelDesc.find("Mixed") != string::npos)
      nomix = false;
  }

  SubstitutionModelSet* modelSet, * modelSet1 = 0;
  modelSet1 = new SubstitutionModelSet(alphabet);
  setSubstitutionModelSet(*modelSet1, alphabet, data, params, suffix, suffixIsOptional, verbose);

  if (modelSet1->hasMixedSubstitutionModel())
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
  const SiteContainer* data,
  map<string, string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose)
{
  modelSet.clear();
  if (!ApplicationTools::parameterExists("nonhomogeneous.number_of_models", params))
    throw Exception("You must specify this parameter: nonhomogeneous.number_of_models .");
  unsigned int nbModels = ApplicationTools::getParameter<unsigned int>("nonhomogeneous.number_of_models", params, 1, suffix, suffixIsOptional, false);
  if (nbModels == 0)
    throw Exception("The number of models can't be 0 !");

  if (verbose)
    ApplicationTools::displayResult("Number of distinct models", TextTools::toString(nbModels));


  // ///////////////////////////////////////////
  // Build a new model set object:

  vector<double> rateFreqs;
  string tmpDesc;
  if (AlphabetTools::isCodonAlphabet(alphabet))
    tmpDesc = ApplicationTools::getStringParameter("model1", params, "CodonRate(model=JC69)", suffix, suffixIsOptional, false);
  else if (AlphabetTools::isWordAlphabet(alphabet))
    tmpDesc = ApplicationTools::getStringParameter("model1", params, "Word(model=JC69)", suffix, suffixIsOptional, false);
  else
    tmpDesc = ApplicationTools::getStringParameter("model1", params, "JC69", suffix, suffixIsOptional, false);

  map<string, string> tmpUnparsedParameterValues;

  auto_ptr<SubstitutionModel> tmp(getSubstitutionModelDefaultInstance(alphabet, tmpDesc, tmpUnparsedParameterValues, true, true, true, false));
  if (tmp->getNumberOfStates() != alphabet->getSize())
  {
    // Markov-Modulated Markov Model...
    unsigned int n = static_cast<unsigned int>(tmp->getNumberOfStates() / alphabet->getSize());
    rateFreqs = vector<double>(n, 1. / static_cast<double>(n)); // Equal rates assumed for now, may be changed later (actually, in the most general case,
  }

  // ////////////////////////////////////
  // Deal with root frequencies

  bool stationarity = ApplicationTools::getBooleanParameter("nonhomogeneous.stationarity", params, false, "", false, false);
  FrequenciesSet* rootFrequencies = 0;
  if (!stationarity)
  {
    rootFrequencies = getRootFrequenciesSet(alphabet, data, params, rateFreqs, suffix, suffixIsOptional, verbose);
    stationarity = !rootFrequencies;
  }
  ApplicationTools::displayBooleanResult("Stationarity assumed", stationarity);

  if (!stationarity)
    modelSet.setRootFrequencies(rootFrequencies);

  // //////////////////////////////////////
  // Now parse all models:

  map<string, double> existingParameters;

  for (unsigned int i = 0; i < nbModels; i++)
  {
    string prefix = "model" + TextTools::toString(i + 1);
    string modelDesc;
    if (AlphabetTools::isCodonAlphabet(alphabet))
      modelDesc = ApplicationTools::getStringParameter(prefix, params, "CodonRate(model=JC69)", suffix, suffixIsOptional, verbose);
    else if (AlphabetTools::isWordAlphabet(alphabet))
      modelDesc = ApplicationTools::getStringParameter(prefix, params, "Word(model=JC69)", suffix, suffixIsOptional, verbose);
    else
      modelDesc = ApplicationTools::getStringParameter(prefix, params, "JC69", suffix, suffixIsOptional, verbose);


    map<string, string> unparsedParameterValues;
    SubstitutionModel* model = getSubstitutionModelDefaultInstance(alphabet, modelDesc, unparsedParameterValues, true, true, true, verbose);
    prefix += ".";

    vector<string> specificParameters, sharedParameters;
    setSubstitutionModelParametersInitialValues(model,
                                                unparsedParameterValues, prefix, data,
                                                existingParameters, specificParameters, sharedParameters,
                                                verbose);
    vector<int> nodesId = ApplicationTools::getVectorParameter<int>(prefix + "nodes_id", params, ',', ':', TextTools::toString(i), suffix, suffixIsOptional, true);
    if (verbose)
      ApplicationTools::displayResult("Model" + TextTools::toString(i + 1) + " is associated to", TextTools::toString(nodesId.size()) + " node(s).");
    // Add model and specific parameters:
    //DEBUG: cout << "Specific parameters:" << endl;
    //DEBUG: VectorTools::print(specificParameters);
    modelSet.addModel(model, nodesId, specificParameters);
    // Now set shared parameters:
    for (unsigned int j = 0; j < sharedParameters.size(); j++)
    {
      string pName = sharedParameters[j];
      //DEBUG: cout << "Shared parameter found: " << pName << endl;
      string::size_type index = pName.find(".");
      if (index == string::npos)
        throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelSet. Bad parameter name: " + pName);
      string name = pName.substr(index + 1) + "_" + pName.substr(5, index - 5);
      // namespace checking:
      vector<unsigned int> models = modelSet.getModelsWithParameter(name);
      if (models.size() == 0)
        throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelSet. Parameter `" + name + "' is not associated to any model.");
      if (model->getNamespace() == modelSet.getModel(models[0])->getNamespace())
        modelSet.setParameterToModel(modelSet.getParameterIndex(name), modelSet.getNumberOfModels() - 1);
      else
      {
        throw Exception("Assigning a value to a parameter with a distinct namespace is not (yet) allowed. Consider using parameter aliasing instead.");
      }
    }
  }
  // Finally check parameter aliasing:
  string aliasDesc = ApplicationTools::getStringParameter("nonhomogeneous.alias", params, "", suffix, suffixIsOptional, verbose);
  StringTokenizer st(aliasDesc, ",");
  while (st.hasMoreToken())
  {
    string alias = st.nextToken();
    string::size_type index = alias.find("->");
    if (index == string::npos)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelSet. Bad alias syntax, should contain `->' symbol: " + alias);
    string p1 = alias.substr(0, index);
    string p2 = alias.substr(index + 2);
    ApplicationTools::displayResult("Parameter alias found", p1 + "->" + p2);
    modelSet.aliasParameters(p1, p2);
  }
}

/******************************************************************************/
void PhylogeneticsApplicationTools::completeMixedSubstitutionModelSet(
  MixedSubstitutionModelSet& mixedModelSet,
  const Alphabet* alphabet,
  const SiteContainer* data,
  map<string, string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose)
{
  // /////////////////////////////////////////
  // Looks for the allowed paths

  unsigned int numd;
  if (!ApplicationTools::parameterExists("site.number_of_paths", params))
    numd = 0;
  else
    numd = ApplicationTools::getParameter<unsigned int>("site.number_of_paths", params, 1, suffix, suffixIsOptional, false);

  if (verbose)
    ApplicationTools::displayResult("Number of distinct paths", TextTools::toString(numd));

  vector<string> vdesc;
  while (numd)
  {
    string desc = ApplicationTools::getStringParameter("site.path" + TextTools::toString(numd), params, "",  suffix, suffixIsOptional, verbose);
    if (desc.size() == 0)
      break;
    else
      vdesc.push_back(desc);
    numd--;
  }

  if (vdesc.size() == 0)
  {
    mixedModelSet.complete();
    mixedModelSet.computeHyperNodesProbabilities();
    return;
  }

  for (vector<string>::iterator it(vdesc.begin()); it != vdesc.end(); it++)
  {
    mixedModelSet.addEmptyHyperNode();
    StringTokenizer st(*it, "&");
    while (st.hasMoreToken())
    {
      string submodel = st.nextToken();
      string::size_type indexo = submodel.find("[");
      string::size_type indexf = submodel.find("]");
      if ((indexo == string::npos) | (indexf == string::npos))
        throw Exception("PhylogeneticsApplicationTools::setMixedSubstitutionModelSet. Bad path syntax, should contain `[]' symbols: " + submodel);
      int num = TextTools::toInt(submodel.substr(5, indexo - 5));
      string p2 = submodel.substr(indexo + 1, indexf - indexo - 1);

      const MixedSubstitutionModel* pSM = dynamic_cast<const MixedSubstitutionModel*>(mixedModelSet.getModel(num - 1));
      if (pSM == NULL)
        throw BadIntegerException("PhylogeneticsApplicationTools::setMixedSubstitutionModelSet: Wron gmodel for number", num - 1);
      Vint submodnb = pSM->getSubmodelNumbers(p2);

      mixedModelSet.addToHyperNode(num - 1, submodnb);
    }

    if (!mixedModelSet.getHyperNode(mixedModelSet.getNumberOfHyperNodes() - 1).isComplete())
      throw Exception("A path should own at least a submodel of each mixed model: " + *it);

    if (verbose)
      ApplicationTools::displayResult("Site Path", *it);
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

DiscreteDistribution* PhylogeneticsApplicationTools::getRateDistributionDefaultInstance(
  const string& distDescription,
  map<string, string>& unparsedParameterValues,
  bool constDistAllowed,
  bool verbose) throw (Exception)
{
  string distName;
  DiscreteDistribution* rDist = 0;
  map<string, string> args;
  KeyvalTools::parseProcedure(distDescription, distName, args);

  if (distName == "Invariant")
  {
    // We have to parse the nested distribution first:
    string nestedDistDescription = args["dist"];
    if (TextTools::isEmpty(nestedDistDescription))
      throw Exception("PhylogeneticsApplicationTools::getRateDistributionDefaultInstance. Missing argument 'dist' for distribution 'Invariant'.");
    if (verbose)
      ApplicationTools::displayResult("Invariant Mixed distribution", distName );
    map<string, string> unparsedParameterValuesNested;
    DiscreteDistribution* nestedDistribution = getRateDistributionDefaultInstance(nestedDistDescription, unparsedParameterValuesNested, constDistAllowed, verbose);

    // Now we create the Invariant rate distribution:
    rDist = new InvariantMixedDiscreteDistribution(nestedDistribution, 0.1, 0.000001); // , "Invariant.");

    // Then we update the parameter set:
    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedParameterValues["InvariantMixed.dist_" + it->first] = it->second;
    }
    if (args.find("p") != args.end())
      unparsedParameterValues["InvariantMixed.p"] = args["p"];
  }
  else if (distName == "Uniform")
  {
    if (!constDistAllowed)
      throw Exception("You can't use a constant distribution here!");
    rDist = new ConstantDistribution(1., true);
  }
  else if (distName == "Gamma")
  {
    if (args.find("n") == args.end())
      throw Exception("Missing argument 'n' (number of classes) in Gamma distribution");
    unsigned int nbClasses = TextTools::to<unsigned int>(args["n"]);
    rDist = new GammaDiscreteDistribution(nbClasses, 1., 1.); // , "Gamma.");
    rDist->aliasParameters("alpha", "beta");
    if (args.find("alpha") != args.end())
      unparsedParameterValues["Gamma.alpha"] = args["alpha"];
  }
  else
  {
    throw Exception("Unknown distribution: " + distName + ".");
  }
  if (verbose)
  {
    ApplicationTools::displayResult("Rate distribution", distName);
    ApplicationTools::displayResult("Number of classes", TextTools::toString(rDist->getNumberOfCategories()));
  }

  return rDist;
}

/******************************************************************************/

DiscreteDistribution* PhylogeneticsApplicationTools::getDistributionDefaultInstance(
  const std::string& distDescription,
  std::map<std::string, std::string>& unparsedParameterValues,
  bool verbose)
throw (Exception)
{
  string distName;
  DiscreteDistribution* rDist = 0;
  map<string, string> args;
  KeyvalTools::parseProcedure(distDescription, distName, args);

  if (distName == "InvariantMixed")
  {
    // We have to parse the nested distribution first:
    string nestedDistDescription = args["dist"];
    if (TextTools::isEmpty(nestedDistDescription))
      throw Exception("PhylogeneticsApplicationTools::getDistributionDefaultInstance. Missing argument 'dist' for distribution 'Invariant'.");
    if (verbose)
      ApplicationTools::displayResult("Invariant Mixed distribution", distName );
    map<string, string> unparsedParameterValuesNested;
    DiscreteDistribution* nestedDistribution = getDistributionDefaultInstance(nestedDistDescription,
                                                                              unparsedParameterValuesNested,
                                                                              verbose);

    // Now we create the Invariant rate distribution:
    rDist = new InvariantMixedDiscreteDistribution(nestedDistribution, 0.1, 0.000001);

    // Then we update the parameter set:
    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin();
         it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedParameterValues["InvarianMixed.dist_" + it->first] = it->second;
    }
    if (args.find("p") != args.end())
      unparsedParameterValues["InvariantMixed.p"] = args["p"];
  }
  else if (distName == "Constant")
  {
    if (args.find("value") == args.end())
      throw Exception("Missing argument 'value' in Constant distribution");
    rDist = new ConstantDistribution(TextTools::to<double>(args["value"]));
    unparsedParameterValues["Constant.value"] = args["value"];
  }
  else if (distName == "Simple")
  {
    if (args.find("values") == args.end())
      throw Exception("Missing argument 'values' in Simple distribution");
    if (args.find("probas") == args.end())
      throw Exception("Missing argument 'probas' in Simple distribution");
    vector<double> probas, values;

    string rf = args["values"];
    StringTokenizer strtok(rf.substr(1, rf.length() - 2), ",");
    while (strtok.hasMoreToken())
      values.push_back(TextTools::toDouble(strtok.nextToken()));

    rf = args["probas"];
    StringTokenizer strtok2(rf.substr(1, rf.length() - 2), ",");
    while (strtok2.hasMoreToken())
      probas.push_back(TextTools::toDouble(strtok2.nextToken()));

    rDist = new SimpleDiscreteDistribution(values, probas);
    vector<string> v = rDist->getParameters().getParameterNames();

    for (unsigned int i = 0; i < v.size(); i++)
    {
      unparsedParameterValues[v[i]] = TextTools::toString(rDist->getParameterValue(rDist->getParameterNameWithoutNamespace(v[i])));
    }
  }
  else if (distName == "Mixture")
  {
    if (args.find("probas") == args.end())
      throw Exception("Missing argument 'probas' in Mixture distribution");
    vector<double> probas;
    vector<DiscreteDistribution*> v_pdd;
    DiscreteDistribution* pdd;
    string rf = args["probas"];
    StringTokenizer strtok2(rf.substr(1, rf.length() - 2), ",");
    while (strtok2.hasMoreToken())
      probas.push_back(TextTools::toDouble(strtok2.nextToken()));

    vector<string> v_nestedDistrDescr;

    unsigned int nbd = 0;
    while (args.find("distribution" + TextTools::toString(++nbd)) != args.end())
      v_nestedDistrDescr.push_back(args["distribution" + TextTools::toString(nbd)]);

    if (v_nestedDistrDescr.size() != probas.size())
      throw Exception("Number of distributions (keyword 'distribution" + TextTools::toString(probas.size()) + "') do not fit the number of probabilities");

    map<string, string> unparsedParameterValuesNested;

    for (unsigned i = 0; i < v_nestedDistrDescr.size(); i++)
    {
      unparsedParameterValuesNested.clear();
      pdd = getDistributionDefaultInstance(v_nestedDistrDescr[i], unparsedParameterValuesNested, false);

      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedParameterValues[distName + "." + TextTools::toString(i + 1) + "_" + it->first] = it->second;
      }
      v_pdd.push_back(pdd);
    }
    rDist = new MixtureOfDiscreteDistributions(v_pdd, probas);
  }
  else
  {
    if (args.find("n") == args.end())
      throw Exception("Missing argument 'n' (number of classes) in " + distName
                      + " distribution");
    unsigned int nbClasses = TextTools::to<unsigned int>(args["n"]);

    if (distName == "Gamma")
    {
      rDist = new GammaDiscreteDistribution(nbClasses, 1, 1);
      if (args.find("alpha") != args.end())
        unparsedParameterValues["Gamma.alpha"] = args["alpha"];
      if (args.find("beta") != args.end())
        unparsedParameterValues["Gamma.beta"] = args["beta"];
    }
    else if (distName == "Gaussian")
    {
      rDist = new GaussianDiscreteDistribution(nbClasses, 0, 1);
      if (args.find("mu") != args.end())
        unparsedParameterValues["Gaussian.mu"] = args["mu"];
      if (args.find("sigma") != args.end())
        unparsedParameterValues["Gaussian.sigma"] = args["sigma"];
    }
    else if (distName == "Beta")
    {
      rDist = new BetaDiscreteDistribution(nbClasses, 1, 1);
      if (args.find("alpha") != args.end())
        unparsedParameterValues["Beta.alpha"] = args["alpha"];
      if (args.find("beta") != args.end())
        unparsedParameterValues["Beta.beta"] = args["beta"];
    }
    else if (distName == "Exponential")
    {
      rDist = new ExponentialDiscreteDistribution(nbClasses, 1);
      if (args.find("lambda") != args.end())
        unparsedParameterValues["Exponential.lambda"] = args["lambda"];
      if (args.find("median") != args.end())
        rDist->setMedian(true);
    }
    else if (distName == "TruncExponential")
    {
      rDist = new TruncatedExponentialDiscreteDistribution(nbClasses, 1, 0);

      if (args.find("median") != args.end())
        rDist->setMedian(true);
      if (args.find("lambda") != args.end())
        unparsedParameterValues["TruncExponential.lambda"] = args["lambda"];
      if (args.find("tp") != args.end())
        unparsedParameterValues["TruncExponential.tp"] = args["tp"];
    }
    else if (distName == "Uniform")
    {
      if (args.find("begin") == args.end())
        throw Exception("Missing argument 'begin' in Uniform distribution");
      if (args.find("end") == args.end())
        throw Exception("Missing argument 'end' in Uniform distribution");
      rDist = new UniformDiscreteDistribution(nbClasses, TextTools::to<double>(args["begin"]),
                                              TextTools::to<double>(args["end"]));
    }
    else
    {
      throw Exception("Unknown distribution: " + distName + ".");
    }
  }
  if (verbose)
  {
    ApplicationTools::displayResult("Distribution", distName);
    ApplicationTools::displayResult("Number of classes", TextTools::toString(rDist->getNumberOfCategories()));
  }

  return rDist;
}

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
    vector<unsigned int> classes;

    string rf = args["alphas"];
    StringTokenizer strtok(rf.substr(1, rf.length() - 2), ",");
    while (strtok.hasMoreToken())
      alphas.push_back(TextTools::toDouble(strtok.nextToken()));

    rf = args["classes"];
    StringTokenizer strtok2(rf.substr(1, rf.length() - 2), ",");
    while (strtok2.hasMoreToken())
      classes.push_back(TextTools::toInt(strtok2.nextToken()));

    pMDD = new DirichletDiscreteDistribution(classes, alphas);
    vector<string> v = pMDD->getParameters().getParameterNames();

    for (unsigned int i = 0; i < v.size(); i++)
    {
      unparsedParameterValues[v[i]] = TextTools::toString(pMDD->getParameterValue(pMDD->getParameterNameWithoutNamespace(v[i])));
    }
  }
  else
    throw Exception("Unknown multiple distribution name: " + distName);

  return pMDD;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::setRateDistributionParametersInitialValues(
  DiscreteDistribution* rDist,
  map<string, string>& unparsedParameterValues,
  bool verbose) throw (Exception)
{
  ParameterList pl = rDist->getIndependentParameters();
  for (unsigned int i = 0; i < pl.size(); i++)
  {
    AutoParameter ap(pl[i]);
    ap.setMessageHandler(ApplicationTools::warning);
    pl.setParameter(i, ap);
  }

  for (unsigned int i = 0; i < pl.size(); i++)
  {
    const string pName = pl[i].getName();
    double value = ApplicationTools::getDoubleParameter(pName, unparsedParameterValues, pl[i].getValue());
    pl[i].setValue(value);
    if (verbose)
      ApplicationTools::displayResult("Parameter found", pName + "=" + TextTools::toString(pl[i].getValue()));
  }
  rDist->matchParametersValues(pl);
  if (verbose)
  {
    for (unsigned int c = 0; c < rDist->getNumberOfCategories(); c++)
    {
      ApplicationTools::displayResult("- Category " + TextTools::toString(c)
                                      + " (Pr = " + TextTools::toString(rDist->getProbability(c)) + ") rate", TextTools::toString(rDist->getCategory(c)));
    }
  }
}


/******************************************************************************/

DiscreteDistribution* PhylogeneticsApplicationTools::getRateDistribution(
  map<string, string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  string distDescription = ApplicationTools::getStringParameter("rate_distribution", params, "Uniform()", suffix, suffixIsOptional);
  map<string, string> unparsedParameterValues;
  DiscreteDistribution* rDist = getRateDistributionDefaultInstance(distDescription, unparsedParameterValues, verbose);
  setRateDistributionParametersInitialValues(rDist, unparsedParameterValues, verbose);
  return rDist;
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
  bool verbose)
throw (Exception)
{
  string optimization = ApplicationTools::getStringParameter("optimization", params, "FullD(derivatives=Newton)", suffix, suffixIsOptional, false);
  if (optimization == "None")
    return tl;
  string optName;
  map<string, string> optArgs;
  KeyvalTools::parseProcedure(optimization, optName, optArgs);

  unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2, suffix, suffixIsOptional);

  string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional);
  OutputStream* messageHandler =
    (mhPath == "none") ? 0 :
    (mhPath == "std") ? ApplicationTools::message :
    new StlOutputStream(auto_ptr<ostream>(new ofstream(mhPath.c_str(), ios::out)));
  if (verbose)
    ApplicationTools::displayResult("Message handler", mhPath);

  string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional);
  OutputStream* profiler =
    (prPath == "none") ? 0 :
    (prPath == "std") ? ApplicationTools::message :
    new StlOutputStream(auto_ptr<ostream>(new ofstream(prPath.c_str(), ios::out)));
  if (profiler)
    profiler->setPrecision(20);
  if (verbose)
    ApplicationTools::displayResult("Profiler", prPath);

  bool scaleFirst = ApplicationTools::getBooleanParameter("optimization.scale_first", params, false, suffix, suffixIsOptional, false);
  if (scaleFirst)
  {
    // We scale the tree before optimizing each branch length separately:
    if (verbose)
      ApplicationTools::displayMessage("Scaling the tree before optimizing each branch length separately.");
    double tolerance = ApplicationTools::getDoubleParameter("optimization.scale_first.tolerance", params, .0001, suffix, suffixIsOptional, true);
    if (verbose)
      ApplicationTools::displayResult("Scaling tolerance", TextTools::toString(tolerance));
    int nbEvalMax = ApplicationTools::getIntParameter("optimization.scale_first.max_number_f_eval", params, 1000000, suffix, suffixIsOptional, true);
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

  size_t i;
  // Should I ignore some parameters?
  ParameterList parametersToEstimate = parameters;
  string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameter", params, "", suffix, suffixIsOptional, false);
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
      else if ((i = param.find("*")) != string::npos)
      {
        string pref = param.substr(0, i);
        vector<string> vs;
        for (unsigned int j = 0; j < parametersToEstimate.size(); j++)
        {
          if (parametersToEstimate[j].getName().find(pref) == 0)
            vs.push_back(parametersToEstimate[j].getName());
        }
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

  unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000, suffix, suffixIsOptional);
  if (verbose)
    ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));

  double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", params, .000001, suffix, suffixIsOptional);
  if (verbose)
    ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));

  // Backing up or restoring?
  auto_ptr<BackupListener> backupListener;
  string backupFile = ApplicationTools::getAFilePath("optimization.backup.file", params, false, false);
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
          unsigned int p = pl.whichParameterHasName(pname);
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

  // There it goes...
  bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", params, false, suffix, suffixIsOptional, false);
  if (verbose)
    ApplicationTools::displayResult("Optimize topology", optimizeTopo ? "yes" : "no");
  string nniMethod = ApplicationTools::getStringParameter("optimization.topology.algorithm_nni.method", params, "phyml", suffix, suffixIsOptional, false);
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


  string order = ApplicationTools::getStringParameter("derivatives", optArgs, "Newton", "", true, false);
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
  bool reparam = ApplicationTools::getBooleanParameter("optimization.reparametrization", params, false);
  if (verbose)
    ApplicationTools::displayResult("Reparametrization", (reparam ? "yes" : "no"));

  // See if we should use a molecular clock constraint:
  string clock = ApplicationTools::getStringParameter("optimization.clock", params, "None", "", true, false);
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

    unsigned int nstep = ApplicationTools::getParameter<unsigned int>("nstep", optArgs, 1, "", true, false);

    if (optimizeTopo)
    {
      bool optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, false);
      unsigned int topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, "", true, false);
      double tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional);
      double tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, suffixIsOptional);
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
      bool optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, false);
      unsigned int topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, "", true, false);
      double tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional);
      double tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, suffixIsOptional);
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

  string finalMethod = ApplicationTools::getStringParameter("optimization.final", params, "none", suffix, suffixIsOptional, true);
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

  if (verbose)
    ApplicationTools::displayResult("Performed", TextTools::toString(n) + " function evaluations.");
  if (backupFile != "none")
  {
    remove(backupFile.c_str());
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
  bool verbose)
throw (Exception)
{
  string optimization = ApplicationTools::getStringParameter("optimization", params, "FullD(derivatives=Newton)", suffix, suffixIsOptional, false);
  if (optimization == "None")
    return;
  string optName;
  map<string, string> optArgs;
  KeyvalTools::parseProcedure(optimization, optName, optArgs);

  unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2, suffix, suffixIsOptional);

  string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional);
  OutputStream* messageHandler =
    (mhPath == "none") ? 0 :
    (mhPath == "std") ? ApplicationTools::message :
    new StlOutputStream(auto_ptr<ostream>(new ofstream(mhPath.c_str(), ios::out)));
  if (verbose)
    ApplicationTools::displayResult("Message handler", mhPath);

  string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional);
  OutputStream* profiler =
    (prPath == "none") ? 0 :
    (prPath == "std") ? ApplicationTools::message :
    new StlOutputStream(auto_ptr<ostream>(new ofstream(prPath.c_str(), ios::out)));
  if (profiler)
    profiler->setPrecision(20);
  if (verbose)
    ApplicationTools::displayResult("Profiler", prPath);

  ParameterList parametersToEstimate = parameters;

  // Should I ignore some parameters?
  string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameter", params, "", suffix, suffixIsOptional, false);
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

  unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000, suffix, suffixIsOptional);
  if (verbose)
    ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));

  double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", params, .000001, suffix, suffixIsOptional);
  if (verbose)
    ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));

  string order  = ApplicationTools::getStringParameter("derivatives", optArgs, "Gradient", "", true, false);
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
  auto_ptr<BackupListener> backupListener;
  string backupFile = ApplicationTools::getAFilePath("optimization.backup.file", params, false, false);
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
          unsigned int p = pl.whichParameterHasName(pname);
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

  unsigned int n = 0;
  if (optName == "D-Brent")
  {
    // Uses Newton-Brent method:
    unsigned int nstep = ApplicationTools::getParameter<unsigned int>("nstep", optArgs, 1, "", true, false);
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

  string finalMethod = ApplicationTools::getStringParameter("optimization.final", params, "none", suffix, suffixIsOptional, false);
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
    remove(backupFile.c_str());
  }
}

/******************************************************************************/

void PhylogeneticsApplicationTools::checkEstimatedParameters(const ParameterList& pl)
{
  for (unsigned int i = 0; i < pl.size(); ++i)
  {
    const Constraint* constraint = pl[i].getConstraint();
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


/*************************************************************/
/**************  OUTPUT **************************************/
/*************************************************************/

/******************************************************************************/

void PhylogeneticsApplicationTools::writeTree(
  const TreeTemplate<Node>& tree,
  map<string, string>& params,
  const string& prefix,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  bool checkOnly) throw (Exception)
{
  string format = ApplicationTools::getStringParameter(prefix + "tree.format", params, "Newick", suffix, suffixIsOptional, false);
  string file = ApplicationTools::getAFilePath(prefix + "tree.file", params, true, false, suffix, suffixIsOptional);
  OTree* treeWriter;
  if (format == "Newick")
    treeWriter = new Newick();
  else if (format == "Nexus")
    treeWriter = new NexusIOTree();
  else if (format == "NHX")
    treeWriter = new Nhx(false);
  else
    throw Exception("Unknown format for tree writing: " + format);
  if (!checkOnly)
    treeWriter->write(tree, file, true);
  delete treeWriter;
  if (verbose)
    ApplicationTools::displayResult("Wrote tree to file ", file);
}

/******************************************************************************/

void PhylogeneticsApplicationTools::writeTrees(
  const vector<Tree*>& trees,
  map<string, string>& params,
  const string& prefix,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  bool checkOnly) throw (Exception)
{
  string format = ApplicationTools::getStringParameter(prefix + "trees.format", params, "Newick", suffix, suffixIsOptional, false);
  string file = ApplicationTools::getAFilePath(prefix + "trees.file", params, true, false, suffix, suffixIsOptional);
  OMultiTree* treeWriter;
  if (format == "Newick")
    treeWriter = new Newick();
  else if (format == "Nexus")
    treeWriter = new NexusIOTree();
  else if (format == "NHX")
    treeWriter = new Nhx();
  else
    throw Exception("Unknow format for tree writing: " + format);
  if (!checkOnly)
    treeWriter->write(trees, file, true);
  delete treeWriter;
  if (verbose)
    ApplicationTools::displayResult("Wrote trees to file ", file);
}

/******************************************************************************/

void PhylogeneticsApplicationTools::describeParameters_(const ParameterAliasable* parametrizable, OutputStream& out, map<string, string>& globalAliases, const vector<string>& names,  std::vector<std::string>& writtenNames, bool printLocalAliases, bool printComma)
{
  ParameterList pl = parametrizable->getIndependentParameters().subList(names);
  unsigned int p = out.getPrecision();
  out.setPrecision(12);
  bool flag(printComma);
  for (unsigned int i = 0; i < pl.size(); i++)
  {
    if (find(writtenNames.begin(),writtenNames.end(),pl[i].getName())==writtenNames.end()){
      if (flag)
        out << ",";
      else
        flag=true;
      string pname = parametrizable->getParameterNameWithoutNamespace(pl[i].getName());
      
      // Check for global aliases:
      if (globalAliases.find(pl[i].getName()) == globalAliases.end())
        {
          (out << pname << "=").enableScientificNotation(false) << pl[i].getValue();
        }
      else
        out << pname << "=" << globalAliases[pl[i].getName()];
      
      // Now check for local aliases:
      if (printLocalAliases)
        {
          vector<string> aliases = parametrizable->getAlias(pname);
          for (unsigned int j = 0; aliases.size(); j++)
            {
              out << ", " << aliases[j] << "=" << pname;
            }
        }
      writtenNames.push_back(pl[i].getName());
    }
  }
  out.setPrecision(p);
}

/******************************************************************************/

void PhylogeneticsApplicationTools::describeParameters_(const Parametrizable* parametrizable, OutputStream& out, const vector<string>& names,  std::vector<std::string>& writtenNames, bool printComma)
{
  ParameterList pl = parametrizable->getParameters();
  unsigned int p = out.getPrecision();
  out.setPrecision(12);
  bool flag(printComma);
  for (unsigned int i = 0; i < pl.size(); i++)
    {
      if (find(writtenNames.begin(),writtenNames.end(),pl[i].getName())==writtenNames.end()){
        if (flag)
          out << ",";
        else
          flag=true;
        string pname = parametrizable->getParameterNameWithoutNamespace(pl[i].getName());
      
        (out << pname << "=").enableScientificNotation(false) << pl[i].getValue();
      
      }
    }
  out.setPrecision(p);
}

/******************************************************************************/

void PhylogeneticsApplicationTools::describeSubstitutionModel_(const SubstitutionModel* model, OutputStream& out, map<string, string>& globalAliases, std::vector<std::string>& writtenNames)
{
  const BppOSubstitutionModelFormat* bIO=new BppOSubstitutionModelFormat();
  
  bIO->write(*model, out, globalAliases, writtenNames);
  delete bIO;
}

/******************************************************************************/


void PhylogeneticsApplicationTools::describeFrequenciesSet_(const FrequenciesSet* pfreqset, OutputStream& out, std::vector<std::string>& writtenNames)
{
  const BppOFrequenciesSetFormat* bIO=new BppOFrequenciesSetFormat();
  
  bIO->write(pfreqset, out, writtenNames);
  delete bIO;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printParameters(const SubstitutionModel* model, OutputStream& out)
{
  out << "model = ";
  map<string, string> globalAliases;
  vector<string> writtenNames;
  describeSubstitutionModel_(model, out, globalAliases, writtenNames);
  out.endLine();
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printParameters(const SubstitutionModelSet* modelSet, OutputStream& out)
{
  (out << "nonhomogeneous=general").endLine();
  (out << "nonhomogeneous.number_of_models=" << modelSet->getNumberOfModels()).endLine();

  // Get the parameter links:
  map< unsigned int, vector<string> > modelLinks; // for each model index, stores the list of global parameters.
  map< string, set<unsigned int> > parameterLinks; // for each parameter name, stores the list of model indices, wich should be sorted.
  vector<string> writtenNames;
  ParameterList pl = modelSet->getParameters();
  ParameterList plroot = modelSet->getRootFrequenciesParameters();
  for (unsigned int i = 0; i < pl.size(); i++)
  {
    if (!plroot.hasParameter(pl[i].getName()))
    {
      string name = pl[i].getName();
      vector<unsigned int> models = modelSet->getModelsWithParameter(name);
      for (size_t j = 0; j < models.size(); ++j)
      {
        modelLinks[models[j]].push_back(name);
        parameterLinks[name].insert(models[j]);
      }
    }
  }

  // Loop over all models:
  for (unsigned int i = 0; i < modelSet->getNumberOfModels(); i++)
  {
    const SubstitutionModel* model = modelSet->getModel(i);

    // First get the global aliases for this model:
    map<string, string> globalAliases;
    vector<string> names = modelLinks[i];
    for (unsigned int j = 0; j < names.size(); j++)
    {
      const string name = names[j];
      if (parameterLinks[name].size() > 1)
      {
        // there is a global alias here
        if (*parameterLinks[name].begin() != i) // Otherwise, this is the 'reference' value
        {
          globalAliases[modelSet->getParameterModelName(name)] = "model" + TextTools::toString((*parameterLinks[name].begin()) + 1) + "." + modelSet->getParameterModelName(name);
        }
      }
    }

    // Now print it:
    writtenNames.clear();
    out.endLine() << "model" << (i + 1) << "=";
    describeSubstitutionModel_(model, out, globalAliases, writtenNames);
    out.endLine();
    vector<int> ids = modelSet->getNodesWithModel(i);
    out << "model" << (i + 1) << ".nodes_id=" << ids[0];
    for (unsigned int j = 1; j < ids.size(); j++)
    {
      out << "," << ids[j];
    }
    out.endLine();
  }

  // Root frequencies:
  out.endLine();
  (out << "# Root frequencies:").endLine();
  out << "nonhomogeneous.root_freq=";
  describeFrequenciesSet_(modelSet->getRootFrequenciesSet(), out, writtenNames);
}

/******************************************************************************/

void PhylogeneticsApplicationTools::describeDiscreteDistribution_(const DiscreteDistribution* rDist, OutputStream& out, map<string, string>& globalAliases, vector<string>& writtenNames)
{
  const InvariantMixedDiscreteDistribution* invar = dynamic_cast<const InvariantMixedDiscreteDistribution*>(rDist);
  const DiscreteDistribution* test = rDist;
  if (invar)
  {
    test = invar->getVariableSubDistribution();
    out << "Invariant(dist=";
    describeDiscreteDistribution_(test, out, globalAliases, writtenNames);
    out << ", ";
    vector<string> names;
    names.push_back(invar->getParameter("p").getName());
    describeParameters_(invar, out, globalAliases, names, writtenNames);
    out << ")";
  }
  else
  {
    test = dynamic_cast<const ConstantDistribution*>(rDist);
    if (test)
      out << "Uniform()";
    else 
    {
      test = dynamic_cast<const GammaDiscreteDistribution*>(rDist);
      if (test)
      {
        out << "Gamma(n=" << rDist->getNumberOfCategories() << ", ";
        describeParameters_(rDist, out, globalAliases, rDist->getIndependentParameters().getParameterNames(), writtenNames, false);
        out << ")";
      }
      else
        throw Exception("PhylogeneticsApplicationTools::printParameters(DiscreteDistribution). Unsupported distribution.");
    }
  }
}

void PhylogeneticsApplicationTools::printParameters(const DiscreteDistribution* rDist, OutputStream& out)
{
  out << "rate_distribution=";
  map<string, string> globalAliases;
  vector<string> writtenNames;
  describeDiscreteDistribution_(rDist, out, globalAliases, writtenNames);
  out.endLine();
}

/******************************************************************************/

