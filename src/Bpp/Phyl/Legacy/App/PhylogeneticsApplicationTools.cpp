// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../../App/PhylogeneticsApplicationTools.h"
#include "../../Io/BppOBranchModelFormat.h"
#include "../../Io/BppOFrequencySetFormat.h"
#include "../../Io/BppOMultiTreeReaderFormat.h"
#include "../../Io/BppOSubstitutionModelFormat.h"
#include "../../Io/BppOTransitionModelFormat.h"
#include "../../Io/BppOTreeReaderFormat.h"
#include "../../Io/Newick.h"
#include "../../Io/NexusIoTree.h"
#include "../../Io/Nhx.h"
#include "../../Model/FrequencySet/MvaFrequencySet.h"
#include "../../Model/MixedTransitionModel.h"
#include "../../Model/Protein/Coala.h"
#include "../../Model/SubstitutionModel.h"
#include "../../Model/WrappedModel.h"
#include "../../Tree/Tree.h"
#include "../../Tree/TreeTools.h"
#include "../OptimizationTools.h"
#include "PhylogeneticsApplicationTools.h"

// From bpp-core
#include <Bpp/Io/BppODiscreteDistributionFormat.h>
#include <Bpp/Io/BppOParametrizableFormat.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Prob/DirichletDiscreteDistribution.h>
#include <Bpp/Numeric/Function/DownhillSimplexMethod.h>
#include <Bpp/Numeric/Function/PowellMultiDimensions.h>
#include <Bpp/Utils/AttributesTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/SymbolListTools.h>

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


map<size_t, shared_ptr<Tree>> LegacyPhylogeneticsApplicationTools::getTrees(
    const map<string, string>& params,
    const map<size_t, std::shared_ptr<AlignmentDataInterface>>& mSeq,
    map<string, string>& unparsedParams,
    const string& prefix,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  vector<string> vTreesName = ApplicationTools::matchingParameters(prefix + "tree*", params);

  map<size_t, shared_ptr<Tree>> mTree;

  for (size_t nT = 0; nT < vTreesName.size(); nT++)
  {
    size_t poseq = vTreesName[nT].find("=");
    size_t num = 0;
    size_t len = (prefix + "tree").size();
    string suff = vTreesName[nT].substr(len, poseq - len);
    bool flag = 0;
    size_t nbTree = 1;

    if (TextTools::isDecimalInteger(suff, '$'))
      num = static_cast<size_t>(TextTools::toInt(suff));
    else
    {
      flag = 1;
      num = 1;
    }

    if (!flag)
    {
      ApplicationTools::displayMessage("");
      ApplicationTools::displayMessage("Tree " + TextTools::toString(num));
    }

    string treeDesc = ApplicationTools::getStringParameter(vTreesName[nT], params, "", suffix, suffixIsOptional);

    string treeName;

    map<string, string> args;

    KeyvalTools::parseProcedure(treeDesc, treeName, args);

    if (treeName == "user")
    {
      string format;

      if (args.find("format") != args.end())
        format = args["format"];
      else
      {
        format = "Newick";
        ApplicationTools::displayWarning("Warning, " + vTreesName[nT] + " format set to Newick");
      }

      string treeFilePath = ApplicationTools::getAFilePath("file", args, true, true, suffix, suffixIsOptional, "none", warn);

      unique_ptr<IMultiTree> treeReader;
      if (format == "Newick")
        treeReader.reset(new Newick(true));
      else if (format == "Nexus")
        treeReader.reset(new NexusIOTree());
      else if (format == "NHX")
        treeReader.reset(new Nhx());
      else
        throw Exception("Unknown format for tree reading: " + format);

      vector<unique_ptr<Tree>> trees;
      treeReader->readTrees(treeFilePath, trees);

      if (verbose)
      {
        if (flag)
        {
          ApplicationTools::displayMessage("");
          ApplicationTools::displayResult("Tree file", treeFilePath);
        }

        ApplicationTools::displayResult("Number of trees in file", trees.size());
      }

      if (flag)
      {
        nbTree = trees.size();

        for (size_t i2 = 0; i2 < trees.size(); i2++)
        {
          if (mTree.find(i2 + 1) != mTree.end())
          {
            ApplicationTools::displayWarning("Tree " + TextTools::toString(i2 + 1) + " already assigned, replaced by new one.");
          }

          mTree[i2 + 1] = std::move(trees[i2]);
          ApplicationTools::displayResult("Number of leaves", trees[i2]->getNumberOfLeaves());
        }
      }
      else
      {
        if (trees.size() > 1)
          throw Exception("Error : Several trees for description of " + vTreesName[nT] + ".");

        if (trees.size() == 1)
        {
          if (mTree.find(num) != mTree.end())
          {
            ApplicationTools::displayWarning("Tree " + TextTools::toString(num) + " already assigned, replaced by new one.");
          }
          mTree[num] = std::move(trees[0]);
          ApplicationTools::displayResult("Number of leaves", trees[0]->getNumberOfLeaves());
        }
      }
    }
    else if (treeName == "random")
    {
      size_t seqNum;

      if (args.find("data") == args.end())
      {
        ApplicationTools::displayWarning("Random tree set from data 1");
        seqNum = 1;
      }
      else
        seqNum = (size_t) TextTools::toInt(args["data"]);


      if (mSeq.find(seqNum) == mSeq.end())
        throw Exception("Error : Wrong number of data " + TextTools::toString(seqNum));

      vector<string> names = mSeq.find(seqNum)->second->getSequenceNames();
      auto tree = TreeTemplateTools::getRandomTree(names);
      tree->setBranchLengths(1.);

      if (mTree.find(num) != mTree.end())
      {
        ApplicationTools::displayWarning("Tree " + TextTools::toString(num) + " already assigned, replaced by new one.");
      }
      mTree[num] = std::move(tree);
      ApplicationTools::displayResult("Number of leaves", tree->getNumberOfLeaves());
    }

    // //////////
    // Setting branch lengths?
    string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", args, "Input", "", true, 1);
    string cmdName;
    map<string, string> cmdArgs;

    KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs);
    if (cmdName == "Input")
    {
      // Is the root has to be moved to the midpoint position along the branch that contains it ? If no, do nothing!
      string midPointRootBrLengths = ApplicationTools::getStringParameter("midPointRootBrLengths", cmdArgs, "no", "", true, 2);
      if (midPointRootBrLengths == "yes")
      {
        if (flag)
        {
          for (size_t i = 0; i < nbTree; i++)
          {
            TreeTools::constrainedMidPointRooting(*mTree[i + 1]);
          }
        }
        else
          TreeTools::constrainedMidPointRooting(*mTree[num]);
      }
    }
    else if (cmdName == "Equal")
    {
      double value = ApplicationTools::getDoubleParameter("value", cmdArgs, 0.1, "", true, 2);
      if (value <= 0)
        throw Exception("Value for branch length must be superior to 0");
      ApplicationTools::displayResult("Branch lengths set to", value);
      if (flag)
      {
        for (size_t i = 0; i < nbTree; i++)
        {
          mTree[i + 1]->setBranchLengths(value);
        }
      }
      else
        mTree[num]->setBranchLengths(value);
    }
    else if (cmdName == "Clock")
    {
      if (flag)
      {
        for (size_t i = 0; i < nbTree; i++)
        {
          TreeTools::convertToClockTree(*mTree[i + 1], mTree[i + 1]->getRootId(), true);
        }
      }
      else
        TreeTools::convertToClockTree(*mTree[num], mTree[num]->getRootId(), true);
    }
    else if (cmdName == "Grafen")
    {
      string grafenHeight = ApplicationTools::getStringParameter("height", cmdArgs, "input", "", true, 2);
      double h;
      if (flag)
      {
        for (size_t i = 0; i < nbTree; i++)
        {
          auto& tree = *mTree[i + 1];
          if (grafenHeight == "input")
          {
            h = TreeTools::getHeight(tree, tree.getRootId());
          }
          else
          {
            h = TextTools::toDouble(grafenHeight);
            if (h <= 0)
              throw Exception("Height must be positive in Grafen's method.");
          }
          ApplicationTools::displayResult("Total height", TextTools::toString(h));

          double rho = ApplicationTools::getDoubleParameter("rho", cmdArgs, 1., "", true, 2);
          ApplicationTools::displayResult("Grafen's rho", rho);
          TreeTools::computeBranchLengthsGrafen(tree, rho);
          double nh = TreeTools::getHeight(tree, tree.getRootId());
          tree.scaleTree(h / nh);
        }
      }
      else
      {
        auto& tree = *mTree[num];
        if (grafenHeight == "input")
        {
          h = TreeTools::getHeight(tree, tree.getRootId());
        }
        else
        {
          h = TextTools::toDouble(grafenHeight);
          if (h <= 0)
            throw Exception("Height must be positive in Grafen's method.");
        }
        ApplicationTools::displayResult("Total height", TextTools::toString(h));

        double rho = ApplicationTools::getDoubleParameter("rho", cmdArgs, 1., "", true, 2);
        ApplicationTools::displayResult("Grafen's rho", rho);
        TreeTools::computeBranchLengthsGrafen(tree, rho);
        double nh = TreeTools::getHeight(tree, tree.getRootId());
        tree.scaleTree(h / nh);
      }
    }
    else
      throw Exception("Method '" + initBrLenMethod + "' unknown for computing branch lengths.");

    // //////////// Setting branch lengths with aliases

    vector<string> vBrNb = ApplicationTools::matchingParameters("BrLen*", args);

    for (size_t ib = 0; ib < vBrNb.size(); ib++)
    {
      string apeq = args[vBrNb[ib]];
      string aveq = vBrNb[ib];

      if (TextTools::isDecimalInteger(apeq))
        mTree[num]->setDistanceToFather(TextTools::toInt(aveq.substr(5, string::npos)), TextTools::toDouble(apeq));
      else
      {
        size_t posun = apeq.find("_");
        size_t posd = aveq.find("_");
        unparsedParams[aveq + (posd != string::npos ? "" : "_" + TextTools::toString(num))] = apeq + (posun != string::npos ? "" : "_" + TextTools::toString(num));
      }
    }

    ApplicationTools::displayResult("Branch lengths", cmdName);
  }

  return mTree;
}


/******************************************************/
/**** SUBSTITUTION MODEL SET **************************/
/******************************************************/

/******************************************************************************/

void LegacyPhylogeneticsApplicationTools::setSubstitutionModelParametersInitialValuesWithAliases(
    BranchModelInterface& model,
    map<string, string>& unparsedParameterValues,
    size_t modelNumber,
    std::shared_ptr<const AlignmentDataInterface> data,
    map<string, string>& sharedParams,
    bool verbose)
{
  string initFreqs = ApplicationTools::getStringParameter(model.getNamespace() + "initFreqs", unparsedParameterValues, "", "", true, 2);

  if (verbose)
    ApplicationTools::displayResult("Frequencies Initialization for model", (initFreqs == "") ? "None" : initFreqs);

  if (initFreqs != "")
  {
    try
    {
      auto& tmodel = dynamic_cast<TransitionModelInterface&>(model);
      if (initFreqs == "observed")
      {
        if (!data)
          throw Exception("Missing data for observed frequencies");
        unsigned int psi = ApplicationTools::getParameter<unsigned int>(model.getNamespace() + "initFreqs.observedPseudoCount", unparsedParameterValues, 0);
        tmodel.setFreqFromData(*data, psi);
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
        tmodel.setFreq(frequencies);
      }
      else
        throw Exception("Unknown initFreqs argument");
    }
    catch (exception& e)
    {
      ApplicationTools::displayMessage("Frequencies initialization not possible for model " + model.getName());
    }
  }

  ParameterList pl = model.getIndependentParameters();
  for (size_t i = 0; i < pl.size(); ++i)
  {
    AutoParameter ap(pl[i]);
    ap.setMessageHandler(ApplicationTools::warning);
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

/******************************************************************************/

unique_ptr<SubstitutionModelSet> LegacyPhylogeneticsApplicationTools::getSubstitutionModelSet(
    shared_ptr<const Alphabet> alphabet,
    shared_ptr<const GeneticCode> gCode,
    shared_ptr<const AlignmentDataInterface> data,
    const map<string, string>& params,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  if (!ApplicationTools::parameterExists("nonhomogeneous.number_of_models", params))
    throw Exception("A value is needed for this parameter: nonhomogeneous.number_of_models .");
  size_t nbModels = ApplicationTools::getParameter<size_t>("nonhomogeneous.number_of_models", params, 1, suffix, suffixIsOptional, warn);
  if (nbModels == 0)
    throw Exception("The number of models can't be 0 !");

  bool nomix = true;
  for (size_t i = 0; nomix& (i < nbModels); i++)
  {
    string prefix = "model" + TextTools::toString(i + 1);
    string modelDesc;
    modelDesc = ApplicationTools::getStringParameter(prefix, params, "", suffix, suffixIsOptional, warn);

    if (modelDesc.find("Mixed") != string::npos)
      nomix = false;
  }

  unique_ptr<SubstitutionModelSet> modelSet = nullptr;
  auto modelSet1 = make_unique<SubstitutionModelSet>(alphabet);
  setSubstitutionModelSet(*modelSet1, alphabet, gCode, data, params, suffix, suffixIsOptional, verbose, warn);

  if (modelSet1->hasMixedTransitionModel())
  {
    modelSet = make_unique<MixedSubstitutionModelSet>(*modelSet1);
    completeMixedSubstitutionModelSet(dynamic_cast<MixedSubstitutionModelSet&>(*modelSet), alphabet, data, params, suffix, suffixIsOptional, verbose);
  }
  else
    modelSet = std::move(modelSet1);

  return modelSet;
}

/******************************************************************************/

void LegacyPhylogeneticsApplicationTools::setSubstitutionModelSet(
    SubstitutionModelSet& modelSet,
    shared_ptr<const Alphabet> alphabet,
    shared_ptr<const GeneticCode> gCode,
    shared_ptr<const AlignmentDataInterface> data,
    const map<string, string>& params,
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

  BppOTransitionModelFormat bIO(BppOTransitionModelFormat::ALL, true, true, true, false, warn);

  // ///////////////////////////////////////////
  // Build a new model set object:

  vector<double> rateFreqs;
  string tmpDesc;
  if (AlphabetTools::isCodonAlphabet(alphabet.get()))
  {
    if (!gCode)
      throw Exception("PhylogeneticsApplicationTools::setSubstitutionModelSet(): a GeneticCode instance is required for instantiating a codon model.");
    bIO.setGeneticCode(gCode);
    tmpDesc = ApplicationTools::getStringParameter("model1", params, "CodonRate(model=JC69)", suffix, suffixIsOptional, warn);
  }
  else if (AlphabetTools::isWordAlphabet(alphabet.get()))
    tmpDesc = ApplicationTools::getStringParameter("model1", params, "Word(model=JC69)", suffix, suffixIsOptional, warn);
  else
    tmpDesc = ApplicationTools::getStringParameter("model1", params, "JC69", suffix, suffixIsOptional, warn);

  shared_ptr<TransitionModelInterface> tmp = bIO.readTransitionModel(alphabet, tmpDesc, *data, false);

  if (tmp->getNumberOfStates() != alphabet->getSize())
  {
    // Markov-Modulated Markov Model...

    size_t n = static_cast<size_t>(tmp->getNumberOfStates() / alphabet->getSize());
    rateFreqs = vector<double>(n, 1. / static_cast<double>(n));
    // Equal rates assumed for now, may be changed later
  }

  // ////////////////////////////////////
  // Deal with root frequencies

  map<string, string> unparsedParameters;

  bool stationarity = ApplicationTools::getBooleanParameter("nonhomogeneous.stationarity", params, false, "", true, warn);
  std::shared_ptr<FrequencySetInterface> rootFrequencies(0);
  if (!stationarity)
  {
    rootFrequencies = PhylogeneticsApplicationTools::getRootFrequencySet(alphabet, gCode, *data, params, unparsedParameters, rateFreqs, suffix, suffixIsOptional, verbose);
    stationarity = !rootFrequencies;
    string freqDescription = ApplicationTools::getStringParameter("nonhomogeneous.root_freq", params, "", suffix, suffixIsOptional, warn);
    if (freqDescription.substr(0, 10) == "MVAprotein")
    {
      auto tmp2 = dynamic_pointer_cast<Coala>(tmp);
      if (tmp2)
        dynamic_pointer_cast<MvaFrequencySet>(rootFrequencies)->initSet(*tmp2);
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

  for (size_t i = 0; i < nbModels; ++i)
  {
    string prefix = "model" + TextTools::toString(i + 1);
    string modelDesc;
    if (AlphabetTools::isCodonAlphabet(alphabet.get()))
      modelDesc = ApplicationTools::getStringParameter(prefix, params, "CodonRate(model=JC69)", suffix, suffixIsOptional, warn);
    else if (AlphabetTools::isWordAlphabet(alphabet.get()))
      modelDesc = ApplicationTools::getStringParameter(prefix, params, "Word(model=JC69)", suffix, suffixIsOptional, warn);
    else
      modelDesc = ApplicationTools::getStringParameter(prefix, params, "JC69", suffix, suffixIsOptional, warn);

    auto model = bIO.readTransitionModel(alphabet, modelDesc, *data, false);

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

    modelSet.addModel(std::move(model), nodesId);
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

void LegacyPhylogeneticsApplicationTools::completeMixedSubstitutionModelSet(
    MixedSubstitutionModelSet& mixedModelSet,
    shared_ptr<const Alphabet> alphabet,
    shared_ptr<const AlignmentDataInterface> data,
    const map<string, string>& params,
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
  while (numd)
  {
    string desc = ApplicationTools::getStringParameter("site.path" + TextTools::toString(numd), params, "",  suffix, suffixIsOptional, warn);
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

      auto pSM = dynamic_pointer_cast<const MixedTransitionModelInterface>(mixedModelSet.getModel(static_cast<size_t>(num - 1)));
      if (pSM == NULL)
        throw BadIntegerException("LegacyPhylogeneticsApplicationTools::setMixedSubstitutionModelSet: Wrong model for number", num - 1);
      Vuint submodnb = pSM->getSubmodelNumbers(p2);

      mixedModelSet.addToHyperNode(static_cast<size_t>(num - 1), submodnb);
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


/*************************************************************/
/*****  OPTIMIZATORS *****************************************/
/*************************************************************/

/******************************************************************************/

shared_ptr<TreeLikelihoodInterface> LegacyPhylogeneticsApplicationTools::optimizeParameters(
    shared_ptr<TreeLikelihoodInterface> tl,
    const ParameterList& parameters,
    const map<string, string>& params,
    const string& suffix,
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

  unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2, suffix, suffixIsOptional, warn);

  string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional);
  shared_ptr<OutputStream> messageHandler =
      (mhPath == "none") ? nullptr :
      (mhPath == "std") ? ApplicationTools::message :
      make_shared<StlOutputStream>(make_unique<ofstream>(mhPath.c_str(), ios::out));
  if (verbose)
    ApplicationTools::displayResult("Message handler", mhPath);

  string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional);
  shared_ptr<OutputStream> profiler =
      (prPath == "none") ? nullptr :
      (prPath == "std") ? ApplicationTools::message :
      make_shared<StlOutputStream>(make_unique<ofstream>(prPath.c_str(), ios::out));
  if (profiler)
    profiler->setPrecision(20);
  if (verbose)
    ApplicationTools::displayResult("Profiler", prPath);

  bool scaleFirst = ApplicationTools::getBooleanParameter("optimization.scale_first", params, false, suffix, suffixIsOptional, warn);
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
    LegacyOptimizationTools::optimizeTreeScale(
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

  string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameter", params, "", suffix, suffixIsOptional, warn);
  if (paramListDesc.length() == 0)
    paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameters", params, "", suffix, suffixIsOptional, warn);
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
        auto nhtl = dynamic_pointer_cast<NonHomogeneousTreeLikelihood>(tl);
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
        auto nhtl = dynamic_pointer_cast<NonHomogeneousTreeLikelihood>(tl);
        if (nhtl)
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
      else if (param == "*")
      {
        parametersToEstimate.reset();
        if (verbose)
          ApplicationTools::displayResult("Parameter ignored", string("All"));
      }
      else if (param.find("*") != string::npos)
      {
        vector<string> vs = ApplicationTools::matchingParameters(param, parNames);

        bool verbhere = verbose;

        if (vs.size() >= 20)
        {
          if (verbose)
            ApplicationTools::displayResult("Number of parameters ignored", vs.size());
          verbhere = false;
        }

        for (auto& it :  vs)
        {
          parametersToEstimate.deleteParameter(it);
          if (verbhere)
            ApplicationTools::displayResult("Parameter ignored", it);
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
      ApplicationTools::displayWarning("Parameter '" + pnfe.parameter() + "' not found, and so can't be ignored!");
    }
  }

  // Should I constrain some parameters?
  vector<string> parToEstNames = parametersToEstimate.getParameterNames();

  paramListDesc = ApplicationTools::getStringParameter("optimization.constrain_parameter", params, "", suffix, suffixIsOptional, warn);
  if (paramListDesc.length() == 0)
    paramListDesc = ApplicationTools::getStringParameter("optimization.constrain_parameters", params, "", suffix, suffixIsOptional, warn);

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
        throw Exception("PhylogeneticsApplicationTools::optimizeParameters. Bad constrain syntax, should contain `=' symbol: " + pc);
      param = pc.substr(0, index);
      constraint = pc.substr(index + 1);
      std::shared_ptr<IntervalConstraint> ic(new IntervalConstraint(constraint));

      vector<string> parNames2;

      if (param == "BrLen")
        parNames2  = tl->getBranchLengthsParameters().getParameterNames();
      else if (param == "Ancient")
      {
        auto nhtl = dynamic_pointer_cast<NonHomogeneousTreeLikelihood>(tl);
        if (!nhtl)
          ApplicationTools::displayWarning("The 'Ancient' parameters do not exist in homogeneous models, and will be ignored.");
        else
          parNames2 = nhtl->getRootFrequenciesParameters().getParameterNames();
      }
      else if (param == "Model")
      {
        vector<string> vs1 = tl->getSubstitutionModelParameters().getParameterNames();
        auto nhtl = dynamic_pointer_cast<NonHomogeneousTreeLikelihood>(tl);
        if (nhtl)
        {
          vector<string> vs2 = nhtl->getRootFrequenciesParameters().getParameterNames();
          VectorTools::diff(vs1, vs2, parNames2);
        }
        else
          parNames2 = vs1;
      }
      else if (param == "*")
        parNames2 = parToEstNames;
      else if (param.find("*") != string::npos)
        parNames2 = ApplicationTools::matchingParameters(param, parToEstNames);
      else
        parNames2.push_back(param);

      for (size_t i = 0; i < parNames2.size(); ++i)
      {
        Parameter& par = parametersToEstimate.parameter(parNames2[i]);
        if (par.hasConstraint())
        {
          par.setConstraint(std::shared_ptr<ConstraintInterface>(*ic & (*par.getConstraint())));
          if (par.constraint().isEmpty())
            throw Exception("Empty interval for parameter " + parNames[i] + par.constraint().getDescription());
        }
        else
          par.setConstraint(ic);

        if (verbose)
          ApplicationTools::displayResult("Parameter constrained " + par.getName(), par.constraint().getDescription());
      }
    }
    catch (ParameterNotFoundException& pnfe)
    {
      ApplicationTools::displayWarning("Parameter '" + pnfe.parameter() + "' not found, and so can't be constrained!");
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
  shared_ptr<BackupListener> backupListener;
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
          if (pl.hasParameter(pname))
          {
            size_t p = pl.whichParameterHasName(pname);
            pl.setParameter(p, AutoParameter(pl[p]));
            pl[p].setValue(TextTools::toDouble(pvalue));
          }
          else
            ApplicationTools::displayMessage("Warning: unknown parameter in backup file : " + pname);
        }
      }
      bck.close();
      tl->setParameters(pl);
      if (convert(abs(tl->getValue() - fval)) > 0.000001)
        ApplicationTools::displayWarning("Warning, incorrect likelihood value after restoring from backup file.");
      ApplicationTools::displayResult("Restoring log-likelihood", -tl->getValue());
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
  bool reparam = ApplicationTools::getBooleanParameter("optimization.reparametrization", params, warn + 1);
  if (verbose)
    ApplicationTools::displayResult("Reparametrization", (reparam ? "yes" : "no"));

  // See if we should use a molecular clock constraint:
  string clock = ApplicationTools::getStringParameter("optimization.clock", params, "None", "", true, warn + 1);
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
      unsigned int topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, "", true, warn + 1);
      double tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional, warn + 1);
      double tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, warn + 1);
      tl = LegacyOptimizationTools::optimizeTreeNNI(
            dynamic_pointer_cast<NNIHomogeneousTreeLikelihood>(tl), parametersToEstimate,
            optNumFirst, tolBefore, tolDuring, nbEvalMax, topoNbStep, messageHandler, profiler,
            reparam, optVerbose, optMethodDeriv, nstep, nniAlgo);
    }

    if (verbose && nstep > 1)
      ApplicationTools::displayResult("# of precision steps", TextTools::toString(nstep));
    parametersToEstimate.matchParametersValues(tl->getParameters());
    n = LegacyOptimizationTools::optimizeNumericalParameters(
          dynamic_pointer_cast<DiscreteRatesAcrossSitesTreeLikelihoodInterface>(tl), parametersToEstimate,
          backupListener, nstep, tolerance, nbEvalMax, messageHandler, profiler, reparam, optVerbose, optMethodDeriv, optMethodModel);
  }
  else if (optName == "FullD")
  {
    // Uses Newton-raphson algorithm with numerical derivatives when required.

    if (optimizeTopo)
    {
      bool optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, warn + 1);
      unsigned int topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, "", true, warn + 1);
      double tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional, warn + 1);
      double tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, suffixIsOptional, warn + 1);
      tl = LegacyOptimizationTools::optimizeTreeNNI2(
            dynamic_pointer_cast<NNIHomogeneousTreeLikelihood>(tl), parametersToEstimate,
            optNumFirst, tolBefore, tolDuring, nbEvalMax, topoNbStep, messageHandler, profiler,
            reparam, optVerbose, optMethodDeriv, nniAlgo);
    }

    parametersToEstimate.matchParametersValues(tl->getParameters());
    n = LegacyOptimizationTools::optimizeNumericalParameters2(
          dynamic_pointer_cast<DiscreteRatesAcrossSitesTreeLikelihoodInterface>(tl), parametersToEstimate,
          backupListener, tolerance, nbEvalMax, messageHandler, profiler, reparam, useClock, optVerbose, optMethodDeriv);
  }
  else
    throw Exception("Unknown optimization method: " + optName);

  string finalMethod = ApplicationTools::getStringParameter("optimization.final", params, "none", suffix, suffixIsOptional, warn);
  shared_ptr<OptimizerInterface> finalOptimizer = nullptr;
  if (finalMethod == "none")
  {}
  else if (finalMethod == "simplex")
  {
    finalOptimizer = make_shared<DownhillSimplexMethod>(tl);
  }
  else if (finalMethod == "powell")
  {
    finalOptimizer = make_shared<PowellMultiDimensions>(tl);
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
  }

  if (verbose)
    ApplicationTools::displayResult("Performed", TextTools::toString(n) + " function evaluations.");
  if (backupFile != "none")
  {
    string bf = backupFile + ".def";
    rename(backupFile.c_str(), bf.c_str());
  }
  return tl;
}


/******************************************************************************/
/**************** Output ************************************/
/******************************************************************************/

void LegacyPhylogeneticsApplicationTools::writeTrees(
    const vector<const Tree*>& trees,
    const map<string, string>& params,
    const string& prefix,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    bool checkOnly,
    int warn)
{
  string format = ApplicationTools::getStringParameter(prefix + "tree.format", params, "Newick", suffix, suffixIsOptional, warn);
  string file = ApplicationTools::getAFilePath(prefix + "tree.file", params, true, false, suffix, suffixIsOptional, "none", warn);
  OMultiTree* treeWriter;
  if (format == "Newick")
    treeWriter = new Newick();
  else if (format == "Nexus")
    treeWriter = new NexusIOTree();
  else if (format == "NHX")
    treeWriter = new Nhx();
  else
    throw Exception("Unknown format for tree writing: " + format);

  if (!checkOnly)
    treeWriter->writeTrees(trees, file, true);

  delete treeWriter;
  if (verbose)
    ApplicationTools::displayResult("Wrote trees to file ", file);
}


/******************************************************************************/

void LegacyPhylogeneticsApplicationTools::printParameters(const SubstitutionModelSet* modelSet, OutputStream& out, int warn, bool withAlias)
{
  (out << "nonhomogeneous=general").endLine();
  (out << "nonhomogeneous.number_of_models=" << modelSet->getNumberOfModels()).endLine();

  if (modelSet->isStationary())
    (out << "nonhomogeneous.stationarity = yes");

  // Get the parameter links:
  map< size_t, vector<string>> modelLinks; // for each model index, stores the list of global parameters.
  map< string, set<size_t>> parameterLinks; // for each parameter name, stores the list of model indices, which should be sorted.
  vector<string> writtenNames;

  // Loop over all models:
  for (size_t i = 0; i < modelSet->getNumberOfModels(); ++i)
  {
    shared_ptr<const BranchModelInterface> model = modelSet->getModel(i);

    // First get the aliases for this model:

    map<string, string> aliases;
    ParameterList pl = model->getParameters();

    if (withAlias)
    {
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
    BppOBranchModelFormat bIOsm(BppOSubstitutionModelFormat::ALL, true, true, true, false, warn);
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
    const auto pFS = modelSet->getRootFrequencySet();

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

    BppOFrequencySetFormat bIO(BppOFrequencySetFormat::ALL, false, warn);
    bIO.writeFrequencySet(*pFS, out, aliases, writtenNames);
  }
}
