// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../Io/BppOBranchModelFormat.h"
#include "../Io/BppOFrequencySetFormat.h"
#include "../Io/BppOMultiTreeReaderFormat.h"
#include "../Io/BppOMultiTreeWriterFormat.h"
#include "../Io/BppORateDistributionFormat.h"
#include "../Io/BppOTreeReaderFormat.h"
#include "../Io/BppOTreeWriterFormat.h"
#include "../Io/Newick.h"
#include "../Io/NexusIoTree.h"
#include "../Io/Nhx.h"
#include "../Likelihood/AutoCorrelationSequenceEvolution.h"
#include "../Likelihood/HmmSequenceEvolution.h"
#include "../Likelihood/MixtureSequenceEvolution.h"
#include "../Likelihood/NonHomogeneousSubstitutionProcess.h"
#include "../Likelihood/OneProcessSequenceEvolution.h"
#include "../Likelihood/ParametrizablePhyloTree.h"
#include "../Likelihood/PartitionSequenceEvolution.h"
#include "../Likelihood/PhyloLikelihoods/AlignedPhyloLikelihoodAutoCorrelation.h"
#include "../Likelihood/PhyloLikelihoods/AutoCorrelationProcessPhyloLikelihood.h"
#include "../Likelihood/PhyloLikelihoods/PhyloLikelihoodFormula.h"
#include "../Likelihood/PhyloLikelihoods/AlignedPhyloLikelihoodHmm.h"
#include "../Likelihood/PhyloLikelihoods/HmmProcessPhyloLikelihood.h"
#include "../Likelihood/PhyloLikelihoods/AlignedPhyloLikelihoodMixture.h"
#include "../Likelihood/PhyloLikelihoods/MixtureProcessPhyloLikelihood.h"
#include "../Likelihood/PhyloLikelihoods/OneProcessSequencePhyloLikelihood.h"
#include "../Likelihood/PhyloLikelihoods/PartitionProcessPhyloLikelihood.h"
#include "../Likelihood/PhyloLikelihoods/AlignedPhyloLikelihoodProduct.h"
#include "../Likelihood/PhyloLikelihoods/SingleDataPhyloLikelihood.h"
#include "../Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h"
#include "../Likelihood/RateAcrossSitesSubstitutionProcess.h"
#include "../Likelihood/SimpleSubstitutionProcess.h"
#include "../Likelihood/SubstitutionProcessCollection.h"
#include "../Likelihood/SubstitutionProcessCollectionMember.h"
#include "../Mapping/DecompositionSubstitutionCount.h"
#include "../Mapping/LaplaceSubstitutionCount.h"
#include "../Mapping/NaiveSubstitutionCount.h"
#include "../Mapping/OneJumpSubstitutionCount.h"
#include "../Mapping/UniformizationSubstitutionCount.h"
#include "../Model/FrequencySet/MvaFrequencySet.h"
#include "../Model/MixedTransitionModel.h"
#include "../Model/Protein/Coala.h"
#include "../Model/SubstitutionModel.h"
#include "../Model/WrappedModel.h"
#include "../Model/RateDistribution/ConstantRateDistribution.h"
#include "../OptimizationTools.h"
#include "../Tree/PhyloTree.h"
#include "../Tree/PhyloTreeTools.h"
#include "PhylogeneticsApplicationTools.h"

// From bpp-core
#include <Bpp/BppString.h>
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
#include <algorithm>

using namespace std;

/*************************************************************/
/*****************  TREES ************************************/
/*************************************************************/


/******************************************************************************/

unique_ptr<Tree> PhylogeneticsApplicationTools::getTree(
    const map<string, string>& params,
    const string& prefix,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  string format = ApplicationTools::getStringParameter(prefix + "tree.format", params, "Newick", suffix, suffixIsOptional, warn);
  string treeFilePath = ApplicationTools::getAFilePath(prefix + "tree.file", params, true, true, suffix, suffixIsOptional, "none", warn);

  BppOTreeReaderFormat bppoReader(warn);
  auto iTree = bppoReader.readITree(format);
  if (verbose)
  {
    ApplicationTools::displayResult("Input tree file " + suffix, treeFilePath);
    ApplicationTools::displayResult("Input tree format " + suffix, iTree->getFormatName());
  }
  return iTree->readTree(treeFilePath);
}

/******************************************************************************/

vector<unique_ptr<Tree>> PhylogeneticsApplicationTools::getTrees(
    const map<string, string>& params,
    const string& prefix,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  string format = ApplicationTools::getStringParameter(prefix + "tree.format", params, "Newick", suffix, suffixIsOptional, warn);
  string treeFilePath = ApplicationTools::getAFilePath(prefix + "tree.file", params, true, true, suffix, suffixIsOptional, "none", warn);

  BppOMultiTreeReaderFormat bppoReader(warn);
  auto iTrees = bppoReader.readIMultiTree(format);
  if (verbose)
  {
    ApplicationTools::displayResult("Input trees file " + suffix, treeFilePath);
    ApplicationTools::displayResult("Input trees format " + suffix, iTrees->getFormatName());
  }
  vector<unique_ptr<Tree>> trees;
  iTrees->readTrees(treeFilePath, trees);

  if (verbose)
  {
    ApplicationTools::displayResult("Number of trees in file", trees.size());
  }
  return trees;
}

/******************************************************************************/

map<size_t, std::shared_ptr<PhyloTree>> PhylogeneticsApplicationTools::getPhyloTrees(
    const map<string, string>& params,
    const map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mSeq,
    map<string, string>& unparsedParams,
    const string& prefix,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  
  vector<string> vTreesName = ApplicationTools::matchingParameters(prefix + "tree*", params);

  map<size_t, shared_ptr<PhyloTree>> mTree;

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

      IMultiPhyloTree* treeReader;
      if (format == "Newick")
        treeReader = new Newick(true);
      else if (format == "Nexus")
        treeReader = new NexusIOTree();
      else if (format == "NHX")
        treeReader = new Nhx();
      else
        throw Exception("Unknown format for tree reading: " + format);

      vector<unique_ptr<PhyloTree>> trees;
      treeReader->readPhyloTrees(treeFilePath, trees);

      if (args.find("data") != args.end())
      {
        auto seqNum = (size_t) TextTools::toInt(args["data"]);
        if (mSeq.find(seqNum) == mSeq.end())
          throw Exception("Error : Wrong number of data " + TextTools::toString(seqNum));
        else
        {
          ApplicationTools::displayMessage("Tree leaves pruned to fit data " + TextTools::toString(seqNum));
          vector<string> names = mSeq.find(seqNum)->second->getSequenceNames();

          for (auto& tree:trees)
          {
            auto nb1 = tree->getNumberOfLeaves();
            tree->pruneTree(names);
            auto nb2 = tree->getNumberOfLeaves();
            if (nb1 != nb2)
            {
              ApplicationTools::displayResult("Number of removed leaves", nb1 - nb2);
            }
          }
        }
      }

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
            mTree.erase(i2 + 1);
          }

          mTree[i2 + 1] = std::move(trees[i2]);
          ApplicationTools::displayResult("Number of leaves", mTree[i2 + 1]->getNumberOfLeaves());
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
            mTree.erase(num);
          }
          mTree[num] = std::move(trees[0]);
          ApplicationTools::displayResult("Number of leaves", mTree[num]->getNumberOfLeaves());
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

      // Not optimal process: make random PhlyoTree directly
      auto treetemp = TreeTemplateTools::getRandomTree(names);
      treetemp->setBranchLengths(1.);
      
      auto tree = PhyloTreeTools::buildFromTreeTemplate(*treetemp);

      if (mTree.find(num) != mTree.end())
      {
        ApplicationTools::displayWarning("Tree " + TextTools::toString(num) + " already assigned, replaced by new one.");
        mTree.erase(num);
      }
      mTree[num] = tree;
      ApplicationTools::displayResult("Number of leaves", tree->getNumberOfLeaves());
    }


    // //////////
    // Setting branch lengths?
    string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", args, "Input", "", true, 1);
    string cmdName;
    map<string, string> cmdArgs;

    KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs);

    ApplicationTools::displayResult("Branch lengths", cmdName);

    if (cmdName == "Input")
    {
      // Is the root has to be moved to the midpoint position along the branch that contains it ? If no, do nothing!
      string midPointRootBrLengths = ApplicationTools::getStringParameter("midPointRootBrLengths", cmdArgs, "no", "", true, 2);
      if (midPointRootBrLengths == "yes")
      {
        ApplicationTools::displayResult(" Mid Point Rooting", midPointRootBrLengths);
        if (flag)
        {
          for (size_t i = 0; i < nbTree; i++)
          {
            PhyloTreeTools::constrainedMidPointRooting(*mTree[i + 1]);
          }
        }
        else
          PhyloTreeTools::constrainedMidPointRooting(*mTree[num]);
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
          PhyloTreeTools::convertToClockTree(*mTree[i + 1], mTree[i + 1]->getRoot());
        }
      }
      else
        PhyloTreeTools::convertToClockTree(*mTree[num], mTree[num]->getRoot());
    }
    else if (cmdName == "Grafen")
    {
      string grafenHeight = ApplicationTools::getStringParameter("height", cmdArgs, "input", "", true, 2);
      double h;
      if (flag)
      {
        for (size_t i = 0; i < nbTree; i++)
        {
          auto tree = mTree[i + 1];
          if (grafenHeight == "input")
          {
            h = PhyloTreeTools::getHeight(*tree, tree->getRoot());
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
          PhyloTreeTools::computeBranchLengthsGrafen(*tree, rho);

          double nh = PhyloTreeTools::getHeight(*tree, tree->getRoot());
          tree->scaleTree(h / nh);
        }
      }
      else
      {
        auto tree = mTree[num];
        if (grafenHeight == "input")
          h = PhyloTreeTools::getHeight(*tree, tree->getRoot());
        else
        {
          h = TextTools::toDouble(grafenHeight);
          if (h <= 0)
            throw Exception("Height must be positive in Grafen's method.");
        }
        ApplicationTools::displayResult("Total height", TextTools::toString(h));

        double rho = ApplicationTools::getDoubleParameter("rho", cmdArgs, 1., "", true, 2);
        ApplicationTools::displayResult("Grafen's rho", rho);

        PhyloTreeTools::computeBranchLengthsGrafen(*tree, rho);
        double nh = PhyloTreeTools::getHeight(*tree, tree->getRoot());

        tree->scaleTree(h / nh);
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
      {
        shared_ptr<PhyloBranch> branch = mTree[num]->getEdgeToFather(mTree[num]->getNode(static_cast<PhyloTree::NodeIndex>(TextTools::toInt(aveq.substr(5, string::npos)))));
        if (branch)
          branch->setLength(TextTools::toDouble(apeq));
      }
      else
      {
        size_t posun = apeq.find("_");
        size_t posd = aveq.find("_");
        unparsedParams[aveq + (posd != string::npos ? "" : "_" + TextTools::toString(num))] = apeq + (posun != string::npos ? "" : "_" + TextTools::toString(num));
      }
    }
  }

  return mTree;
}

/******************************************************/
/**** SUBSTITUTION RATES *******************************/
/******************************************************/

std::unique_ptr<SubstitutionModelInterface> PhylogeneticsApplicationTools::getSubstitutionModel(
    std::shared_ptr<const Alphabet> alphabet,
    std::shared_ptr<const GeneticCode> gCode,
    std::shared_ptr<const AlignmentDataInterface> data,
    const map<string, string>& params,
    map<string, string>& unparsedParams,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  BppOSubstitutionModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, verbose, warn + 1);
  string modelDescription;
  auto ca = dynamic_pointer_cast<const CodonAlphabet>(alphabet);
  if (ca)
  {
    modelDescription = ApplicationTools::getStringParameter("model", params, "CodonRate(model=JC69)", suffix, suffixIsOptional, warn);
    if (!gCode)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModel(): a GeneticCode instance is required for instantiating a codon model.");
    bIO.setGeneticCode(gCode);
  }
  else if (AlphabetTools::isWordAlphabet(*alphabet))
    modelDescription = ApplicationTools::getStringParameter("model", params, "Word(model=JC69)", suffix, suffixIsOptional, warn);
  else
    modelDescription = ApplicationTools::getStringParameter("model", params, "JC69", suffix, suffixIsOptional, warn);

  std::map<size_t, std::shared_ptr<const AlignmentDataInterface>> mData;
  mData[1]=data;
  
  auto model = bIO.readSubstitutionModel(alphabet, modelDescription, mData, 1, true);

  unparsedParams.insert(bIO.getUnparsedArguments().begin(), bIO.getUnparsedArguments().end());

  return model;
}

std::unique_ptr<BranchModelInterface> PhylogeneticsApplicationTools::getBranchModel(
    std::shared_ptr<const Alphabet> alphabet,
    std::shared_ptr<const GeneticCode> gCode,
    std::shared_ptr<const AlignmentDataInterface> data,
    const map<string, string>& params,
    map<string, string>& unparsedParams,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  BppOBranchModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, verbose, warn + 1);
  string modelDescription;
  auto ca = dynamic_pointer_cast<const CodonAlphabet>(alphabet);
  if (ca)
  {
    modelDescription = ApplicationTools::getStringParameter("model", params, "CodonRate(model=JC69)", suffix, suffixIsOptional, warn);
    if (!gCode)
      throw Exception("PhylogeneticsApplicationTools::getBranchModel(): a GeneticCode instance is required for instantiating a codon model.");
    bIO.setGeneticCode(gCode);
  }
  else if (AlphabetTools::isWordAlphabet(*alphabet))
    modelDescription = ApplicationTools::getStringParameter("model", params, "Word(model=JC69)", suffix, suffixIsOptional, warn);
  else
    modelDescription = ApplicationTools::getStringParameter("model", params, "JC69", suffix, suffixIsOptional, warn);

  std::map<size_t, std::shared_ptr<const AlignmentDataInterface>> mData;
  mData[1]=data;

  auto model = bIO.readBranchModel(alphabet, modelDescription, mData, 1, true);
  map<string, string> tmpUnparsedParameterValues(bIO.getUnparsedArguments());

  unparsedParams.insert(tmpUnparsedParameterValues.begin(), tmpUnparsedParameterValues.end());

  return model;
}

/******************************************************************************/

map<size_t, std::shared_ptr<DiscreteDistributionInterface>> PhylogeneticsApplicationTools::getRateDistributions(
    const map<string, string>& params,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose)
{
  string DistFilePath = ApplicationTools::getAFilePath("rate_distribution.file", params, false, false, suffix, suffixIsOptional, "none", 1);

  map<string, string> paramDist;

  if (DistFilePath != "none")
    paramDist = AttributesTools::getAttributesMapFromFile(DistFilePath, "=");

  paramDist.insert(params.begin(), params.end());

  vector<string> vratesName = ApplicationTools::matchingParameters("rate_distribution*", paramDist);

  BppORateDistributionFormat bIO(true);
  map<size_t, std::shared_ptr<DiscreteDistributionInterface>> mDist;


  for (size_t i = 0; i < vratesName.size(); ++i)
  {
    size_t poseq = vratesName[i].find("=");
    size_t num = 0;
    string suff = vratesName[i].substr(17, poseq - 17);
    bool flag = 0;


    if (TextTools::isDecimalInteger(suff, '$'))
      num = static_cast<size_t>(TextTools::toInt(suff));
    else
    {
      flag = 1;
      num = 0;
    }

    if (verbose)
      if (!flag)
      {
        ApplicationTools::displayMessage("");
        ApplicationTools::displayMessage("Rate " + TextTools::toString(num));
      }

    string distDescription = ApplicationTools::getStringParameter(vratesName[i], paramDist, "", suffix, suffixIsOptional);

    if (num!=0)
      mDist[num] = std::shared_ptr<DiscreteDistributionInterface>(bIO.readDiscreteDistribution(distDescription, true));
  }

  if (mDist.size() == 0)
  {
    ApplicationTools::displayMessage("");
    ApplicationTools::displayMessage("Rate 1");
    string distDescription = ApplicationTools::getStringParameter("rate_distribution", paramDist, "Constant()", suffix, suffixIsOptional);
    mDist[1] = std::shared_ptr<DiscreteDistributionInterface>(bIO.readDiscreteDistribution(distDescription, true));
  }

  return mDist;
}


/*************************************************************/
/******* MODELS **********************************************/
/*************************************************************/

map<size_t, std::unique_ptr<BranchModelInterface>> PhylogeneticsApplicationTools::getBranchModels(
    std::shared_ptr<const Alphabet> alphabet,
    std::shared_ptr<const GeneticCode> gCode,
    const map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mData,
    const map<string, string>& params,
    map<string, string>& unparsedParams,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  if (dynamic_pointer_cast<const CodonAlphabet>(alphabet) && !gCode)
    throw Exception("PhylogeneticsApplicationTools::getBranchModels(): a GeneticCode instance is required for instantiating codon models.");

  string ModelFilePath = ApplicationTools::getAFilePath("models.file", params, false, false, suffix, suffixIsOptional,  "none", 1);

  map<string, string> paramModel;

  if (ModelFilePath != "none")
    paramModel = AttributesTools::getAttributesMapFromFile(ModelFilePath, "=");

  paramModel.insert(params.begin(), params.end());

  vector<string> modelsName = ApplicationTools::matchingParameters("model*", paramModel);

  vector<size_t> modelsNum;
  for (const auto& name : modelsName)
  {
    size_t poseq = name.find("=");
    if (name.find("nodes_id") == string::npos)
    {
      modelsNum.push_back(TextTools::to<size_t>(name.substr(5, poseq - 5)));
    }
  }

  map<size_t, std::unique_ptr<BranchModelInterface>> mModel;

  BppOBranchModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, verbose, warn);
  bIO.setGeneticCode(gCode);

  for (size_t i = 0; i < modelsNum.size(); ++i)
  {
    if (i >= 10)
    {
      bIO.setVerbose(false);
      warn = 10;
      if (i == 10)
        ApplicationTools::displayMessage("");
      ApplicationTools::displayResult("Model " + TextTools::toString(modelsNum[i]), string("..."));
    }
    else
    {
      ApplicationTools::displayMessage("");
      ApplicationTools::displayMessage("Model " + TextTools::toString(modelsNum[i]));
    }

    string modelDescription = ApplicationTools::getStringParameter("model" + TextTools::toString(modelsNum[i]), paramModel, "", suffix, suffixIsOptional, warn);

    unique_ptr<BranchModelInterface> model;
    model = bIO.readBranchModel(alphabet, modelDescription, mData, 0, true);

    map<string, string> tmpUnparsedParameterValues(bIO.getUnparsedArguments());

    for (auto& it : tmpUnparsedParameterValues)
    {
      unparsedParams[it.first + "_" + TextTools::toString(modelsNum[i])] = it.second;
    }

    mModel[modelsNum[i]] = std::move(model);
  }

  return mModel;
}


/******************************************************/
/**** FREQUENCIES SET *********************************/
/******************************************************/

std::unique_ptr<FrequencySetInterface> PhylogeneticsApplicationTools::getFrequencySet(
    std::shared_ptr<const Alphabet> alphabet,
    std::shared_ptr<const GeneticCode> gCode,
    const string& freqDescription,
    const std::map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mData,
    size_t nData,
    map<string, string>& sharedparams,
    const vector<double>& rateFreqs,
    bool verbose,
    int warn)
{
  map<string, string> unparsedParameterValues;
  BppOFrequencySetFormat bIO(BppOFrequencySetFormat::ALL, verbose, warn);
  if (AlphabetTools::isCodonAlphabet(*alphabet))
  {
    if (!gCode)
      throw Exception("PhylogeneticsApplicationTools::getFrequencySet(): a GeneticCode instance is required for instantiating a codon frequencies set.");
    bIO.setGeneticCode(gCode);
  }

  auto pFS = bIO.readFrequencySet(alphabet, freqDescription, mData, nData, true);

  map<string, string> unparsedparam = bIO.getUnparsedArguments();

  sharedparams.insert(unparsedparam.begin(), unparsedparam.end());

  // /////// To be changed for input normalization
  if (rateFreqs.size() > 0)
  {
    pFS = std::make_unique<MarkovModulatedFrequencySet>(std::move(pFS), rateFreqs);
  }

  return pFS;
}


std::unique_ptr<FrequencySetInterface> PhylogeneticsApplicationTools::getRootFrequencySet(
    std::shared_ptr<const Alphabet> alphabet,
    std::shared_ptr<const GeneticCode> gCode,
    const std::map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mData,
    size_t nData,
    const map<string, string>& params,
    map<string, string>& sharedparams,
    const vector<double>& rateFreqs,
    const string& suffix,
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

    auto freq = getFrequencySet(alphabet, gCode, freqDescription, mData, nData, unparams, rateFreqs, verbose, warn + 1);
    freq->setNamespace("root." + freq->getNamespace());

    for (auto& it : unparams)
    {
      sharedparams["root." + it.first] = it.second;
    }

    if (verbose)
      ApplicationTools::displayResult("Root frequencies ", freq->getName());
    return freq;
  }
}


map<size_t, std::unique_ptr<FrequencySetInterface>> PhylogeneticsApplicationTools::getRootFrequencySets(
    std::shared_ptr<const Alphabet> alphabet,
    std::shared_ptr<const GeneticCode> gCode,
    const map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mData,
    const map<string, string>& params,
    map<string, string>& sharedparams,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  if (dynamic_pointer_cast<const CodonAlphabet>(alphabet) && !gCode)
    throw Exception("PhylogeneticsApplicationTools::getRootFrequencySets(): a GeneticCode instance is required for instantiating codon frequencies sets.");

  string RootFilePath = ApplicationTools::getAFilePath("root_freq.file", params, false, false, suffix, suffixIsOptional,  "none", 1);
  map<string, string> paramRF;

  if (RootFilePath != "none")
    paramRF = AttributesTools::getAttributesMapFromFile(RootFilePath, "=");

  paramRF.insert(params.begin(), params.end());

  vector<string> vrfName = ApplicationTools::matchingParameters("root_freq*", paramRF);

  vector<size_t> rfNum;
  for (const auto& rfName : vrfName)
  {
    size_t poseq = rfName.find("=");
    try
    {
      rfNum.push_back(TextTools::to<size_t>(rfName.substr(9, poseq - 9)));
    }
    catch (Exception& e)
    {}
  }

  BppOFrequencySetFormat bIO(BppOFrequencySetFormat::ALL, verbose, warn);
  bIO.setGeneticCode(gCode);

  map<size_t, std::unique_ptr<FrequencySetInterface>> mFS;

  for (size_t i = 0; i < rfNum.size(); ++i)
  {
    ApplicationTools::displayMessage("");
    ApplicationTools::displayMessage("Root Frequencies Set " + TextTools::toString(rfNum[i]));

    string freqDescription = ApplicationTools::getStringParameter("root_freq" + TextTools::toString(rfNum[i]), paramRF, "", suffix, suffixIsOptional, warn);

    map<string, string> args;
    string freqName;

    KeyvalTools::parseProcedure(freqDescription, freqName, args);

    size_t nData = 0;

    if (args.find("data") != args.end())
      nData = TextTools::to<size_t>(args["data"]);

    unique_ptr<FrequencySetInterface> rFS;

    rFS = bIO.readFrequencySet(alphabet, freqDescription, mData, nData, true);

    rFS->setNamespace("root." + rFS->getNamespace());

    map<string, string> unparsedparam = bIO.getUnparsedArguments();

    for (auto& it : unparsedparam)
    {
      sharedparams["root." + it.first + "_" + TextTools::toString(rfNum[i])] = it.second;
    }

    if (verbose)
    {
      // ApplicationTools::displayResult("Root Frequencies Set " + TextTools::toString(rfNum[i]), rFS->getName());
      if (nData != 0)
        ApplicationTools::displayResult("Data used ", TextTools::toString(nData));
    }

    mFS[rfNum[i]] = std::move(rFS);
  }

  return mFS;
}

/******************************************************/
/**** SETOFMODELPATH **********************************/
/******************************************************/

map<size_t, std::unique_ptr<ModelPath>> PhylogeneticsApplicationTools::getModelPaths(
    const std::map<std::string, std::string>& params,
    const map<size_t, std::shared_ptr<BranchModelInterface>>& mModel,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  string ModelPathsPath = ApplicationTools::getAFilePath("path.file", params, false, false, suffix, suffixIsOptional,  "none", warn);
  map<string, string> paramMP;

  if (ModelPathsPath != "none")
    paramMP = AttributesTools::getAttributesMapFromFile(ModelPathsPath, "=");

  paramMP.insert(params.begin(), params.end());

  vector<string> vmpName = ApplicationTools::matchingParameters("path*", paramMP);

  map<size_t, std::unique_ptr<ModelPath>> modelPaths;

  for (size_t i = 0; i < vmpName.size(); ++i)
  {
    const auto& name = vmpName[i];

    string desc = ApplicationTools::getStringParameter(name, paramMP, "", "", true);

    size_t num;
    try
    {
      num = TextTools::to<size_t>(name.substr(4));
    }
    catch (const Exception& e)
    {
      throw Exception("PhylogeneticsApplicationTools::getModelPaths: bad path number in line " + name);
    }

    modelPaths[num] = std::make_unique<ModelPath>();

    if (verbose)
    {
      if (i >= 10)
      {
        if (i == 10)
          ApplicationTools::displayMessage("");
        ApplicationTools::displayResult("Path " + TextTools::toString(num), string("..."));
      }
      else
      {
        ApplicationTools::displayMessage("");
        ApplicationTools::displayMessage("Path " + TextTools::toString(num));
      }
    }

    StringTokenizer st(desc, "&");
    while (st.hasMoreToken())
    {
      string submodel = st.nextToken();
      Vuint submodelNb;
      auto indexo = submodel.find("[");
      auto indexf = submodel.find("]");
      if ((indexo == string::npos) | (indexf == string::npos))
        throw Exception("PhylogeneticsApplicationTools::getModelPaths. Bad path syntax, should contain `[]' symbols: " + submodel);

      auto pos = submodel.find("model");
      if (pos == string::npos)
        throw Exception("PhylogeneticsApplicationTools::getModelPaths. Missing identifier 'model' in description: " + submodel);

      size_t num2 = TextTools::to<size_t>(submodel.substr(pos + 5, indexo - 5 - pos));
      if (mModel.find(num2) == mModel.end())
        throw BadIntegerException("PhylogeneticsApplicationTools::getModelPaths: Wrong model number", static_cast<int>(num2));

      auto pSM = std::dynamic_pointer_cast<MixedTransitionModelInterface>(mModel.at(num2));
      if (!pSM)
        throw Exception("PhylogeneticsApplicationTools::getModelPaths: Model number " + TextTools::toString(num2) + " ( " + mModel.at(num2)->getName() + " ) is not Mixed.");

      string lp2 = submodel.substr(indexo + 1, indexf - indexo - 1);
      StringTokenizer stp2(lp2, ",");
      while (stp2.hasMoreToken())
      {
        string p2 = stp2.nextToken();

        unsigned int n2;
        bool n2ok = true;
        try
        {
          n2 = TextTools::to<unsigned int>(p2);
          if (n2 <= 0 || n2 > pSM->getNumberOfModels())
            n2ok = false;
          else
            submodelNb.push_back(n2 - 1);
        }
        catch (Exception& e)
        {
          Vuint submodnb = pSM->getSubmodelNumbers(p2);
          if (submodelNb.size() == 0)
            submodelNb = submodnb;
          else
            submodelNb = VectorTools::vectorIntersection(submodelNb, submodnb);
        }

        if (!n2ok)
          throw BadIntegerException("PhylogeneticsApplicationTools::getModelPaths: Wrong model number for model " + TextTools::toString(num2), int(n2));
      }

      modelPaths[num]->setModel(pSM, submodelNb);
      if (!modelPaths[num]->getLeadModel())
        modelPaths[num]->setLeadModel(pSM);
    }

    if (verbose &&  (i < 10))
      ApplicationTools::displayResult("Model Path", desc);
  }

  return modelPaths;
}


map<size_t, std::unique_ptr<ModelScenario>> PhylogeneticsApplicationTools::getModelScenarios(
    const std::map<std::string, std::string>& params,
    const map<size_t, std::shared_ptr<ModelPath>>& mModelPath,
    const map<size_t, std::shared_ptr<BranchModelInterface>>& mModel,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  string ModelPathsPath = ApplicationTools::getAFilePath("scenario.file", params, false, false, suffix, suffixIsOptional,  "none", warn);
  map<string, string> paramMS;

  if (ModelPathsPath != "none")
    paramMS = AttributesTools::getAttributesMapFromFile(ModelPathsPath, "=");

  paramMS.insert(params.begin(), params.end());

  vector<string> vmsName = ApplicationTools::matchingParameters("scenario*", paramMS);

  map<size_t, std::unique_ptr<ModelScenario>> somp;

  for (size_t i = 0; i < vmsName.size(); ++i)
  {
    const auto& name = vmsName[i];

    string desc = ApplicationTools::getStringParameter(name, paramMS, "", suffix, suffixIsOptional, warn);

    size_t num;
    try
    {
      num = TextTools::to<size_t>(name.substr(8));
    }
    catch (const Exception& e)
    {
      throw Exception("PhylogeneticsApplicationTools::getModelScenarios: bad scenario number in line " + name);
    }

    somp[num] = std::make_unique<ModelScenario>();

    if (verbose)
    {
      if (i >= 10)
      {
        if (i == 10)
          ApplicationTools::displayMessage("");
        ApplicationTools::displayResult("Scenario " + TextTools::toString(num), string("..."));
      }
      else
      {
        ApplicationTools::displayMessage("");
        ApplicationTools::displayMessage("Scenario " + TextTools::toString(num));
      }
    }

    bool complete = false;
    size_t numpath;

    StringTokenizer st(desc, "&");
    while (st.hasMoreToken())
    {
      string path = st.nextToken();
      bool numok = true;
      try
      {
        if (path == "complete")
          complete = true;
        else if (path.substr(0, 5) == "split")
        {
          auto pos = path.find("model");
          if (pos == string::npos)
            throw Exception("PhylogeneticsApplicationTools::getModelScenarios. Missing identifier 'model' in scenario description: " + path);

          auto poseq = path.find("=", pos);
          size_t num2 = TextTools::to<size_t>(path.substr(poseq + 1));

          if (mModel.find(num2) == mModel.end())
            throw BadIntegerException("PhylogeneticsApplicationTools::getModelScenarios: Wrong model number", static_cast<int>(num2));

          auto pSM = std::dynamic_pointer_cast<MixedTransitionModelInterface>(mModel.at(num2));
          if (!pSM)
            throw Exception("PhylogeneticsApplicationTools::getModelScenarios: Model number " + TextTools::toString(num2) + " ( " + mModel.at(num2)->getName() + " ) is not Mixed.");

          std::vector<std::shared_ptr<ModelPath>> modelPaths;

          auto nmod = pSM->getNumberOfModels();

          for (unsigned int nm = 0; nm < static_cast<unsigned int>(nmod); ++nm)
          {
            auto mp = std::make_shared<ModelPath>();
            mp->setModel(pSM, Vuint({nm}));
            mp->setLeadModel(pSM);
            somp[num]->addModelPath(mp);
          }
        }
        else
        {
          numpath = TextTools::to<size_t>(path.substr(4));
          if (mModelPath.find(numpath) == mModelPath.end())
            numok = false;
          else
            somp[num]->addModelPath(mModelPath.at(numpath));
        }
      }
      catch (Exception& e)
      {
        Exception("PhylogeneticsApplicationTools::getModelScenarios: wrong path description " + path);
      }

      if (!numok)
        throw BadIntegerException("PhylogeneticsApplicationTools::getModelScenarios: Wrong path number", static_cast<int>(numpath));
    }


    if (verbose &&  (i < 10))
      ApplicationTools::displayResult("Model Scenario", desc);

    if (complete)
    {
      if (somp[num]->getNumberOfModelPaths() == 0)
        throw Exception("PhylogeneticsApplicationTools::getModelScenarios: 'complete' is not possible on empty scenarios");
      somp[num]->complete();
    }

    somp[num]->computeModelPathsProbabilities();
  }

  return somp;
}

/******************************************************/
/********** SUBSTITUTION PROCESS   ********************/
/******************************************************/

unique_ptr<AutonomousSubstitutionProcessInterface> PhylogeneticsApplicationTools::getSubstitutionProcess(
    std::shared_ptr<const Alphabet> alphabet,
    std::shared_ptr<const GeneticCode> gCode,
    std::shared_ptr<const AlignmentDataInterface> pData,
    const vector< shared_ptr<PhyloTree>>& vTree,
    const map<string, string>& params,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  // Read files with same process as SubstitutionCollection

  map<string, string> unparsedParams;

  map<size_t, std::shared_ptr<const AlignmentDataInterface>> mData;
  mData[1] = pData;

  map<size_t, std::shared_ptr<PhyloTree>> mTree;
  size_t i = 1;
  for (auto it : vTree)
  {
    mTree[i++] = it;
  }

  auto mModU = getBranchModels(alphabet, gCode, mData, params, unparsedParams, suffix, suffixIsOptional, verbose, warn);
  auto mMod = uniqueToSharedMap<BranchModelInterface>(mModU);

  auto mRootFreqU = getRootFrequencySets(alphabet, gCode, mData, params, unparsedParams, suffix, suffixIsOptional, verbose, warn);
  auto mRootFreq = uniqueToSharedMap<FrequencySetInterface>(mRootFreqU);

  auto mDist = getRateDistributions(params, suffix, suffixIsOptional, verbose);

  auto mPathU = getModelPaths(params, mMod, suffix, suffixIsOptional, verbose, warn);
  auto mPath = uniqueToSharedMap<ModelPath>(mPathU);

  auto mScenU = getModelScenarios(params, mPath, mMod, suffix, suffixIsOptional, verbose, warn);
  auto mScen = uniqueToSharedMap<ModelScenario>(mScenU);

  auto SPC = getSubstitutionProcessCollection(alphabet, gCode, mTree,
        mMod, mRootFreq, mDist, mScen,
        params, unparsedParams, suffix, suffixIsOptional, verbose, warn);

  // Get relevant objects from Collection to build an AutonomousSubstitutionProcess
  unique_ptr<AutonomousSubstitutionProcessInterface> ASP;

  auto psNum = SPC->getSubstitutionProcessNumbers();
  if (psNum.size() == 0)
    throw Exception("PhylogeneticsApplicationTools::getSubstitutionProcess : missing process in parameters.");

  size_t maxps = *max_element(psNum.begin(), psNum.end());

  SubstitutionProcessCollectionMember& procm = SPC->substitutionProcess(maxps);

  auto distproc = procm.getRateDistribution();

  auto rootproc = procm.getRootFrequencySet();

  auto scen = procm.getModelScenario();

  auto vmodnb = procm.getModelNumbers();

  if (vmodnb.size() == 1)
  {
    if (!distproc)
      ASP = make_unique<SimpleSubstitutionProcess>(procm.getModel(1), procm.getParametrizablePhyloTree(), rootproc);
    else
      ASP = make_unique<RateAcrossSitesSubstitutionProcess>(procm.getModel(1), procm.getRateDistribution(), procm.getParametrizablePhyloTree(), rootproc);
  }
  else
  {
    auto NHSP = make_unique<NonHomogeneousSubstitutionProcess>(procm.getRateDistribution(), procm.getParametrizablePhyloTree(), rootproc);

    for (auto nb:vmodnb)
    {
      NHSP->addModel(procm.getModel(nb), procm.getNodesWithModel(nb));
    }

    if (!NHSP->isFullySetUp(false))
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionProcess: process not fully set up.");

    ASP = std::move(NHSP);
  }

  if (procm.getModelScenario())
    ASP->setModelScenario(procm.getModelScenario());

  return ASP;
}


/************************************************************/

bool PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember(
    SubstitutionProcessCollection& SubProColl,
    size_t procNum,
    const map<string, string>& params,
    bool verbose,
    int warn)
{
  string procName = "";
  map<string, string> args;

  string procDesc = ApplicationTools::getStringParameter("process", params, "", TextTools::toString(procNum), warn);

  KeyvalTools::parseProcedure(procDesc, procName, args);

  if ((procName != "OnePerBranch") && (procName != "Homogeneous") && (procName != "Nonhomogeneous") &&  (procName != "NonHomogeneous"))
  {
    if (warn >= 2)
      ApplicationTools::displayWarning("Warning, unknown process name: " + procName);

    return 0;
  }

  // ///
  // tree number

  size_t numTree;

  if (args.find("tree") == args.end())
  {
    if (warn)
      ApplicationTools::displayWarning("Warning, missing tree for  process name: " + procName);
    numTree = 0;
  }
  else
  {
    numTree = (size_t) ApplicationTools::getIntParameter("tree", args, 1, "", true, warn);

    if (numTree != 0 && !SubProColl.hasTreeNumber(numTree))
      throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : unknown tree number", (int)numTree);
  }

  // /////
  // rate number

  size_t numRate = 0;
  if (args.find("rate") == args.end())
  {
    const auto& vrdn = SubProColl.getRateDistributionNumbers();
    numRate = 0;
    for (auto rdn:vrdn)
    {
      if (SubProColl.rateDistribution(rdn).getName() == "Constant")
      {
        numRate = rdn;
        break;
      }
    }
    if (numRate == 0)
    {
      for (uint i = 1; i <= *std::max_element(vrdn.begin(), vrdn.end()) + 1; i++)
      {
        if (std::find(vrdn.begin(), vrdn.end(), i) == vrdn.end())
        {
          numRate = i;
          break;
        }
      }
      SubProColl.addDistribution(std::make_shared<ConstantRateDistribution>(), numRate);
    }
  }
  else
  {
    string sRate = ApplicationTools::getStringParameter("rate", args, "1", "", true, warn);

    size_t pp = sRate.find(".");

    numRate = TextTools::to<size_t>(sRate.substr(0, pp));

    if (!SubProColl.hasDistributionNumber(numRate))
      throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : unknown rate number", (int)numRate);

    if (pp != string::npos)
    {
      size_t numSRate = TextTools::to<size_t>(sRate.substr(pp + 1));
      SubProColl.addDistribution(
          std::make_shared<ConstantDistribution>(
          SubProColl.rateDistribution(numRate).getCategory(numSRate)),
          10000 * (numRate + 1) + numSRate);

      numRate = 10000 * (numRate + 1) + numSRate;
    }
  }

  // ////////
  // root freq number

  bool stationarity = (args.find("root_freq") == args.end());
  size_t numFreq = 0;

  if (!stationarity)
  {
    numFreq = (size_t) ApplicationTools::getIntParameter("root_freq", args, 1, "", true, warn);
    if (!SubProColl.hasFrequenciesNumber(numFreq))
      throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : unknown root frequencies number", (int)numFreq);
  }

  // ///
  // scenario number

  size_t numScen = 0;

  if (args.find("scenario") != args.end())
  {
    numScen = (size_t) ApplicationTools::getIntParameter("scenario", args, 1, "", true, warn);

    if (!SubProColl.hasModelScenario(numScen))
      throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : unknown scenario number", (int)numScen);
  }

  // ////////////////
  // / models

  if (verbose)
  {
    ApplicationTools::displayMessage("");
    ApplicationTools::displayMessage("Process " + TextTools::toString(procNum));
  }

  if (procName == "Homogeneous")
  {
    if (args.find("model") == args.end())
      throw Exception("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember. A model number is compulsory.");

    size_t numModel = (size_t) ApplicationTools::getIntParameter("model", args, 1, "", true, warn);

    if (!SubProColl.hasModelNumber(numModel))
      throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : unknown model number", static_cast<int>(numModel));

    map<size_t, vector<unsigned int>> mModBr;

    vector<uint> vNodes;
    if (numTree != 0)
      vNodes = SubProColl.tree(numTree).getAllEdgesIndexes();
    else
      vNodes = {0}
    ;
    mModBr[numModel] = vNodes;

    if (verbose)
    {
      ApplicationTools::displayResult("Process type", string("Homogeneous"));
      ApplicationTools::displayResult (" Model number", TextTools::toString(numModel));
      ApplicationTools::displayResult (" Tree number", TextTools::toString(numTree));
      if (numRate < 10000)
        ApplicationTools::displayResult (" Rate number", TextTools::toString(numRate));
      else
        ApplicationTools::displayResult (" Rate number", TextTools::toString(numRate / 10000 - 1) + "." + TextTools::toString(numRate % 10000));

      if (numScen != 0)
        ApplicationTools::displayResult (" Scenario number", TextTools::toString(numScen));

      if (!stationarity)
        ApplicationTools::displayResult (" Root frequencies number", TextTools::toString(numFreq));
      else
        ApplicationTools::displayMessage(" Stationarity assumed.");
    }

    if (stationarity)
      SubProColl.addSubstitutionProcess(procNum, mModBr, numTree, numRate);
    else
      SubProColl.addSubstitutionProcess(procNum, mModBr, numTree, numRate, numFreq);
  }

  else if ((procName == "Nonhomogeneous") ||  (procName == "NonHomogeneous"))
  {
    if (numTree == 0)
      throw Exception("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : missing tree number for process " + TextTools::toString(procName));

    size_t indModel = 1;
    map<size_t, vector<unsigned int>> mModBr;

    while (args.find("model" + TextTools::toString(indModel)) != args.end())
    {
      size_t numModel = (size_t) ApplicationTools::getIntParameter("model" + TextTools::toString(indModel), args, 1, "", true, warn);

      if (mModBr.find(numModel) != mModBr.end())
        throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : model number seen twice.", (int)numModel);

      vector<unsigned int> nodesId;

      auto snodesid = "model" + TextTools::toString(indModel)  + ".nodes_id";
      auto descnodes = ApplicationTools::getStringParameter(snodesid, args, "", "", true, warn);


      auto tree = SubProColl.getTree(numTree);
      if (descnodes == "All")
      {
        nodesId = tree->getEdgeIndexes(tree->getSubtreeEdges(tree->getRoot()));
      }
      else if (descnodes == "Leaves")
      {
        nodesId = tree->getNodeIndexes(tree->getLeavesUnderNode(tree->getRoot()));
      }
      else if (descnodes == "NoLeaves")
      {
        auto allIds = tree->getEdgeIndexes(tree->getSubtreeEdges(tree->getRoot()));
        auto leavesId = tree->getNodeIndexes(tree->getLeavesUnderNode(tree->getRoot()));
        VectorTools::diff(allIds, leavesId, nodesId);
      }
      else
        nodesId = ApplicationTools::getVectorParameter<unsigned int>(snodesid, args, ',', ':', "", "", true, warn);

      mModBr[numModel] = nodesId;
      indModel++;
    }

    if (verbose)
    {
      ApplicationTools::displayResult("Process type", string("NonHomogeneous"));

      for (auto& it : mModBr)
      {
        ApplicationTools::displayResult (" Model number" + TextTools::toString(it.first) + " associated to", TextTools::toString(it.second.size()) + " node(s).");
      }
      ApplicationTools::displayResult (" Tree number", TextTools::toString(numTree));
      ApplicationTools::displayResult (" Rate number", TextTools::toString(numRate));
      if (!stationarity)
        ApplicationTools::displayResult (" Root frequencies number", TextTools::toString(numFreq));
      else
        ApplicationTools::displayMessage(" Stationarity assumed.");
    }

    if (stationarity)
      SubProColl.addSubstitutionProcess(procNum, mModBr, numTree, numRate);
    else
      SubProColl.addSubstitutionProcess(procNum, mModBr, numTree, numRate, numFreq);
  }
  else if (procName == "OnePerBranch")
  {
    if (numTree == 0)
      throw Exception("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : missing tree number for process " + TextTools::toString(procName));

    if (args.find("model") == args.end())
      throw Exception("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember. A model number is compulsory.");

    size_t numModel = (size_t) ApplicationTools::getIntParameter("model", args, 1, "", true, warn);

    if (!SubProColl.hasModelNumber(numModel))
      throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : unknown model number", (int)numModel);

    vector<string> sharedParameters = ApplicationTools::getVectorParameter<string>("shared_parameters", args, ',', "", "", true, 1);

    if (stationarity)
      SubProColl.addOnePerBranchSubstitutionProcess(procNum, numModel, numTree, numRate, sharedParameters);
    else
      SubProColl.addOnePerBranchSubstitutionProcess(procNum, numModel, numTree, numRate, numFreq, sharedParameters);

    if (verbose)
    {
      ApplicationTools::displayResult("Process type", string("OnePerBranch"));

      ApplicationTools::displayResult (" Model number", TextTools::toString(numModel));
      ApplicationTools::displayResult (" Tree number", TextTools::toString(numTree));
      ApplicationTools::displayResult (" Rate number", TextTools::toString(numRate));
      if (!stationarity)
        ApplicationTools::displayResult (" Root frequencies number", TextTools::toString(numFreq));
      else
        ApplicationTools::displayMessage(" Stationarity assumed.");

      for (const auto& sP : sharedParameters)
      {
        ApplicationTools::displayResult(" Shared parameter", sP);
      }
    }
  }

  if (numScen != 0)
    SubProColl.substitutionProcess(procNum).setModelScenario(numScen);

  return true;
}


/******************************************************************************/

unique_ptr<SubstitutionProcessCollection> PhylogeneticsApplicationTools::getSubstitutionProcessCollection(
    std::shared_ptr<const Alphabet> alphabet,
    std::shared_ptr<const GeneticCode> gCode,
    const map<size_t, std::shared_ptr<PhyloTree>>& mTree,
    const map<size_t, std::shared_ptr<BranchModelInterface>>& mMod,
    const map<size_t, std::shared_ptr<FrequencySetInterface>>& mRootFreq,
    const map<size_t, std::shared_ptr<DiscreteDistributionInterface>>& mDist,
    const map<size_t, std::shared_ptr<ModelScenario>>& mScen,
    const map<string, string>& params,
    map<string, string>& unparsedParams,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  auto SPC = make_unique<SubstitutionProcessCollection>();

  map<string, double> existingParameters;

  // ///////////////////////
  // Trees

  for (const auto& itt : mTree)
  {
    if (itt.second)
    {
      SPC->addTree(std::make_shared<ParametrizablePhyloTree>(*(itt.second)), itt.first);
    }
  }

  // ///////////////////////
  // Rates

  for (const auto& itd : mDist)
  {
    SPC->addDistribution(itd.second, itd.first);
  }

  // ////////////////////////
  // Models

  for (const auto& itm : mMod)
  {
    SPC->addModel(itm.second, itm.first);
  }

  // ///////////////////////////
  // Root Frequencies

  for (const auto& itr : mRootFreq)
  {
    SPC->addFrequencies(itr.second, itr.first);
  }

  // ///////////////////////
  // Scenarios

  for (const auto& itt : mScen)
  {
    SPC->addScenario(itt.second, itt.first);
  }

  // //////////////////////////////
  // Now processes

  vector<string> vProcName = ApplicationTools::matchingParameters("process*", params);

  if (vProcName.size() == 0)
    throw Exception("Missing process in construction of SubstitutionProcessCollection.");

  for (size_t nT = 0; nT < vProcName.size(); nT++)
  {
    size_t poseq = vProcName[nT].find("=");
    size_t num;
    size_t len = 7;

    string suff = vProcName[nT].substr(len, poseq - len);

    if (TextTools::isDecimalInteger(suff, '$'))
      num = TextTools::to<size_t>(suff);
    else
      num = 1;

    bool addok = addSubstitutionProcessCollectionMember(*SPC, num, params, (nT < 10 ? verbose : false), warn);

    if (addok)
    {
      if (nT == 10)
        ApplicationTools::displayMessage("");

      if (nT >= 10)
        ApplicationTools::displayResult("Process" + TextTools::toString(num), string("..."));
    }
  }


  // string ProcessFilePath = ApplicationTools::getAFilePath("processes.file", params, false, false, suffix, suffixIsOptional);

  // map<string, string> paramProcess;

  // if (ProcessFilePath!="none")
  //   paramModel=AttributesTools::getAttributesMapFromFile(ProcessFilePath,"=");

  // paramProcess.insert(params.begin(), params.end());

  // vector<string> processName=ApplicationTools::matchingParameters("process*", paramProcess);

  // vector<size_t> processNum;
  // for (size_t i=0; i< processName.size(); i++)
  // {
  //   size_t poseq=processName[i].find("=");
  //   processNum.push_back((size_t)TextTools::toInt(processName[i].substr(7,poseq-7)));
  // }

  // if (processNum.size()==0)
  //   throw Exception("Missing process in construction of SubstitutionProcessCollection.");

  // for (size_t i=0; i<processNum.size(); i++)
  //   addSubstitutionProcessCollectionMember(*SPC, params, processNum[i]);


  // /////////////////////////
  // Now set shared parameters:

  // ////// Aliasing
  // Finally check parameter aliasing:

  for (const auto& param : params)
  {
    try
    {
      auto v2 = TextTools::toDouble(param.second);
      SPC->setParameterValue(param.first, v2);
    }
    catch (Exception& e)
    {
      if (SPC->hasParameter(param.first))
        unparsedParams[param.first] = param.second;
    }
  }

  SPC->aliasParameters(unparsedParams, verbose);

  return SPC;
}

/******************************************************/
/**** SEQUENCE EVOLUTIONS *****************************/
/******************************************************/

map<size_t, unique_ptr<SequenceEvolution>> PhylogeneticsApplicationTools::getSequenceEvolutions(
    shared_ptr<SubstitutionProcessCollection> SPC,
    const map<string, string>& params,
    map<string, string>& unparsedParams,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  map<string, string> paramEvol;

  paramEvol.insert(params.begin(), params.end());

  vector<string> evolsName = ApplicationTools::matchingParameters("process*", paramEvol);

  vector<size_t> evolsNum;
  for (size_t i = 0; i < evolsName.size(); ++i)
  {
    size_t poseq = evolsName[i].find("=");
    evolsNum.push_back(TextTools::to<size_t>(evolsName[i].substr(7, poseq - 7)));
  }

  map<size_t, unique_ptr<SequenceEvolution>> mEvol;

  for (size_t mPi = 0; mPi < evolsNum.size(); ++mPi)
  {
    if (SPC->hasSubstitutionProcessNumber(evolsNum[mPi]))
      continue;

    unique_ptr<SequenceEvolution> nEvol;

    string evolName = "";
    map<string, string> args;

    string evolDesc = ApplicationTools::getStringParameter("process", params, "", TextTools::toString(evolsNum[mPi]), warn);

    KeyvalTools::parseProcedure(evolDesc, evolName, args);

    // Process
    if (verbose)
    {
      ApplicationTools::displayMessage("");
      ApplicationTools::displayMessage("Process " + TextTools::toString(evolsNum[mPi]));
    }

    if (evolName == "Simple")
    {
      size_t nproc = (size_t) ApplicationTools::getIntParameter("process", args, ',', "");
      if (!SPC->hasSubstitutionProcessNumber(nproc))
        throw BadIntegerException("PhylogeneticsApplicationTools::getEvolutions. Unknown process number:", (int)nproc);

      nEvol = make_unique<OneProcessSequenceEvolution>(SPC->getSubstitutionProcess(nproc), nproc);
      if (verbose)
      {
        ApplicationTools::displayResult("Process type", string("Simple"));
        ApplicationTools::displayResult (" Process number", TextTools::toString(nproc));
      }
    }
    else
    {
      size_t indProc = 1;
      vector<size_t> vproc;

      while (args.find("process" + TextTools::toString(indProc)) != args.end())
      {
        size_t numProc = (size_t) ApplicationTools::getIntParameter("process" + TextTools::toString(indProc), args, 1, "", true, warn);

        vproc.push_back(numProc);
        indProc++;
      }

      if (vproc.size() == 0)
        throw BadIntegerException("PhylogeneticsApplicationTools::getEvolutions. A process number is compulsory for process", (int)indProc);

      for (size_t i = 0; i < vproc.size(); ++i)
      {
        if (!SPC->hasSubstitutionProcessNumber(vproc[i]))
          throw BadIntegerException("PhylogeneticsApplicationTools::getEvolutions. Unknown process number:", (int)vproc[i]);
      }

      if (evolName == "Partition")
      {
        // parse all processes sites

        vector<size_t> vMap;

        map<size_t, size_t> posProc;

        for (size_t i = 0; i < vproc.size(); ++i)
        {
          string prefix = "process" + TextTools::toString(i + 1);

          vector<size_t> procPos = ApplicationTools::getVectorParameter<size_t>(prefix + ".sites", args, ',', ':', TextTools::toString(i), "", true, true);

          for (size_t j = 0; j < procPos.size(); ++j)
          {
            if (posProc.find(procPos[j]) != posProc.end())
              throw BadIntegerException("A process position is defined twice ", (int)procPos[j]);
            else
              posProc[procPos[j]] = vproc[i];
          }
        }

        size_t pos = posProc.begin()->first;

        while (posProc.find(pos) != posProc.end())
        {
          vMap.push_back(posProc[pos]);
          pos++;
        }

        if (vMap.size() != posProc.size())
          throw Exception("Error : there are gaps in the process sites");

        if (verbose)
          ApplicationTools::displayResult("Process type", string("Partition"));

        auto pMP = make_unique<PartitionSequenceEvolution>(SPC, vMap);

        nEvol = std::move(pMP);
      }
      else if (evolName == "Mixture")
      {
        auto pMP = make_unique<MixtureSequenceEvolution>(SPC, vproc);

        if (verbose)
          ApplicationTools::displayResult("Process type", string("Mixture"));

        size_t nbP = pMP->getNumberOfSubstitutionProcess();

        vector<double> vprob = ApplicationTools::getVectorParameter<double>("probas", args, ',', "(" + VectorTools::paste(vector<double>(nbP, 1. / (double)nbP)) + ")");
        if (vprob.size() != 1)
        {
          if (vprob.size() != nbP)
            throw BadSizeException("Wrong size of probas description in Mixture", vprob.size(), nbP);
          Simplex si(vprob);
          pMP->setSubProcessProb(si);
        }

        nEvol = std::move(pMP);
      }
      else if (evolName == "HMM")
      {
        auto pMP = make_unique<HmmSequenceEvolution>(SPC, vproc);

        if (verbose)
          ApplicationTools::displayResult("Process type", string("HMM"));

        size_t nbP = pMP->getNumberOfSubstitutionProcess();

        string vs = "(" + VectorTools::paste(vector<double>(nbP, 1. / (double)nbP), ",") + ")";
        string vvs = "(";
        for (size_t i = 0; i < nbP; ++i)
        {
          vvs += (i == 0 ? "" : ",") + vs;
        }
        vvs += ")";

        RowMatrix<double> mat = ApplicationTools::getMatrixParameter<double>("probas", args, ',', vvs);

        FullHmmTransitionMatrix fhtm(pMP->hmmTransitionMatrix().getHmmStateAlphabet(), pMP->getNamespace());
        fhtm.setTransitionProbabilities(mat);

        pMP->matchParametersValues(fhtm.getParameters());

        nEvol = std::move(pMP);
      }
      else if (evolName == "AutoCorr")
      {
        auto pMP = make_unique<AutoCorrelationSequenceEvolution>(SPC, vproc);

        size_t nbP = pMP->getNumberOfSubstitutionProcess();

        if (verbose)
          ApplicationTools::displayResult("Process type", string("AutoCorr"));

        string vs = "(" + VectorTools::paste(vector<double>(nbP, 1. / static_cast<int>(nbP)), ",") + ")";

        vector<double> v = ApplicationTools::getVectorParameter<double>("lambdas", args, ',', vs);

        ParameterList pl;

        for (size_t i = 0; i < v.size(); ++i)
        {
          pl.addParameter(Parameter("AutoCorr.lambda" + TextTools::toString(i + 1), v[i]));
        }

        pMP->matchParametersValues(pl);

        nEvol = std::move(pMP);
      }
      else
        throw Exception("Unknown Process description : " + evolName);

      if (verbose)
      {
        ApplicationTools::displayResult (" Process numbers", VectorTools::paste(vproc, ","));
        ApplicationTools::displayMessage("");
      }
    }

    mEvol[evolsNum[mPi]] = std::move(nEvol);
  }

  return mEvol;
}


/******************************************************/
/**** PHYLO LIKELIHOODS *********************************/
/******************************************************/

std::shared_ptr<PhyloLikelihoodContainer> PhylogeneticsApplicationTools::getPhyloLikelihoodContainer(
    Context& context,
    shared_ptr<SubstitutionProcessCollection> SPC,
    map<size_t, std::shared_ptr<SequenceEvolution>>& mSeqEvol,
    const map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mData,
    const map<string, string>& params,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  auto mPhylo = std::make_shared<PhyloLikelihoodContainer>(context, SPC);

  // get all members of the collection and link then to Configured Objects
  auto collNodes = mPhylo->getCollectionNodes();

  // the phylo members
  map<string, string> paramPhyl;
  paramPhyl.insert(params.begin(), params.end());
  vector<string> phylosName = ApplicationTools::matchingParameters("phylo*", paramPhyl);

  // map of dependencies between phylolikelihoods

  map<size_t, vector<size_t>> phylosMap;

  for (size_t i = 0; i < phylosName.size(); ++i)
  {
    size_t poseq = phylosName[i].find("=");
    size_t phyln = TextTools::to<size_t>(phylosName[i].substr(5, poseq - 5));

    if (phyln == 0)
      throw BadIntegerException("PhylogeneticsApplicationTools::getPhyloLikelihoodContainer : Forbidden Phylo Number", 0);

    string phyloDesc = ApplicationTools::getStringParameter("phylo", params, "Single", TextTools::toString(phyln), 2);

    map<string, string> args;

    string phyloName = "";

    KeyvalTools::parseProcedure(phyloDesc, phyloName, args);

    size_t indPhyl = 1;
    vector<size_t> vphyl;

    while (args.find("phylo" + TextTools::toString(indPhyl)) != args.end())
    {
      size_t numPhyl = (size_t) ApplicationTools::getIntParameter("phylo" + TextTools::toString(indPhyl), args, 1, "", true, warn);
      vphyl.push_back(numPhyl);
      indPhyl++;
    }

    phylosMap[phyln] = vphyl;
  }

  vector<size_t> usedPhylo;

  // //////////////////////////////////////////
  // First the phylos that do not depend on other phylos

  uint nbPh(0);
  bool verbhere(verbose);

  for (const auto& it : phylosMap)
  {
    nbPh++;

    if (it.second.size() != 0)
      continue;

    size_t phylonum = it.first;

    std::shared_ptr<PhyloLikelihoodInterface> nPL;
    string phyloName = "";

    map<string, string> args;

    string phyloDesc = ApplicationTools::getStringParameter("phylo", params, "Single", TextTools::toString(phylonum), warn);

    if (verbose)
    {
      if (nbPh <= 20)
        ApplicationTools::displayMessage("");
      else
        verbhere = false;

      ApplicationTools::displayMessage("Phylolikelihood " + TextTools::toString(phylonum));
    }

    KeyvalTools::parseProcedure(phyloDesc, phyloName, args);

    // Data

    size_t nData = (args.find("data") == args.end() ? 1 : TextTools::to<size_t>(args["data"]));

    if (mData.find(nData) == mData.end())
    {
      ApplicationTools::displayWarning("PhylogeneticsApplicationTools::getPhyloLikelihoodContainer. Data number is wrong:" + TextTools::toString(nData) + ". Not built.");
      continue;
    }

    auto data = mData.find(nData)->second;

    if (!data)
    {
      ApplicationTools::displayWarning("PhylogeneticsApplicationTools::getPhyloLikelihoodContainer. Data " + TextTools::toString(nData) + " does not match with aligned sequences. Not built.");
      continue;
    }

    if (verbhere)
      ApplicationTools::displayResult(" Data used ", TextTools::toString(nData));

    // Sequence Evolution or process

    size_t nProcess = (args.find("process") == args.end() ? 1 : TextTools::to<size_t>(args["process"]));
    if (verbhere)
      ApplicationTools::displayResult(" Process ", TextTools::toString(nProcess));


    // Construction

    if (SPC->hasSubstitutionProcessNumber(nProcess))
    {
      auto l = std::make_shared<LikelihoodCalculationSingleProcess>(collNodes, data, nProcess);
      nPL = make_unique<SingleProcessPhyloLikelihood>(context, l, nProcess, nData);
    }
    else if (mSeqEvol.find(nProcess) != mSeqEvol.end())
    {
      // ////////////////
      // / from sequence evolutions to phyloLikelihoods

      auto opse = dynamic_pointer_cast<OneProcessSequenceEvolution>(mSeqEvol[nProcess]);

      if (opse)
        nPL = make_unique<OneProcessSequencePhyloLikelihood>(data, opse, collNodes, nProcess, nData);
      else
      {
        auto mse = dynamic_pointer_cast<MixtureSequenceEvolution>(mSeqEvol[nProcess]);

        if (mse)
          nPL = make_unique<MixtureProcessPhyloLikelihood>(data, mse, collNodes, nProcess, nData);

        else
        {
          auto hse = dynamic_pointer_cast<HmmSequenceEvolution>(mSeqEvol[nProcess]);

          if (hse)
            nPL = make_unique<HmmProcessPhyloLikelihood>(data, hse, collNodes, nProcess, nData);

          else
          {
            auto ase = dynamic_pointer_cast<AutoCorrelationSequenceEvolution>(mSeqEvol[nProcess]);

            if (ase)
              nPL = make_unique<AutoCorrelationProcessPhyloLikelihood>(data, ase, collNodes, nProcess, nData);
            else
            {
              auto pse = dynamic_pointer_cast<PartitionSequenceEvolution>(mSeqEvol[nProcess]);

              if (pse)
                nPL = make_unique<PartitionProcessPhyloLikelihood>(data, pse, collNodes, nProcess, nData);

              else
                throw Exception("PhylogeneticsApplicationTools::getPhyloLikelihoodContainer : Unknown Sequence Evolution.");
            }
          }
        }
      }
    }
    else
      throw BadIntegerException("PhylogeneticsApplicationTools::getPhyloLikelihoodContainer : Unknown Process number.", (int)nProcess);

    mPhylo->addPhyloLikelihood(phylonum, std::move(nPL));
    usedPhylo.push_back(phylonum);
  }

  // Now clean the map
  for (auto it = phylosMap.begin(); it != phylosMap.end();)
  {
    if (it->second.size() == 0)
    {
      phylosMap.erase(it++);
      continue;
    }

    vector<size_t>& vphyl = it->second;
    for (size_t i = vphyl.size(); i > 0; i--)
    {
      vector<size_t>::iterator posp = find(usedPhylo.begin(), usedPhylo.end(), vphyl[i - 1]);
      if (posp != usedPhylo.end())
        vphyl.erase(vphyl.begin() + static_cast<ptrdiff_t>(i - 1));
    }
    ++it;
  }

  // Proceed the other phylos

  while (phylosMap.size() != 0)  // there is still phylos to be treated
  {
    if (usedPhylo.size() == 0)
    {
      ApplicationTools::displayWarning("Warning, some phylolikelihoods are not used.");
      break;
    }

    usedPhylo.clear();

    for (map<size_t, vector<size_t>>::iterator it = phylosMap.begin(); it != phylosMap.end(); it++)
    {
      nbPh++;

      if (it->second.size() == 0)
      {
        size_t phylonum = it->first;

        unique_ptr<PhyloLikelihoodInterface> nPL;
        string phyloName = "";

        map<string, string> args;

        string phyloDesc = ApplicationTools::getStringParameter("phylo", params, "Single", TextTools::toString(phylonum), warn);
        KeyvalTools::parseProcedure(phyloDesc, phyloName, args);

        if (verbose)
        {
          if (nbPh <= 20)
            ApplicationTools::displayMessage("");
          else
            verbhere = false;

          ApplicationTools::displayMessage("Phylolikelihood " + TextTools::toString(phylonum));
        }

        KeyvalTools::parseProcedure(phyloDesc, phyloName, args);

        size_t indPhylo = 1;
        vector<size_t> vPhylo;

        while (args.find("phylo" + TextTools::toString(indPhylo)) != args.end())
        {
          size_t numPhylo = (size_t) ApplicationTools::getIntParameter("phylo" + TextTools::toString(indPhylo), args, 1, "", true, warn);

          vPhylo.push_back(numPhylo);
          indPhylo++;
        }

        if (phyloName == "Mixture")
        {
          auto pMA = make_unique<AlignedPhyloLikelihoodMixture>(context, std::move(mPhylo), vPhylo);
          vector<double> vprob = ApplicationTools::getVectorParameter<double>("probas", args, ',', "(" + VectorTools::paste(vector<double>(vPhylo.size(), 1. / (double)vPhylo.size())) + ")");
          if (vprob.size() != 1)
          {
            if (vprob.size() != vPhylo.size())
              throw BadSizeException("Wrong size of probas description in Mixture", vprob.size(), vPhylo.size());
            Simplex si(vprob);
            pMA->setPhyloProb(si);
          }

          nPL = std::move(pMA);
        }
        else if (phyloName == "HMM")
        {
          auto pMA = make_unique<AlignedPhyloLikelihoodHmm>(context, std::move(mPhylo), vPhylo);

          size_t nbP = pMA->getNumbersOfPhyloLikelihoods().size();

          string vs = "(" + VectorTools::paste(vector<double>(nbP, 1. / static_cast<double>(nbP)), ",") + ")";
          string vvs = "(";
          for (size_t i = 0; i < nbP; ++i)
          {
            vvs += (i == 0 ? "" : ",") + vs;
          }
          vvs += ")";

          RowMatrix<double> mat = ApplicationTools::getMatrixParameter<double>("probas", args, ',', vvs);

          FullHmmTransitionMatrix fhtm(pMA->getHmmStateAlphabet(), pMA->getNamespace());
          fhtm.setTransitionProbabilities(mat);

          pMA->matchParametersValues(fhtm.getParameters());

          nPL = std::move(pMA);
        }
        else if (phyloName == "AutoCorr")
        {
          auto pMA = make_unique<AlignedPhyloLikelihoodAutoCorrelation>(context, std::move(mPhylo), vPhylo);

          size_t nbP = pMA->getNumbersOfPhyloLikelihoods().size();

          string vs = "(" + VectorTools::paste(vector<double>(nbP, 1. / static_cast<double>(nbP)), ",") + ")";

          vector<double> v = ApplicationTools::getVectorParameter<double>("lambdas", args, ',', vs);

          ParameterList pl;

          for (size_t i = 0; i < v.size(); ++i)
          {
            pl.addParameter(Parameter("AutoCorr.lambda" + TextTools::toString(i + 1), v[i]));
          }

          pMA->matchParametersValues(pl);

          nPL = std::move(pMA);
        }
        else if (phyloName == "Product")
        {
          auto pAP = make_unique<AlignedPhyloLikelihoodProduct>(context, std::move(mPhylo), vPhylo);

          nPL = std::move(pAP);
        }
        else
          throw Exception("PhylogeneticsApplicationTools::getPhyloLikelihoodContainer : Unknown Phylo name " + phyloName);

        if (verbhere)
        {
          ApplicationTools::displayResult(" Phylolikelihood type", phyloName);
          ApplicationTools::displayResult(" Phylo numbers", VectorTools::paste(vPhylo, ","));
        }

        mPhylo->addPhyloLikelihood(phylonum, std::move(nPL));
        usedPhylo.push_back(phylonum);
      }
    }

    // Now clean the map
    for (auto it = phylosMap.begin(); it != phylosMap.end();)
    {
      if (it->second.size() == 0)
      {
        phylosMap.erase(it++);
        continue;
      }

      vector<size_t>& vphyl = it->second;
      for (size_t i = vphyl.size(); i > 0; i--)
      {
        vector<size_t>::iterator posp = find(usedPhylo.begin(), usedPhylo.end(), vphyl[i - 1]);
        if (posp != usedPhylo.end())
          vphyl.erase(vphyl.begin() + static_cast<ptrdiff_t>(i - 1));
      }
      ++it;
    }
  }


  if (mPhylo->getNumbersOfPhyloLikelihoods().size() == 0)
  {
    if (warn)
      ApplicationTools::displayMessage("PhylogeneticsApplicationTools::getPhyloLikelihoodContainer : No phyloLikelihoods described");
    return mPhylo;
  }

  // get the result phylogeny => with number 0 in the
  // PhyloLikelihoodContainer

  ApplicationTools::displayMessage("");
  ApplicationTools::displayMessage("Result Phylolikelihood ");

  string sumAll;
  const vector<size_t>& nPhyl = mPhylo->getNumbersOfPhyloLikelihoods();

  for (size_t i = 0; i < nPhyl.size(); i++)
  {
    if (i != 0)
      sumAll += " + ";

    sumAll += "phylo" + TextTools::toString(nPhyl[i]);
  }

  string resultDesc = ApplicationTools::getStringParameter("result", params, sumAll);

  // check if really formula, or previous phylo

  std::shared_ptr<PhyloLikelihoodInterface> nPL;
  size_t nP(0);
  bool flag(resultDesc.substr(0, 5) == "phylo");

  if (flag)
  {
    try
    {
      nP = TextTools::to<size_t>(resultDesc.substr(5));
    }
    catch (Exception& e)
    {
      flag = false;
    }
  }

  if (!flag)
  {
    nPL = make_shared<PhyloLikelihoodFormula>(context, mPhylo, resultDesc);
    if (verbose)
      ApplicationTools::displayResult(" Result", dynamic_cast<PhyloLikelihoodFormula*>(nPL.get())->output());
  }
  else
  {
    if (!mPhylo->hasPhyloLikelihood(nP))
      throw BadIntegerException("Unknown Phylolikelihood number for result", (int)nP);
    else
      nPL = mPhylo->getPhyloLikelihood(nP);
    if (verbose)
      ApplicationTools::displayResult(" Result", resultDesc);
  }

  mPhylo->addPhyloLikelihood(0, nPL);

  return mPhylo;
}


/******************************************************/
/*** DISTRIBUTIONS ********************************/
/******************************************************/


MultipleDiscreteDistribution* PhylogeneticsApplicationTools::getMultipleDistributionDefaultInstance(
    const string& distDescription,
    map<string, string>& unparsedParameterValues,
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

unique_ptr<DiscreteDistributionInterface> PhylogeneticsApplicationTools::getRateDistribution(
    const map<string, string>& params,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose)
{
  string distDescription = ApplicationTools::getStringParameter("rate_distribution", params, "Constant()", suffix, suffixIsOptional);

  string distName;
  map<string, string> args;
  KeyvalTools::parseProcedure(distDescription, distName, args);

  BppORateDistributionFormat bIO(true);
  unique_ptr<DiscreteDistributionInterface> rDist(bIO.readDiscreteDistribution(distDescription, true));

  if (verbose)
  {
    ApplicationTools::displayResult("Rate distribution", distName);
    ApplicationTools::displayResult("Number of classes", TextTools::toString(rDist->getNumberOfCategories()));
  }

  return rDist;
}


/*************************************************************/
/*****  OPTIMIZATORS *****************************************/
/*************************************************************/

std::shared_ptr<PhyloLikelihoodInterface> PhylogeneticsApplicationTools::optimizeParameters(
    std::shared_ptr<PhyloLikelihoodInterface> lik,
    const map<string, string>& params,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    int warn)
{
  OptimizationTools::OptimizationOptions optopt(lik, params,  suffix, suffixIsOptional, verbose, warn);
  
  if (optopt.optMethodModel == "None")
    return lik;


  unsigned int n = 0;

  if ((optopt.optMethodModel == OptimizationTools::OPTIMIZATION_BRENT) || (optopt.optMethodModel == OptimizationTools::OPTIMIZATION_BFGS))
  {
    if (verbose && optopt.nstep > 1)
      ApplicationTools::displayResult("# of precision steps", TextTools::toString(optopt.nstep));
    
    optopt.parameters.matchParametersValues(lik->getParameters());
    n = OptimizationTools::optimizeNumericalParameters(lik, optopt);
  }
  else if (optopt.optMethodModel == "FullD")
  {
    // Uses Newton-raphson algorithm with numerical derivatives when required.
    optopt.parameters.matchParametersValues(lik->getParameters());
    if (dynamic_pointer_cast<SingleProcessPhyloLikelihood>(lik))
      n = OptimizationTools::optimizeNumericalParameters2(
            dynamic_pointer_cast<SingleProcessPhyloLikelihood>(lik),
            optopt);
    else
      n = OptimizationTools::optimizeNumericalParameters2(lik, optopt);
  }
  else
    throw Exception("Unknown optimization method: " + optopt.optMethodModel);

  string finalMethod = ApplicationTools::getStringParameter("optimization.final", params, "none", suffix, suffixIsOptional, warn + 1);
  unique_ptr<OptimizerInterface> finalOptimizer = nullptr;
  if (finalMethod == "none")
  {}
  else if (finalMethod == "simplex")
  {
    finalOptimizer = make_unique<DownhillSimplexMethod>(lik);
  }
  else if (finalMethod == "powell")
  {
    finalOptimizer = make_unique<PowellMultiDimensions>(lik);
  }
  else
    throw Exception("Unknown final optimization method: " + finalMethod);

  if (finalOptimizer)
  {
    optopt.parameters.matchParametersValues(lik->getParameters());
    if (verbose)
      ApplicationTools::displayResult("Final optimization step", finalMethod);
    finalOptimizer->setProfiler(optopt.profiler);
    finalOptimizer->setMessageHandler(optopt.messenger);
    finalOptimizer->setMaximumNumberOfEvaluations(optopt.nbEvalMax);
    finalOptimizer->getStopCondition()->setTolerance(optopt.tolerance);
    finalOptimizer->setVerbose(verbose);
    finalOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    finalOptimizer->init(optopt.parameters);
    finalOptimizer->optimize();
    n += finalOptimizer->getNumberOfEvaluations();
  }

  if (verbose)
    ApplicationTools::displayResult("Performed", TextTools::toString(n) + " function evaluations.");
  if (optopt.backupFile != "none")
  {
    string bf = optopt.backupFile + ".def";
    rename(optopt.backupFile.c_str(), bf.c_str());
  }
  return lik;
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
/**************** Output ************************************/
/******************************************************************************/

void PhylogeneticsApplicationTools::writeTree(
    const TreeTemplate<Node>& tree,
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
    treeWriter->writeTree(tree, file, true);
  delete treeWriter;
  if (verbose)
    ApplicationTools::displayResult("Wrote tree to file ", file);
}

/******************************************************************************/

void PhylogeneticsApplicationTools::writePhyloTrees(
    const vector<const PhyloTree*>& trees,
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
  OMultiPhyloTree* treeWriter;
  if (format == "Newick")
    treeWriter = new Newick();
  else if (format == "Nexus")
    treeWriter = new NexusIOTree();
  else if (format == "NHX")
    treeWriter = new Nhx();
  else
    throw Exception("Unknown format for tree writing: " + format);

  if (!checkOnly)
  {
    treeWriter->writePhyloTrees(trees, file, true);

    if (verbose)
      ApplicationTools::displayResult("Wrote trees to file ", file);
  }

  delete treeWriter;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::writePhyloTrees(
    const SubstitutionProcessCollection& spc,
    const map<string, string>& params,
    const string& prefix,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose,
    bool checkOnly,
    bool withIds,
    int warn)
{
  string format = ApplicationTools::getStringParameter(prefix + "tree.format", params, "Newick", suffix, suffixIsOptional, warn + 1);
  string file = ApplicationTools::getAFilePath(prefix + "tree.file", params, true, false, suffix, suffixIsOptional);

  OPhyloTree* treeWriter;
  if (format == "Newick")
    treeWriter = new Newick();
  else if (format == "Nexus")
    treeWriter = new NexusIOTree();
  else if (format == "NHX")
    treeWriter = new Nhx();
  else
    throw Exception("Unknown format for tree writing: " + format);

  if (!checkOnly)
  {
    vector<size_t> vTN = spc.getTreeNumbers();

    for (size_t i = 0; i < vTN.size(); i++)
    {
      auto tree = spc.getTree(vTN[i]);

      std::vector<shared_ptr<PhyloNode>> nodes = tree->getAllNodes();

      for (auto& node : nodes)
      {
        if (tree->isLeaf(node) && withIds)
          node->setName(TextTools::toString(tree->getNodeIndex(node)) + "_" + node->getName());
        else
          node->setProperty("NodeId", BppString(TextTools::toString(tree->getNodeIndex(node))));
      }

      Newick* nt = dynamic_cast<Newick*>(treeWriter);
      if (nt)
        nt->enableExtendedBootstrapProperty("NodeId");

      treeWriter->writePhyloTree(*tree, file + "_" + TextTools::toString(vTN[i]), true);
    }
    if (verbose)
      ApplicationTools::displayResult("Wrote trees to files : ", file + "_...");
  }

  delete treeWriter;
}

void PhylogeneticsApplicationTools::printParameters(const BranchModelInterface& model, OutputStream& out, int warn)
{
  out << "model=";
  map<string, string> globalAliases;
  vector<string> writtenNames;
  BppOBranchModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, false, warn);
  bIO.write(model, out, globalAliases, writtenNames);
  out.endLine();
}

void PhylogeneticsApplicationTools::printParameters(
    const SubstitutionProcessInterface& process,
    OutputStream& out,
    int warn)
{
  try
  {
    auto& sp = dynamic_cast<const SimpleSubstitutionProcess&>(process);
    (out << "nonhomogeneous=no").endLine();

    out << "model=";
    map<string, string> globalAliases;
    vector<string> writtenNames;
    BppOBranchModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, false, warn);
    bIO.write(sp.model(0, 0), out, globalAliases, writtenNames);
    out.endLine();
    return;
  }
  catch (bad_cast& e)
  {}

  try
  {
    const RateAcrossSitesSubstitutionProcess& pRA = dynamic_cast<const RateAcrossSitesSubstitutionProcess&>(process);

    (out << "nonhomogeneous=no").endLine();

    out << "model=";
    map<string, string> globalAliases;
    vector<string> writtenNames;
    const BppOBranchModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, false, warn);
    bIO.write(process.model(0, 0), out, globalAliases, writtenNames);
    out.endLine();
    out.endLine();

    // Rate distribution

    out << "rate_distribution=";
    const BppORateDistributionFormat bIOR(true);
    bIOR.writeDiscreteDistribution(*pRA.getRateDistribution(), out, globalAliases, writtenNames);
    out.endLine();
    return;
  }
  catch (bad_cast& e)
  {}

  try
  {
    const NonHomogeneousSubstitutionProcess& pNH = dynamic_cast<const NonHomogeneousSubstitutionProcess&>(process);

    (out << "nonhomogeneous=general").endLine();
    (out << "nonhomogeneous.number_of_models=" << pNH.getNumberOfModels()).endLine();

    vector<string> writtenNames;

    // Loop over all models:
    for (size_t i = 0; i < pNH.getNumberOfModels(); ++i)
    {
      const auto model = pNH.getModel(i);

      // First get the aliases for this model:
      map<string, string> aliases;

      ParameterList pl = model->getParameters();

      for (size_t np = 0; np < pl.size(); ++np)
      {
        string nfrom = pNH.getFrom(pl[np].getName() + "_" + TextTools::toString(i + 1));
        if (nfrom != "")
          aliases[pl[np].getName()] = nfrom;
      }

      // Now print it:
      writtenNames.clear();
      out.endLine() << "model" << (i + 1) << "=";
      BppOBranchModelFormat bIOsm(BppOSubstitutionModelFormat::ALL, true, true, true, false, warn);
      bIOsm.write(*model, out, aliases, writtenNames);
      out.endLine();
      vector<unsigned int> ids = pNH.getNodesWithModel(i);
      out << "model" << (i + 1) << ".nodes_id=" << ids[0];
      for (size_t j = 1; j < ids.size(); ++j)
      {
        out << "," << ids[j];
      }
      out.endLine();
    }

    // Root frequencies:
    out.endLine();
    if (pNH.getRootFrequencySet())
    {
      out << "nonhomogeneous.root_freq=";

      map<string, string> aliases;

      ParameterList pl = pNH.rootFrequencySet().getParameters();

      for (size_t np = 0; np < pl.size(); ++np)
      {
        string nfrom = pNH.getFrom(pl[np].getName());
        if (nfrom != "")
          aliases[pl[np].getName()] = nfrom;
      }

      BppOFrequencySetFormat bIO(BppOFrequencySetFormat::ALL, false, warn);
      bIO.writeFrequencySet(pNH.rootFrequencySet(), out, aliases, writtenNames);
    }
    else
      out << "nonhomogeneous.stationarity=true";
    out.endLine();

    // Rate distribution

    map<string, string> aliases;
    auto& pdd = pNH.rateDistribution();

    ParameterList pl = pdd.getParameters();
    for (size_t np = 0; np < pl.size(); ++np)
    {
      string nfrom = pNH.getFrom(pl[np].getName());
      if (nfrom != "")
        aliases[pl[np].getName()] = nfrom;
    }
    out.endLine();
    out << "rate_distribution=";
    const BppORateDistributionFormat bIO(true);
    bIO.writeDiscreteDistribution(pdd, out, aliases, writtenNames);
    out.endLine();
    return;
  }
  catch (bad_cast& e)
  {}

  return;
}

void PhylogeneticsApplicationTools::printParameters(const SubstitutionProcessCollection& collection, OutputStream& out, int warn, bool withAlias)
{
  vector<string> writtenNames;

  // The models
  vector<size_t> vModN = collection.getModelNumbers();

  for (auto modn : vModN)
  {
    const auto& model = *collection.getModel(modn);

    // First get the aliases for this model:
    map<string, string> aliases;

    if (withAlias)
    {
      ParameterList pl = model.getParameters();

      for (size_t np = 0; np < pl.size(); ++np)
      {
        string nfrom = collection.getFrom(pl[np].getName() + "_" + TextTools::toString(modn));
        if (nfrom != "")
          aliases[pl[np].getName()] = nfrom;
      }
    }

    // Now print it:
    writtenNames.clear();
    out << "model" << modn << "=";
    BppOBranchModelFormat bIOsm(BppOSubstitutionModelFormat::ALL, true, true, true, false, warn);
    bIOsm.write(model, out, aliases, writtenNames);
    out.endLine();
    out.endLine();
  }

  // Root frequencies:
  vector<size_t> rootFreqN = collection.getFrequenciesNumbers();

  for (size_t i = 0; i < rootFreqN.size(); ++i)
  {
    auto rootFreq = collection.getFrequencySet(rootFreqN[i]);

    // Now print it:
    writtenNames.clear();
    out.endLine() << "root_freq" << rootFreqN[i] << "=";
    BppOFrequencySetFormat bIOf(BppOFrequencySetFormat::ALL, true, warn);

    map<string, string> aliases;

    if (withAlias)
    {
      ParameterList pl = rootFreq->getParameters();

      for (size_t np = 0; np < pl.size(); ++np)
      {
        string nfrom = collection.getFrom(pl[np].getName() + "_" + TextTools::toString(rootFreqN[i]));
        if (nfrom != "")
          aliases[pl[np].getName()] = nfrom;
      }
    }

    bIOf.writeFrequencySet(*rootFreq, out, aliases, writtenNames);
    out.endLine();
  }

  // Rate distribution

  vector<size_t> vDistN = collection.getRateDistributionNumbers();

  for (auto distn : vDistN)
  {
    if (distn < 10000)
    {
      auto dist = collection.getRateDistribution(distn);

      // First get the aliases for this model:
      map<string, string> aliases;

      if (withAlias)
      {
        ParameterList pl = dist->getParameters();

        for (size_t np = 0; np < pl.size(); ++np)
        {
          string nfrom = collection.getFrom(pl[np].getName() + "_" + TextTools::toString(distn));
          if (nfrom != "")
            aliases[pl[np].getName()] = nfrom;
        }
      }

      // Now print it:
      writtenNames.clear();
      out.endLine() << "rate_distribution" << distn << "=";
      BppORateDistributionFormat bIOd(true);
      bIOd.writeDiscreteDistribution(*dist, out, aliases, writtenNames);
      out.endLine();
    }
  }

  // scenarios

  vector<size_t> vSce = collection.getScenarioNumbers();

  if (vSce.size() > 0)
    out.endLine();

  vector<const ModelPath*> vMP;

  // first output the scenarios
  for (const auto& scennum : vSce)
  {
    const auto scen = collection.getModelScenario(scennum);

    out.endLine();

    out << "scenario" << scennum << "=";

    size_t nbMP = scen->getNumberOfModelPaths();

    for (size_t mpn = 0; mpn < nbMP; mpn++)
    {
      const auto& mp = scen->getModelPath(mpn);

      auto itmp = find(vMP.begin(), vMP.end(), mp.get());
      auto inmp = std::distance(vMP.begin(), itmp);
      if (itmp == vMP.end())
        vMP.push_back(mp.get());

      if (mpn != 0)
        out << "&";
      out << "path" << TextTools::toString(inmp + 1);
    }
    out.endLine();
  }

  // then the model path
  for (size_t inmp = 0; inmp < vMP.size(); inmp++)
  {
    out.endLine();
    out << "path" << inmp + 1 << "=";

    const ModelPath& mp = *vMP[inmp];

    auto vMod = mp.getModels();

    bool dem = true;
    for (const auto& mod:vMod)
    {
      // look for model number in collection
      size_t modN = collection.getModelIndex(mod);

      if (!dem)
        out << "&";

      out << "model" << modN;
      out << "[" << mp.getPathNode(mod).to_string() <<  "]";
      dem = false;
    }
    out.endLine();
  }

  // processes
  out.endLine();

  vector<size_t> vprocN = collection.getSubstitutionProcessNumbers();

  for (size_t i = 0; i < vprocN.size(); ++i)
  {
    const auto& spcm = collection.substitutionProcess(vprocN[i]);

    out << "process" << vprocN[i] << "=";

    if (spcm.getNumberOfModels() == 1)
      out << "Homogeneous(model=" << spcm.getModelNumbers()[0];
    else
    {
      out << "Nonhomogeneous(";
      vector<size_t> vMN = spcm.getModelNumbers();
      for (size_t j = 0; j < vMN.size(); ++j)
      {
        if (j != 0)
          out << ",";

        out << "model" << (j + 1) << "=" << vMN[j];
        out << ",";

        vector<unsigned int> ids = spcm.getNodesWithModel(vMN[j]);
        out << "model" << (j + 1) << ".nodes_id=(" << ids[0];
        for (size_t k = 1; k < ids.size(); ++k)
        {
          out << "," << ids[k];
        }
        out << ")";
      }
    }

    out << ", tree=" << spcm.getTreeNumber();

    out << ", rate=";
    size_t dN = spcm.getRateDistributionNumber();

    if (dN < 10000)
      out << dN;
    else
      out << size_t(dN / 10000 - 1) << "." << dN % 10000;

    if (spcm.getRootFrequencySet())
      out << ", root_freq=" << spcm.getRootFrequenciesNumber();

    if (spcm.getModelScenario())
      out << ", scenario=" << spcm.getModelScenarioNumber();

    out << ")";
    out.endLine();
    out.endLine();
  }
}


void PhylogeneticsApplicationTools::printParameters(const PhyloLikelihoodContainer& phylocont, OutputStream& out, int warn)
{
  out << "# Log likelihood = ";

  std::shared_ptr<const PhyloLikelihoodInterface> result = phylocont[0];

  if (!result)
  {
    out << "Nan";
    out.endLine();
    return;
  }

  out.setPrecision(20) << (-result->getValue());
  out.endLine();
  out.endLine();


  // First output result
  out << "result=";

  auto pop = dynamic_pointer_cast<const PhyloLikelihoodFormula>(result);

  vector<size_t> phyldep;

  if (!pop)
  {
    out << "phylo1";
    phyldep.push_back(1);
  }
  else
  {
    string popout = pop->output();

    out << popout;

    StringTokenizer st(popout, "phylo", true, true);
    st.nextToken();


    while (st.hasMoreToken())
    {
      string ex = st.nextToken();
      phyldep.push_back((size_t)(atoi(ex.c_str())));
    }
  }

  out.endLine();
  out.endLine();

  // Then the other phylolikelihoods

  while (phyldep.size() != 0)
  {
    size_t num = phyldep[0];
    std::shared_ptr<const PhyloLikelihoodInterface> phyloLike = phylocont[num];

    // remove phylolikelihoods with this number
    auto itf = find(phyldep.begin(), phyldep.end(), num);
    while (itf != phyldep.end())
    {
      phyldep.erase(itf);
      itf = find(itf, phyldep.end(), num);
    }


    // then output

    if (dynamic_pointer_cast<const SingleDataPhyloLikelihoodInterface>(phyloLike))
      printParameters(*dynamic_pointer_cast<const SingleDataPhyloLikelihoodInterface>(phyloLike), out, num, warn);
    else
    {
      out << "phylo" << num << "=";

      auto mDP = dynamic_pointer_cast<const PhyloLikelihoodSetInterface>(phyloLike);
      if (mDP)
      {
        if (dynamic_pointer_cast<const AlignedPhyloLikelihoodMixture>(phyloLike))
        {
          auto pM = dynamic_pointer_cast<const AlignedPhyloLikelihoodMixture>(phyloLike);

          out << "Mixture(probas=(" << VectorTools::paste(pM->getPhyloProbabilities(), ",");

          out << "),";
        }

        else if (dynamic_pointer_cast<const AlignedPhyloLikelihoodHmm>(phyloLike))
        {
          auto pM = dynamic_pointer_cast<const AlignedPhyloLikelihoodHmm>(phyloLike);
          out << "HMM(probas=";

          RowMatrix<double> tMt;
          copyEigenToBpp(pM->getHmmTransitionMatrix(), tMt);
          MatrixTools::print(tMt, out);

          out << ",";
        }

        else if (dynamic_pointer_cast<const AlignedPhyloLikelihoodAutoCorrelation>(phyloLike))
        {
          auto pM = dynamic_pointer_cast<const AlignedPhyloLikelihoodAutoCorrelation>(phyloLike);

          out << "AutoCorr(lambdas=(";

          Vdouble vP;
          for (unsigned int i = 0; i < pM->getHmmTransitionMatrix().cols(); ++i)
          {
            vP.push_back(pM->getHmmTransitionMatrix()(i, i));
          }

          out << VectorTools::paste(vP, ",");

          out << "),";
        }
        else if (dynamic_pointer_cast<const AlignedPhyloLikelihoodProduct>(phyloLike))
        {
          out << "Product(";
        }
        else
          throw Exception("PhylogeneticsApplicationTools::printParameters - unknown phylolikelihood type : phylo " + TextTools::toString(num));


        vector<size_t> vPN = mDP->getNumbersOfPhyloLikelihoods();

        for (size_t i = 0; i < vPN.size(); i++)
        {
          out << "phylo" << i + 1 << "=" << vPN[i];
          if (i != vPN.size() - 1)
            out << ",";
        }

        out << ")";

        // update phyldep
        for (size_t i = 0; i < vPN.size(); i++)
        {
          if (find(phyldep.begin(), phyldep.end(), vPN[i]) == phyldep.end())
            phyldep.push_back(vPN[i]);
        }
      }
      out.endLine();
    }
    out.endLine();
  }
  out.endLine();
}


void PhylogeneticsApplicationTools::printParameters(const SingleDataPhyloLikelihoodInterface& phyloLike, OutputStream& out, size_t nPhylo, int warn)
{
  out << "phylo" << TextTools::toString(nPhylo) << "=";

  out << "Single(";

  try
  {
    auto& pMP = dynamic_cast<const SequencePhyloLikelihoodInterface&>(phyloLike);
    out << "process=" << pMP.getSequenceEvolutionNumber();
    goto finish;
  }
  catch (bad_cast& e)
  {}

  try
  {
    auto& pS = dynamic_cast<const SingleProcessPhyloLikelihood&>(phyloLike);
    out << "process=" << pS.getSubstitutionProcessNumber();
  }
  catch (bad_cast& e)
  {}

finish:
  out << ", data=" << TextTools::toString(phyloLike.getNData()) << ")";
  out.endLine();
}

void PhylogeneticsApplicationTools::printParameters(const SequenceEvolution& evol, OutputStream& out, size_t nEvol, int warn)
{
  out << "process" << TextTools::toString(nEvol) << "=";

  if (dynamic_cast<const OneProcessSequenceEvolution*>(&evol))
  {
    const OneProcessSequenceEvolution* pOP = dynamic_cast<const OneProcessSequenceEvolution*>(&evol);

    out << "Simple(process=" <<  pOP->getSubstitutionProcessNumber() << ")";
  }
  else if (dynamic_cast<const MultiProcessSequenceEvolution*>(&evol))
  {
    const MultiProcessSequenceEvolution* pMP = dynamic_cast<const MultiProcessSequenceEvolution*>(&evol);

    if (dynamic_cast<const MixtureSequenceEvolution*>(&evol))
    {
      const MixtureSequenceEvolution* pM = dynamic_cast<const MixtureSequenceEvolution*>(&evol);

      out << "Mixture(probas=(" << VectorTools::paste(pM->getSubProcessProbabilities(), ",");
      out << "),";
    }

    else if (dynamic_cast<const HmmSequenceEvolution*>(&evol))
    {
      const HmmSequenceEvolution* pM = dynamic_cast<const HmmSequenceEvolution*>(&evol);
      out << "HMM(probas=";

      const Matrix<double>& tMt = pM->hmmTransitionMatrix().getPij();
      MatrixTools::print(tMt, out);

      out << ",";
    }
    else if (dynamic_cast<const AutoCorrelationSequenceEvolution*>(&evol))
    {
      const AutoCorrelationSequenceEvolution* pM = dynamic_cast<const AutoCorrelationSequenceEvolution*>(&evol);

      out << "AutoCorr(lambdas=(";

      Vdouble vP;
      for (unsigned int i = 0; i < pM->getNumberOfSubstitutionProcess(); i++)
      {
        vP.push_back(pM->hmmTransitionMatrix().Pij(i, i));
      }

      out << VectorTools::paste(vP, ",");

      out << "),";
    }
    else if (dynamic_cast<const PartitionSequenceEvolution*>(&evol))
    {
      const PartitionSequenceEvolution* pM = dynamic_cast<const PartitionSequenceEvolution*>(&evol);

      out << "Partition(";

      const map<size_t, vector<size_t>>& mProcPos = pM->mapOfProcessSites();

      vector<size_t> vP = pMP->getSubstitutionProcessNumbers();

      for (unsigned int i = 0; i < vP.size(); i++)
      {
        out << "process" << TextTools::toString(i + 1) << ".sites=";

        vector<size_t> v = mProcPos.find(vP[i])->second + 1;

        if (v.size() > 1)
          out << "(";

        VectorTools::printRange(v, out, ",", ":");

        if (v.size() > 1)
          out << ")";

        out << ",";
      }
    }

    vector<size_t> vPN = pMP->getSubstitutionProcessNumbers();

    for (size_t i = 0; i < vPN.size(); i++)
    {
      out << "process" << i + 1 << "=" << vPN[i];
      if (i != vPN.size() - 1)
        out << ",";
    }

    out << ")";
  }

  out.endLine();
}


// ///////////////////////////////////////////////////////
// Analysis Information
// //////////////////////////////////////////////////////


void PhylogeneticsApplicationTools::printAnalysisInformation(const PhyloLikelihoodContainer& phylocont, const string& infosFile, int warn)
{
  std::shared_ptr<const PhyloLikelihoodInterface> result = phylocont[0];

  if (!result)
    return;

  vector<size_t> phyldep = phylocont.getNumbersOfPhyloLikelihoods();

  while (phyldep.size() != 0)
  {
    size_t num = phyldep[0];

    phyldep.erase(phyldep.begin());

    std::shared_ptr<const PhyloLikelihoodInterface> phyloLike = phylocont[num];
    // output

    string info_out = infosFile + "_" + TextTools::toString(num);

    if (dynamic_pointer_cast<const SingleDataPhyloLikelihoodInterface>(phyloLike) && num != 0)
      printAnalysisInformation(*dynamic_pointer_cast<const SingleDataPhyloLikelihoodInterface>(phyloLike), info_out, warn);
    else
    {
      auto sOAP = dynamic_pointer_cast<const AlignedPhyloLikelihoodSetInterface>(phyloLike);
      if (sOAP)
      {
        if (num != 0)
          printAnalysisInformation(*sOAP, info_out, warn);

        vector<size_t> vPN = sOAP->getNumbersOfPhyloLikelihoods();

        // update phyldep
        phyldep.assign(vPN.begin(), vPN.end());
        phyldep = VectorTools::unique(phyldep);
      }
      else
      {
        auto sOAB = dynamic_pointer_cast<const PhyloLikelihoodSetInterface>(phyloLike);
        if (sOAB)
        {
          vector<size_t> vPN = sOAB->getNumbersOfPhyloLikelihoods();

          // update phyldep
          phyldep.assign(vPN.begin(), vPN.end());
          phyldep = VectorTools::unique(phyldep);
        }
      }
    }
  }
}


void PhylogeneticsApplicationTools::printAnalysisInformation(const AlignedPhyloLikelihoodSetInterface& sOAP, const string& infosFile, int warn)
{
  const AlignedPhyloLikelihoodMixture* mOAP = nullptr;
  const AlignedPhyloLikelihoodHmm* hOAP = nullptr;
  const AlignedPhyloLikelihoodAutoCorrelation* aCOAP = nullptr;

  vector<size_t> phyloNum = sOAP.getNumbersOfPhyloLikelihoods();
  size_t nbP = phyloNum.size();

  if (dynamic_cast<const AlignedPhyloLikelihoodProduct*>(&sOAP) == nullptr)
  {
    StlOutputStream out(make_unique<ofstream>(infosFile.c_str(), ios::out));

    if (dynamic_cast<const AlignedPhyloLikelihoodMixture*>(&sOAP) != nullptr)
      mOAP = dynamic_cast<const AlignedPhyloLikelihoodMixture*>(&sOAP);
    else if (dynamic_cast<const AlignedPhyloLikelihoodHmm*>(&sOAP) != nullptr)
      hOAP = dynamic_cast<const AlignedPhyloLikelihoodHmm*>(&sOAP);
    else if (dynamic_cast<const AlignedPhyloLikelihoodAutoCorrelation*>(&sOAP) != nullptr)
      aCOAP = dynamic_cast<const AlignedPhyloLikelihoodAutoCorrelation*>(&sOAP);

    vector<string> colNames;
    colNames.push_back("Sites");
    colNames.push_back("lnL");

    for (size_t i = 0; i < nbP; i++)
    {
      colNames.push_back("lnL_phylo_" + TextTools::toString(phyloNum[i]));
    }
    for (size_t i = 0; i < nbP; i++)
    {
      colNames.push_back("prob_phylo_" + TextTools::toString(phyloNum[i]));
    }


    vector<string> row(2 + (nbP > 1 ? 2 * nbP : 0));
    DataTable* infos = new DataTable(colNames);

    Vdouble vap(0);
    Vdouble vlog(nbP);

    if (mOAP)
      vap = mOAP->getPhyloProbabilities();

    size_t nSites = sOAP.getNumberOfSites();

    for (size_t i = 0; i < nSites; i++)
    {
      row[0] = (string("[" + TextTools::toString(i) + "]"));
      row[1] = TextTools::toString(sOAP.getLogLikelihoodForASite(i));

      if (nbP > 1)
      {
        for (size_t j = 0; j < nbP; j++)
        {
          vlog[j] = sOAP.getLogLikelihoodForASiteForAPhyloLikelihood(i, phyloNum[j]);
        }

        for (size_t j = 0; j < nbP; j++)
        {
          row[2 + j] = TextTools::toString(vlog[j]);
        }

        if (mOAP)
        {
          double sum = VectorTools::sumExp(vlog, vap);
          for (size_t j = 0; j < nbP; j++)
          {
            row[2 + nbP + j] = TextTools::toString(exp(vlog[j]) * vap[j] / sum);
          }
        }
        else
        {
          if (hOAP)
            vap = hOAP->getPosteriorProbabilitiesForASitePerAligned(i);
          else if (aCOAP)
            vap = aCOAP->getPosteriorProbabilitiesForASitePerAligned(i);

          for (size_t j = 0; j < vap.size(); j++)
          {
            row[2 + nbP + j] = TextTools::toString(vap[j]);
          }
        }
      }

      infos->addRow(row);
    }

    DataTable::write(*infos, out, "\t");
    delete infos;
  }
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printAnalysisInformation(
    const SingleDataPhyloLikelihoodInterface& phyloLike,
    const string& infosFile,
    int warn)
{
  try
  {
    auto& pSPL = dynamic_cast<const SingleProcessPhyloLikelihood&>(phyloLike);

    StlOutputStream out(make_unique<ofstream>(infosFile.c_str(), ios::out));

    std::shared_ptr<const SubstitutionProcessInterface> pSP = pSPL.getSubstitutionProcess();

    vector<string> colNames;
    colNames.push_back("Sites");
    colNames.push_back("is.complete");
    colNames.push_back("is.constant");
    colNames.push_back("lnL");

    auto pDD = pSP->getRateDistribution();
    size_t nbR = 0;

    if (pDD)
    {
      nbR = pDD->getNumberOfCategories();

      if (nbR > 1)
        for (size_t i = 0; i < nbR; ++i)
        {
          colNames.push_back("Pr_rate=" + TextTools::toString(pDD->getCategory(i)));
        }
    }
    colNames.push_back("rc");
    colNames.push_back("pr");

    std::shared_ptr<const AlignmentDataInterface> sites = phyloLike.getData();

    vector<string> row(6 + (nbR > 1 ? nbR : 0));
    auto infos = make_unique<DataTable>(colNames);

    VVdouble vvPP(pSPL.getPosteriorProbabilitiesPerSitePerClass());

    for (size_t i = 0; i < sites->getNumberOfSites(); ++i)
    {
      double lnL = phyloLike.getLogLikelihoodForASite(i);

      const CoreSiteInterface& currentSite = sites->site(i);
      int currentSiteCoordinate = currentSite.getCoordinate();
      string isCompl = "NA";
      string isConst = "NA";
      try
      {
        isCompl = (SymbolListTools::isComplete(currentSite) ? "1" : "0");
      }
      catch (EmptySiteException& ex)
      {}
      try
      {
        isConst = (SymbolListTools::isConstant(currentSite) ? "1" : "0");
      }
      catch (EmptySiteException& ex)
      {}
      row[0] = (string("[" + TextTools::toString(currentSiteCoordinate) + "]"));
      row[1] = isCompl;
      row[2] = isConst;
      row[3] = TextTools::toString(lnL);

      if (nbR > 1)
      {
        double pr = 0;
        for (size_t j = 0; j < nbR; ++j)
        {
          row[4 + j] = TextTools::toString(vvPP[i][j]);
          pr += vvPP[i][j] * pDD->getCategory(j);
        }

        row[4 + nbR] = TextTools::toString(VectorTools::whichMax(vvPP[i]) + 1);
        row[5 + nbR] = TextTools::toString(pr);
      }
      else
      {
        row[4] = "1";
        row[5] = "1";
      }

      infos->addRow(row);
    }

    DataTable::write(*infos, out, "\t");
    return;
  }
  catch (bad_cast& e)
  {}

  try
  {
    auto& pPPL = dynamic_cast<const PartitionProcessPhyloLikelihood&>(phyloLike);

    auto& pSE = dynamic_cast<const PartitionSequenceEvolution&>(pPPL.sequenceEvolution());

    const map<size_t, vector<size_t>>& mProcPos = pSE.mapOfProcessSites();

    vector<size_t> nbProc = pSE.getSubstitutionProcessNumbers();

    map<size_t, size_t> mNbr;

    for (auto nP : nbProc)
    {
      auto& sp = pSE.substitutionProcess(nP);
      auto pDD = sp.getRateDistribution();
      mNbr[nP] = (pDD ? pDD->getNumberOfCategories() : 1);
    }

    size_t maxR = max_element(mNbr.begin(), mNbr.end(), [](const std::pair<size_t, size_t>& p1, const std::pair<size_t, size_t>& p2){
      return p1.second < p2.second;
    })->second;

    StlOutputStream out(make_unique<ofstream>(infosFile.c_str(), ios::out));

    vector<string> colNames;
    colNames.push_back("Sites");
    colNames.push_back("is.complete");
    colNames.push_back("is.constant");
    colNames.push_back("lnL");

    if (maxR > 1)
      for (size_t i = 0; i < maxR; i++)
      {
        colNames.push_back("prob" + TextTools::toString(i + 1));
      }

    size_t nbSites = pSE.getNumberOfSites();

    vector<string> row(4 + (maxR > 1 ? maxR : 0));
    auto infos = make_unique<DataTable>(nbSites, colNames);
    for (auto nP : nbProc)
    {
      auto pSPPL = dynamic_pointer_cast<const SingleProcessPhyloLikelihood>(pPPL.getPhyloLikelihood(nP));

      if (!pSPPL)
        throw Exception("PhylogeneticsApplicationTools::printAnalysisInformation : no SingleProcessPhyloLikelihood in PartitionProcessPhyloLikelihood.");

      size_t nbr = mNbr[pSPPL->getSubstitutionProcessNumber()];

      const vector<size_t>& mPos = mProcPos.at(nP);

      auto sites = pSPPL->getData();

      for (size_t i = 0; i < sites->getNumberOfSites(); ++i)
      {
        double lnL = pSPPL->getLogLikelihoodForASite(i);

        const CoreSiteInterface& currentSite = sites->site(i);
        int currentSiteCoordinate = currentSite.getCoordinate();
        string isCompl = "NA";
        string isConst = "NA";
        try
        {
          isCompl = (SymbolListTools::isComplete(currentSite) ? "1" : "0");
        }
        catch (EmptySiteException& ex)
        {}
        try
        {
          isConst = (SymbolListTools::isConstant(currentSite) ? "1" : "0");
        }
        catch (EmptySiteException& ex)
        {}
        row[0] = (string("[" + TextTools::toString(currentSiteCoordinate) + "]"));
        row[1] = isCompl;
        row[2] = isConst;
        row[3] = TextTools::toString(lnL);

        if (nbr > 1)
        {
          Vdouble vPP = pSPPL->getPosteriorProbabilitiesForSitePerClass(i);

          for (size_t j = 0; j < nbr; ++j)
          {
            row[4 + j] = TextTools::toString(vPP[j]);
          }
        }

        for (size_t j = nbr; j < maxR; j++)
        {
          row[4 + j] = "NA";
        }

        infos->setRow(mPos[i], row);
      }
    }

    DataTable::write(*infos, out, "\t");
    return;
  }
  catch (bad_cast& e)
  {}

  try
  {
    auto& pMPL = dynamic_cast<const MultiProcessSequencePhyloLikelihood&>(phyloLike);

    StlOutputStream out(make_unique<ofstream>(infosFile.c_str(), ios::out));

    vector<string> colNames;
    colNames.push_back("Sites");
    colNames.push_back("is.complete");
    colNames.push_back("is.constant");
    colNames.push_back("lnL");

    size_t nbP = pMPL.getNumberOfSubstitutionProcess();

    if (nbP > 1)
    {
      for (size_t i = 0; i < nbP; ++i)
      {
        colNames.push_back("lnL" + TextTools::toString(i + 1));
      }
      for (size_t i = 0; i < nbP; ++i)
      {
        colNames.push_back("prob" + TextTools::toString(i + 1));
      }
    }

    auto sites = phyloLike.getData();

    vector<string> row(4 + (nbP > 1 ? 2 * nbP : 0));
    DataTable* infos = new DataTable(colNames);

    VVdouble vvPP = pMPL.getPosteriorProbabilitiesPerSitePerProcess();
    VVdouble vvL = pMPL.getLikelihoodPerSitePerProcess();

    for (size_t i = 0; i < sites->getNumberOfSites(); ++i)
    {
      double lnL = phyloLike.getLogLikelihoodForASite(i);
      const CoreSiteInterface& currentSite = sites->site(i);
      int currentSiteCoordinate = currentSite.getCoordinate();
      string isCompl = "NA";
      string isConst = "NA";
      try
      {
        isCompl = (SymbolListTools::isComplete(currentSite) ? "1" : "0");
      }
      catch (EmptySiteException& ex)
      {}
      try
      {
        isConst = (SymbolListTools::isConstant(currentSite) ? "1" : "0");
      }
      catch (EmptySiteException& ex)
      {}
      row[0] = (string("[" + TextTools::toString(currentSiteCoordinate) + "]"));
      row[1] = isCompl;
      row[2] = isConst;
      row[3] = TextTools::toString(lnL);

      if (nbP > 1)
      {
        for (size_t j = 0; j < nbP; j++)
        {
          row[4 + j] = TextTools::toString(std::log(vvL[i][j]));
        }
        for (size_t j = 0; j < nbP; j++)
        {
          row[4 + nbP + j] = TextTools::toString(vvPP[i][j]);
        }
      }
      infos->addRow(row);
    }

    DataTable::write(*infos, out, "\t");
    delete infos;
    return;
  }
  catch (bad_cast& e)
  {}

  return;
}


/******************************************************************************/

void PhylogeneticsApplicationTools::printParameters(const DiscreteDistributionInterface& rDist, OutputStream& out, bool withAlias)
{
  out << "rate_distribution=";
  map<string, string> globalAliases;
  vector<string> writtenNames;
  const BppORateDistributionFormat* bIO = new BppORateDistributionFormat(true);
  bIO->writeDiscreteDistribution(rDist, out, globalAliases, writtenNames);
  delete bIO;
  out.endLine();
}

/************************
* Substitution Mapping *
************************/
unique_ptr<SubstitutionCountInterface> PhylogeneticsApplicationTools::getSubstitutionCount(
    std::shared_ptr<const Alphabet> alphabet,
    std::shared_ptr<const SubstitutionModelInterface> model,
    const map<string, string>& params,
    string suffix,
    bool verbose,
    int warn)
{
  auto stateMap = model->getStateMap();

  unique_ptr<SubstitutionCountInterface> substitutionCount = nullptr;
  string nijtOption;
  map<string, string> nijtParams;
  string nijtText = ApplicationTools::getStringParameter("nijt", params, "Uniformization", suffix, true, warn);
  KeyvalTools::parseProcedure(nijtText, nijtOption, nijtParams);

  if (nijtOption == "Laplace")
  {
    size_t trunc = ApplicationTools::getParameter<size_t>("trunc", nijtParams, 10, suffix, true, warn + 1);
    substitutionCount = make_unique<LaplaceSubstitutionCount>(model, trunc);
  }
  else if (nijtOption == "Uniformization")
  {
    string weightOption = ApplicationTools::getStringParameter("weight", nijtParams, "None", "", true, warn + 1);
    shared_ptr<const AlphabetIndex2> weights(SequenceApplicationTools::getAlphabetIndex2(alphabet, weightOption, "Substitution weight scheme:"));
    string distanceOption = ApplicationTools::getStringParameter("distance", nijtParams, "None", "", true, warn + 1);
    shared_ptr<const AlphabetIndex2> distances(SequenceApplicationTools::getAlphabetIndex2(alphabet, distanceOption, "Substitution distances:"));
    substitutionCount = make_unique<UniformizationSubstitutionCount>(model, make_shared<TotalSubstitutionRegister>(stateMap), weights, distances);
  }
  else if (nijtOption == "Decomposition")
  {
    string weightOption = ApplicationTools::getStringParameter("weight", nijtParams, "None", "", true, warn + 1);
    shared_ptr<const AlphabetIndex2> weights(SequenceApplicationTools::getAlphabetIndex2(alphabet, weightOption, "Substitution weight scheme:"));
    string distanceOption = ApplicationTools::getStringParameter("distance", nijtParams, "None", "", true, warn + 1);
    shared_ptr<const AlphabetIndex2> distances(SequenceApplicationTools::getAlphabetIndex2(alphabet, distanceOption, "Substitution distances:"));
    auto revModel = dynamic_pointer_cast<const ReversibleSubstitutionModelInterface>(model);
    if (revModel)
      substitutionCount = make_unique<DecompositionSubstitutionCount>(revModel, make_shared<TotalSubstitutionRegister>(stateMap), weights, distances);
    else
      throw Exception("Decomposition method can only be used with reversible substitution models.");
  }
  else if (nijtOption == "Naive")
  {
    string weightOption = ApplicationTools::getStringParameter("weight", nijtParams, "None", "", true, warn + 1);
    std::shared_ptr<const AlphabetIndex2> weights(SequenceApplicationTools::getAlphabetIndex2(alphabet, weightOption, "Substitution weight scheme:"));
    string distanceOption = ApplicationTools::getStringParameter("distance", nijtParams, "", "", true, warn + 1);
    if (distanceOption != "")
      ApplicationTools::displayMessage("Naive substitution count: distances not handled");

    substitutionCount = make_unique<NaiveSubstitutionCount>(model, make_shared<TotalSubstitutionRegister>(stateMap), false, weights);
  }
  else if (nijtOption == "Label")
  {
    substitutionCount = make_unique<LabelSubstitutionCount>(model);
  }
  else if (nijtOption == "ProbOneJump")
  {
    substitutionCount = make_unique<OneJumpSubstitutionCount>(model);
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


unique_ptr<SubstitutionRegisterInterface> PhylogeneticsApplicationTools::getSubstitutionRegister(
    const string& regTypeDesc,
    std::shared_ptr<const StateMapInterface> stateMap,
    std::shared_ptr<const GeneticCode> genCode,
    std::shared_ptr<AlphabetIndex2>& weights,
    std::shared_ptr<AlphabetIndex2>& distances,
    bool verbose)
{
  string regType = "";
  map<string, string> regArgs;
  KeyvalTools::parseProcedure(regTypeDesc, regType, regArgs);

  auto alphabet = stateMap->getAlphabet();

  unique_ptr<SubstitutionRegisterInterface> reg = nullptr;
  weights = nullptr;
  distances = nullptr;

  string weightOption = ApplicationTools::getStringParameter("weight", regArgs, "None", "", true, 1);
  string distanceOption = ApplicationTools::getStringParameter("distance", regArgs, "None", "", true, 1);

  if (AlphabetTools::isCodonAlphabet(*alphabet))
  {
    weights = SequenceApplicationTools::getAlphabetIndex2(dynamic_pointer_cast<const CodonAlphabet>(alphabet), genCode, weightOption, "Substitution weight scheme:");
    distances = SequenceApplicationTools::getAlphabetIndex2(dynamic_pointer_cast<const CodonAlphabet>(alphabet), genCode, distanceOption, "Substitution distances:");
  }
  else
  {
    weights = SequenceApplicationTools::getAlphabetIndex2(alphabet, weightOption, "Substitution weight scheme:");
    distances = SequenceApplicationTools::getAlphabetIndex2(alphabet, distanceOption, "Substitution distances:");
  }

  if (regType == "Combination")
  {
    shared_ptr<AlphabetIndex2> w2;
    shared_ptr<AlphabetIndex2> d2;

    auto vreg = make_unique<VectorOfSubstitutionRegisters>(stateMap);

    size_t i = 0;
    while (++i)
    {
      string regDesc = ApplicationTools::getStringParameter("reg" + TextTools::toString(i), regArgs, "", "", false, 1);
      if (regDesc == "")
        break;

      auto sreg = getSubstitutionRegister(regDesc, stateMap, genCode, w2, d2);

      vreg->addRegister(std::move(sreg));
    }

    if (vreg->getNumberOfSubstitutionTypes()==0)
      throw Exception("Missing registers reg1, reg2, ... in description of Combination");

    reg = std::move(vreg);
  }
  else if (regType == "All")
  {
    reg = make_unique<ComprehensiveSubstitutionRegister>(stateMap, false);
  }
  else if (regType == "Total")
  {
    reg = make_unique<TotalSubstitutionRegister>(stateMap);
  }
  else if (regType == "Selected")
  {
    string subsList = ApplicationTools::getStringParameter("substitution.list", regArgs, "All", "", true, false);
    reg = make_unique<SelectedSubstitutionRegister>(stateMap, subsList);
  }


  // Alphabet dependent registers
  else if (AlphabetTools::isNucleicAlphabet(*alphabet))
  {
    if (regType == "GC")
      reg = make_unique<GCSubstitutionRegister>(stateMap, false);
    else if (regType == "GCw")
      reg = make_unique<GCSubstitutionRegister>(stateMap, true);
    else if (regType == "TsTv")
      reg = make_unique<TsTvSubstitutionRegister>(stateMap);
    else if (regType == "SW")
      reg = make_unique<SWSubstitutionRegister>(stateMap);
    else
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionRegister: unsupported substitution categorization:" + regType + " for alphabet " + alphabet->getAlphabetType());
  }
  else if (AlphabetTools::isCodonAlphabet(*alphabet))
  {
    if (regType == "IntraAA")
      reg = make_unique<AAInteriorSubstitutionRegister>(stateMap, genCode);
    else if (regType == "InterAA")
      reg = make_unique<AAExteriorSubstitutionRegister>(stateMap, genCode);
    else if (regType == "GC")
      reg = make_unique<GCSynonymousSubstitutionRegister>(stateMap, genCode);
    else if (regType == "GC1")
      reg = make_unique<GCPositionSubstitutionRegister>(stateMap, genCode, 0);
    else if (regType == "GC2")
      reg = make_unique<GCPositionSubstitutionRegister>(stateMap, genCode, 1);
    else if (regType == "GC3")
      reg = make_unique<GCPositionSubstitutionRegister>(stateMap, genCode, 2);
    else if (regType == "GCw")
      reg = make_unique<GCSynonymousSubstitutionRegister>(stateMap, genCode, true);
    else if (regType == "GC1w")
      reg = make_unique<GCPositionSubstitutionRegister>(stateMap, genCode, 0, true);
    else if (regType == "GC2w")
      reg = make_unique<GCPositionSubstitutionRegister>(stateMap, genCode, 1, true);
    else if (regType == "GC3w")
      reg = make_unique<GCPositionSubstitutionRegister>(stateMap, genCode, 2, true);
    else if (regType == "TsTv")
      reg = make_unique<TsTvSubstitutionRegister>(stateMap, genCode);
    else if (regType == "SW")
      reg = make_unique<SWSubstitutionRegister>(stateMap, genCode);
    else if (regType == "KrKc")
      reg = make_unique<KrKcSubstitutionRegister>(stateMap, genCode);
    else if (regType == "DnDs")
      reg = make_unique<DnDsSubstitutionRegister>(stateMap, genCode, false);
    else
      throw Exception("Unsupported substitution categorization: " + regType + " for alphabet " + alphabet->getAlphabetType());
  }

  else if (AlphabetTools::isProteicAlphabet(*alphabet))
  {
    if (regType == "KrKc")
      reg = make_unique<KrKcSubstitutionRegister>(stateMap);
    else
      throw Exception("Unsupported substitution categorization: " + regType + " for alphabet " + alphabet->getAlphabetType());
  }

  auto csr = dynamic_cast<CategorySubstitutionRegister*>(reg.get());
  if (csr)
    csr->setStationarity(ApplicationTools::getBooleanParameter("stationarity", regArgs, true));

  if (verbose)
    ApplicationTools::displayResult("Substitution Register", regType);

  return reg;
}
