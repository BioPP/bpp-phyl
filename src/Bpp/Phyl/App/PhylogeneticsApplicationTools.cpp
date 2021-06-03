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
#include "../Model/FrequencySet/MvaFrequencySet.h"
#include "../Likelihood/TreeLikelihood.h"
#include "../Mapping/LaplaceSubstitutionCount.h"
#include "../Mapping/UniformizationSubstitutionCount.h"
#include "../Mapping/DecompositionSubstitutionCount.h"
#include "../Mapping/NaiveSubstitutionCount.h"
#include "../Mapping/OneJumpSubstitutionCount.h"
#include "../OptimizationTools.h"

#include "../Tree/Tree.h"
#include "../Tree/PhyloTreeTools.h"
#include "../Tree/TreeTools.h"

#include "../Tree/PhyloTree.h"
#include "../Io/Newick.h"
#include "../Io/NexusIoTree.h"
#include "../Io/Nhx.h"
#include "../Io/BppOTreeReaderFormat.h"
#include "../Io/BppOMultiTreeReaderFormat.h"
#include "../Io/BppOTreeWriterFormat.h"
#include "../Io/BppOMultiTreeWriterFormat.h"
#include "../Io/BppOBranchModelFormat.h"
#include "../Io/BppOFrequencySetFormat.h"
#include "../Io/BppORateDistributionFormat.h"

#include "../NewLikelihood/OneProcessSequenceEvolution.h"
#include "../NewLikelihood/MixtureSequenceEvolution.h"
#include "../NewLikelihood/PartitionSequenceEvolution.h"
#include "../NewLikelihood/AutoCorrelationSequenceEvolution.h"
#include "../NewLikelihood/HmmSequenceEvolution.h"

#include "../NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h"
#include "../NewLikelihood/PhyloLikelihoods/SingleDataPhyloLikelihood.h"
#include "../NewLikelihood/PhyloLikelihoods/OneProcessSequencePhyloLikelihood.h"
#include "../NewLikelihood/PhyloLikelihoods/PartitionProcessPhyloLikelihood.h"
#include "../NewLikelihood/PhyloLikelihoods/MixtureProcessPhyloLikelihood.h"
#include "../NewLikelihood/PhyloLikelihoods/AutoCorrelationProcessPhyloLikelihood.h"
#include "../NewLikelihood/PhyloLikelihoods/HmmProcessPhyloLikelihood.h"
#include "../NewLikelihood/PhyloLikelihoods/MixtureOfAlignedPhyloLikelihood.h"
#include "../NewLikelihood/PhyloLikelihoods/HmmOfAlignedPhyloLikelihood.h"
#include "../NewLikelihood/PhyloLikelihoods/AutoCorrelationOfAlignedPhyloLikelihood.h"
#include "../NewLikelihood/PhyloLikelihoods/FormulaOfPhyloLikelihood.h"
#include "../NewLikelihood/PhyloLikelihoods/ProductOfAlignedPhyloLikelihood.h"
#include "../NewLikelihood/NonHomogeneousSubstitutionProcess.h"
#include "../NewLikelihood/SimpleSubstitutionProcess.h"
#include "../NewLikelihood/SubstitutionProcessCollection.h"
#include "../NewLikelihood/RateAcrossSitesSubstitutionProcess.h"

#include "../NewLikelihood/ParametrizablePhyloTree.h"

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

Tree* PhylogeneticsApplicationTools::getTree(
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
  unique_ptr<ITree> iTree(bppoReader.readITree(format));
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
  unique_ptr<IMultiTree> iTrees(bppoReader.readIMultiTree(format));
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

/******************************************************************************/

 
map<size_t, Tree*> PhylogeneticsApplicationTools::getTrees(
  const map<string, string>& params,
  const map<size_t, AlignedValuesContainer*>& mSeq,
  map<string, string>& unparsedParams,
  const string& prefix,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  vector<string> vTreesName = ApplicationTools::matchingParameters(prefix + "tree*", params);

  map<size_t, Tree*> mTree;

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
      treeReader->readTrees(treeFilePath, trees);
      delete treeReader;

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
            delete mTree[i2 + 1];
          }

          mTree[i2 + 1] = trees[i2];
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
            delete mTree[num];
          }
          mTree[num] = trees[0];
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

      vector<string> names = mSeq.find(seqNum)->second->getSequencesNames();
      Tree* tree = TreeTemplateTools::getRandomTree(names);
      tree->setBranchLengths(1.);

      if (mTree.find(num) != mTree.end())
      {
        ApplicationTools::displayWarning("Tree " + TextTools::toString(num) + " already assigned, replaced by new one.");
        delete mTree[num];
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
          Tree* tree = mTree[i + 1];
          if (grafenHeight == "input")
          {
            h = TreeTools::getHeight(*tree, tree->getRootId());
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
          TreeTools::computeBranchLengthsGrafen(*tree, rho);
          double nh = TreeTools::getHeight(*tree, tree->getRootId());
          tree->scaleTree(h / nh);
        }
      }
      else
      {
        Tree* tree = mTree[num];
        if (grafenHeight == "input")
        {
          h = TreeTools::getHeight(*tree, tree->getRootId());
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
        TreeTools::computeBranchLengthsGrafen(*tree, rho);
        double nh = TreeTools::getHeight(*tree, tree->getRootId());
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

map<size_t, std::shared_ptr<PhyloTree>> PhylogeneticsApplicationTools::getPhyloTrees(
  const map<string, string>& params,
  const map<size_t, AlignedValuesContainer*>& mSeq,
  map<string, string>& unparsedParams,
  const string& prefix,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  vector<string> vTreesName = ApplicationTools::matchingParameters(prefix + "tree*", params);

  map<size_t, std::shared_ptr<PhyloTree>> mTree;

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

      IMultiTree* treeReader;
      if (format == "Newick")
        treeReader = new Newick(true);
      else if (format == "Nexus")
        treeReader = new NexusIOTree();
      else if (format == "NHX")
        treeReader = new Nhx();
      else
        throw Exception("Unknow format for tree reading: " + format);

      vector<PhyloTree*> trees;
      treeReader->readTrees(treeFilePath, trees);
      delete treeReader;

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

          mTree[i2 + 1] = std::shared_ptr<PhyloTree>(trees[i2]);
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
            mTree.erase(num);
          }
          mTree[num] = std::shared_ptr<PhyloTree>(trees[0]);
          ApplicationTools::displayResult("Number of leaves", trees[0]->getNumberOfLeaves());
        }
      }
    }
    else if (treeName == "random")
    {
      throw Exception("Random Phylotrees not defined yet. Ask developpers.");
      
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

      vector<string> names = mSeq.find(seqNum)->second->getSequencesNames();
//      PhyloTree* tree = TreeTemplateTools::getRandomTree(names);
      PhyloTree* tree = 0;
      tree->setBranchLengths(1.);

      if (mTree.find(num) != mTree.end())
      {
        ApplicationTools::displayWarning("Tree " + TextTools::toString(num) + " already assigned, replaced by new one.");
        mTree.erase(num);
      }
      mTree[num] = std::shared_ptr<PhyloTree>(tree);
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
        shared_ptr<PhyloBranch> branch=mTree[num]->getEdgeToFather(mTree[num]->getNode(static_cast<PhyloTree::NodeIndex> (TextTools::toInt(aveq.substr(5, string::npos)))));
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

SubstitutionModel* PhylogeneticsApplicationTools::getSubstitutionModel(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const AlignedValuesContainer* data,
  const map<string, string>& params,
  map<string, string>& unparsedParams,
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

  unparsedParams.insert(bIO.getUnparsedArguments().begin(), bIO.getUnparsedArguments().end());

  return model;
}

BranchModel* PhylogeneticsApplicationTools::getBranchModel(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const AlignedValuesContainer* data,
  const map<string, string>& params,
  map<string, string>& unparsedParams,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  BppOBranchModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, verbose, warn + 1);
  string modelDescription;
  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(alphabet);
  if (ca)
  {
    modelDescription = ApplicationTools::getStringParameter("model", params, "CodonRate(model=JC69)", suffix, suffixIsOptional, warn);
    if (!gCode)
      throw Exception("PhylogeneticsApplicationTools::getBranchModel(): a GeneticCode instance is required for instanciating a codon model.");
    bIO.setGeneticCode(gCode);
  }
  else if (AlphabetTools::isWordAlphabet(alphabet))
    modelDescription = ApplicationTools::getStringParameter("model", params, "Word(model=JC69)", suffix, suffixIsOptional, warn);
  else
    modelDescription = ApplicationTools::getStringParameter("model", params, "JC69", suffix, suffixIsOptional, warn);

  BranchModel* model = bIO.readBranchModel(alphabet, modelDescription, data, true);
  map<string, string> tmpUnparsedParameterValues(bIO.getUnparsedArguments());

  unparsedParams.insert(tmpUnparsedParameterValues.begin(), tmpUnparsedParameterValues.end());

  return model;
}      

/******************************************************************************/

map<size_t, std::shared_ptr<DiscreteDistribution>> PhylogeneticsApplicationTools::getRateDistributions(
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
  map<size_t, std::shared_ptr<DiscreteDistribution>> mDist;


  for (size_t i = 0; i < vratesName.size(); i++)
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

    mDist[num] = std::shared_ptr<DiscreteDistribution>(bIO.readDiscreteDistribution(distDescription, true));
  }

  if (mDist.size() == 0)
  {
    string distDescription = ApplicationTools::getStringParameter("rate_distribution", paramDist, "Constant()", suffix, suffixIsOptional);
    mDist[0]= std::shared_ptr<DiscreteDistribution>(bIO.readDiscreteDistribution(distDescription, true));
  }

  return mDist;
}


/*************************************************************/
/******* MODELS **********************************************/
/*************************************************************/

map<size_t, std::shared_ptr<BranchModel>> PhylogeneticsApplicationTools::getBranchModels(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const map<size_t, AlignedValuesContainer*>& mData,
  const map<string, string>& params,
  map<string, string>& unparsedParams,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  if (dynamic_cast<const CodonAlphabet*>(alphabet) && !gCode)
    throw Exception("PhylogeneticsApplicationTools::getBranchModels(): a GeneticCode instance is required for instanciating codon models.");

  string ModelFilePath = ApplicationTools::getAFilePath("models.file", params, false, false, suffix, suffixIsOptional,  "none", 1);

  map<string, string> paramModel;

  if (ModelFilePath != "none")
    paramModel = AttributesTools::getAttributesMapFromFile(ModelFilePath, "=");

  paramModel.insert(params.begin(), params.end());

  vector<string> modelsName = ApplicationTools::matchingParameters("model*", paramModel);

  vector<size_t> modelsNum;
  for (const auto& name:modelsName)
  {
    size_t poseq = name.find("=");
    if (name.find("nodes_id") == string::npos)
      modelsNum.push_back((size_t) TextTools::toInt(name.substr(5, poseq - 5)));
  }

  map<size_t, std::shared_ptr<BranchModel>> mModel;

  BppOBranchModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, verbose, warn);
  bIO.setGeneticCode(gCode);

  for (size_t i = 0; i < modelsNum.size(); i++)
  {
    if (i>=10)
    {
      bIO.setVerbose(false);
      warn=10;
      if (i==10)
        ApplicationTools::displayMessage("");
      ApplicationTools::displayResult("Model " + TextTools::toString(modelsNum[i]), string("..."));
    }
    else
    {
      ApplicationTools::displayMessage("");
      ApplicationTools::displayMessage("Model " + TextTools::toString(modelsNum[i]));
    }
    
    string modelDescription = ApplicationTools::getStringParameter("model" + TextTools::toString(modelsNum[i]), paramModel, "", suffix, suffixIsOptional, warn);

    map<string, string> args;
    string modelName;

    KeyvalTools::parseProcedure(modelDescription, modelName, args);

    size_t nData = 0;

    if (args.find("data") != args.end())
      nData = (size_t) TextTools::toInt(args["data"]);

    
    shared_ptr<BranchModel> model;
    if (args.find("data") != args.end() && mData.find(nData)!= mData.end())
      model = shared_ptr<BranchModel>(bIO.readBranchModel(alphabet, modelDescription, mData.find(nData)->second, true));
    else
      model = shared_ptr<BranchModel>(bIO.readBranchModel(alphabet, modelDescription, 0, true));
    
    map<string, string> tmpUnparsedParameterValues(bIO.getUnparsedArguments());

    for (auto& it : tmpUnparsedParameterValues)
      unparsedParams[it.first + "_" + TextTools::toString(modelsNum[i])] = it.second;

    if (verbose)
    {
      //   ApplicationTools::displayResult("substitution Model " + TextTools::toString(modelsNum[i]), model->getName());
      if (nData != 0)
        ApplicationTools::displayResult("Data used ", TextTools::toString(nData));
    }

    mModel[modelsNum[i]] = model;
  }

  return mModel;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::setSubstitutionModelParametersInitialValuesWithAliases(
  BranchModel& model,
  map<string, string>& unparsedParameterValues,
  size_t modelNumber,
  const AlignedValuesContainer* data,
  map<string, string>& sharedParams,
  bool verbose)
{
  string initFreqs = ApplicationTools::getStringParameter(model.getNamespace() + "initFreqs", unparsedParameterValues, "", "", true, 2);

  if (verbose)
    ApplicationTools::displayResult("Frequencies Initialization for model", (initFreqs == "") ? "None" : initFreqs);

  if (initFreqs != "")
  {
    auto tmodel=dynamic_cast<TransitionModel*>(&model);
    if (!tmodel)
      ApplicationTools::displayMessage("Frequencies initialization not possible for model " + model.getName());
    else
    {
      if (initFreqs == "observed")
      {
        if (!data)
          throw Exception("Missing data for observed frequencies");
        unsigned int psi = ApplicationTools::getParameter<unsigned int>(model.getNamespace() + "initFreqs.observedPseudoCount", unparsedParameterValues, 0);
        tmodel->setFreqFromData(*data, psi);
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
        tmodel->setFreq(frequencies);
      }
      else
        throw Exception("Unknown initFreqs argument");
    }
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

std::shared_ptr<FrequencySet> PhylogeneticsApplicationTools::getFrequencySet(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const string& freqDescription,
  const AlignedValuesContainer* data,
  map<string, string>& sharedparams,
  const vector<double>& rateFreqs,
  bool verbose,
  int warn)
{
  map<string, string> unparsedParameterValues;
  BppOFrequencySetFormat bIO(BppOFrequencySetFormat::ALL, verbose, warn);
  if (AlphabetTools::isCodonAlphabet(alphabet))
  {
    if (!gCode)
      throw Exception("PhylogeneticsApplicationTools::getFrequencySet(): a GeneticCode instance is required for instanciating a codon frequencies set.");
    bIO.setGeneticCode(gCode);
  }
  auto pFS=bIO.readFrequencySet(alphabet, freqDescription, data, true);

  map<string, string> unparsedparam = bIO.getUnparsedArguments();

  sharedparams.insert(unparsedparam.begin(), unparsedparam.end());

  // /////// To be changed for input normalization
  if (rateFreqs.size() > 0)
  {
    pFS=std::make_shared<MarkovModulatedFrequencySet>(pFS, rateFreqs);
  }

  return pFS;
}


std::shared_ptr<FrequencySet> PhylogeneticsApplicationTools::getRootFrequencySet(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const AlignedValuesContainer* data,
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

    auto freq = getFrequencySet(alphabet, gCode, freqDescription, data, unparams, rateFreqs, verbose, warn + 1);
    freq->setNamespace("root." + freq->getNamespace());

    for (auto& it : unparams)
      sharedparams["root." + it.first] = it.second;

    if (verbose)
      ApplicationTools::displayResult("Root frequencies ", freq->getName());
    return freq;
  }
}


map<size_t, std::shared_ptr<FrequencySet>> PhylogeneticsApplicationTools::getRootFrequencySets(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const map<size_t, AlignedValuesContainer*>& mData,
  const map<string, string>& params,
  map<string, string>& sharedparams,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  if (dynamic_cast<const CodonAlphabet*>(alphabet) && !gCode)
    throw Exception("PhylogeneticsApplicationTools::getRootFrequencySets(): a GeneticCode instance is required for instanciating codon frequencies sets.");

  string RootFilePath = ApplicationTools::getAFilePath("root_freq.file", params, false, false, suffix, suffixIsOptional,  "none", 1);
  map<string, string> paramRF;

  if (RootFilePath != "none")
    paramRF = AttributesTools::getAttributesMapFromFile(RootFilePath, "=");

  paramRF.insert(params.begin(), params.end());

  vector<string> vrfName = ApplicationTools::matchingParameters("root_freq*", paramRF);

  vector<size_t> rfNum;
  for (const auto& rfName: vrfName)
  {
    size_t poseq = rfName.find("=");
    try
    {
      rfNum.push_back((size_t) TextTools::toInt(rfName.substr(9, poseq - 9)));
    }
    catch (Exception& e)
    {}
  }

  BppOFrequencySetFormat bIO(BppOFrequencySetFormat::ALL, verbose, warn);
  bIO.setGeneticCode(gCode);

  map<size_t, std::shared_ptr<FrequencySet>> mFS;

  for (size_t i = 0; i < rfNum.size(); i++)
  {
    ApplicationTools::displayMessage("");
    ApplicationTools::displayMessage("Root Frequencies Set " + TextTools::toString(rfNum[i]));

    string freqDescription = ApplicationTools::getStringParameter("root_freq" + TextTools::toString(rfNum[i]), paramRF, "", suffix, suffixIsOptional, warn);

    map<string, string> args;
    string freqName;

    KeyvalTools::parseProcedure(freqDescription, freqName, args);

    size_t nData = 0;

    if (args.find("data") != args.end())
      nData = (size_t) TextTools::toInt(args["data"]);

    shared_ptr<FrequencySet> rFS(bIO.readFrequencySet(alphabet, freqDescription, (args.find("data") != args.end()) ? mData.find(nData)->second : 0, true));
    rFS->setNamespace("root." + rFS->getNamespace());
    map<string, string> unparsedparam = bIO.getUnparsedArguments();

    for (auto& it : unparsedparam)
      sharedparams["root." + it.first + "_" + TextTools::toString(rfNum[i])] = it.second;

    if (verbose)
    {
      // ApplicationTools::displayResult("Root Frequencies Set " + TextTools::toString(rfNum[i]), rFS->getName());
      if (nData != 0)
        ApplicationTools::displayResult("Data used ", TextTools::toString(nData));
    }

    mFS[rfNum[i]] = rFS;
  }

  return mFS;
}

/******************************************************/
/**** SETOFMODELPATH **********************************/
/******************************************************/

map<size_t, std::shared_ptr<ModelPath>> PhylogeneticsApplicationTools::getModelPaths(
  const std::map<std::string, std::string>& params,
  const map<size_t, std::shared_ptr<BranchModel>>& mModel,
  bool verbose)
{
  string ModelPathsPath = ApplicationTools::getAFilePath("path.file", params, false, false, "", true,  "none", 1);
  map<string, string> paramMP;

  if (ModelPathsPath != "none")
    paramMP = AttributesTools::getAttributesMapFromFile(ModelPathsPath, "=");

  paramMP.insert(params.begin(), params.end());

  vector<string> vmpName = ApplicationTools::matchingParameters("path*", paramMP);

  map<size_t, std::shared_ptr<ModelPath>> modelPaths;

  for (size_t i = 0; i < vmpName.size(); i++)
  {
    const auto& name=vmpName[i];

    string desc = ApplicationTools::getStringParameter(name, paramMP, "", "", true);

    size_t num;
    try{
      num=size_t(TextTools::toInt(name.substr(4)));
    }
    catch (const Exception& e)
    {
      throw Exception("PhylogeneticsApplicationTools::getModelPaths: bad path number in line " + name);
    }

    modelPaths[num]=std::make_shared<ModelPath>();

    if (verbose)
    {
      if (i >= 10)
      {
        if (i==10)
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
      string::size_type indexo = submodel.find("[");
      string::size_type indexf = submodel.find("]");
      if ((indexo == string::npos) | (indexf == string::npos))
        throw Exception("PhylogeneticsApplicationTools::getModelPaths. Bad path syntax, should contain `[]' symbols: " + submodel);

      auto pos=submodel.find("model");
      if (pos==string::npos)
        throw Exception("PhylogeneticsApplicationTools::getModelPaths. Missing identifier 'model' in description: " + submodel);
        
      size_t num2 = TextTools::to<size_t>(submodel.substr(pos+5, indexo - 5 - pos));
      if (mModel.find(num2)==mModel.end())
        throw BadIntegerException("PhylogeneticsApplicationTools::getModelPaths: Wrong model number", static_cast<int>(num2));
      
      auto pSM = std::dynamic_pointer_cast<MixedTransitionModel>(mModel.at(num2));
      if (!pSM)
        throw Exception("PhylogeneticsApplicationTools::getModelPaths: Model number "+ TextTools::toString(num2) + " ( " + mModel.at(num2)->getName() + " ) is not Mixed.");
      
      string lp2 = submodel.substr(indexo + 1, indexf - indexo - 1);      
      StringTokenizer stp2(lp2, ",");
      while (stp2.hasMoreToken())
      {
        string p2=stp2.nextToken();

        uint n2;
        bool n2ok=true;
        try  {
          n2=TextTools::to<uint>(p2);
          if (n2<=0 || n2>pSM->getNumberOfModels())
            n2ok=false;
          else
            submodelNb.push_back(n2-1);
        }
        catch (Exception& e)
        {
          Vuint submodnb = pSM->getSubmodelNumbers(p2);
          if (submodelNb.size()==0)
            submodelNb=submodnb;
          else
            submodelNb=VectorTools::vectorIntersection(submodelNb,submodnb);
        }

        if (!n2ok)
          throw BadIntegerException("PhylogeneticsApplicationTools::getModelPaths: Wrong model number for model " + TextTools::toString(num2), int(n2));
      }

      modelPaths[num]->setModel(pSM,submodelNb);
      if (!modelPaths[num]->getLeadModel())
        modelPaths[num]->setLeadModel(pSM);
    }

    if (verbose &&  (i < 10))
      ApplicationTools::displayResult("Model Path", desc);
  }
  
  return modelPaths;
}


map<size_t, std::shared_ptr<ModelScenario>> PhylogeneticsApplicationTools::getModelScenarios(
  const std::map<std::string, std::string>& params,
  const map<size_t, std::shared_ptr<ModelPath>>& mModelPath,
  const map<size_t, std::shared_ptr<BranchModel>>& mModel,
  bool verbose)
{
  string ModelPathsPath = ApplicationTools::getAFilePath("scenario.file", params, false, false, "", true,  "none", 1);
  map<string, string> paramMS;

  if (ModelPathsPath != "none")
    paramMS = AttributesTools::getAttributesMapFromFile(ModelPathsPath, "=");

  paramMS.insert(params.begin(), params.end());

  vector<string> vmsName = ApplicationTools::matchingParameters("scenario*", paramMS);

  map<size_t, std::shared_ptr<ModelScenario>> somp;

  for (size_t i = 0; i < vmsName.size(); i++)
  {
    const auto& name=vmsName[i];

    string desc = ApplicationTools::getStringParameter(name, paramMS, "", "", true);

    size_t num;
    try{
      num=size_t(TextTools::toInt(name.substr(8)));
    }
    catch (const Exception& e)
    {
      throw Exception("PhylogeneticsApplicationTools::getModelScenarios: bad scenario number in line " + name);
    }

    somp[num]=std::make_shared<ModelScenario>();

    if (verbose)
    {
      if (i >= 10)
      {
        if (i==10)
          ApplicationTools::displayMessage("");
        ApplicationTools::displayResult("Scenario " + TextTools::toString(num), string("..."));
      }
      else
      {
        ApplicationTools::displayMessage("");
        ApplicationTools::displayMessage("Scenario " + TextTools::toString(num));
      }
    }

    bool complete=false;
    size_t numpath;

    StringTokenizer st(desc, "&");
    while (st.hasMoreToken())
    {      
      string path = st.nextToken();
      bool numok=true;
      try {
        if (path=="complete")
          complete=true;
        else if (path.substr(0,5)=="split")
        {
          auto pos=path.find("model");
          if (pos==string::npos)
            throw Exception("PhylogeneticsApplicationTools::getModelScenarios. Missing identifier 'model' in scenarion description: " + path);

          auto poseq = path.find("=",pos);
          size_t num2 = TextTools::to<size_t>(path.substr(poseq+1));
          
          if (mModel.find(num2)==mModel.end())
            throw BadIntegerException("PhylogeneticsApplicationTools::getModelScenarios: Wrong model number", static_cast<int>(num2));
      
          auto pSM = std::dynamic_pointer_cast<MixedTransitionModel>(mModel.at(num2));
          if (!pSM)
            throw Exception("PhylogeneticsApplicationTools::getModelScenarios: Model number "+ TextTools::toString(num2) + " ( " + mModel.at(num2)->getName() + " ) is not Mixed.");

          std::vector<std::shared_ptr<ModelPath>> modelPaths;

          auto nmod = pSM->getNumberOfModels();
          
          for (uint nm = 0; nm < (uint)nmod; nm++)
          {
            auto mp = std::make_shared<ModelPath>();
            mp->setModel(pSM,Vuint({nm}));
            mp->setLeadModel(pSM);
            somp[num]->addModelPath(mp);
          }
        }
        else
        {
          numpath = TextTools::to<size_t>(path.substr(4));
          if (mModelPath.find(numpath)==mModelPath.end())
            numok=false;
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
      if (somp[num]->getNumberOfModelPaths()==0)
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


SubstitutionProcess* PhylogeneticsApplicationTools::getSubstitutionProcess(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const AlignedValuesContainer* pData,
  const vector<PhyloTree*>& vTree,
  const map<string, string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  SubstitutionProcess* SP = 0;

  map<string, string> unparsedParams;

  string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", params, "no", "", true, warn);
  ApplicationTools::displayResult("Heterogeneous process", nhOpt);

  // ///////////////////////
  // Tree

  unique_ptr<ParametrizablePhyloTree> pTree(new ParametrizablePhyloTree(*vTree[0]));

  // ////////////////////////
  // Rates

  unique_ptr<DiscreteDistribution> rDist(getRateDistribution(params));

  BppOBranchModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, false, warn);
  bIO.setGeneticCode(gCode);


  // /////////////////////////
  // / Models

  string tmpDesc;

  if (nhOpt == "no")
  {
    // Homogeneous & stationary models

    shared_ptr<BranchModel> tmp(getBranchModel(alphabet, gCode, pData, params, unparsedParams));

    if (tmp->getNumberOfStates() >= 2 * tmp->getAlphabet()->getSize() || (rDist->getName() == "Constant")) // first test is for Markov-modulated Markov model!
      SP = new SimpleSubstitutionProcess(tmp, pTree.release());
    else
      SP = new RateAcrossSitesSubstitutionProcess(tmp, rDist.release(), pTree.release());
  }

  // Non-homogeneous models
  else
  {
    string fName = (nhOpt == "one_per_branch" ? "model" : "model1");

    tmpDesc = ApplicationTools::getStringParameter(fName, params, "", suffix, suffixIsOptional, warn);
    shared_ptr<BranchModel> tmp(bIO.readBranchModel(alphabet, tmpDesc, pData, true));


    // ////////////////////////////////////
    // Root frequencies

    bool stationarity = ApplicationTools::getBooleanParameter("nonhomogeneous.stationarity", params, false, "", false, warn);

    shared_ptr<FrequencySet> rootFrequencies;

    if (!stationarity)
    {
      // Markov Modulated  models
      vector<double> rateFreqs;
      if (tmp->getNumberOfStates() != alphabet->getSize())
      {
        // Markov-Modulated Markov Model...
        size_t n = static_cast<size_t>(tmp->getNumberOfStates() / alphabet->getSize());
        rateFreqs = vector<double>(n, 1. / static_cast<double>(n)); // Equal rates assumed for now, may be changed later (actually, in the most general case,
      }

      // MVA models

      string freqDescription = ApplicationTools::getStringParameter("nonhomogeneous.root_freq", params, "", suffix, suffixIsOptional, warn);
      if (freqDescription.substr(0, 10) == "MVAprotein")
      {
        if (dynamic_cast<Coala*>(tmp.get()))
          dynamic_pointer_cast<MvaFrequencySet>(rootFrequencies)->initSet(dynamic_cast<CoalaCore*>(tmp.get()));
        else
          throw Exception("The MVAprotein frequencies set at the root can only be used if a Coala model is used on branches.");
      }
      else
        rootFrequencies=getRootFrequencySet(alphabet, gCode, pData, params, unparsedParams, rateFreqs, suffix, suffixIsOptional, warn);

      stationarity = !rootFrequencies.get();
    }

    ApplicationTools::displayBooleanResult("Stationarity assumed", stationarity);

    // /////////////////////////////////////
    // One_per_branch

    if (nhOpt == "one_per_branch")
    {
      vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", params, ',', "");

      for (unsigned int i = 0; i < globalParameters.size(); i++)
      {
        ApplicationTools::displayResult("Global parameter", globalParameters[i]);
      }

      SP = NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(
        tmp,
        rDist.release(),
        pTree.release(),
        rootFrequencies,
        globalParameters);
    }
    else
    {
      // //////////////////////////////
      // General

      size_t nbModels = ApplicationTools::getParameter<size_t>("nonhomogeneous.number_of_models", params, 1, suffix, suffixIsOptional, warn);

      if (nbModels == 0)
        throw Exception("The number of models can't be 0 !");

      if (verbose)
        ApplicationTools::displayResult("Number of distinct models", TextTools::toString(nbModels));

      // //////////////////////////////////////
      // Now parse all models:

      bIO.setVerbose(true);

      SP = new NonHomogeneousSubstitutionProcess(rDist.release(), pTree.release(), rootFrequencies->clone());

      NonHomogeneousSubstitutionProcess* nhSP = dynamic_cast<NonHomogeneousSubstitutionProcess*>(SP);

      for (size_t i = 0; i < nbModels; i++)
      {
        string prefix = "model" + TextTools::toString(i + 1);
        string modelDesc;
        modelDesc = ApplicationTools::getStringParameter(prefix, params, "", suffix, suffixIsOptional, warn);

        shared_ptr<BranchModel> model(bIO.readBranchModel(alphabet, modelDesc, pData, true));
        map<string, string> tmpUnparsedParameterValues(bIO.getUnparsedArguments());

        for (auto& it : tmpUnparsedParameterValues)
          unparsedParams[it.first + "_" + TextTools::toString(i + 1)] = it.second;

        vector<unsigned int> nodesId;

        auto snodesid = prefix + ".nodes_id";
        auto descnodes = ApplicationTools::getStringParameter(snodesid, params, "", suffix, suffixIsOptional, warn);

        const auto& tree= SP->getParametrizablePhyloTree();
        if (descnodes == "All")
        {
          nodesId = pTree->getEdgeIndexes(pTree->getSubtreeEdges(tree.getRoot()));
        }
        else if (descnodes == "Leaves")
        {
          nodesId = pTree->getNodeIndexes(pTree->getLeavesUnderNode(tree.getRoot()));
        }
        else if (descnodes == "NoLeaves")
        {
          auto allIds= pTree->getEdgeIndexes(pTree->getSubtreeEdges(tree.getRoot()));
          auto leavesId = pTree->getNodeIndexes(pTree->getLeavesUnderNode(tree.getRoot()));
          VectorTools::diff(allIds, leavesId, nodesId);
        }
        else
          nodesId = ApplicationTools::getVectorParameter<unsigned int>(snodesid, params, ',', ':', "", suffix, suffixIsOptional, warn);

        if (verbose)
          ApplicationTools::displayResult("Model" + TextTools::toString(i + 1) + " is associated to", TextTools::toString(nodesId.size()) + " node(s).");

        nhSP->addModel(model, nodesId);
      }

      nhSP->isFullySetUp();
    }
  }


  // ////// Aliasing
  // Finally check parameter aliasing:

  string aliasDesc = ApplicationTools::getStringParameter("nonhomogeneous.alias", params, "", suffix, suffixIsOptional, warn);

  StringTokenizer st(aliasDesc, ",");
  while (st.hasMoreToken())
  {
    string alias = st.nextToken();
    string::size_type index = alias.find("->");
    if (index == string::npos)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionProcess. Bad alias syntax, should contain `->' symbol: " + alias);
    string p1 = alias.substr(0, index);
    string p2 = alias.substr(index + 2);
    unparsedParams[p1] = p2;
  }

  SP->aliasParameters(unparsedParams, verbose);

  return SP;
}

/************************************************************/

bool PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember(
  SubstitutionProcessCollection* SubProColl,
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
    if (warn)
      ApplicationTools::displayWarning("Warning, unknown process name: " + procName);

    return 0;
  }

  // ///
  // tree number

  if (args.find("tree") == args.end())
    throw Exception("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember. A tree number is compulsory.");

  size_t numTree = (size_t) ApplicationTools::getIntParameter("tree", args, 1, "", true, warn);

  if (!SubProColl->hasTreeNumber(numTree))
    throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : unknown tree number", (int)numTree);


  // /////
  // rate number

  if (args.find("rate") == args.end())
    throw Exception("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember. A rate number is compulsory.");

  size_t numRate;

  string sRate = ApplicationTools::getStringParameter("rate", args, "1", "", true, warn);

  size_t pp = sRate.find(".");

  numRate = static_cast<size_t> (TextTools::toInt(sRate.substr(0, pp)));
  if (!SubProColl->hasDistributionNumber(numRate))
    throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : unknown rate number", (int)numRate);

  if (pp != string::npos)
  {
    size_t numSRate = static_cast<size_t> (TextTools::toInt(sRate.substr(pp + 1)));
    SubProColl->addDistribution(std::make_shared<ConstantDistribution>(SubProColl->getRateDistribution(numRate).getCategory(numSRate)), 10000 * (numRate + 1) + numSRate);

    numRate = 10000 * (numRate + 1) + numSRate;
  }


  // ////////
  // root freq number

  bool stationarity = (args.find("root_freq") == args.end());
  size_t numFreq = 0;

  if (!stationarity)
  {
    numFreq = (size_t) ApplicationTools::getIntParameter("root_freq", args, 1, "", true, warn);
    if (!SubProColl->hasFrequenciesNumber(numFreq))
      throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : unknown root frequencies number", (int)numFreq);
  }

  // ///
  // scenario number

  size_t numScen = 0;
  
  if (args.find("scenario") != args.end())
  {
    numScen = (size_t) ApplicationTools::getIntParameter("scenario", args, 1, "", true, warn);

    if (!SubProColl->hasModelScenario(numScen))
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
    
    if (!SubProColl->hasModelNumber(numModel))
      throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : unknown model number", static_cast<int> (numModel));

    vector<uint> vNodes = SubProColl->getTree(numTree).getAllEdgesIndexes();

    map<size_t, vector<unsigned int> > mModBr;
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

      if (numScen!=0)
        ApplicationTools::displayResult (" Scenario number", TextTools::toString(numScen));
      
      if (!stationarity)
        ApplicationTools::displayResult (" Root frequencies number", TextTools::toString(numFreq));
      else
        ApplicationTools::displayMessage(" Stationarity assumed.");
    }

    if (stationarity)
      SubProColl->addSubstitutionProcess(procNum, mModBr, numTree, numRate);
    else
      SubProColl->addSubstitutionProcess(procNum, mModBr, numTree, numRate, numFreq);
  }

  else if ((procName == "Nonhomogeneous") ||  (procName == "NonHomogeneous"))
  {
    size_t indModel = 1;
    map<size_t, vector<unsigned int> > mModBr;

    while (args.find("model" + TextTools::toString(indModel)) != args.end())
    {
      size_t numModel = (size_t) ApplicationTools::getIntParameter("model" + TextTools::toString(indModel), args, 1, "", true, warn);

      if (mModBr.find(numModel) != mModBr.end())
        throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : model number seen twice.", (int)numModel);

      vector<unsigned int> nodesId;

      auto snodesid = "model" + TextTools::toString(indModel)  + ".nodes_id";
      auto descnodes = ApplicationTools::getStringParameter(snodesid, args, "", "", true, warn);

      auto& tree= SubProColl->getTree(numTree);
      if (descnodes == "All")
      {
        nodesId = tree.getEdgeIndexes(tree.getSubtreeEdges(tree.getRoot()));
      }
      else if (descnodes == "Leaves")
      {
        nodesId = tree.getNodeIndexes(tree.getLeavesUnderNode(tree.getRoot()));
      }
      else if (descnodes == "NoLeaves")
      {
        auto allIds= tree.getEdgeIndexes(tree.getSubtreeEdges(tree.getRoot()));
        auto leavesId = tree.getNodeIndexes(tree.getLeavesUnderNode(tree.getRoot()));
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
      SubProColl->addSubstitutionProcess(procNum, mModBr, numTree, numRate);
    else
      SubProColl->addSubstitutionProcess(procNum, mModBr, numTree, numRate, numFreq);
  }
  else if (procName == "OnePerBranch")
  {
    if (args.find("model") == args.end())
      throw Exception("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember. A model number is compulsory.");

    size_t numModel = (size_t) ApplicationTools::getIntParameter("model", args, 1, "", true, warn);

    if (!SubProColl->hasModelNumber(numModel))
      throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : unknown model number", (int)numModel);

    vector<string> sharedParameters = ApplicationTools::getVectorParameter<string>("shared_parameters", args, ',', "", "", true, 1);

    if (stationarity)
      SubProColl->addOnePerBranchSubstitutionProcess(procNum, numModel, numTree, numRate, sharedParameters);
    else
      SubProColl->addOnePerBranchSubstitutionProcess(procNum, numModel, numTree, numRate, numFreq, sharedParameters);

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
        ApplicationTools::displayResult(" Shared parameter", sP);
    }
  }

  if (numScen!=0)
    SubProColl->getSubstitutionProcess(procNum).setModelScenario(numScen);

  return true;
}



/******************************************************************************/

SubstitutionProcessCollection* PhylogeneticsApplicationTools::getSubstitutionProcessCollection(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const map<size_t, std::shared_ptr<PhyloTree>>& mTree,
  const map<size_t, std::shared_ptr<BranchModel>>& mMod,
  const map<size_t, std::shared_ptr<FrequencySet>>& mRootFreq,
  const map<size_t, std::shared_ptr<DiscreteDistribution>>& mDist,
  const map<size_t, std::shared_ptr<ModelScenario>>& mScen,
  const map<string, string>& params,
  map<string, string>& unparsedParams,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  SubstitutionProcessCollection*  SPC = new SubstitutionProcessCollection();

  map<string, double> existingParameters;

  // ///////////////////////
  // Trees

  if (mTree.size() == 0)
    throw Exception("Missing tree in construction of SubstitutionProcessCollection.");
  for (const auto& itt : mTree)
    SPC->addTree(std::make_shared<ParametrizablePhyloTree>(*(itt.second)), itt.first);
  
  // ///////////////////////
  // Rates

  if (mDist.size() == 0)
    throw Exception("Missing rate distribution in construction of SubstitutionProcessCollection.");

  for (const auto& itd : mDist)
    SPC->addDistribution(itd.second, itd.first);

  // ////////////////////////
  // Models

  if (mMod.size() == 0)
    throw Exception("Missing model in construction of SubstitutionProcessCollection.");

  for (const auto& itm : mMod)
    SPC->addModel(itm.second, itm.first);

  // ///////////////////////////
  // Root Frequencies

  for (const auto& itr : mRootFreq)
    SPC->addFrequencies(itr.second, itr.first);

  // ///////////////////////
  // Scenarios

  for (const auto& itt : mScen)
    SPC->addScenario(itt.second, itt.first);
  
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
      num = static_cast<size_t>(TextTools::toInt(suff));
    else
      num = 1;

    bool addok=addSubstitutionProcessCollectionMember(SPC, num, params, (nT<10?verbose:false), warn);

    if (addok)
    {
      if (nT==10)
        ApplicationTools::displayMessage("");
      
      if (nT>=10)
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
  //   addSubstitutionProcessCollectionMember(SPC, params, processNum[i]);


  // /////////////////////////
  // Now set shared parameters:

  // ////// Aliasing
  // Finally check parameter aliasing:

  string aliasDesc = ApplicationTools::getStringParameter("likelihood.alias", params, "", suffix, suffixIsOptional, warn);

  StringTokenizer st(aliasDesc, ",");
  while (st.hasMoreToken())
  {
    string alias = st.nextToken();
    string::size_type index = alias.find("->");
    if (index == string::npos)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionProcessCollection. Bad alias syntax, should contain `->' symbol: " + alias);
    string p1 = alias.substr(0, index);
    string p2 = alias.substr(index + 2);
    unparsedParams[p1] = p2;
  }

  SPC->aliasParameters(unparsedParams, verbose);

  return SPC;
}

/******************************************************/
/**** SEQUENCE EVOLUTIONS *****************************/
/******************************************************/

map<size_t, SequenceEvolution*> PhylogeneticsApplicationTools::getSequenceEvolutions(
  SubstitutionProcessCollection& SPC,
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
  for (size_t i = 0; i < evolsName.size(); i++)
  {
    size_t poseq = evolsName[i].find("=");
    evolsNum.push_back((size_t) TextTools::toInt(evolsName[i].substr(7, poseq - 7)));
  }

  map<size_t, SequenceEvolution*> mEvol;

  for (size_t mPi = 0; mPi < evolsNum.size(); mPi++)
  {
    if (SPC.hasSubstitutionProcessNumber(evolsNum[mPi]))
      continue;

    SequenceEvolution* nEvol;

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
      if (!SPC.hasSubstitutionProcessNumber(nproc))
        throw BadIntegerException("PhylogeneticsApplicationTools::getEvolutions. Unknown process number:", (int)nproc);

      nEvol = new OneProcessSequenceEvolution(SPC.getSubstitutionProcess(nproc), nproc);
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

      for (size_t i = 0; i < vproc.size(); i++)
      {
        if (!SPC.hasSubstitutionProcessNumber(vproc[i]))
          throw BadIntegerException("PhylogeneticsApplicationTools::getEvolutions. Unknown process number:", (int)vproc[i]);
      }

      if (evolName == "Partition")
      {
        // parse all processes sites

        vector<size_t> vMap;

        map<size_t, size_t> posProc;

        for (size_t i = 0; i < vproc.size(); i++)
        {
          string prefix = "process" + TextTools::toString(i + 1);

          vector<size_t> procPos = ApplicationTools::getVectorParameter<size_t>(prefix + ".sites", args, ',', ':', TextTools::toString(i), "", true, true);

          for (size_t j = 0; j < procPos.size(); j++)
          {
            if (posProc.find(procPos[j]) != posProc.end())
              throw BadIntegerException("A process position is defined twice ", (int)procPos[j]);
            else
              posProc[procPos[j]] = vproc[i];
          }
        }

        size_t pos = 1;

        while (posProc.find(pos) != posProc.end())
        {
          vMap.push_back(posProc[pos]);
          pos++;
        }

        if (vMap.size() != posProc.size())
          throw Exception("Error : there are gaps in the process sites");

        if (verbose)
          ApplicationTools::displayResult("Process type", string("Partition"));

        PartitionSequenceEvolution* pMP = new PartitionSequenceEvolution(&SPC, vMap);

        nEvol = pMP;
      }
      else if (evolName == "Mixture")
      {
        MixtureSequenceEvolution* pMP = new MixtureSequenceEvolution(&SPC, vproc);

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

        nEvol = pMP;
      }
      else if (evolName == "HMM")
      {
        HmmSequenceEvolution* pMP = new HmmSequenceEvolution(&SPC, vproc);

        if (verbose)
          ApplicationTools::displayResult("Process type", string("HMM"));

        size_t nbP = pMP->getNumberOfSubstitutionProcess();

        string vs = "(" + VectorTools::paste(vector<double>(nbP, 1. / (double)nbP), ",") + ")";
        string vvs = "(";
        for (size_t i = 0; i < nbP; i++)
        {
          vvs += (i == 0 ? "" : ",") + vs;
        }
        vvs += ")";

        RowMatrix<double> mat = ApplicationTools::getMatrixParameter<double>("probas", args, ',', vvs);

        FullHmmTransitionMatrix fhtm(pMP->getHmmTransitionMatrix().getHmmStateAlphabet(), pMP->getNamespace());
        fhtm.setTransitionProbabilities(mat);

        pMP->matchParametersValues(fhtm.getParameters());

        nEvol = pMP;
      }
      else if (evolName == "AutoCorr")
      {
        AutoCorrelationSequenceEvolution* pMP = new AutoCorrelationSequenceEvolution(&SPC, vproc);

        size_t nbP = pMP->getNumberOfSubstitutionProcess();

        if (verbose)
          ApplicationTools::displayResult("Process type", string("AutoCorr"));

        string vs = "(" + VectorTools::paste(vector<double>(nbP, 1. / (int)nbP), ",") + ")";

        vector<double> v = ApplicationTools::getVectorParameter<double>("lambdas", args, ',', vs);

        ParameterList pl;
        
        for (size_t i = 0; i < v.size(); i++)
        {
          pl.addParameter(Parameter("AutoCorr.lambda" + TextTools::toString(i + 1), v[i]));
        }

        pMP->matchParametersValues(pl);

        nEvol = pMP;
      }
      else
        throw Exception("Unknown Process description : " + evolName);

      if (verbose){
        ApplicationTools::displayResult (" Process numbers", VectorTools::paste(vproc, ","));
        ApplicationTools::displayMessage("");
      }
    }
    
    mEvol[evolsNum[mPi]] = nEvol;
  }

  return mEvol;
}


/******************************************************/
/**** PHYLO LIKELIHOODS *********************************/
/******************************************************/

std::shared_ptr<PhyloLikelihoodContainer> PhylogeneticsApplicationTools::getPhyloLikelihoodContainer(
  Context& context,
  SubstitutionProcessCollection& SPC,
  map<size_t, SequenceEvolution*>& mSeqEvol,
  const map<size_t, AlignedValuesContainer*>& mData,
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

  map<size_t, vector<size_t> > phylosMap;

  for (size_t i = 0; i < phylosName.size(); i++)
  {
    size_t poseq = phylosName[i].find("=");
    size_t phyln = (size_t) TextTools::toInt(phylosName[i].substr(5, poseq - 5));

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

    PhyloLikelihood* nPL;
    string phyloName = "";

    map<string, string> args;

    string phyloDesc = ApplicationTools::getStringParameter("phylo", params, "Single", TextTools::toString(phylonum), warn);

    if (verbose)
    {
      if (nbPh<=20)
        ApplicationTools::displayMessage("");
      else
        verbhere=false;
      
      ApplicationTools::displayMessage("Phylolikelihood " + TextTools::toString(phylonum));
    }

    KeyvalTools::parseProcedure(phyloDesc, phyloName, args);

    // Data

    size_t nData = (args.find("data") == args.end()? 1: (size_t)TextTools::toInt(args["data"]));
    
    if (mData.find(nData) == mData.end())
    {
      ApplicationTools::displayWarning("PhylogeneticsApplicationTools::getPhyloLikelihoodContainer. Data number is wrong:" + TextTools::toString(nData) + ". Not built.");
      continue;
    }

    const AlignedValuesContainer* data = dynamic_cast<const AlignedValuesContainer*>(mData.find(nData)->second);

    if (!data)
    {
      ApplicationTools::displayWarning("PhylogeneticsApplicationTools::getPhyloLikelihoodContainer. Data " + TextTools::toString(nData) + " does not match with aligned sequences. Not built.");
      continue;
    }
    
    if (verbhere)
      ApplicationTools::displayResult(" Data used ", TextTools::toString(nData));

    // Sequence Evolution or process

    size_t nProcess = (args.find("process") == args.end()? 1: (size_t) TextTools::toInt(args["process"]));
    if (verbhere)
      ApplicationTools::displayResult(" Process ", TextTools::toString(nProcess));


    // Compression

    char compression = (args.find("compression") != args.end() && args["compression"] == "recursive")?'R':'S';
    if (verbhere)
      ApplicationTools::displayResult(" Compression ", (compression == 'R') ? "recursive" : "simple");

    // Construction

    if (SPC.hasSubstitutionProcessNumber(nProcess))
    {
      auto l = std::make_shared<LikelihoodCalculationSingleProcess>(*collNodes, *data, nProcess);
      nPL = new SingleProcessPhyloLikelihood(context, l, nProcess, nData);
    }
    else if (mSeqEvol.find(nProcess) != mSeqEvol.end())
    {
      // ////////////////
      // / from sequence evolutions to phyloLikelihoods

      OneProcessSequenceEvolution* opse = dynamic_cast<OneProcessSequenceEvolution*>(mSeqEvol[nProcess]);

      if (opse != NULL)
        nPL = new OneProcessSequencePhyloLikelihood(*data, *opse, *collNodes, nProcess, nData, true, compression == 'R');
      else
      {
        MixtureSequenceEvolution* mse = dynamic_cast<MixtureSequenceEvolution*>(mSeqEvol[nProcess]);

        if (mse != NULL)
          nPL = new MixtureProcessPhyloLikelihood(*data, *mse, *collNodes, nProcess, nData, true, compression == 'R');

        else
        {
          HmmSequenceEvolution* hse = dynamic_cast<HmmSequenceEvolution*>(mSeqEvol[nProcess]);

          if (hse != NULL)
            nPL = new HmmProcessPhyloLikelihood(*data, *hse, *collNodes, nProcess, nData, true, compression == 'R');

          else
          {
            AutoCorrelationSequenceEvolution* ase = dynamic_cast<AutoCorrelationSequenceEvolution*>(mSeqEvol[nProcess]);

            if (ase != NULL)
              nPL = new AutoCorrelationProcessPhyloLikelihood(*data, *ase, *collNodes, nProcess, nData, true, compression == 'R');
            else
            {
              PartitionSequenceEvolution* pse = dynamic_cast<PartitionSequenceEvolution*>(mSeqEvol[nProcess]);

              if (pse != NULL)
                nPL = new PartitionProcessPhyloLikelihood(*data, *pse, collNodes, nProcess, nData, true, compression == 'R');

              else
                throw Exception("PhylogeneticsApplicationTools::getPhyloLikelihoodContainer : Unknown Sequence Evolution.");
            }
          }
        }
      }
    }
    else
      throw BadIntegerException("PhylogeneticsApplicationTools::getPhyloLikelihoodContainer : Unknown Process number.", (int)nProcess);

    mPhylo->addPhyloLikelihood(phylonum, nPL);
    usedPhylo.push_back(phylonum);
  }

  // Now clean the map
  for (map<size_t, vector<size_t> >::iterator it = phylosMap.begin(); it != phylosMap.end(); )
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
        vphyl.erase(vphyl.begin() + static_cast<ptrdiff_t> (i - 1));
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

    for (map<size_t, vector<size_t> >::iterator it = phylosMap.begin(); it != phylosMap.end(); it++)
    {
      nbPh++;
      
      if (it->second.size() == 0)
      {
        size_t phylonum = it->first;

        PhyloLikelihood* nPL;
        string phyloName = "";

        map<string, string> args;

        string phyloDesc = ApplicationTools::getStringParameter("phylo", params, "Single", TextTools::toString(phylonum), warn);
        KeyvalTools::parseProcedure(phyloDesc, phyloName, args);

        if (verbose)
        {
          if (nbPh<=20)
            ApplicationTools::displayMessage("");
          else
            verbhere=false;
      
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
          MixtureOfAlignedPhyloLikelihood* pMA = new MixtureOfAlignedPhyloLikelihood(context, mPhylo, vPhylo);
          vector<double> vprob = ApplicationTools::getVectorParameter<double>("probas", args, ',', "(" + VectorTools::paste(vector<double>(vPhylo.size(), 1. / (double)vPhylo.size())) + ")");
          if (vprob.size() != 1)
          {
            if (vprob.size() != vPhylo.size())
              throw BadSizeException("Wrong size of probas description in Mixture", vprob.size(), vPhylo.size());
            Simplex si(vprob);
            pMA->setPhyloProb(si);
          }

          nPL = dynamic_cast<PhyloLikelihood*>(pMA);
        }
        else if (phyloName == "HMM")
        {
          HmmOfAlignedPhyloLikelihood* pMA = new HmmOfAlignedPhyloLikelihood(context, mPhylo, vPhylo);

          size_t nbP = pMA->getNumbersOfPhyloLikelihoods().size();

          string vs = "(" + VectorTools::paste(vector<double>(nbP, 1. / (double)nbP), ",") + ")";
          string vvs = "(";
          for (size_t i = 0; i < nbP; i++)
          {
            vvs += (i == 0 ? "" : ",") + vs;
          }
          vvs += ")";

          RowMatrix<double> mat = ApplicationTools::getMatrixParameter<double>("probas", args, ',', vvs);

          FullHmmTransitionMatrix fhtm(pMA->getHmmStateAlphabet(), pMA->getNamespace());
          fhtm.setTransitionProbabilities(mat);

          pMA->matchParametersValues(fhtm.getParameters());

          nPL = pMA;
        }
        else if (phyloName == "AutoCorr")
        {
          AutoCorrelationOfAlignedPhyloLikelihood* pMA = new AutoCorrelationOfAlignedPhyloLikelihood(context, mPhylo, vPhylo);

          size_t nbP = pMA->getNumbersOfPhyloLikelihoods().size();

          string vs = "(" + VectorTools::paste(vector<double>(nbP, 1. / (double)nbP), ",") + ")";

          vector<double> v = ApplicationTools::getVectorParameter<double>("lambdas", args, ',', vs);

          ParameterList pl;

          for (size_t i = 0; i < v.size(); i++)
          {
            pl.addParameter(Parameter("AutoCorr.lambda" + TextTools::toString(i + 1), v[i]));
          }

          pMA->matchParametersValues(pl);

          nPL = pMA;
        }
        else if (phyloName == "Product")
        {
          ProductOfAlignedPhyloLikelihood* pAP = new ProductOfAlignedPhyloLikelihood(context, mPhylo, vPhylo);

          nPL = pAP;
        }
        else
          throw Exception("PhylogeneticsApplicationTools::getPhyloLikelihoodContainer : Unknown Phylo name " + phyloName);

        if (verbhere)
        {
          ApplicationTools::displayResult(" Phylolikelihood type", phyloName);
          ApplicationTools::displayResult(" Phylo numbers", VectorTools::paste(vPhylo, ","));
        }

        mPhylo->addPhyloLikelihood(phylonum, nPL);
        usedPhylo.push_back(phylonum);
      }
    }

    // Now clean the map
    for (map<size_t, vector<size_t> >::iterator it = phylosMap.begin(); it != phylosMap.end(); )
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
          vphyl.erase(vphyl.begin() + static_cast<ptrdiff_t> (i - 1));
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
    if (i!=0)
      sumAll+=" + ";
    
    sumAll += "phylo"+TextTools::toString(nPhyl[i]);
  }
  
  string resultDesc = ApplicationTools::getStringParameter("result", params, sumAll);

  // check if really formula, or previous phylo

  std::shared_ptr<PhyloLikelihood> nPL(0);
  size_t nP(0);
  bool flag(resultDesc.substr(0,5)=="phylo");

  if (flag)
  {
    try {
      nP=(size_t)TextTools::toInt(resultDesc.substr(5));
    }
    catch (Exception& e)
    {
      flag=false;
    }
  }

  if (!flag)
  {
    nPL = shared_ptr<PhyloLikelihood>(new FormulaOfPhyloLikelihood(context, mPhylo, resultDesc));
    if (verbose)
      ApplicationTools::displayResult(" Result", dynamic_cast<FormulaOfPhyloLikelihood*>(nPL.get())->output());
  }
  else
  {
    if (!mPhylo->hasPhyloLikelihood(nP))
      throw BadIntegerException("Unknown Phylolikelihood number for result",(int)nP);
    else
      nPL=mPhylo->getPhyloLikelihood(nP);
    if (verbose)
      ApplicationTools::displayResult(" Result", resultDesc);
  }

  mPhylo->sharePhyloLikelihood(0, nPL);
  return mPhylo;
}


/******************************************************/
/**** SUBSTITUTION MODEL SET **************************/
/******************************************************/

/******************************************************************************/

SubstitutionModelSet* PhylogeneticsApplicationTools::getSubstitutionModelSet(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const AlignedValuesContainer* data,
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
  for (size_t i = 0; nomix &(i < nbModels); i++)
  {
    string prefix = "model" + TextTools::toString(i + 1);
    string modelDesc;
    modelDesc = ApplicationTools::getStringParameter(prefix, params, "", suffix, suffixIsOptional, warn);

    if (modelDesc.find("Mixed") != string::npos)
      nomix = false;
  }

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
  const AlignedValuesContainer* data,
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
    rateFreqs = vector<double>(n, 1. / static_cast<double>(n));
    // Equal rates assumed for now, may be changed later
  }

  // ////////////////////////////////////
  // Deal with root frequencies

  map<string, string> unparsedParameters;

  bool stationarity = ApplicationTools::getBooleanParameter("nonhomogeneous.stationarity", params, false, "", true, warn);
  std::shared_ptr<FrequencySet> rootFrequencies(0);
  if (!stationarity)
  {
    rootFrequencies = getRootFrequencySet(alphabet, gCode, data, params, unparsedParameters, rateFreqs, suffix, suffixIsOptional, verbose);
    stationarity = !rootFrequencies;
    string freqDescription = ApplicationTools::getStringParameter("nonhomogeneous.root_freq", params, "", suffix, suffixIsOptional, warn);
    if (freqDescription.substr(0, 10) == "MVAprotein")
    {
      if (dynamic_cast<Coala*>(tmp.get()))
        dynamic_pointer_cast<MvaFrequencySet>(rootFrequencies)->initSet(dynamic_cast<CoalaCore*>(tmp.get()));
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
  const AlignedValuesContainer* data,
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

      const MixedTransitionModel* pSM = dynamic_cast<const MixedTransitionModel*>(mixedModelSet.getModel(static_cast<size_t> (num - 1)));
      if (pSM == NULL)
        throw BadIntegerException("PhylogeneticsApplicationTools::setMixedSubstitutionModelSet: Wrong model for number", num - 1);
      Vuint submodnb = pSM->getSubmodelNumbers(p2);

      mixedModelSet.addToHyperNode(static_cast<size_t> (num - 1), submodnb);
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
      classes.push_back(static_cast<size_t> (TextTools::toInt(strtok2.nextToken())));

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
  unique_ptr<DiscreteDistribution> rDist(bIO.readDiscreteDistribution(distDescription, true));

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
  OutputStream* messageHandler =
    (mhPath == "none") ? 0 :
    (mhPath == "std") ? ApplicationTools::message.get() :
    new StlOutputStream(new ofstream(mhPath.c_str(), ios::out));
  if (verbose)
    ApplicationTools::displayResult("Message handler", mhPath);

  string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional);
  OutputStream* profiler =
    (prPath == "none") ? 0 :
    (prPath == "std") ? ApplicationTools::message.get() :
    new StlOutputStream(new ofstream(prPath.c_str(), ios::out));
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
      else if (param == "*")
      {
        parametersToEstimate.reset();
        if (verbose)
          ApplicationTools::displayResult("Parameter ignored", string("All"));
      }
      else if (param.find("*") != string::npos)
      {
        vector<string> vs = ApplicationTools::matchingParameters(param, parNames);

        bool verbhere=verbose;
        
        if (vs.size()>=20)
        {
          if (verbose)
            ApplicationTools::displayResult("Number of parameters ignored", vs.size());
          verbhere=false;
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
      ApplicationTools::displayWarning("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
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
        NonHomogeneousTreeLikelihood* nhtl = dynamic_cast<NonHomogeneousTreeLikelihood*>(tl);
        if (!nhtl)
          ApplicationTools::displayWarning("The 'Ancient' parameters do not exist in homogeneous models, and will be ignored.");
        else
           parNames2 = nhtl->getRootFrequenciesParameters().getParameterNames();
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
      else if (param == "*")
        parNames2 = parToEstNames;
      else if (param.find("*") != string::npos)
        parNames2 = ApplicationTools::matchingParameters(param, parToEstNames);
      else
        parNames2.push_back(param);

      for (size_t i = 0; i < parNames2.size(); i++)
      {
        Parameter& par = parametersToEstimate.getParameter(parNames2[i]);
        if (par.hasConstraint())
        {
          par.setConstraint(std::shared_ptr<Constraint>(*ic & (*par.getConstraint())));
          if (par.getConstraint()->isEmpty())
            throw Exception("Empty interval for parameter " + parNames[i] + par.getConstraint()->getDescription());
        }
        else
          par.setConstraint(ic);

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
      if (abs(tl->getValue() - fval) > 0.000001)
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
      unsigned int topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, "", true, warn + 1);
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

  string finalMethod = ApplicationTools::getStringParameter("optimization.final", params, "none", suffix, suffixIsOptional, warn);
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
    string bf=backupFile+".def";
    rename(backupFile.c_str(),bf.c_str());
  }
  return tl;
}

/******************************************************************************/

PhyloLikelihood* PhylogeneticsApplicationTools::optimizeParameters(
  PhyloLikelihood* lik,
  const ParameterList& parameters,
  const map<string, string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  int warn)
{
  string optimization = ApplicationTools::getStringParameter("optimization", params, "FullD(derivatives=Newton)", suffix, suffixIsOptional, warn);
  if (optimization == "None")
    return lik;
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
    ApplicationTools::displayError("Sorry, optimization.scale_first not implemented yet for process.");
    exit(-1);
  }

  //     // We scale the tree before optimizing each branch length separately:
  //     if (verbose)
  //       ApplicationTools::displayMessage("Scaling the tree before optimizing each branch length separately.");
  //     double tolerance = ApplicationTools::getDoubleParameter("optimization.scale_first.tolerance", params, .0001, suffix, suffixIsOptional, true);
  //     if (verbose)
  //       ApplicationTools::displayResult("Scaling tolerance", TextTools::toString(tolerance));
  //     int nbEvalMax = ApplicationTools::getIntParameter("optimization.scale_first.max_number_f_eval", params, 1000000, suffix, suffixIsOptional, true);
  //     if (verbose)
  //       ApplicationTools::displayResult("Scaling max # f eval", TextTools::toString(nbEvalMax));

  //     OptimizationTools::optimizeTreeScale(
  //                                          tl,
  //                                          tolerance,
  //                                          nbEvalMax,
  //                                          messageHandler,
  //                                          profiler);
  //     if (verbose)
  //       ApplicationTools::displayResult("New tree likelihood", -tl->getValue());
  //   }

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
        vector<string> vs = lik->getBranchLengthParameters().getParameterNames();
        parametersToEstimate.deleteParameters(vs);
        if (verbose)
          ApplicationTools::displayResult("Parameter ignored", string("Branch lengths"));
      }
      else if (param == "Ancient")
      {
        vector<string> vs = lik->getRootFrequenciesParameters().getParameterNames();
        parametersToEstimate.deleteParameters(vs);
        if (verbose)
          ApplicationTools::displayResult("Parameter ignored", string("Root frequencies"));
      }
      else if (param == "Model")
      {
        vector<string> vs = lik->getSubstitutionModelParameters().getParameterNames();
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

        bool verbhere=verbose;
        
        if (vs.size()>=20)
        {
          if (verbose)
            ApplicationTools::displayResult("Number of parameters ignored", vs.size());
          verbhere=false;
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
      std::shared_ptr<IntervalConstraint> ic(new IntervalConstraint(constraint));

      vector<string> parNames2;

      if (param == "BrLen")
        parNames2  = lik->getBranchLengthParameters().getParameterNames();
      else if (param == "Ancient")
        parNames2 = lik->getRootFrequenciesParameters().getParameterNames();
      else if (param == "Model")
      {
        vector<string> vs = lik->getSubstitutionModelParameters().getParameterNames();
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
          par.setConstraint(std::shared_ptr<Constraint>(*ic & (*par.getConstraint())));
          if (par.getConstraint()->isEmpty())
            throw Exception("Empty interval for parameter " + parNames[i] + par.getConstraint()->getDescription());
        }
        else
          par.setConstraint(ic);

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
      ParameterList pl = lik->getParameters();
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
      lik->setParameters(pl);
      if (abs(lik->getValue() - fval) > 0.000001)
        ApplicationTools::displayMessage("Changed likelihood from backup file.");
      ApplicationTools::displayResult("Restoring log-likelihood", -lik->getValue());
    }
  }

  // There it goes...
  bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", params, false, suffix, suffixIsOptional, warn + 1);
  if (optimizeTopo)
    throw Exception("Topology opmitization not implemented yet for processes");

  // if (verbose)
  //   ApplicationTools::displayResult("Optimize topology", optimizeTopo ? "yes" : "no");
  // string nniMethod = ApplicationTools::getStringParameter("optimization.topology.algorithm_nni.method", params, "phyml", suffix, suffixIsOptional, warn + 1);
  // string nniAlgo;
  // if (nniMethod == "fast")
  //   {
  //     nniAlgo = NNITopologySearch::FAST;
  //   }
  // else if (nniMethod == "better")
  //   {
  //     nniAlgo = NNITopologySearch::BETTER;
  //   }
  // else if (nniMethod == "phyml")
  //   {
  //     nniAlgo = NNITopologySearch::PHYML;
  //   }
  // else
  //   throw Exception("Unknown NNI algorithm: '" + nniMethod + "'.");


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

    // if (optimizeTopo)
    //   {
    //     bool optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, warn + 1);
    //     unsigned int topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, "", true, warn + 1);
    //     double tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional);
    //     double tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, suffixIsOptional);
    //     tl = OptimizationTools::optimizeTreeNNI(
    //                                             dynamic_cast<NNIHomogeneousTreeLikelihood*>(tl), parametersToEstimate,
    //                                             optNumFirst, tolBefore, tolDuring, nbEvalMax, topoNbStep, messageHandler, profiler,
    //                                             reparam, optVerbose, optMethodDeriv, nstep, nniAlgo);
    //   }

    if (verbose && nstep > 1)
      ApplicationTools::displayResult("# of precision steps", TextTools::toString(nstep));
    parametersToEstimate.matchParametersValues(lik->getParameters());
    n = OptimizationTools::optimizeNumericalParameters(
      lik, parametersToEstimate,
      backupListener.get(), nstep, tolerance, nbEvalMax, messageHandler, profiler, reparam, optVerbose, optMethodDeriv, optMethodModel);
  }
  else if (optName == "FullD")
  {
    // Uses Newton-raphson algorithm with numerical derivatives when required.
    parametersToEstimate.matchParametersValues(lik->getParameters());
    if (dynamic_cast<SingleProcessPhyloLikelihood*>(lik))
      n = OptimizationTools::optimizeNumericalParameters2(
        *dynamic_cast<SingleProcessPhyloLikelihood*>(lik), parametersToEstimate,
        backupListener.get(), tolerance, nbEvalMax, messageHandler, profiler, reparam, useClock, optVerbose, optMethodDeriv);
    else
      n = OptimizationTools::optimizeNumericalParameters2(
        lik, parametersToEstimate,
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
    finalOptimizer = new DownhillSimplexMethod(lik);
  }
  else if (finalMethod == "powell")
  {
    finalOptimizer = new PowellMultiDimensions(lik);
  }
  else
    throw Exception("Unknown final optimization method: " + finalMethod);

  if (finalOptimizer)
  {
    parametersToEstimate.matchParametersValues(lik->getParameters());
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
    string bf=backupFile+".def";
    rename(backupFile.c_str(),bf.c_str());
  }
  return lik;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::optimizeParameters(
  DiscreteRatesAcrossSitesClockTreeLikelihood* tl,
  const ParameterList& parameters,
  const map<string, string>& params,
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
        ApplicationTools::displayMessage("Changed likelihood from backup file.");
      ApplicationTools::displayResult("Restoring log-likelihood", -tl->getValue());
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

  if (prPath != "none" && prPath != "std")
    delete profiler;
  if (mhPath != "none" && mhPath != "std")
    delete messageHandler;

  if (verbose)
    ApplicationTools::displayResult("Performed", TextTools::toString(n) + " function evaluations.");
  if (backupFile != "none")
  {
    string bf = backupFile + ".def";
    rename(backupFile.c_str(), bf.c_str());
  }
}

/******************************************************************************/

void PhylogeneticsApplicationTools::checkEstimatedParameters(const ParameterList& pl)
{
  for (size_t i = 0; i < pl.size(); ++i)
  {
    std::shared_ptr<Constraint> constraint = pl[i].getConstraint();
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

void PhylogeneticsApplicationTools::writeTrees(
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
    throw Exception("Unknow format for tree writing: " + format);

  if (!checkOnly)
    treeWriter->writeTrees(trees, file, true);

  delete treeWriter;
  if (verbose)
    ApplicationTools::displayResult("Wrote trees to file ", file);
}

void PhylogeneticsApplicationTools::writeTrees(
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
  {
    treeWriter->writeTrees(trees, file, true);
    
    if (verbose)
      ApplicationTools::displayResult("Wrote trees to file ", file);
  }
  
  delete treeWriter;
}

void PhylogeneticsApplicationTools::writeTrees(
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

  OTree* treeWriter;
  if (format == "Newick")
    treeWriter = new Newick();
  else if (format == "Nexus")
    treeWriter = new NexusIOTree();
  else if (format == "NHX")
    treeWriter = new Nhx();
  else
    throw Exception("Unknow format for tree writing: " + format);

  if (!checkOnly)
  {
    vector<size_t> vTN = spc.getTreeNumbers();

    for (size_t i = 0; i < vTN.size(); i++)
    {
      PhyloTree tree(spc.getTree(vTN[i]));

      std::vector<shared_ptr<PhyloNode> > nodes = tree.getAllNodes();

      for (auto& node : nodes)
      {
        if (tree.isLeaf(node) && withIds)
          node->setName(TextTools::toString(tree.getNodeIndex(node)) + "_" + node->getName());
        else
          node->setProperty("NodeId", BppString(TextTools::toString(tree.getNodeIndex(node))));
      }

      Newick* nt=dynamic_cast<Newick*>(treeWriter);
      if (nt)
        nt->enableExtendedBootstrapProperty("NodeId");

      treeWriter->writeTree(tree, file + "_" + TextTools::toString(vTN[i]), true);
    }
    if (verbose)
      ApplicationTools::displayResult("Wrote trees to files : ", file + "_...");
  }

  delete treeWriter;
}

void PhylogeneticsApplicationTools::printParameters(const BranchModel* model, OutputStream& out, int warn)
{
  out << "model=";
  map<string, string> globalAliases;
  vector<string> writtenNames;
  BppOBranchModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, false, warn);
  bIO.write(*model, out, globalAliases, writtenNames);
  out.endLine();
}

void PhylogeneticsApplicationTools::printParameters(const SubstitutionProcess* process, OutputStream& out, int warn)
{
  if (dynamic_cast<const SimpleSubstitutionProcess*>(process) != NULL)
  {
    (out << "nonhomogeneous=no").endLine();

    out << "model=";
    map<string, string> globalAliases;
    vector<string> writtenNames;
    BppOBranchModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, false, warn);
    bIO.write(*process->getModel(0, 0), out, globalAliases, writtenNames);
    out.endLine();
  }

  else if (dynamic_cast<const RateAcrossSitesSubstitutionProcess*>(process) != NULL)
  {
    const RateAcrossSitesSubstitutionProcess* pRA = dynamic_cast<const RateAcrossSitesSubstitutionProcess*>(process);

    (out << "nonhomogeneous=no").endLine();

    out << "model=";
    map<string, string> globalAliases;
    vector<string> writtenNames;
    BppOBranchModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, false, warn);
    bIO.write(*process->getModel(0, 0), out, globalAliases, writtenNames);
    out.endLine();
    out.endLine();

    // Rate distribution

    out << "rate_distribution=";
    const BppORateDistributionFormat* bIOR = new BppORateDistributionFormat(true);
    bIOR->writeDiscreteDistribution(*pRA->getRateDistribution(), out, globalAliases, writtenNames);
    delete bIOR;
    out.endLine();
  }

  else if (dynamic_cast<const NonHomogeneousSubstitutionProcess*>(process) != NULL)
  {
    const NonHomogeneousSubstitutionProcess* pNH = dynamic_cast<const NonHomogeneousSubstitutionProcess*>(process);

    (out << "nonhomogeneous=general").endLine();
    (out << "nonhomogeneous.number_of_models=" << pNH->getNumberOfModels()).endLine();

    vector<string> writtenNames;

    // Loop over all models:
    for (size_t i = 0; i < pNH->getNumberOfModels(); i++)
    {
      const auto model = pNH->getModel(i);

      // First get the aliases for this model:
      map<string, string> aliases;

      ParameterList pl = model->getParameters();

      for (size_t np = 0; np < pl.size(); np++)
      {
        string nfrom = pNH->getFrom(pl[np].getName() + "_" + TextTools::toString(i + 1));
        if (nfrom != "")
          aliases[pl[np].getName()] = nfrom;
      }

      // Now print it:
      writtenNames.clear();
      out.endLine() << "model" << (i + 1) << "=";
      BppOBranchModelFormat bIOsm(BppOSubstitutionModelFormat::ALL, true, true, true, false, warn);
      bIOsm.write(*model, out, aliases, writtenNames);
      out.endLine();
      vector<unsigned int> ids = pNH->getNodesWithModel(i);
      out << "model" << (i + 1) << ".nodes_id=" << ids[0];
      for (size_t j = 1; j < ids.size(); ++j)
      {
        out << "," << ids[j];
      }
      out.endLine();
    }

    // Root frequencies:
    out.endLine();
    if (pNH->getRootFrequencySet())
    {
      out << "nonhomogeneous.root_freq=";

      map<string, string> aliases;

      ParameterList pl = pNH->getRootFrequencySet()->getParameters();

      for (size_t np = 0; np < pl.size(); np++)
      {
        string nfrom = pNH->getFrom(pl[np].getName());
        if (nfrom != "")
          aliases[pl[np].getName()] = nfrom;
      }

      BppOFrequencySetFormat bIO(BppOFrequencySetFormat::ALL, false, warn);
      bIO.writeFrequencySet(pNH->getRootFrequencySet().get(), out, aliases, writtenNames);
    }
    else
      out << "nonhomogeneous.stationarity=true";
    out.endLine();

    // Rate distribution

    map<string, string> aliases;
    const DiscreteDistribution* pdd = pNH->getRateDistribution();

    ParameterList pl = pdd->getParameters();
    for (size_t np = 0; np < pl.size(); np++)
    {
      string nfrom = pNH->getFrom(pl[np].getName());
      if (nfrom != "")
        aliases[pl[np].getName()] = nfrom;
    }
    out.endLine();
    out << "rate_distribution=";
    const BppORateDistributionFormat* bIO = new BppORateDistributionFormat(true);
    bIO->writeDiscreteDistribution(*pdd, out, aliases, writtenNames);
    delete bIO;
    out.endLine();
  }
}

void PhylogeneticsApplicationTools::printParameters(const SubstitutionProcessCollection* collection, OutputStream& out, int warn, bool withAlias)
{
  vector<string> writtenNames;

  // The models
  vector<size_t> vModN = collection->getModelNumbers();

  for (auto modn : vModN)
  {
    const auto& model = *collection->getModel(modn);

    // First get the aliases for this model:
    map<string, string> aliases;

    if (withAlias)
    {
      ParameterList pl = model.getParameters();

      for (size_t np = 0; np < pl.size(); np++)
      {
        string nfrom = collection->getFrom(pl[np].getName() + "_" + TextTools::toString(modn));
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
  }

  // Root frequencies:
  vector<size_t> rootFreqN = collection->getFrequenciesNumbers();

  for (size_t i = 0; i < rootFreqN.size(); i++)
  {
    auto rootFreq = collection->shareFrequencies(rootFreqN[i]);

    // Now print it:
    writtenNames.clear();
    out.endLine() << "root_freq" << rootFreqN[i] << "=";
    BppOFrequencySetFormat bIOf(BppOFrequencySetFormat::ALL, true, warn);

    map<string, string> aliases;

    if (withAlias)
    {
      ParameterList pl = rootFreq->getParameters();

      for (size_t np = 0; np < pl.size(); np++)
      {
        string nfrom = collection->getFrom(pl[np].getName() + "_" + TextTools::toString(rootFreqN[i]));
        if (nfrom != "")
          aliases[pl[np].getName()] = nfrom;
      }
    }

    bIOf.writeFrequencySet(rootFreq.get(), out, aliases, writtenNames);
    out.endLine();
  }

  // Rate distribution

  vector<size_t> vDistN = collection->getRateDistributionNumbers();

  for (auto distn : vDistN)
  {
    if (distn < 10000)
    {
      const DiscreteDistribution& dist = collection->getRateDistribution(distn);

      // First get the aliases for this model:
      map<string, string> aliases;

      if (withAlias)
      {
        ParameterList pl = dist.getParameters();

        for (size_t np = 0; np < pl.size(); np++)
        {
          string nfrom = collection->getFrom(pl[np].getName() + "_" + TextTools::toString(distn));
          if (nfrom != "")
            aliases[pl[np].getName()] = nfrom;
        }
      }

      // Now print it:
      writtenNames.clear();
      out.endLine() << "rate_distribution" << distn << "=";
      BppORateDistributionFormat bIOd(true);
      bIOd.writeDiscreteDistribution(dist, out, aliases, writtenNames);
      out.endLine();
    }
  }

  // scenarios

  vector<size_t> vSce = collection->getScenarioNumbers();

  if (vSce.size()>0)
    out.endLine();

  vector<const ModelPath*> vMP;

  // first output the scenarios 
  for (const auto& scennum : vSce)
  {
    const auto& scen = collection->getModelScenario(scennum);

    out.endLine();

    out << "scenario" << scennum << "=";

    size_t nbMP=scen.getNumberOfModelPaths();
    
    for (size_t mpn = 0; mpn < nbMP; mpn++)
    {
      const auto& mp = scen.getModelPath(mpn);

      auto itmp=find(vMP.begin(),vMP.end(),mp.get());
      auto inmp=std::distance(vMP.begin(), itmp);
      if (itmp==vMP.end())
        vMP.push_back(mp.get());

      if (mpn!=0)
        out << "&";
      out << "path" << TextTools::toString(inmp+1);
    }
    out.endLine();
  }

  // then the model path 
  for (size_t inmp = 0; inmp < vMP.size(); inmp++)
  {
    out.endLine();
    out << "path" << inmp+1 << "=";

    const ModelPath& mp = *vMP[inmp];

    auto vMod = mp.getModels();

    bool dem=true;
    for (const auto& mod:vMod)
    {
      // look for model number in collection
      size_t modN=collection->getModelIndex(mod);

      if (!dem)
      {
        out << "&";
        dem=false;
      }
      
      out << "model" << modN;
      out << "[" << mp.getPathNode(mod).to_string() <<  "]";
    }
    out.endLine();
  }

  // processes
  out.endLine();

  vector<size_t> vprocN = collection->getSubstitutionProcessNumbers();

  for (size_t i = 0; i < vprocN.size(); i++)
  {
    const auto& spcm = collection->getSubstitutionProcess(vprocN[i]);

    out << "process" << vprocN[i] << "=";

    if (spcm.getNumberOfModels() == 1)
      out << "Homogeneous(model=" << spcm.getModelNumbers()[0];
    else
    {
      out << "Nonhomogeneous(";
      vector<size_t> vMN = spcm.getModelNumbers();
      for (size_t j = 0; j < vMN.size(); j++)
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

    if (spcm.hasModelScenario())
      out << ", scenario=" << spcm.getModelScenarioNumber();
    
    out << ")";
    out.endLine();
    out.endLine();
  }
}


void PhylogeneticsApplicationTools::printParameters(const PhyloLikelihoodContainer& phylocont, OutputStream& out, int warn)
{
  out << "# Log likelihood = ";

  const PhyloLikelihood* result = phylocont[0];

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

  const FormulaOfPhyloLikelihood* pop = dynamic_cast<const FormulaOfPhyloLikelihood*>(result);

  vector<size_t> phyldep;

  if (!pop)
  {  
    out << "phylo1";
    phyldep.push_back(1);
  }
  else
  {
    string popout=pop->output();

    out << popout;
  
    StringTokenizer st(popout,"phylo",true, true);
    st.nextToken();
    
    
    while (st.hasMoreToken())
    {
      string ex=st.nextToken();
      phyldep.push_back((size_t)(atoi(ex.c_str())));
    }
  }
  
  out.endLine();
  out.endLine();
  
  // Then the other phylolikelihoods

  while (phyldep.size() != 0)
  {
    size_t num = phyldep[0];
    const PhyloLikelihood* phyloLike = phylocont[num];

    // remove phylolikelihoods with this number
    vector<size_t>::iterator itf = find(phyldep.begin(), phyldep.end(), num);
    while (itf != phyldep.end())
    {
      phyldep.erase(itf);
      itf = find(itf, phyldep.end(), num);
    }


    // then output

    if (dynamic_cast<const SingleDataPhyloLikelihood*>(phyloLike) != NULL)
      printParameters(dynamic_cast<const SingleDataPhyloLikelihood&>(*phyloLike), out, num, warn);
    else
    {
      out << "phylo" << num << "=";

      const SetOfAbstractPhyloLikelihood* mDP = dynamic_cast<const SetOfAbstractPhyloLikelihood*>(phyloLike);
      if (mDP)
      {
        if (dynamic_cast<const MixtureOfAlignedPhyloLikelihood*>(phyloLike) != NULL)
        {
          const MixtureOfAlignedPhyloLikelihood* pM = dynamic_cast<const MixtureOfAlignedPhyloLikelihood*>(phyloLike);

          out << "Mixture(probas=(" << VectorTools::paste(pM->getPhyloProbabilities(), ",");

          out << "),";
        }

        else if (dynamic_cast<const HmmOfAlignedPhyloLikelihood*>(phyloLike) != NULL)
        {
          const HmmOfAlignedPhyloLikelihood* pM = dynamic_cast<const HmmOfAlignedPhyloLikelihood*>(phyloLike);
          out << "HMM(probas=";

          RowMatrix<double> tMt;
          copyEigenToBpp(pM->getHmmTransitionMatrix(), tMt);
          MatrixTools::print(tMt, out);

          out << ",";
        }

        else if (dynamic_cast<const AutoCorrelationOfAlignedPhyloLikelihood*>(phyloLike) != NULL)
        {
          const AutoCorrelationOfAlignedPhyloLikelihood* pM = dynamic_cast<const AutoCorrelationOfAlignedPhyloLikelihood*>(phyloLike);

          out << "AutoCorr(lambdas=(";

          Vdouble vP;
          for (unsigned int i = 0; i < pM->getHmmTransitionMatrix().cols(); i++)
          {
            vP.push_back(pM->getHmmTransitionMatrix()(i, i));
          }

          out << VectorTools::paste(vP, ",");

          out << "),";
        }
        else if (dynamic_cast<const ProductOfAlignedPhyloLikelihood*>(phyloLike) != NULL)
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


void PhylogeneticsApplicationTools::printParameters(const SingleDataPhyloLikelihood& phyloLike, OutputStream& out, size_t nPhylo, int warn)
{
  out << "phylo" << TextTools::toString(nPhylo) << "=";

  out << "Single(";

  if (dynamic_cast<const SequencePhyloLikelihood*>(&phyloLike) != NULL)
  {
    const SequencePhyloLikelihood* pMP = dynamic_cast<const SequencePhyloLikelihood*>(&phyloLike);

    out << "process=" << pMP->getSequenceEvolutionNumber();
  }
  else
  {
    const SingleProcessPhyloLikelihood* pS = dynamic_cast<const SingleProcessPhyloLikelihood*>(&phyloLike);

    if (pS)
      out << "process=" << pS->getSubstitutionProcessNumber();
  }

  out << ",data=" << TextTools::toString(phyloLike.getNData()) << ")";
  out.endLine();
}

void PhylogeneticsApplicationTools::printParameters(const SequenceEvolution* evol, OutputStream& out, size_t nEvol, int warn)
{
  out << "process" << TextTools::toString(nEvol) << "=";

  if (dynamic_cast<const OneProcessSequenceEvolution*>(evol) != NULL)
  {
    const OneProcessSequenceEvolution* pOP = dynamic_cast<const OneProcessSequenceEvolution*>(evol);

    out << "Simple(process=" <<  pOP->getSubstitutionProcessNumber() << ")";
  }
  else if (dynamic_cast<const MultiProcessSequenceEvolution*>(evol) != NULL)
  {
    const MultiProcessSequenceEvolution* pMP = dynamic_cast<const MultiProcessSequenceEvolution*>(evol);

    if (dynamic_cast<const MixtureSequenceEvolution*>(evol) != NULL)
    {
      const MixtureSequenceEvolution* pM = dynamic_cast<const MixtureSequenceEvolution*>(evol);

      out << "Mixture(probas=(" << VectorTools::paste(pM->getSubProcessProbabilities(), ",");
      out << "),";
    }

    else if (dynamic_cast<const HmmSequenceEvolution*>(evol) != NULL)
    {
      const HmmSequenceEvolution* pM = dynamic_cast<const HmmSequenceEvolution*>(evol);
      out << "HMM(probas=";

      const Matrix<double>& tMt = pM->getHmmTransitionMatrix().getPij();
      MatrixTools::print(tMt, out);

      out << ",";
    }
    else if (dynamic_cast<const AutoCorrelationSequenceEvolution*>(evol) != NULL)
    {
      const AutoCorrelationSequenceEvolution* pM = dynamic_cast<const AutoCorrelationSequenceEvolution*>(evol);

      out << "AutoCorr(lambdas=(";

      Vdouble vP;
      for (unsigned int i = 0; i < pM->getNumberOfSubstitutionProcess(); i++)
      {
        vP.push_back(pM->getHmmTransitionMatrix().Pij(i, i));
      }

      out << VectorTools::paste(vP, ",");

      out << "),";
    }
    else if (dynamic_cast<const PartitionSequenceEvolution*>(evol) != NULL)
    {
      const PartitionSequenceEvolution* pM = dynamic_cast<const PartitionSequenceEvolution*>(evol);

      out << "Partition(";

      const map<size_t, vector<size_t> >& mProcPos = pM->getMapOfProcessSites();

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
  const PhyloLikelihood* result = phylocont[0];

  if (!result)
    return;

  vector<size_t> phyldep = phylocont.getNumbersOfPhyloLikelihoods();
  
  while (phyldep.size() != 0)
  {
    size_t num = phyldep[0];

    phyldep.erase(phyldep.begin());

    const PhyloLikelihood* phyloLike = phylocont[num];
    // output

    string info_out = infosFile + "_" + TextTools::toString(num);

    if (dynamic_cast<const SingleDataPhyloLikelihood*>(phyloLike) != NULL
        && num!=0 )
      printAnalysisInformation(dynamic_cast<const SingleDataPhyloLikelihood&>(*phyloLike), info_out, warn);
    else
    {
      const SetOfAlignedPhyloLikelihood* sOAP = dynamic_cast<const SetOfAlignedPhyloLikelihood*>(phyloLike);
      if (sOAP != NULL)
      {
        if (num!=0)
          printAnalysisInformation(*sOAP, info_out, warn);

        vector<size_t> vPN = sOAP->getNumbersOfPhyloLikelihoods();

        // update phyldep
        phyldep.assign(vPN.begin(),vPN.end());
        phyldep=VectorTools::unique(phyldep);
      }
      else
      {
        const SetOfAbstractPhyloLikelihood* sOAB = dynamic_cast<const SetOfAbstractPhyloLikelihood*>(phyloLike);
        if (sOAB != NULL)
        {
          vector<size_t> vPN = sOAB->getNumbersOfPhyloLikelihoods();

          // update phyldep
          phyldep.assign(vPN.begin(),vPN.end());
          phyldep=VectorTools::unique(phyldep);
        }
      }
    }
  }
}


void PhylogeneticsApplicationTools::printAnalysisInformation(const SetOfAlignedPhyloLikelihood& sOAP, const string& infosFile, int warn)
{
  const MixtureOfAlignedPhyloLikelihood* mOAP = NULL;
  const HmmOfAlignedPhyloLikelihood* hOAP = NULL;
  const AutoCorrelationOfAlignedPhyloLikelihood* aCOAP = NULL;

  vector<size_t> phyloNum = sOAP.getNumbersOfPhyloLikelihoods();
  size_t nbP = phyloNum.size();

  if (dynamic_cast<const ProductOfAlignedPhyloLikelihood*>(&sOAP) == NULL)
  {
    StlOutputStream out(new ofstream(infosFile.c_str(), ios::out));

    if (dynamic_cast<const MixtureOfAlignedPhyloLikelihood*>(&sOAP) != NULL)
      mOAP = dynamic_cast<const MixtureOfAlignedPhyloLikelihood*>(&sOAP);
    else if (dynamic_cast<const HmmOfAlignedPhyloLikelihood*>(&sOAP) != NULL)
      hOAP = dynamic_cast<const HmmOfAlignedPhyloLikelihood*>(&sOAP);
    else if (dynamic_cast<const AutoCorrelationOfAlignedPhyloLikelihood*>(&sOAP) != NULL)
      aCOAP = dynamic_cast<const AutoCorrelationOfAlignedPhyloLikelihood*>(&sOAP);

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

void PhylogeneticsApplicationTools::printAnalysisInformation(const SingleDataPhyloLikelihood& phyloLike, const string& infosFile, int warn)
{
  if (dynamic_cast<const SingleProcessPhyloLikelihood*>(&phyloLike) != NULL)
  {
    auto pSPL = dynamic_cast<const SingleProcessPhyloLikelihood*>(&phyloLike);
    
    StlOutputStream out(new ofstream(infosFile.c_str(), ios::out));
    
    const SubstitutionProcess* pSP = &pSPL->getSubstitutionProcess();
    
    vector<string> colNames;
    colNames.push_back("Sites");
    colNames.push_back("is.complete");
    colNames.push_back("is.constant");
    colNames.push_back("lnL");
    
    const DiscreteDistribution* pDD = pSP->getRateDistribution();
    size_t nbR = 0;
    
    if (pDD != NULL)
    {
      nbR = pDD->getNumberOfCategories();

      if (nbR > 1)
        for (size_t i = 0; i < nbR; i++)
          colNames.push_back("prob" + TextTools::toString(i + 1));
    }

    const AlignedValuesContainer* sites = phyloLike.getData();
    
    vector<string> row(4 + (nbR > 1 ? nbR : 0));
    DataTable* infos = new DataTable(colNames);

    VVdouble vvPP(pSPL->getPosteriorProbabilitiesPerSitePerClass());

    for (size_t i = 0; i < sites->getNumberOfSites(); i++)
    {
      double lnL = phyloLike.getLogLikelihoodForASite(i);
      
      const CruxSymbolListSite& currentSite = sites->getSymbolListSite(i);
      int currentSitePosition = currentSite.getPosition();
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
      row[0] = (string("[" + TextTools::toString(currentSitePosition) + "]"));
      row[1] = isCompl;
      row[2] = isConst;
      row[3] = TextTools::toString(lnL);
      
      if (nbR > 1)
        for (size_t j = 0; j < nbR; j++)
        {
          row[4 + j] = TextTools::toString(vvPP[i][j]);
        }
      
      infos->addRow(row);
    }

    DataTable::write(*infos, out, "\t");
    delete infos;
  }
  else if (dynamic_cast<const PartitionProcessPhyloLikelihood*>(&phyloLike) != NULL)
  {
    const PartitionProcessPhyloLikelihood* pPPL = dynamic_cast<const PartitionProcessPhyloLikelihood*>(&phyloLike);
    
    const PartitionSequenceEvolution& pSE=dynamic_cast<const PartitionSequenceEvolution&>(pPPL->getSequenceEvolution());

    const map<size_t, vector<size_t> >& mProcPos=pSE.getMapOfProcessSites();
    
    vector<size_t> nbProc=pSE.getSubstitutionProcessNumbers();
    
    map<size_t, size_t> mNbr;

    for (auto nP : nbProc)
    {
      const SubstitutionProcess& sp=pSE.getSubstitutionProcess(nP);
      const DiscreteDistribution* pDD = sp.getRateDistribution();
      mNbr[nP]=(pDD?pDD->getNumberOfCategories():1);
    }

    size_t maxR=max_element(mNbr.begin(), mNbr.end(), [](const std::pair<size_t, size_t>& p1, const std::pair<size_t, size_t>& p2){ return p1.second < p2.second;})->second;

    StlOutputStream out(new ofstream(infosFile.c_str(), ios::out));
    
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
    DataTable* infos = new DataTable(nbSites,colNames);
    
    for (auto nP : nbProc)
    {
      auto pSPPL = dynamic_cast<const SingleProcessPhyloLikelihood*>(pPPL->getAbstractPhyloLikelihood(nP));

      if (!pSPPL)
        throw Exception("PhylogeneticsApplicationTools::printAnalysisInformation : no SingleProcessPhyloLikelihood in PartitionProcessPhyloLikelihood.");
      
      size_t nbr=mNbr[pSPPL->getSubstitutionProcessNumber()];

      const vector<size_t>& mPos=mProcPos.at(nP);
      
      const AlignedValuesContainer* sites = pSPPL->getData();

      for (size_t i = 0; i < sites->getNumberOfSites(); i++)
      {
        double lnL = pSPPL->getLogLikelihoodForASite(i);

        const CruxSymbolListSite& currentSite = sites->getSymbolListSite(i);
        int currentSitePosition = currentSite.getPosition();
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
        row[0] = (string("[" + TextTools::toString(currentSitePosition) + "]"));
        row[1] = isCompl;
        row[2] = isConst;
        row[3] = TextTools::toString(lnL);
      
        if (nbr>1)
        {
          Vdouble vPP = pSPPL->getPosteriorProbabilitiesForSitePerClass(i);

          for (size_t j = 0; j < nbr; j++)
            row[4 + j] = TextTools::toString(vPP[j]);
        }

        for (size_t j = nbr; j<maxR; j++)
          row[4 + j] = "NA";

        infos->setRow(mPos[i],row);
      }
    }
    
    DataTable::write(*infos, out, "\t");
    delete infos;
  }
  else if (dynamic_cast<const MultiProcessSequencePhyloLikelihood*>(&phyloLike) != NULL)
  {
    const MultiProcessSequencePhyloLikelihood* pMPL = dynamic_cast<const MultiProcessSequencePhyloLikelihood*>(&phyloLike);

    StlOutputStream out(new ofstream(infosFile.c_str(), ios::out));

    vector<string> colNames;
    colNames.push_back("Sites");
    colNames.push_back("is.complete");
    colNames.push_back("is.constant");
    colNames.push_back("lnL");

    size_t nbP = pMPL->getNumberOfSubstitutionProcess();

    if (nbP > 1)
    {
      for (size_t i = 0; i < nbP; i++)
      {
        colNames.push_back("lnL" + TextTools::toString(i + 1));
      }
      for (size_t i = 0; i < nbP; i++)
      {
        colNames.push_back("prob" + TextTools::toString(i + 1));
      }
    }

    const AlignedValuesContainer* sites = phyloLike.getData();
    
    vector<string> row(4 + (nbP > 1 ? 2 * nbP : 0));
    DataTable* infos = new DataTable(colNames);

    VVdouble vvPP = pMPL->getPosteriorProbabilitiesPerSitePerProcess();
    VVdouble vvL = pMPL->getLikelihoodPerSitePerProcess();

    for (size_t i = 0; i < sites->getNumberOfSites(); i++)
    {
      double lnL = phyloLike.getLogLikelihoodForASite(i);
      const CruxSymbolListSite& currentSite = sites->getSymbolListSite(i);
      int currentSitePosition = currentSite.getPosition();
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
      row[0] = (string("[" + TextTools::toString(currentSitePosition) + "]"));
      row[1] = isCompl;
      row[2] = isConst;
      row[3] = TextTools::toString(lnL);

      if (nbP > 1)
      {
        for (size_t j = 0; j < nbP; j++)
        {
          row[4 + j] = TextTools::toString(log(vvL[i][j]));
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
  }
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
    const BranchModel* model = modelSet->getModel(i);

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
    bIO.writeFrequencySet(pFS.get(), out, aliases, writtenNames);
  }
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printParameters(const DiscreteDistribution* rDist, OutputStream& out, bool withAlias)
{
  out << "rate_distribution=";
  map<string, string> globalAliases;
  vector<string> writtenNames;
  const BppORateDistributionFormat* bIO = new BppORateDistributionFormat(true);
  bIO->writeDiscreteDistribution(*rDist, out, globalAliases, writtenNames);
  delete bIO;
  out.endLine();
}

/************************
* Substitution Mapping *
************************/
SubstitutionCount* PhylogeneticsApplicationTools::getSubstitutionCount(
  const Alphabet* alphabet,
  const SubstitutionModel* model,
  const map<string, string>& params,
  string suffix,
  bool verbose,
  int warn)
{
  const StateMap& stateMap=model->getStateMap();
  
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
    shared_ptr<const AlphabetIndex2> weights(SequenceApplicationTools::getAlphabetIndex2(alphabet, weightOption, "Substitution weight scheme:"));
    string distanceOption = ApplicationTools::getStringParameter("distance", nijtParams, "None", "", true, warn + 1);
    shared_ptr<const AlphabetIndex2> distances(SequenceApplicationTools::getAlphabetIndex2(alphabet, distanceOption, "Substitution distances:"));
    substitutionCount = new UniformizationSubstitutionCount(model, new TotalSubstitutionRegister(stateMap), weights, distances);
  }
  else if (nijtOption == "Decomposition")
  {
    string weightOption = ApplicationTools::getStringParameter("weight", nijtParams, "None", "", true, warn + 1);
    shared_ptr<const AlphabetIndex2> weights(SequenceApplicationTools::getAlphabetIndex2(alphabet, weightOption, "Substitution weight scheme:"));
    string distanceOption = ApplicationTools::getStringParameter("distance", nijtParams, "None", "", true, warn + 1);
    shared_ptr<const AlphabetIndex2> distances(SequenceApplicationTools::getAlphabetIndex2(alphabet, distanceOption, "Substitution distances:"));
    const ReversibleSubstitutionModel* revModel = dynamic_cast<const ReversibleSubstitutionModel*>(model);
    if (revModel)
      substitutionCount = new DecompositionSubstitutionCount(revModel, new TotalSubstitutionRegister(stateMap), weights, distances);
    else
      throw Exception("Decomposition method can only be used with reversible substitution models.");
  }
  else if (nijtOption == "Naive")
  {
    string weightOption = ApplicationTools::getStringParameter("weight", nijtParams, "None", "", true, warn + 1);
    std::shared_ptr<const AlphabetIndex2> weights(SequenceApplicationTools::getAlphabetIndex2(alphabet, weightOption, "Substitution weight scheme:"));
    string distanceOption = ApplicationTools::getStringParameter("distance", nijtParams, "", "", true, warn + 1);
    if (distanceOption!="")
      ApplicationTools::displayMessage("Naive substitution count: distances not handled");
    
    substitutionCount = new NaiveSubstitutionCount(model, new TotalSubstitutionRegister(stateMap), false, weights);
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


SubstitutionRegister* PhylogeneticsApplicationTools::getSubstitutionRegister(const string& regTypeDesc, const StateMap& stateMap, const GeneticCode* genCode, AlphabetIndex2*& weights, AlphabetIndex2*& distances, bool verbose)
{
  string regType = "";
  map<string, string> regArgs;
  KeyvalTools::parseProcedure(regTypeDesc, regType, regArgs);

  const Alphabet* alphabet=stateMap.getAlphabet();
  
  SubstitutionRegister* reg = 0;
  weights = 0;
  distances = 0;
  
  string weightOption = ApplicationTools::getStringParameter("weight", regArgs, "None", "", true, 1);
  string distanceOption = ApplicationTools::getStringParameter("distance", regArgs, "None", "", true, 1);

  if (AlphabetTools::isCodonAlphabet(alphabet))
  {
    weights = SequenceApplicationTools::getAlphabetIndex2(dynamic_cast<const CodonAlphabet*>(alphabet), genCode, weightOption, "Substitution weight scheme:");
    distances = SequenceApplicationTools::getAlphabetIndex2(dynamic_cast<const CodonAlphabet*>(alphabet), genCode, distanceOption, "Substitution distances:");
  }
  else
  {
    weights = SequenceApplicationTools::getAlphabetIndex2(alphabet, weightOption, "Substitution weight scheme:");
    distances = SequenceApplicationTools::getAlphabetIndex2(alphabet, distanceOption, "Substitution distances:");
  }
  
  if (regType=="Combination")
  {
    AlphabetIndex2* w2=0;
    AlphabetIndex2* d2=0;
    
    VectorOfSubstitionRegisters* vreg= new VectorOfSubstitionRegisters(stateMap);

    size_t i = 0;
    while (++i)
    {
      string regDesc = ApplicationTools::getStringParameter("reg" + TextTools::toString(i), regArgs, "", "", false, 1);
      if (regDesc == "")
        break;
      
      SubstitutionRegister* sreg=getSubstitutionRegister(regDesc, stateMap, genCode, w2, d2);

      vreg->addRegister(sreg);
    }

    reg = vreg;
  }
  else if (regType == "All")
  {
    reg = new ComprehensiveSubstitutionRegister(stateMap, false);
  }
  else if (regType == "Total")
  {
    reg = new TotalSubstitutionRegister(stateMap);
  }    
  else if (regType == "Selected"){  
    string subsList = ApplicationTools::getStringParameter("substitution.list", regArgs, "All", "", true, false);
    reg = new SelectedSubstitutionRegister(stateMap, subsList);  
  }


  // Alphabet dependent registers

  else if (AlphabetTools::isNucleicAlphabet(alphabet))
  {    
    if (regType == "GC")
      reg = new GCSubstitutionRegister(stateMap, false);
    else if (regType == "TsTv")
      reg = new TsTvSubstitutionRegister(stateMap);
    else if (regType == "SW")
      reg = new SWSubstitutionRegister(stateMap);
    else
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionRegister: unsupported substitution categorization:" + regType + " for alphabet " + alphabet->getAlphabetType());
  }  
  else if (AlphabetTools::isCodonAlphabet(alphabet))
  {    
    if (regType == "IntraAA")
      reg = new AAInteriorSubstitutionRegister(stateMap, *genCode);
    else if (regType == "InterAA")
      reg = new AAExteriorSubstitutionRegister(stateMap, *genCode);
    else if (regType == "GC")
      reg = new GCSynonymousSubstitutionRegister(stateMap, *genCode);
    else if (regType == "TsTv")
      reg = new TsTvSubstitutionRegister(stateMap, *genCode);
    else if (regType == "SW")
      reg = new SWSubstitutionRegister(stateMap, *genCode);
    else if (regType == "KrKc")
      reg = new KrKcSubstitutionRegister(stateMap, *genCode);
    else if (regType == "DnDs")
      reg = new DnDsSubstitutionRegister(stateMap, *genCode, false);
    else
      throw Exception("Unsupported substitution categorization: " + regType + " for alphabet " + alphabet->getAlphabetType());
  }
  
  else if (AlphabetTools::isProteicAlphabet(alphabet))
  {  
    if (regType == "KrKc")
      reg = new KrKcSubstitutionRegister(stateMap);
    else
      throw Exception("Unsupported substitution categorization: " + regType + " for alphabet " + alphabet->getAlphabetType());
  }

  CategorySubstitutionRegister* csr = dynamic_cast<CategorySubstitutionRegister*>(reg);
  if (csr)
    csr->setStationarity(ApplicationTools::getBooleanParameter("stationarity", regArgs, true));

  if (verbose)
    ApplicationTools::displayResult("Substitution Register", regType);

  return reg;
}
