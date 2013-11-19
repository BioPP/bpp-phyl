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
#include "../Model/Protein/Coala.h"
#include "../Model/FrequenciesSet/MvaFrequenciesSet.h"
#include "../Likelihood/TreeLikelihood.h"
#include "../Mapping/LaplaceSubstitutionCount.h"
#include "../Mapping/UniformizationSubstitutionCount.h"
#include "../Mapping/DecompositionSubstitutionCount.h"
#include "../Mapping/NaiveSubstitutionCount.h"
#include "../Mapping/OneJumpSubstitutionCount.h"
#include "../OptimizationTools.h"
#include "../Tree/Tree.h"
#include "../Io/Newick.h"
#include "../Io/NexusIoTree.h"
#include "../Io/Nhx.h"
#include "../Io/BppOSubstitutionModelFormat.h"
#include "../Io/BppOFrequenciesSetFormat.h"
#include "../Io/BppORateDistributionFormat.h"

#include "../NewLikelihood/SinglePhyloLikelihood.h"
#include "../NewLikelihood/MixturePhyloLikelihood.h"
#include "../NewLikelihood/HmmPhyloLikelihood.h"
#include "../NewLikelihood/ParametrizableTree.h"
#include "../NewLikelihood/NonHomogeneousSubstitutionProcess.h"
#include "../NewLikelihood/SimpleSubstitutionProcess.h"
#include "../NewLikelihood/SubstitutionProcessCollection.h"
#include "../NewLikelihood/RateAcrossSitesSubstitutionProcess.h"

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

// From bpp-seq:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/SiteTools.h>

using namespace bpp;

// From the STL:
#include <fstream>
#include <memory>
#include <set>
#include <vector>

using namespace std;

using namespace bpp::newlik;

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
  string format = ApplicationTools::getStringParameter(prefix + "tree.format", params, "Newick", suffix, suffixIsOptional, true);
  string treeFilePath = ApplicationTools::getAFilePath(prefix + "tree.file", params, true, true, suffix, suffixIsOptional);

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

/******************************************************************************/

map<size_t, DiscreteDistribution*> PhylogeneticsApplicationTools::getRateDistributions(
      map<string, string>& params,
      const string& suffix,
      bool suffixIsOptional,
      bool verbose) throw (Exception)
{
  string DistFilePath = ApplicationTools::getAFilePath("rate_distribution.file", params, false, false, suffix, suffixIsOptional);

  map<string, string> paramDist;
  
  if (DistFilePath!="none")
    paramDist=AttributesTools::getAttributesMapFromFile(DistFilePath,"=");

  paramDist.insert(params.begin(), params.end());

  vector<string> vratesName=ApplicationTools::matchingParameters("rate_distribution*", paramDist);

  vector<size_t> vratesNum;
  for (size_t i=0; i< vratesName.size(); i++)
  {
    size_t poseq=vratesName[i].find("=");
    vratesNum.push_back((size_t)TextTools::toInt(vratesName[i].substr(17,poseq-17)));
  }

  BppORateDistributionFormat bIO(true);

  map<size_t, DiscreteDistribution*> mDist;
  
  for (size_t i=0; i<vratesNum.size(); i++)
  {
    string distDescription = ApplicationTools::getStringParameter("rate_distribution"+TextTools::toString(vratesNum[i]), paramDist, "", suffix, suffixIsOptional);

    auto_ptr<DiscreteDistribution> rDist(bIO.read(distDescription, true));
    
    if (verbose)
      {
        ApplicationTools::displayResult("Rate distribution " + TextTools::toString(vratesNum[i]), rDist->getName());
        ApplicationTools::displayResult("Number of classes", TextTools::toString(rDist->getNumberOfCategories()));
      }
    
    mDist[vratesNum[i]]=rDist.release();    
  }
  
  return mDist;
}

std::map<size_t, FrequenciesSet*> PhylogeneticsApplicationTools::getRootFrequenciesSets(
        const Alphabet* alphabet,
        const GeneticCode* gCode,
        const SiteContainer* data,
        std::map<std::string, std::string>& params,
        const std::string& suffix,
        bool suffixIsOptional,
        bool verbose) throw (Exception)
{
  if (dynamic_cast<const CodonAlphabet*>(alphabet) && !gCode)
    throw Exception("PhylogeneticsApplicationTools::getRootFrequenciesSets(): a GeneticCode instance is required for instanciating codon frequencies sets.");

  string RootFilePath = ApplicationTools::getAFilePath("root_freq.file", params, false, false, suffix, suffixIsOptional);
  map<string, string> paramRF;
  
  if (RootFilePath!="none")
    paramRF=AttributesTools::getAttributesMapFromFile(RootFilePath,"=");

  paramRF.insert(params.begin(), params.end());

  
  vector<string> vrfName=ApplicationTools::matchingParameters("root_freq*", paramRF);

  vector<size_t> rfNum;
  for (size_t i=0; i< vrfName.size(); i++)
    {
      size_t poseq=vrfName[i].find("=");
      try {
        rfNum.push_back((size_t)TextTools::toInt(vrfName[i].substr(9,poseq-9)));
      }
      catch (Exception& e) {}
    }

  BppOFrequenciesSetFormat bIO(BppOFrequenciesSetFormat::ALL, verbose);

  map<size_t, FrequenciesSet*> mFS;

  for (size_t i=0; i<rfNum.size(); i++)
  {
    string freqDescription = ApplicationTools::getStringParameter("root_freq"+TextTools::toString(rfNum[i]), paramRF, "", suffix, suffixIsOptional);

    auto_ptr<FrequenciesSet> rFS(bIO.read(alphabet, freqDescription, data, true));

    if (verbose)
    {
      ApplicationTools::displayResult("Root Frequencies Set " + TextTools::toString(rfNum[i]), rFS->getName());
    }
    
    mFS[rfNum[i]]=rFS.release();
  }

  return mFS;
}

/*************************************************************/
/******* MODELS **********************************************/
/*************************************************************/

/******************************************************************************/

map<size_t, SubstitutionModel*> PhylogeneticsApplicationTools::getSubstitutionModels(
     const Alphabet* alphabet,
     const GeneticCode* gCode,
     const SiteContainer* data,
     map<string, string>& params,
     map<string, string>& unparsedParams,
     const string& suffix,
     bool suffixIsOptional,
     bool verbose) throw (Exception)
{
  if (dynamic_cast<const CodonAlphabet*>(alphabet) && !gCode)
    throw Exception("PhylogeneticsApplicationTools::getSubstitutionModels(): a GeneticCode instance is required for instanciating codon models.");

  string ModelFilePath = ApplicationTools::getAFilePath("models.file", params, false, false, suffix, suffixIsOptional);

  map<string, string> paramModel;
  
  if (ModelFilePath!="none")
    paramModel=AttributesTools::getAttributesMapFromFile(ModelFilePath,"=");

  paramModel.insert(params.begin(), params.end());

  vector<string> modelsName=ApplicationTools::matchingParameters("model*", paramModel);

  vector<size_t> modelsNum;
  for (size_t i=0; i< modelsName.size(); i++)
    {
      size_t poseq=modelsName[i].find("=");
      modelsNum.push_back((size_t)TextTools::toInt(modelsName[i].substr(5,poseq-5)));
    }

  map<size_t, SubstitutionModel*> mModel;

  BppOSubstitutionModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, verbose);
  bIO.setGeneticCode(gCode);

  for (size_t i=0; i<modelsNum.size(); i++)
    {
      string modelDescription = ApplicationTools::getStringParameter("model"+TextTools::toString(modelsNum[i]), paramModel, "", suffix, suffixIsOptional);

      auto_ptr<SubstitutionModel> model(bIO.read(alphabet, modelDescription, data, false));
      map<string, string> tmpUnparsedParameterValues(bIO.getUnparsedArguments());
      
      map<string, string>::iterator it;
      for (it=tmpUnparsedParameterValues.begin(); it != tmpUnparsedParameterValues.end(); it++)
        unparsedParams[it->first+"_"+TextTools::toString(modelsNum[i])]=it->second;

      if (verbose)
        {
          ApplicationTools::displayResult("Substitution Model " + TextTools::toString(modelsNum[i]), model->getName());
        }
    
      mModel[modelsNum[i]]=model.release();
    }

  return mModel;
}

/******************************************************************************/

SubstitutionModel* PhylogeneticsApplicationTools::getSubstitutionModel(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const SiteContainer* data,
  std::map<std::string, std::string>& params,
  map<string, string>& unparsedParams,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  BppOSubstitutionModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, verbose);
  string modelDescription;
  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(alphabet);
  if (ca) {
    modelDescription = ApplicationTools::getStringParameter("model", params, "CodonRate(model=JC69)", suffix, suffixIsOptional, verbose);
    if (!gCode)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModel(): a GeneticCode instance is required for instanciating a codon model.");
    bIO.setGeneticCode(gCode);
  } else if (AlphabetTools::isWordAlphabet(alphabet))
    modelDescription = ApplicationTools::getStringParameter("model", params, "Word(model=JC69)", suffix, suffixIsOptional, verbose);
  else
    modelDescription = ApplicationTools::getStringParameter("model", params, "JC69", suffix, suffixIsOptional, verbose);

  SubstitutionModel* model = bIO.read(alphabet, modelDescription, data, true);
  map<string, string> tmpUnparsedParameterValues(bIO.getUnparsedArguments());
      
  map<string, string>::iterator it;
  for (it=tmpUnparsedParameterValues.begin(); it != tmpUnparsedParameterValues.end(); it++)
    unparsedParams[it->first]=it->second;

  return model;
}

// /******************************************************************************/

void PhylogeneticsApplicationTools::setSubstitutionModelParametersInitialValuesWithAliases(
                                                                                           SubstitutionModel& model,
                                                                                           std::map<std::string, std::string>& unparsedParameterValues,
                                                                                           size_t modelNumber,
                                                                                           const SiteContainer* data,
                                                                                           std::map<std::string, double>& existingParams,
                                                                                           std::map<std::string, std::string>& sharedParams,
                                                                                           bool verbose) throw (Exception)
{
  string initFreqs = ApplicationTools::getStringParameter(model.getNamespace() + "initFreqs", unparsedParameterValues, "", "", true, false);

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
      if (value.rfind("_")!=string::npos)
      {
        if (existingParams.find(value) != existingParams.end())
          {
            pl[i].setValue(existingParams[value]);
            sharedParams[pl[i].getName()+"_"+TextTools::toString(modelNumber)]=value;
          }
        else
          throw Exception("Error, unknown parameter " + value);
      }
      else
        pl[i].setValue(TextTools::toDouble(value));
    }

    existingParams[pName+"_"+TextTools::toString(modelNumber)] = pl[i].getValue();
    
    if (verbose)
      ApplicationTools::displayResult("Parameter found", pName + +"_"+TextTools::toString(modelNumber) + "=" + TextTools::toString(pl[i].getValue()));
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
    FrequenciesSet* freq = getFrequenciesSet(alphabet, gCode, freqDescription, data, rateFreqs, verbose);
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
  const std::vector<double>& rateFreqs,
  bool verbose) throw (Exception)
{
  map<string, string> unparsedParameterValues;
  BppOFrequenciesSetFormat bIO(BppOFrequenciesSetFormat::ALL, verbose);
  if (AlphabetTools::isCodonAlphabet(alphabet)) {
    if (!gCode)
      throw Exception("PhylogeneticsApplicationTools::getFrequenciesSet(): a GeneticCode instance is required for instanciating a codon frequencies set.");
    bIO.setGeneticCode(gCode);
  }
  auto_ptr<FrequenciesSet> pFS(bIO.read(alphabet, freqDescription, data, true));
  
  // /////// To be changed for input normalization
  if (rateFreqs.size() > 0)
  {
    pFS.reset(new MarkovModulatedFrequenciesSet(pFS.release(), rateFreqs));
  }

  return pFS.release();
}

/******************************************************/
/**** SUBSTITUTION PROCESS   **************************/
/******************************************************/

/******************************************************************************/


SubstitutionProcess* PhylogeneticsApplicationTools::getSubstitutionProcess(
  const Alphabet* alphabet,
  const GeneticCode* gCode,
  const SiteContainer* data,
  const vector<Tree*> vTree, 
  map<string, string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose)
{
  SubstitutionProcess* SP=0;

  map<string, string> unparsedParams;

  string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", params, "no", "", true, false);
  ApplicationTools::displayResult("Heterogeneous process", nhOpt);

  /////////////////////////
  // Tree
  
  auto_ptr<ParametrizableTree> pTree(new ParametrizableTree(*vTree[0]));

  //////////////////////////
  // Rates
  
  auto_ptr<DiscreteDistribution> rDist(getRateDistribution(params));

  BppOSubstitutionModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, false);
  bIO.setGeneticCode(gCode);


  ///////////////////////////
  /// Models
  
  string tmpDesc;

  if (nhOpt=="no")
  {
    // Homogeneous & stationary models
  
    auto_ptr<SubstitutionModel> tmp(getSubstitutionModel(alphabet, gCode, data, params, unparsedParams));

    if (tmp->getNumberOfStates() >= 2 * tmp->getAlphabet()->getSize() || (rDist->getName()=="Constant"))// first test is for Markov-modulated Markov model!
      SP = new SimpleSubstitutionProcess(tmp.release(), pTree.release(), true);
    else
      SP = new RateAcrossSitesSubstitutionProcess(tmp.release(), rDist.release(), pTree.release());
  }

  // Non-homogeneous models
  else
  {
    string fName=(nhOpt=="one_per_branch"?"model":"model1");
    
    tmpDesc = ApplicationTools::getStringParameter(fName, params, "", suffix, suffixIsOptional, false);
    auto_ptr<SubstitutionModel> tmp(bIO.read(alphabet, tmpDesc, data, true));
    

    
    // ////////////////////////////////////
    // Root frequencies
    
    bool stationarity = ApplicationTools::getBooleanParameter("nonhomogeneous.stationarity", params, false, "", false, false);
    
    auto_ptr<FrequenciesSet> rootFrequencies(0);
  
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
    
      string freqDescription = ApplicationTools::getStringParameter("nonhomogeneous.root_freq", params, "", suffix, suffixIsOptional);
      if (freqDescription.substr(0, 10) == "MVAprotein")
      {
        if (dynamic_cast<Coala*>(tmp.get()))
          dynamic_cast<MvaFrequenciesSet*>(rootFrequencies.get())->initSet(dynamic_cast<CoalaCore*>(tmp.get()));
        else
          throw Exception("The MVAprotein frequencies set at the root can only be used if a Coala model is used on branches.");
      }
      else
        rootFrequencies.reset(getRootFrequenciesSet(alphabet, gCode, data, params, rateFreqs, suffix, suffixIsOptional, verbose));
    
      stationarity = !rootFrequencies.get();
    }

    ApplicationTools::displayBooleanResult("Stationarity assumed", stationarity);

    ///////////////////////////////////////
    // One_per_branch

    if (nhOpt=="one_per_branch"){
      vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", params, ',', "");

      for (unsigned int i = 0; i < globalParameters.size(); i++)
        ApplicationTools::displayResult("Global parameter", globalParameters[i]);

      SP = NonHomogeneousSubstitutionProcess::createNonHomogeneousSubstitutionProcess(
                            tmp.release(),
                            rDist.release(),
                            rootFrequencies.release(),
                            pTree.release(),
                            globalParameters);
      
    }
    else
    {

      // //////////////////////////////
      // General
    
      size_t nbModels = ApplicationTools::getParameter<size_t>("nonhomogeneous.number_of_models", params, 1, suffix, suffixIsOptional, false);

      if (nbModels == 0)
        throw Exception("The number of models can't be 0 !");

      if (verbose)
        ApplicationTools::displayResult("Number of distinct models", TextTools::toString(nbModels));


      // //////////////////////////////////////
      // Now parse all models:

      bIO.setVerbose(true);
      
      SP = new NonHomogeneousSubstitutionProcess(rDist.release(), pTree.release(),rootFrequencies.release());

      NonHomogeneousSubstitutionProcess* nhSP=dynamic_cast<NonHomogeneousSubstitutionProcess*>(SP);
      
      for (size_t i = 0; i < nbModels; i++)
      {
        string prefix = "model" + TextTools::toString(i + 1);
        string modelDesc;
        modelDesc = ApplicationTools::getStringParameter(prefix, params, "", suffix, suffixIsOptional, verbose);
        
        auto_ptr<SubstitutionModel> model(bIO.read(alphabet, modelDesc, data, true));
        map<string, string> tmpUnparsedParameterValues(bIO.getUnparsedArguments());

        map<string, string>::iterator it;
        for (it=tmpUnparsedParameterValues.begin(); it != tmpUnparsedParameterValues.end(); it++)
          unparsedParams[it->first+"_"+TextTools::toString(i+1)]=it->second;

        vector<int> nodesId = ApplicationTools::getVectorParameter<int>(prefix + ".nodes_id", params, ',', ':', TextTools::toString(i), suffix, suffixIsOptional, true);
        
        if (verbose)
          ApplicationTools::displayResult("Model" + TextTools::toString(i + 1) + " is associated to", TextTools::toString(nodesId.size()) + " node(s).");

        nhSP->addModel(model.release(), nodesId);
      }
      
      nhSP->isFullySetUp();

    }
  }


  //////// Aliasing
  // Finally check parameter aliasing:

  string aliasDesc = ApplicationTools::getStringParameter("nonhomogeneous.alias", params, "", suffix, suffixIsOptional, false);

  StringTokenizer st(aliasDesc, ",");
  while (st.hasMoreToken())
  {
    string alias = st.nextToken();
    string::size_type index = alias.find("->");
    if (index == string::npos)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionProcess. Bad alias syntax, should contain `->' symbol: " + alias);
    string p1 = alias.substr(0, index);
    string p2 = alias.substr(index + 2);
    unparsedParams[p1]=p2;
  }

  SP->aliasParameters(unparsedParams, verbose);
  
  return SP;
}


/******************************************************************************/

void PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember(
          SubstitutionProcessCollection* SubProColl, 
          map<string, string>& params,
          const string& suffix,
          bool verbose)
{
  string procName = "";
  map<string, string> args;

  string procDesc = ApplicationTools::getStringParameter("process", params, "Hxomogeneous", suffix, false);

  KeyvalTools::parseProcedure(procDesc, procName, args);


  /////
  // tree number

  if (args.find("tree")==args.end())
    throw Exception("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember. A tree number is compulsory.");

  size_t numTree=(size_t)ApplicationTools::getIntParameter("tree", args, 1, "", true, verbose);

  if (! SubProColl->hasTreeNumber(numTree))
    throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : unknown tree number", (int)numTree);

  ///////
  // rate number
      
  if (args.find("rate")==args.end())
    throw Exception("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember. A rate number is compulsory.");

  size_t numRate=(size_t)ApplicationTools::getIntParameter("rate", args, 1, "", true, verbose);

  if (! SubProColl->hasDistributionNumber(numRate))
    throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : unknown rate number", (int)numRate);


  //////////
  // root freq number

  bool stationarity=(args.find("root_freq")==args.end());
  size_t numFreq=0;
  
  if (stationarity)
    ApplicationTools::displayMessage("Stationarity assumed.");
  else
  {
    numFreq=(size_t)ApplicationTools::getIntParameter("root_freq", args, 1, "", true, verbose);
    if (! SubProColl->hasFrequenciesNumber(numFreq))
      throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : unknown root frequencies number", (int)numFreq);
  }

  //////////////////
  /// models

  if (procName=="Homogeneous")
  {
    if (args.find("model")==args.end())
      throw Exception("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember. A model number is compulsory.");

    size_t numModel=(size_t)ApplicationTools::getIntParameter("model", args, 1, "", true, verbose);

    if (! SubProColl->hasModelNumber(numModel))
      throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : unknown model number", (int)numModel);
    
    size_t nNodes=SubProColl->getTree(numTree)->getNumberOfBranches();

    vector<int> vNodes;
    for (int i=0; i<(int)nNodes; i++)
      vNodes.push_back(i);

    map<size_t, vector<int> > mModBr;
    mModBr[numModel]=vNodes;

    if (verbose){
      ApplicationTools::displayMessage("Homogeneous process : ");
      ApplicationTools::displayResult (" Model number",TextTools::toString(numModel));
      ApplicationTools::displayResult (" Tree number",TextTools::toString(numTree));
      ApplicationTools::displayResult (" Rate number",TextTools::toString(numRate));
      if (!stationarity)
        ApplicationTools::displayResult (" Root frequencies number",TextTools::toString(numFreq));
    }

    if (stationarity)
      SubProColl->addSubstitutionProcess(mModBr, numTree, numRate);
    else
      SubProColl->addSubstitutionProcess(mModBr, numTree, numRate, numFreq);
  }
  
  else if ((procName=="Nonhomogeneous") ||  (procName=="Nonhomogeneous"))
  {
    size_t indModel=1;
    map<size_t, vector<int> > mModBr;

    while (args.find("model"+TextTools::toString(indModel))!=args.end())
    {
      size_t numModel=(size_t)ApplicationTools::getIntParameter("model"+TextTools::toString(indModel), args, 1, "", true, verbose);

      if (mModBr.find(numModel)!=mModBr.end())
        throw BadIntegerException("PhylogeneticsApplicationTools::addSubstitutionProcessCollectionMember : model number seen twice.", (int)numModel);

      vector<int> nodesId = ApplicationTools::getVectorParameter<int>("model"+TextTools::toString(indModel)  + ".nodes_id", args, ',', ':', "0");
      
      mModBr[numModel]=nodesId;
      
      indModel++;
    }

    if (verbose){
      ApplicationTools::displayMessage("Nonhomogeneous process : ");
      map<size_t, vector<int> >::const_iterator it;
      for (it=mModBr.begin(); it!=mModBr.end(); it++)
        ApplicationTools::displayResult (" Model number" + TextTools::toString(it->first) + " associated to", TextTools::toString(it->second.size()) + " node(s).");
      ApplicationTools::displayResult (" Tree number",TextTools::toString(numTree));
      ApplicationTools::displayResult (" Rate number",TextTools::toString(numRate));
      if (!stationarity)
        ApplicationTools::displayResult (" Root frequencies number",TextTools::toString(numFreq));
    }
    
    if (stationarity)
      SubProColl->addSubstitutionProcess(mModBr, numTree, numRate);
    else
      SubProColl->addSubstitutionProcess(mModBr, numTree, numRate, numFreq);
  }
    
}


/******************************************************************************/


SubstitutionProcessCollection* PhylogeneticsApplicationTools::getSubstitutionProcessCollection(
       const Alphabet* alphabet,
       const GeneticCode* gCode,
       const SiteContainer* data,
       const vector<Tree*> vTree, 
       map<string, string>& params,
       const string& suffix,
       bool suffixIsOptional,
       bool verbose)
{
  SubstitutionProcessCollection*  SPC=new SubstitutionProcessCollection();

  map<string, double> existingParameters;
 
  /////////////////////////
  // Trees
  
  vector<ParametrizableTree*> vpTree;
  
  for (size_t i=0; i<vTree.size(); i++)
    vpTree.push_back(new ParametrizableTree(*(vTree[i])));

  if (vpTree.size()==0)
    throw Exception("Missing tree in construction of SubstitutionProcessCollection.");
  
  for (size_t i=1;i<=vpTree.size(); i++)
    SPC->addTree(vpTree[i-1], i);

  /////////////////////////
  // Rates
  
  map<size_t, DiscreteDistribution*> mDist=getRateDistributions(params);

  if (mDist.size()==0)
    throw Exception("Missing rate distribution in construction of SubstitutionProcessCollection.");
  
  map<size_t, DiscreteDistribution*>::iterator itd;
  
  for (itd=mDist.begin();itd!=mDist.end(); itd++)
    SPC->addDistribution(itd->second, itd->first);

  //////////////////////////
  // Models

  map<string, string> unparsedParams;

  map<size_t, SubstitutionModel*> mModel=getSubstitutionModels(alphabet, gCode, data, params, unparsedParams, suffix, suffixIsOptional);

  if (mModel.size()==0)
    throw Exception("Missing model in construction of SubstitutionProcessCollection.");
  
  map<size_t, SubstitutionModel*>::iterator itm;
  
  for (itm=mModel.begin();itm!=mModel.end(); itm++)
    SPC->addModel(itm->second, itm->first);


  /////////////////////////////
  // Root Frequencies

  map<size_t, FrequenciesSet*> mFreq=getRootFrequenciesSets(alphabet, gCode, data, params, suffix, suffixIsOptional);

  map<size_t, FrequenciesSet*>::iterator itr;
  
  for (itr=mFreq.begin();itr!=mFreq.end(); itr++)
    SPC->addFrequencies(itr->second, itr->first);

  
  ////////////////////////////////
  // Now processes

  size_t indProc=1;

  while (params.find("process"+TextTools::toString(indProc))!=params.end())
  {
    addSubstitutionProcessCollectionMember(SPC, params, TextTools::toString(indProc));

    indProc++;
  }

  if (indProc==1)
    throw Exception("Missing process in construction of SubstitutionProcessCollection.");


  ///////////////////////////
  // Now set shared parameters:

  //////// Aliasing
  // Finally check parameter aliasing:

  string aliasDesc = ApplicationTools::getStringParameter("collection.alias", params, "", suffix, suffixIsOptional, false);

  StringTokenizer st(aliasDesc, ",");
  while (st.hasMoreToken())
  {
    string alias = st.nextToken();
    string::size_type index = alias.find("->");
    if (index == string::npos)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionProcessCollection. Bad alias syntax, should contain `->' symbol: " + alias);
    string p1 = alias.substr(0, index);
    string p2 = alias.substr(index + 2);
    unparsedParams[p1]=p2;
  }

  SPC->aliasParameters(unparsedParams, verbose);
  
  return SPC;
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
                                                                             bool verbose)
{
  if (!ApplicationTools::parameterExists("nonhomogeneous.number_of_models", params))
    throw Exception("You must specify this parameter: nonhomogeneous.number_of_models .");
  size_t nbModels = ApplicationTools::getParameter<size_t>("nonhomogeneous.number_of_models", params, 1, suffix, suffixIsOptional, false);
  if (nbModels == 0)
    throw Exception("The number of models can't be 0 !");

  bool nomix = true;
  for (size_t i = 0; nomix &(i < nbModels); i++)
    {
      string prefix = "model" + TextTools::toString(i + 1);
      string modelDesc;
      modelDesc = ApplicationTools::getStringParameter(prefix, params, "", suffix, suffixIsOptional, verbose);

      if (modelDesc.find("Mixed") != string::npos)
        nomix = false;
    }

  SubstitutionModelSet* modelSet, * modelSet1 = 0;
  modelSet1 = new SubstitutionModelSet(alphabet);
  setSubstitutionModelSet(*modelSet1, alphabet, gCode, data, params, suffix, suffixIsOptional, verbose);

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
  const GeneticCode* gCode,
  const SiteContainer* data,
  map<string, string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose)
{
  modelSet.clear();
  if (!ApplicationTools::parameterExists("nonhomogeneous.number_of_models", params))
    throw Exception("You must specify this parameter: nonhomogeneous.number_of_models .");
  size_t nbModels = ApplicationTools::getParameter<size_t>("nonhomogeneous.number_of_models", params, 1, suffix, suffixIsOptional, false);
  if (nbModels == 0)
    throw Exception("The number of models can't be 0 !");

  if (verbose)
    ApplicationTools::displayResult("Number of distinct models", TextTools::toString(nbModels));

  BppOSubstitutionModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, false);

  // ///////////////////////////////////////////
  // Build a new model set object:

  vector<double> rateFreqs;
  string tmpDesc;
  if (AlphabetTools::isCodonAlphabet(alphabet)) {
    if (!gCode)
      throw Exception("PhylogeneticsApplicationTools::setSubstitutionModelSet(): a GeneticCode instance is required for instanciating a codon model.");
    bIO.setGeneticCode(gCode);
    tmpDesc = ApplicationTools::getStringParameter("model1", params, "CodonRate(model=JC69)", suffix, suffixIsOptional, false);
  } else if (AlphabetTools::isWordAlphabet(alphabet))
    tmpDesc = ApplicationTools::getStringParameter("model1", params, "Word(model=JC69)", suffix, suffixIsOptional, false);
  else
    tmpDesc = ApplicationTools::getStringParameter("model1", params, "JC69", suffix, suffixIsOptional, false);

  auto_ptr<SubstitutionModel> tmp(bIO.read(alphabet, tmpDesc, data, true));
  //  map<string, string> tmpUnparsedParameterValues(bIO.getUnparsedArguments());

  if (tmp->getNumberOfStates() != alphabet->getSize())
    {
      // Markov-Modulated Markov Model...
      size_t n = static_cast<size_t>(tmp->getNumberOfStates() / alphabet->getSize());
      rateFreqs = vector<double>(n, 1. / static_cast<double>(n)); // Equal rates assumed for now, may be changed later (actually, in the most general case,
    }

  // ////////////////////////////////////
  // Deal with root frequencies

  bool stationarity = ApplicationTools::getBooleanParameter("nonhomogeneous.stationarity", params, false, "", false, false);
  FrequenciesSet* rootFrequencies = 0;
  if (!stationarity)
  {
    rootFrequencies = getRootFrequenciesSet(alphabet, gCode, data, params, rateFreqs, suffix, suffixIsOptional, verbose);
    stationarity = !rootFrequencies;
    string freqDescription = ApplicationTools::getStringParameter("nonhomogeneous.root_freq", params, "", suffix, suffixIsOptional);
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
  
  map<string, double> existingParameters;

  for (size_t i = 0; i < nbModels; i++)
    {
      string prefix = "model" + TextTools::toString(i + 1);
      string modelDesc;
      if (AlphabetTools::isCodonAlphabet(alphabet))
        modelDesc = ApplicationTools::getStringParameter(prefix, params, "CodonRate(model=JC69)", suffix, suffixIsOptional, verbose);
      else if (AlphabetTools::isWordAlphabet(alphabet))
        modelDesc = ApplicationTools::getStringParameter(prefix, params, "Word(model=JC69)", suffix, suffixIsOptional, verbose);
      else
        modelDesc = ApplicationTools::getStringParameter(prefix, params, "JC69", suffix, suffixIsOptional, verbose);

      auto_ptr<SubstitutionModel> model(bIO.read(alphabet, modelDesc, data, false));
      map<string, string> unparsedParameterValues(bIO.getUnparsedArguments());

      map<string, string> sharedParameters;
      setSubstitutionModelParametersInitialValuesWithAliases(
                                                             *model,
                                                             unparsedParameterValues, i+1, data,
                                                             existingParameters, sharedParameters,
                                                             verbose);

      vector<int> nodesId = ApplicationTools::getVectorParameter<int>(prefix + ".nodes_id", params, ',', ':', TextTools::toString(i), suffix, suffixIsOptional, true);

      if (verbose)
        ApplicationTools::displayResult("Model" + TextTools::toString(i + 1) + " is associated to", TextTools::toString(nodesId.size()) + " node(s).");

      modelSet.addModel(model.get(), nodesId);

      // Now set shared parameters:
      map<string, string>::const_iterator it;
      for (it=sharedParameters.begin(); it!=sharedParameters.end(); it++)
        modelSet.aliasParameters(it->second, it->first);
    
      model.release();
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

  size_t numd;
  if (!ApplicationTools::parameterExists("site.number_of_paths", params))
    numd = 0;
  else
    numd = ApplicationTools::getParameter<size_t>("site.number_of_paths", params, 1, suffix, suffixIsOptional, false);

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
      classes.push_back(TextTools::toInt(strtok2.nextToken()));

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
  bool verbose) throw (Exception)
{
  string distDescription = ApplicationTools::getStringParameter("rate_distribution", params, "Constant()", suffix, suffixIsOptional);

  string distName;
  map<string, string> args;
  KeyvalTools::parseProcedure(distDescription, distName, args);

  BppORateDistributionFormat bIO(true);
  auto_ptr<DiscreteDistribution> rDist(bIO.read(distDescription, true));

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

bpp::TreeLikelihood* PhylogeneticsApplicationTools::optimizeParameters(
  bpp::TreeLikelihood* tl,
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
    new StlOutputStream(new ofstream(mhPath.c_str(), ios::out));
  if (verbose)
    ApplicationTools::displayResult("Message handler", mhPath);

  string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional);
  OutputStream* profiler =
    (prPath == "none") ? 0 :
    (prPath == "std") ? ApplicationTools::message :
    new StlOutputStream(new ofstream(prPath.c_str(), ios::out));
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

  // Should I ignore some parameters?
  ParameterList parametersToEstimate = parameters;
  string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameter", params, "", suffix, suffixIsOptional, false);
  if (paramListDesc.length() == 0)
    paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameters", params, "", suffix, suffixIsOptional, false);
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
          if (nhtl!=NULL){
            vector<string> vs2 = nhtl->getRootFrequenciesParameters().getParameterNames();
            VectorTools::diff(vs1,vs2,vs);
            }
          else
            vs=vs1;

          parametersToEstimate.deleteParameters(vs);
          if (verbose)
            ApplicationTools::displayResult("Parameter ignored", string("Model"));          
        }
      else if (param.find("*") != string::npos)
      {
        vector<string> vs;
        for (size_t j = 0; j < parametersToEstimate.size(); j++)
        {
          StringTokenizer stj(param, "*", true, false);
          size_t pos1, pos2;
          string parn = parametersToEstimate[j].getName();
          bool flag(true);
          string g = stj.nextToken();
          pos1 = parn.find(g);
          if (pos1 != 0)
            flag = false;
          pos1 += g.length();
          while (flag && stj.hasMoreToken())
          {
            g = stj.nextToken();
            pos2 = parn.find(g, pos1);
            if (pos2 == string::npos)
            {
              flag = false;
              break;
            }
            pos1 = pos2 + g.length();
          }
          if (flag &&
              ((g.length() == 0) || (pos1 == parn.length()) || (parn.rfind(g) == parn.length() - g.length())))
            vs.push_back(parn);
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


PhyloLikelihood* PhylogeneticsApplicationTools::optimizeParameters(
    PhyloLikelihood* lik,
    const ParameterList& parameters,
    std::map<std::string, std::string>& params,
    const std::string& suffix,
    bool suffixIsOptional,
    bool verbose)
  throw (Exception)
{
  string optimization = ApplicationTools::getStringParameter("optimization", params, "FullD(derivatives=Newton)", suffix, suffixIsOptional, false);
  if (optimization == "None")
    return lik;
  string optName;
  map<string, string> optArgs;
  KeyvalTools::parseProcedure(optimization, optName, optArgs);

  unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2, suffix, suffixIsOptional);

  string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional);
  OutputStream* messageHandler =
    (mhPath == "none") ? 0 :
    (mhPath == "std") ? ApplicationTools::message :
    new StlOutputStream(new ofstream(mhPath.c_str(), ios::out));
  if (verbose)
    ApplicationTools::displayResult("Message handler", mhPath);

  string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional);
  OutputStream* profiler =
    (prPath == "none") ? 0 :
    (prPath == "std") ? ApplicationTools::message :
    new StlOutputStream(new ofstream(prPath.c_str(), ios::out));
  if (profiler)
    profiler->setPrecision(20);
  if (verbose)
    ApplicationTools::displayResult("Profiler", prPath);

  bool scaleFirst = ApplicationTools::getBooleanParameter("optimization.scale_first", params, false, suffix, suffixIsOptional, false);
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
  string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameter", params, "", suffix, suffixIsOptional, false);
  if (paramListDesc.length() == 0)
    paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameters", params, "", suffix, suffixIsOptional, false);
  StringTokenizer st(paramListDesc, ",");
  while (st.hasMoreToken())
    {
      try
        {
          string param = st.nextToken();
          if (param == "BrLen")
            {
              vector<string> vs = lik->getBranchLengthsParameters().getParameterNames();
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
              vector<string> vs= lik->getSubstitutionModelParameters().getParameterNames();
              parametersToEstimate.deleteParameters(vs);
              if (verbose)
                ApplicationTools::displayResult("Parameter ignored", string("Model"));          
            }
          else if (param.find("*") != string::npos)
            {
              vector<string> vs;
              for (size_t j = 0; j < parametersToEstimate.size(); j++)
                {
                  StringTokenizer stj(param, "*", true, false);
                  size_t pos1, pos2;
                  string parn = parametersToEstimate[j].getName();
                  bool flag(true);
                  string g = stj.nextToken();
                  pos1 = parn.find(g);
                  if (pos1 != 0)
                    flag = false;
                  pos1 += g.length();
                  while (flag && stj.hasMoreToken())
                    {
                      g = stj.nextToken();
                      pos2 = parn.find(g, pos1);
                      if (pos2 == string::npos)
                        {
                          flag = false;
                          break;
                        }
                      pos1 = pos2 + g.length();
                    }
                  if (flag &&
                      ((g.length() == 0) || (pos1 == parn.length()) || (parn.rfind(g) == parn.length() - g.length())))
                    vs.push_back(parn);
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
                  size_t p = pl.whichParameterHasName(pname);
                  pl.setParameter(p, AutoParameter(pl[p]));
                  pl[p].setValue(TextTools::toDouble(pvalue));
                }
            }
          bck.close();
          lik->setParameters(pl);
          if (abs(lik->getValue() - fval) > 0.000001)
            throw Exception("Incorrect likelihood value after restoring, from backup file. Remove backup file and start from scratch :s");
          ApplicationTools::displayResult("Restoring log-likelihood", -fval);
        }
    }

  // There it goes...
  bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", params, false, suffix, suffixIsOptional, false);
  if (optimizeTopo)
    throw Exception("Topology opmitization not implemented yet for processes");
  
  // if (verbose)
  //   ApplicationTools::displayResult("Optimize topology", optimizeTopo ? "yes" : "no");
  // string nniMethod = ApplicationTools::getStringParameter("optimization.topology.algorithm_nni.method", params, "phyml", suffix, suffixIsOptional, false);
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

      // if (optimizeTopo)
      //   {
      //     bool optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, false);
      //     unsigned int topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, "", true, false);
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

      // if (optimizeTopo)
      //   {
      //     bool optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, false);
      //     unsigned int topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, "", true, false);
      //     double tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional);
      //     double tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, suffixIsOptional);
      //     tl = OptimizationTools::optimizeTreeNNI2(
      //                                              dynamic_cast<NNIHomogeneousTreeLikelihood*>(tl), parametersToEstimate,
      //                                              optNumFirst, tolBefore, tolDuring, nbEvalMax, topoNbStep, messageHandler, profiler,
      //                                              reparam, optVerbose, optMethodDeriv, nniAlgo);
      //   }

      parametersToEstimate.matchParametersValues(lik->getParameters());
      n = OptimizationTools::optimizeNumericalParameters2(
                                                          lik, parametersToEstimate,
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
      remove(backupFile.c_str());
    }
  return lik;
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
    new StlOutputStream(new ofstream(mhPath.c_str(), ios::out));
  if (verbose)
    ApplicationTools::displayResult("Message handler", mhPath);

  string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional);
  OutputStream* profiler =
    (prPath == "none") ? 0 :
    (prPath == "std") ? ApplicationTools::message :
    new StlOutputStream(new ofstream(prPath.c_str(), ios::out));
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
  for (size_t i = 0; i < pl.size(); ++i)
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
  const vector<const Tree*>& trees,
  map<string, string>& params,
  const string& prefix,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose,
  bool checkOnly) throw (Exception)
{
  string format = ApplicationTools::getStringParameter(prefix + "tree.format", params, "Newick", suffix, suffixIsOptional, false);
  string file = ApplicationTools::getAFilePath(prefix + "tree.file", params, true, false, suffix, suffixIsOptional);
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

void PhylogeneticsApplicationTools::writeTrees(
   const vector<const TreeTemplate<Node>*>& trees,
   map<string, string>& params,
   const string& prefix,
   const string& suffix,
   bool suffixIsOptional,
   bool verbose,
   bool checkOnly) throw (Exception)
{
  string format = ApplicationTools::getStringParameter(prefix + "tree.format", params, "Newick", suffix, suffixIsOptional, false);
  string file = ApplicationTools::getAFilePath(prefix + "tree.file", params, true, false, suffix, suffixIsOptional);
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
    vector<const Tree*> vT;
    for (size_t i=0; i< trees.size(); i++){
      vT.push_back(dynamic_cast<const Tree*>(trees[i]));
    }
    treeWriter->write(vT, file, true);
  }

  delete treeWriter;
  if (verbose)
    ApplicationTools::displayResult("Wrote trees to file ", file);
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printParameters(const SubstitutionModel* model, OutputStream& out)
{
  out << "model=";
  map<string, string> globalAliases;
  vector<string> writtenNames;
  BppOSubstitutionModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, false);
  bIO.write(*model, out, globalAliases, writtenNames);
  out.endLine();
}

void PhylogeneticsApplicationTools::printParameters(const SubstitutionProcess* process, OutputStream& out)
{
  if (dynamic_cast<const SimpleSubstitutionProcess*>(process)!=NULL)
  {
    (out << "nonhomogeneous=no").endLine();
    
    out << "model=";
    map<string, string> globalAliases;
    vector<string> writtenNames;
    BppOSubstitutionModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, false);
    bIO.write(process->getSubstitutionModel(0,0), out, globalAliases, writtenNames);
    out.endLine();
  }

  else if (dynamic_cast<const RateAcrossSitesSubstitutionProcess*>(process)!=NULL)
  {
    const RateAcrossSitesSubstitutionProcess* pRA = dynamic_cast<const RateAcrossSitesSubstitutionProcess*>(process);
      
    (out << "nonhomogeneous=no").endLine();
    
    out << "model=";
    map<string, string> globalAliases;
    vector<string> writtenNames;
    BppOSubstitutionModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, false);
    bIO.write(process->getSubstitutionModel(0,0), out, globalAliases, writtenNames);
    out.endLine();
    out.endLine();
    
    // Rate distribution
    
    out << "rate_distribution=";
    const BppORateDistributionFormat* bIOR = new BppORateDistributionFormat(true);
    bIOR->write(pRA->getRateDistribution(), out, globalAliases, writtenNames);
    delete bIOR;
    out.endLine();
  }

  else if (dynamic_cast<const NonHomogeneousSubstitutionProcess*>(process)!=NULL)
  {
    const NonHomogeneousSubstitutionProcess* pNH = dynamic_cast<const NonHomogeneousSubstitutionProcess*>(process);

    (out << "nonhomogeneous=general").endLine();
    (out << "nonhomogeneous.number_of_models=" << pNH->getNumberOfModels()).endLine();

    vector<string> writtenNames;

    // Loop over all models:
    for (size_t i = 0; i < pNH->getNumberOfModels(); i++)
    {
      const SubstitutionModel* model = pNH->getModel(i);

      // First get the aliases for this model:
      map<string, string> aliases;
      
      ParameterList pl=model->getParameters();
        
      for (size_t np = 0 ; np< pl.size() ; np++)
      {
        string nfrom=pNH->getFrom(pl[np].getName()+"_"+TextTools::toString(i+1));
        if (nfrom!="")
          aliases[pl[np].getName()]=nfrom;
      }

      // Now print it:
      writtenNames.clear();
      out.endLine() << "model" << (i + 1) << "=";
      BppOSubstitutionModelFormat bIOsm(BppOSubstitutionModelFormat::ALL, true, true, true, false);
      map<string, string>::iterator it;
      bIOsm.write(*model, out, aliases, writtenNames);
      out.endLine();
      vector<int> ids = pNH->getNodesWithModel(i);
      out << "model" << (i + 1) << ".nodes_id=" << ids[0];
      for (size_t j = 1; j < ids.size(); ++j)
      {
        out << "," << ids[j];
      }
      out.endLine();
    }

    // Root frequencies:
    out.endLine();
    if (pNH->getRootFrequenciesSet())
    {
      out << "nonhomogeneous.root_freq=";

      BppOFrequenciesSetFormat bIO(BppOFrequenciesSetFormat::ALL, false);
      bIO.write(pNH->getRootFrequenciesSet(), out, writtenNames);
    }
    else
      out << "nonhomogeneous.stationarity=true";
    out.endLine();

    // Rate distribution
    
    map<string, string> aliases;
    const DiscreteDistribution* pdd=&pNH->getRateDistribution();
    
    ParameterList pl=pdd->getParameters();
    for (size_t np = 0 ; np< pl.size() ; np++)
    {
      string nfrom=pNH->getFrom(pl[np].getName());
      if (nfrom!="")
        aliases[pl[np].getName()]=nfrom;
    }
    out.endLine();
    out << "rate_distribution=";
    const BppORateDistributionFormat* bIO = new BppORateDistributionFormat(true);
    bIO->write(*pdd, out, aliases, writtenNames);
    delete bIO;
    out.endLine();
  }
}

void PhylogeneticsApplicationTools::printParameters(const SubstitutionProcessCollection* collection, OutputStream& out)
{
  vector<string> writtenNames;

  // The models
  vector<size_t> modN=collection->getModelNumbers();

  for (size_t i = 0; i < modN.size(); i++)
  {
    const SubstitutionModel* model =collection->getModel(modN[i]);

    // First get the aliases for this model:
    map<string, string> aliases;
    
    ParameterList pl=model->getParameters();
    
    for (size_t np = 0 ; np< pl.size() ; np++)
    {
      string nfrom=collection->getFrom(pl[np].getName()+"_"+TextTools::toString(modN[i]));
      if (nfrom!="")
        aliases[pl[np].getName()]=nfrom;
    }

    // Now print it:
    writtenNames.clear();
    out.endLine() << "model" << modN[i] << "=";
    BppOSubstitutionModelFormat bIOsm(BppOSubstitutionModelFormat::ALL, true, true, true, false);
    map<string, string>::iterator it;
    bIOsm.write(*model, out, aliases, writtenNames);
    out.endLine();
    // vector<int> ids = pNH->getNodesWithModel(i);
    // out << "model" << (i + 1) << ".nodes_id=" << ids[0];
    // for (size_t j = 1; j < ids.size(); ++j)
    //   {
    //     out << "," << ids[j];
    //   }
    // out.endLine();
  }

  // Root frequencies:
  vector<size_t> rootFreqN=collection->getFrequenciesNumbers();

  for (size_t i = 0; i < rootFreqN.size(); i++)
    {
      const FrequenciesSet* rootFreq =collection->getFrequencies(rootFreqN[i]);

      // Now print it:
      writtenNames.clear();
      out.endLine() << "root_freq" << rootFreqN[i] << "=";
      BppOFrequenciesSetFormat bIOf(BppOFrequenciesSetFormat::ALL, true);
      map<string, string>::iterator it;
      bIOf.write(rootFreq, out, writtenNames);
      out.endLine();
      // vector<int> ids = pNH->getNodesWithModel(i);
      // out << "model" << (i + 1) << ".nodes_id=" << ids[0];
      // for (size_t j = 1; j < ids.size(); ++j)
      //   {
      //     out << "," << ids[j];
      //   }
      // out.endLine();
    }

  // Rate distribution

  out.endLine();
  vector<size_t> distN=collection->getDistributionNumbers();

  for (size_t i = 0; i < distN.size(); i++)
  {
    const DiscreteDistribution* dist =collection->getDistribution(distN[i]);

    // First get the aliases for this model:
    map<string, string> aliases;
    
    ParameterList pl=dist->getParameters();
    
    for (size_t np = 0 ; np< pl.size() ; np++)
    {
      string nfrom=collection->getFrom(pl[np].getName()+"_"+TextTools::toString(distN[i]));
      if (nfrom!="")
        aliases[pl[np].getName()]=nfrom;
    }

    // Now print it:
    writtenNames.clear();
    out.endLine() << "rate_distribution" << modN[i] << "=";
    BppORateDistributionFormat bIOd(true);
    map<string, string>::iterator it;
    bIOd.write(*dist, out, aliases, writtenNames);
    out.endLine();
  }

  // processes
  out.endLine();

  size_t procN=collection->getNumberOfSubstitutionProcess();

  for (size_t i=0;i<procN;i++)
  {
    const SubstitutionProcessCollectionMember* spcm=dynamic_cast<const SubstitutionProcessCollectionMember*>(collection->getSubstitutionProcess(i));
    
    out << "process" << i+1 << "=";
      
    if (spcm->getNumberOfModels()==1)
      out << "Homogeneous(model=" << spcm->getModelNumbers()[0];
    else
    {
      out << "Nonhomogeneous(";
      vector<size_t> vMN=spcm->getModelNumbers();
      for (size_t j=0;j<vMN.size();j++)
      {
        out << "model" << (j+1) << "=" << vMN[j];
        vector<int> ids = spcm->getNodesWithModel(vMN[j]);
        out << "model" << (j+1) << ".nodes_id=(" << ids[0];
        for (size_t k = 1; k < ids.size(); ++k)
        {
          out << "," << ids[k];
        }
        out << ")";
      }
    }

    out << ", tree=" << spcm->getTreeNumber();
    out << ", rate=" << spcm->getDistributionNumber();
    if (spcm->getRootFrequenciesSet())
      out << ", root_freq=" << spcm->getRootFrequenciesNumber();
    out << ")";
    out.endLine();
  }
  
}

void PhylogeneticsApplicationTools::printParameters(const PhyloLikelihood* phylolike, OutputStream& out)
{
  out << "# Log likelihood = ";
  out.setPrecision(20) << (-phylolike->getValue());
  out.endLine();
  out << "# Number of sites = ";
  out.setPrecision(20) << phylolike->getNumberOfSites();
  out.endLine();
  out.endLine();
  out << "# Substitution model parameters:";
  out.endLine();

  if (dynamic_cast<const MultiPhyloLikelihood*>(phylolike)!=NULL)
  {
    if (dynamic_cast<const MixturePhyloLikelihood*>(phylolike)!=NULL)
    {
      const MixturePhyloLikelihood* pM=dynamic_cast<const MixturePhyloLikelihood*>(phylolike);
    
      PhylogeneticsApplicationTools::printParameters(pM->getCollection(), out);
    
      out.endLine();
      out << "collection=Mixture(probas=(" << pM->getSubProcessProb(0);
    
      for (size_t i=1; i< pM->getCollection()->getNumberOfSubstitutionProcess(); i++)
        out << "," << pM->getSubProcessProb(i);
    
      out << "))";
      out.endLine();
    }
    else if (dynamic_cast<const HmmPhyloLikelihood*>(phylolike)!=NULL)
    {
      const HmmPhyloLikelihood* pM=dynamic_cast<const HmmPhyloLikelihood*>(phylolike);
      
      PhylogeneticsApplicationTools::printParameters(pM->getCollection(), out);
      out.endLine();
      out << "collection=HMM(probas=";

      const Matrix<double>& tMt = pM->getHmmTransitionMatrix().getPij();
      MatrixTools::print(tMt, out);

      out.endLine();
    }
  }
  else
  {
    const SinglePhyloLikelihood* pS=dynamic_cast<const SinglePhyloLikelihood*>(phylolike);
    
    PhylogeneticsApplicationTools::printParameters(&pS->getSubstitutionProcess(), out);
  }
}

void PhylogeneticsApplicationTools::printAnalysisInformation(const PhyloLikelihood* phylolike, OutputStream& out)
{
  if (dynamic_cast<const SinglePhyloLikelihood*>(phylolike) != NULL)
  {
    const SinglePhyloLikelihood* pSPL = dynamic_cast<const SinglePhyloLikelihood*>(phylolike);
    const SubstitutionProcess* pSP = &pSPL->getSubstitutionProcess();

    vector<string> colNames;
    colNames.push_back("Sites");
    colNames.push_back("is.complete");
    colNames.push_back("is.constant");
    colNames.push_back("lnL");

    const DiscreteDistribution* pDD = 0;
    size_t nbR = 0;
    
    if (dynamic_cast<const RateAcrossSitesSubstitutionProcess*>(pSP) != NULL)
      pDD = &dynamic_cast<const RateAcrossSitesSubstitutionProcess*>(pSP)->getRateDistribution();
    else if (dynamic_cast<const NonHomogeneousSubstitutionProcess*>(pSP) != NULL)
      pDD = &dynamic_cast<const NonHomogeneousSubstitutionProcess*>(pSP)->getRateDistribution();

    if (pDD != NULL) {
      nbR = pDD->getNumberOfCategories();
      
      pDD->print(out);

      out.endLine();
      out.endLine();

      if (nbR>1)
        for (size_t i = 0; i < nbR; i++)
          colNames.push_back("prob"+ TextTools::toString(i+1));
    }
    
    const SiteContainer* sites = phylolike -> getData();

    vector<string> row(4+(nbR>1?nbR:0));
    DataTable* infos = new DataTable(colNames);

    VVdouble vvPP = pSPL->getPosteriorProbabilitiesOfEachClass();
    
    for (size_t i = 0; i < sites->getNumberOfSites(); i++)
    {
      double lnL = phylolike->getLogLikelihoodForASite(i);
      const Site* currentSite = &sites->getSite(i);
      int currentSitePosition = currentSite->getPosition();
      string isCompl = "NA";
      string isConst = "NA";
      try { isCompl = (SiteTools::isComplete(*currentSite) ? "1" : "0"); }
      catch(EmptySiteException& ex) {}
      try { isConst = (SiteTools::isConstant(*currentSite) ? "1" : "0"); }
      catch(EmptySiteException& ex) {}
      row[0] = (string("[" + TextTools::toString(currentSitePosition) + "]"));
      row[1] = isCompl;
      row[2] = isConst;
      row[3] = TextTools::toString(lnL);

      if (nbR>1)
        for (size_t j=0; j<nbR; j++)
          row[4+j] = TextTools::toString(vvPP[i][j]);
      
      infos->addRow(row);
    }
          
    DataTable::write(*infos, out, "\t");
    delete infos;
  }
  else if (dynamic_cast<const MultiPhyloLikelihood*>(phylolike) != NULL)
  {
    const MultiPhyloLikelihood* pMPL = dynamic_cast<const MultiPhyloLikelihood*>(phylolike);

    vector<string> colNames;
    colNames.push_back("Sites");
    colNames.push_back("is.complete");
    colNames.push_back("is.constant");
    colNames.push_back("lnL");
    
    size_t nbP = pMPL->getNumberOfSubstitutionProcess();
      
    if (nbP>1){
      for (size_t i = 0; i < nbP; i++)
        colNames.push_back("lnL"+ TextTools::toString(i+1));
      for (size_t i = 0; i < nbP; i++)
        colNames.push_back("prob"+ TextTools::toString(i+1));
    }
      
    const SiteContainer* sites = phylolike -> getData();

    vector<string> row(4+(nbP>1?2*nbP:0));
    DataTable* infos = new DataTable(colNames);

    VVdouble vvPP = pMPL->getPosteriorProbabilitiesForEachSiteForEachProcess();
    VVdouble vvL = pMPL->getLikelihoodForEachSiteForEachProcess();
    
    for (size_t i = 0; i < sites->getNumberOfSites(); i++)
    {
      double lnL = phylolike->getLogLikelihoodForASite(i);
      const Site* currentSite = &sites->getSite(i);
      int currentSitePosition = currentSite->getPosition();
      string isCompl = "NA";
      string isConst = "NA";
      try { isCompl = (SiteTools::isComplete(*currentSite) ? "1" : "0"); }
      catch(EmptySiteException& ex) {}
      try { isConst = (SiteTools::isConstant(*currentSite) ? "1" : "0"); }
      catch(EmptySiteException& ex) {}
      row[0] = (string("[" + TextTools::toString(currentSitePosition) + "]"));
      row[1] = isCompl;
      row[2] = isConst;
      row[3] = TextTools::toString(lnL);
      
      if (nbP>1){
        for (size_t j=0; j<nbP; j++)
          row[4+j] = TextTools::toString(log(vvL[i][j]));
        for (size_t j=0; j<nbP; j++)
          row[4+nbP+j] = TextTools::toString(vvPP[i][j]);
      }
      infos->addRow(row);
    }
    
    DataTable::write(*infos, out, "\t");
    delete infos;
  }
  
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printParameters(const SubstitutionModelSet* modelSet, OutputStream& out)
{
  (out << "nonhomogeneous=general").endLine();
  (out << "nonhomogeneous.number_of_models=" << modelSet->getNumberOfModels()).endLine();

  // Get the parameter links:
  map< size_t, vector<string> > modelLinks; // for each model index, stores the list of global parameters.
  map< string, set<size_t> > parameterLinks; // for each parameter name, stores the list of model indices, wich should be sorted.
  vector<string> writtenNames;

  // Loop over all models:
  for (size_t i = 0; i < modelSet->getNumberOfModels(); i++)
  {
    const SubstitutionModel* model = modelSet->getModel(i);

    // First get the aliases for this model:
    map<string, string> aliases;

    ParameterList pl=model->getParameters();

    for (size_t np = 0 ; np< pl.size() ; np++)
      {
        string nfrom=modelSet->getFrom(pl[np].getName()+"_"+TextTools::toString(i+1));
        if (nfrom!="")
          aliases[pl[np].getName()]=nfrom;
      }

    // Now print it:
    writtenNames.clear();
    out.endLine() << "model" << (i + 1) << "=";
    BppOSubstitutionModelFormat bIOsm(BppOSubstitutionModelFormat::ALL, true, true, true, false);
    map<string, string>::iterator it;
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

  // Root frequencies:
  out.endLine();
  (out << "# Root frequencies:").endLine();
  out << "nonhomogeneous.root_freq=";

  BppOFrequenciesSetFormat bIO(BppOFrequenciesSetFormat::ALL, false);
  bIO.write(modelSet->getRootFrequenciesSet(), out, writtenNames);
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printParameters(const DiscreteDistribution* rDist, OutputStream& out)
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
  string suffix)
{
  SubstitutionCount* substitutionCount = 0;
  string nijtOption;
  map<string, string> nijtParams;
  string nijtText = ApplicationTools::getStringParameter("nijt", params, "Uniformization", suffix, true);
  KeyvalTools::parseProcedure(nijtText, nijtOption, nijtParams);

  if (nijtOption == "Laplace")
  {
    int trunc = ApplicationTools::getIntParameter("trunc", nijtParams, 10, suffix, true);
    substitutionCount = new LaplaceSubstitutionCount(model, trunc);
  }
  else if (nijtOption == "Uniformization")
  {
    string weightOption = ApplicationTools::getStringParameter("weight", nijtParams, "None", "", true, false);
    AlphabetIndex2* weights = SequenceApplicationTools::getAlphabetIndex2(alphabet, weightOption, "Substitution weight scheme:");
    substitutionCount = new UniformizationSubstitutionCount(model, new TotalSubstitutionRegister(alphabet), weights);
  }
  else if (nijtOption == "Decomposition")
  {
    string weightOption = ApplicationTools::getStringParameter("weight", nijtParams, "None", "", true, false);
    AlphabetIndex2* weights = SequenceApplicationTools::getAlphabetIndex2(alphabet, weightOption, "Substitution weight scheme:");
    const ReversibleSubstitutionModel* revModel = dynamic_cast<const ReversibleSubstitutionModel*>(model);
    if (revModel)
      substitutionCount = new DecompositionSubstitutionCount(revModel, new TotalSubstitutionRegister(alphabet), weights);
    else
      throw Exception("Decomposition method can only be used with reversible substitution models.");
  }
  else if (nijtOption == "Naive")
  {
    string weightOption = ApplicationTools::getStringParameter("weight", nijtParams, "None", "", true, false);
    AlphabetIndex2* weights = SequenceApplicationTools::getAlphabetIndex2(alphabet, weightOption, "Substitution weight scheme:");
    substitutionCount = new NaiveSubstitutionCount(new TotalSubstitutionRegister(alphabet), false, weights);
  }
  else if (nijtOption == "Label")
  {
    substitutionCount = reinterpret_cast<SubstitutionCount*>(new LabelSubstitutionCount(alphabet));
  }
  else if (nijtOption == "ProbOneJump")
  {
    substitutionCount = reinterpret_cast<SubstitutionCount*>(new OneJumpSubstitutionCount(model));
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

/******************************************************************************/

