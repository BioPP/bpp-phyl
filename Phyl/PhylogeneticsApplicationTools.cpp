//
// File: PhylogeneticsApplicationTools.cpp
// Created by: Julien Dutheil
// Created on: Fri Oct 21 16:49 2005
// from old file ApplicationTools.cpp created on Sun Dec 14 09:36:26 2003
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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
#include "models"
#include "OptimizationTools.h"
#include "Tree.h"
#include "Newick.h"
#include "RHomogeneousTreeLikelihood.h"
#include "RNonHomogeneousTreeLikelihood.h"
#include "DRHomogeneousTreeLikelihood.h"
#include "NNIHomogeneousTreeLikelihood.h"
#include "RHomogeneousClockTreeLikelihood.h"

// From Utils:
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>
#include <Utils/ApplicationTools.h>
#include <Utils/StringTokenizer.h>
#include <Utils/KeyvalTools.h>

// From NumCalc:
#include <NumCalc/ConstantDistribution.h>
#include <NumCalc/GammaDiscreteDistribution.h>
#include <NumCalc/InvariantMixedDiscreteDistribution.h>
#include <NumCalc/optimizers>

// From SeqLib:
#include <Seq/AlphabetTools.h>
#include <Seq/SequenceContainerTools.h>

using namespace bpp;

// From the STL:
#include <fstream>
#include <iomanip>

using namespace std;

/******************************************************************************/

TreeTemplate<Node> * PhylogeneticsApplicationTools::getTree(
  map<string, string> & params,
  const string & suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  string treeFilePath = ApplicationTools::getAFilePath("tree.file", params, true, true, suffix, suffixIsOptional);
  
  //Read the tree file:
  Newick newick(true);
  TreeTemplate<Node> * tree = dynamic_cast<TreeTemplate<Node> *>(newick.read(treeFilePath));
  if(verbose) ApplicationTools::displayResult("Tree file", treeFilePath);
  return tree;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printInputTreeHelp()
{
  if(!ApplicationTools::message) return;
  *ApplicationTools::message << "Input tree parameters:" << endl;
  *ApplicationTools::message << "tree.file                     | file from where to read the tree" << endl;
  *ApplicationTools::message << "______________________________|___________________________________________" << endl;
}

/******************************************************************************/

SubstitutionModel* PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance(
    const Alphabet* alphabet,
    const string& modelDescription,
    map<string, string>& unparsedParameterValues,
    bool allowCovarions,
    bool allowGaps,
    bool verbose) throw (Exception)
{
  SubstitutionModel* model = 0;
  string modelName = "", left = "";
  map<string, string> args;
  KeyvalTools::parseProcedure(modelDescription, modelName, args);

  if (modelName == "RE08")
  {
    if (!allowGaps)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance. No Gap model allowed here.");
    
    //We have to parse the nested model first:
    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance. Missing argument 'model' for model 'RE08'.");
    if (verbose)
      ApplicationTools::displayResult("Gap model" , modelName);
    map<string, string> unparsedParameterValuesNested;
    SubstitutionModel* nestedModel = getSubstitutionModelDefaultInstance(alphabet, nestedModelDescription, unparsedParameterValuesNested, allowCovarions, false, verbose);
    
    //Now we create the RE08 substitution model:
    ReversibleSubstitutionModel * tmp = dynamic_cast<ReversibleSubstitutionModel *>(nestedModel);
    model = new RE08(tmp);

    //Then we update the parameter set:
    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      unparsedParameterValues["RE08." + it->first] = it->second;
    if (args.find("lambda") == args.end())
      unparsedParameterValues["RE08.lambda"] = args["lambda"];
    if (args.find("mu") == args.end())
      unparsedParameterValues["RE08.mu"] = args["mu"];
  }  
  else if (modelName == "TS98")
  {
    if (!allowCovarions)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance. No Covarion model allowed here.");

    //We have to parse the nested model first:
    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance. Missing argument 'model' for model 'TS98'.");
    if (verbose)
      ApplicationTools::displayResult("Covarion model" , modelName);
    map<string, string> unparsedParameterValuesNested;
    SubstitutionModel* nestedModel = getSubstitutionModelDefaultInstance(alphabet, nestedModelDescription, unparsedParameterValuesNested, false, allowGaps, verbose);
    
    //Now we create the TS98 substitution model:
    ReversibleSubstitutionModel * tmp = dynamic_cast<ReversibleSubstitutionModel *>(nestedModel);
    model = new TS98(tmp);

    //Then we update the parameter set:
    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      unparsedParameterValues["TS98." + it->first] = it->second;
    if (args.find("s1") == args.end())
      unparsedParameterValues["TS98.s1"] = args["s1"];
    if (args.find("s2") == args.end())
      unparsedParameterValues["TS98.s2"] = args["s2"];
  }
  else if (modelName == "G01")
  {
    if(!allowCovarions)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance. No Covarion model allowed here.");

    //We have to parse the nested model first:
    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance. Missing argument 'model' for model 'G01'.");
    string nestedRateDistDescription = args["rdist"];
    if (TextTools::isEmpty(nestedRateDistDescription))
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance. Missing argument 'rdist' for model 'G01'.");
    if (verbose)
      ApplicationTools::displayResult("Covarion model" , modelName);
   
    map<string, string> unparsedParameterValuesNestedModel;
    SubstitutionModel* nestedModel = getSubstitutionModelDefaultInstance(alphabet, nestedModelDescription, unparsedParameterValuesNestedModel, false, allowGaps, verbose);
    map<string, string> unparsedParameterValuesNestedDist;
    DiscreteDistribution* nestedRDist = getRateDistributionDefaultInstance(nestedRateDistDescription, unparsedParameterValuesNestedDist, false, verbose);

    //Now we create the G01 substitution model:
    ReversibleSubstitutionModel * tmp = dynamic_cast<ReversibleSubstitutionModel *>(nestedModel);
    model = new G2001(tmp, nestedRDist);
    
    //Then we update the parameter set:
    for (map<string, string>::iterator it = unparsedParameterValuesNestedModel.begin(); it != unparsedParameterValuesNestedModel.end(); it++)
      unparsedParameterValues["G01." + it->first] = it->second;
    for (map<string, string>::iterator it = unparsedParameterValuesNestedDist.begin(); it != unparsedParameterValuesNestedDist.end(); it++)
      unparsedParameterValues["G01." + it->first] = it->second;
    if (args.find("nu") == args.end())
      unparsedParameterValues["G01.nu"] = args["nu"];
  }
  else
  {
    //This is a 'simple' model...
    if (AlphabetTools::isNucleicAlphabet(alphabet))
    {
      const NucleicAlphabet * alpha = dynamic_cast<const NucleicAlphabet *>(alphabet);
    
      if (modelName == "GTR")
      {
        model = new GTR(alpha);
        if (args.find("a") != args.end())
          unparsedParameterValues["GTR.a"] = args["a"];
        if (args.find("b") != args.end())
           unparsedParameterValues["GTR.b"] = args["b"];
        if (args.find("c") != args.end())
          unparsedParameterValues["GTR.c"] = args["c"];
        if (args.find("d") != args.end())
          unparsedParameterValues["GTR.d"] = args["d"];
        if (args.find("e") != args.end())
          unparsedParameterValues["GTR.e"] = args["e"];
        if (args.find("theta") != args.end())
          unparsedParameterValues["GTR.theta"] = args["theta"];
        if (args.find("theta1") != args.end())
          unparsedParameterValues["GTR.theta1"] = args["theta1"];
        if (args.find("theta2") != args.end())
          unparsedParameterValues["GTR.theta2"] = args["theta2"];
      }
      else if (modelName == "L95")
      {
        model = new L95(alpha);
        if (args.find("beta") != args.end())
          unparsedParameterValues["L95.beta"] = args["beta"];
        if (args.find("gamma") != args.end())
          unparsedParameterValues["L95.gamma"] = args["gamma"];
        if (args.find("delta") != args.end())
          unparsedParameterValues["L95.delta"] = args["delta"];
        if (args.find("theta") != args.end())
          unparsedParameterValues["L95.theta"] = args["theta"];
      }
      else if (modelName == "TN93")
      {
        model = new TN93(alpha);
        if (args.find("kappa1") != args.end())
          unparsedParameterValues["TN93.kappa1"] = args["kappa1"];
        if (args.find("kappa2") != args.end())
          unparsedParameterValues["TN93.kappa2"] = args["kappa2"];
        if (args.find("theta") != args.end())
          unparsedParameterValues["TN93.theta"] = args["theta"];
        if (args.find("theta1") != args.end())
          unparsedParameterValues["TN93.theta1"] = args["theta1"];
        if (args.find("theta2") != args.end())
          unparsedParameterValues["TN93.theta2"] = args["theta2"];
      }
      else if (modelName == "HKY85")
      {
        model = new HKY85(alpha);
        if (args.find("kappa") != args.end())
          unparsedParameterValues["HKY85.kappa"] = args["kappa"];
        if (args.find("theta") != args.end())
          unparsedParameterValues["HKY85.theta"] = args["theta"];
        if (args.find("theta1") != args.end())
          unparsedParameterValues["HKY85.theta1"] = args["theta1"];
        if (args.find("theta2") != args.end())
          unparsedParameterValues["HKY85.theta2"] = args["theta2"];
      }
      else if (modelName == "F84")
      {
        model = new F84(alpha);
        if (args.find("kappa") != args.end())
          unparsedParameterValues["F84.kappa"] = args["kappa"];
        if (args.find("theta") != args.end())
          unparsedParameterValues["F84.theta"] = args["theta"];
        if (args.find("theta1") != args.end())
          unparsedParameterValues["F84.theta1"] = args["theta1"];
        if (args.find("theta2") != args.end())
          unparsedParameterValues["F84.theta2"] = args["theta2"];
      }
      else if (modelName == "T92")
      {
        model = new T92(alpha);
        if (args.find("kappa") != args.end())
          unparsedParameterValues["T92.kappa"] = args["kappa"];
        if (args.find("theta") != args.end())
          unparsedParameterValues["T92.theta"] = args["theta"];
      }
      else if (modelName == "K80")
      {
        model = new K80(alpha);
        if (args.find("kappa") != args.end())
          unparsedParameterValues["K80.kappa"] = args["kappa"];
      }
      else if (modelName == "JC69")
      {
        model = new JCnuc(alpha);
      }
      else
        throw Exception("Model '" + modelName + "' unknown.");
    }
    else
    { 
      const ProteicAlphabet * alpha = dynamic_cast<const ProteicAlphabet *>(alphabet);
    
      if (modelName == "JC69+F")
        model = new JCprotF(alpha);
      else if (modelName == "DSO78+F")
        model = new DSO78F(alpha);
      else if (modelName == "JTT92+F")
        model = new JTT92F(alpha);
      else if (modelName == "Empirical+F")
      {
        string prefix = args["name"];
        if (TextTools::isEmpty(prefix))
          throw Exception("'name' argument missing for user-defined substitution model.");
        model = new UserProteinSubstitutionModelF(alpha, args["file"], prefix + ".");
      }
      else if (modelName == "JC69")
        model = new JCprot(alpha);
      else if (modelName == "DSO78")
        model = new DSO78(alpha);
      else if (modelName == "JTT92")
        model = new JTT92(alpha);
      else if (modelName == "Empirical")
      {
        string prefix = args["name"];
        if( TextTools::isEmpty(prefix))
          throw Exception("'name' argument missing for user-defined substitution model.");
        model = new UserProteinSubstitutionModel(alpha, args["file"], prefix);
      }
      else
        throw Exception("Model '" + modelName + "' unknown.");
    }
  }
  if (verbose)
    ApplicationTools::displayResult("Substitution model", modelName);

  //Now look if some parameters are aliased:
  ParameterList pl = model->getIndependentParameters();
  string pname, pval, pname2;
  for (unsigned int i = 0; i < pl.size(); i++)
  {
    pname = model->getParameterNameWithoutNamespace(pl[i]->getName());
    if (args.find(pname) == args.end()) continue;
    pval = args[pname];
    if(pval.length() >= 5 && pval.substr(0, 5) == "model")
      continue;
    bool found = false;
    for (unsigned int j = 0; j < pl.size() && !found; j++)
    {
      pname2 = model->getParameterNameWithoutNamespace(pl[j]->getName());
      if (j == i || args.find(pname2) == args.end()) continue;
      if (pval == pname2)
      {
        //This is an alias...
        //NB: this may throw an exception if uncorrect! We leave it as is for now :s
        model->aliasParameters(pname2, pname);
        if (verbose)
          ApplicationTools::displayResult("Parameter alias found", pname + "->" + pname2);
        found = true;
      }
    }
    if (!TextTools::isDecimalNumber(pval) && !found)
      throw Exception("Incorrect parameter syntax: parameter " + pval + " was not found and can't be used as a value for parameter " + pname + ".");
  }
        
  if (args.find("useObservedFreqs") != args.end())
    unparsedParameterValues[model->getNamespace() + "useObservedFreqs"] = args["useObservedFreqs"];
  if (args.find("useObservedFreqs.pseudoCount") != args.end())
    unparsedParameterValues[model->getNamespace() + "useObservedFreqs.pseudoCount"] = args["useObservedFreqs.pseudoCount"];
    
  return model;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::setSubstitutionModelParametersInitialValues(
    SubstitutionModel* model,
    map<string, string>& unparsedParameterValues,
    const SiteContainer* data,
    bool verbose) throw (Exception)
{
  bool useObsFreq = ApplicationTools::getBooleanParameter(model->getNamespace() + "useObservedFreqs", unparsedParameterValues, false, "", true, false);
  if(verbose) ApplicationTools::displayResult("Use observed frequencies", useObsFreq ? "yes" : "no");
  if(useObsFreq && data != NULL) 
  {
    unsigned int psi = ApplicationTools::getParameter<unsigned int>(model->getNamespace() + "useObservedFreqs.pseudoCount", unparsedParameterValues, 0);
    model->setFreqFromData(*data, psi);
  }
  ParameterList pl = model->getIndependentParameters();
	for(unsigned int i = 0; i < pl.size(); i++)
  {
		Parameter* p = pl[i];
		AutoParameter* ap = new AutoParameter(* p);
		ap->setMessageHandler(ApplicationTools::warning);
		pl[i] = ap;
		delete p;
	}
  for(unsigned int i = 0; i < pl.size(); i++)
  {
    const string pName = pl[i]->getName();
    if(!useObsFreq || (model->getParameterNameWithoutNamespace(pName).substr(0,5) != "theta"))
    {
      double value = ApplicationTools::getDoubleParameter(pName, unparsedParameterValues, pl[i]->getValue()); 
      pl[i]->setValue(value);
    }
    if(verbose)
      ApplicationTools::displayResult("Parameter found", pName + "=" + TextTools::toString(pl[i]->getValue()));
  }
  model->matchParametersValues(pl);
}

/******************************************************************************/
 
SubstitutionModel * PhylogeneticsApplicationTools::getSubstitutionModel(
    const Alphabet* alphabet,
    const SiteContainer* data,
    map<string, string>& params,
    const string& suffix,
    bool suffixIsOptional,
    bool verbose) throw (Exception)
{
  string modelDescription = ApplicationTools::getStringParameter("model", params, "JC69()", suffix, suffixIsOptional, verbose);
  map<string, string> unparsedParameterValues;
  SubstitutionModel* model = getSubstitutionModelDefaultInstance(alphabet, modelDescription, unparsedParameterValues, true, true, verbose);
  setSubstitutionModelParametersInitialValues(model, unparsedParameterValues, data, verbose);
  return model;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::setSubstitutionModelParametersInitialValues(
    SubstitutionModel* model,
    map<string, string>& unparsedParameterValues,
    const string& modelPrefix,
    const SiteContainer* data,
    map<string, double> & existingParams,
    vector<string> & specificParams,
    vector<string> & sharedParams,
    bool verbose) throw (Exception)
{
  bool useObsFreq = ApplicationTools::getBooleanParameter(model->getNamespace() + "useObservedFreqs", unparsedParameterValues, false);
  if(verbose) ApplicationTools::displayResult("Use observed frequencies", useObsFreq ? "yes" : "no");
  if(useObsFreq && data != NULL) 
  {
    unsigned int psi = ApplicationTools::getParameter<unsigned int>(model->getNamespace() + "useObservedFreqs.pseudoCount", unparsedParameterValues, 0);
    model->setFreqFromData(*data, psi);
  }

  ParameterList pl = model->getIndependentParameters();
 	for(unsigned int i = 0; i < pl.size(); i++)
  {
		Parameter* p = pl[i];
		AutoParameter* ap = new AutoParameter(* p);
		ap->setMessageHandler(ApplicationTools::warning);
		pl[i] = ap;
		delete p;
	}
  for(unsigned int i = 0; i < pl.size(); i++)
  {
    const string pName = pl[i]->getName();
    string value;
    if(!useObsFreq || (model->getParameterNameWithoutNamespace(pName).substr(0,5) != "theta"))
    {
      value = ApplicationTools::getStringParameter(pName, unparsedParameterValues, TextTools::toString(pl[i]->getValue()));
      if(value.size() > 5 && value.substr(0, 5) == "model")
      {
        if(existingParams.find(value) != existingParams.end())
        {
          pl[i]->setValue(existingParams[value]);
          sharedParams.push_back(value);
        }
        else
          throw Exception("Error, unknown parameter" + modelPrefix + pName);
      }
      else
      {
        double value2 = TextTools::toDouble(value);
        existingParams[modelPrefix + pName] = value2;
        specificParams.push_back(pName);
        pl[i]->setValue(value2);
      }
    }
    else
    {
      existingParams[modelPrefix + pName] = pl[i]->getValue();
      specificParams.push_back(pName);
    }
    if(verbose)
      ApplicationTools::displayResult("Parameter found", modelPrefix + pName + "=" + TextTools::toString(pl[i]->getValue()));
  }
  model->matchParametersValues(pl);
}

/******************************************************************************/
 
FrequenciesSet * PhylogeneticsApplicationTools::getFrequenciesSet(
  const Alphabet * alphabet,
  const SiteContainer * data,
  map<string, string> & params,
  const vector<double> & rateFreqs,
  const string & suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  FrequenciesSet * rootFrequencies = NULL;
  string freqDescription = ApplicationTools::getStringParameter("nonhomogeneous.root_freq", params, "observed", suffix, suffixIsOptional);
  string freqName;
  map<string, string> args;
  KeyvalTools::parseProcedure(freqDescription, freqName, args);
  if(verbose) ApplicationTools::displayResult("Ancestral frequences method", freqName);
  if(freqName == "Observed" && data)
  {
    map<int, double> freqs = SequenceContainerTools::getFrequencies(*data);
    double t = 0;
    vector<double> rootFreq(alphabet->getSize());
    for(unsigned int i = 0; i < alphabet->getSize(); i++) t += freqs[i];
    for(unsigned int i = 0; i < alphabet->getSize(); i++) rootFreq[i] = freqs[i] / t;
    if(AlphabetTools::isNucleicAlphabet(alphabet))
    {
      double theta  = (rootFreq[1] + rootFreq[2]) / (rootFreq[0] + rootFreq[1] + rootFreq[2] + rootFreq[3]);
      double theta1 = rootFreq[0] / (rootFreq[0] + rootFreq[3]);
      double theta2 = rootFreq[2] / (rootFreq[1] + rootFreq[2]);
      rootFrequencies = new FullNAFrequenciesSet(dynamic_cast<const NucleicAlphabet *>(alphabet), theta, theta1, theta2, "RootFreq");
    }
    else if(AlphabetTools::isProteicAlphabet(alphabet))
      rootFrequencies = new FullProteinFrequenciesSet(dynamic_cast<const ProteicAlphabet *>(alphabet), rootFreq, "RootFreq");
    else
      rootFrequencies = new FullFrequenciesSet(alphabet, rootFreq, "RootFreq");
  }
  else if(freqName == "Balanced")
  {
    if(AlphabetTools::isNucleicAlphabet(alphabet))
      rootFrequencies = new FullNAFrequenciesSet(dynamic_cast<const NucleicAlphabet*>(alphabet), "RootFreq");
    else if(AlphabetTools::isProteicAlphabet(alphabet))
      rootFrequencies = new FullProteinFrequenciesSet(dynamic_cast<const ProteicAlphabet*>(alphabet), "RootFreq");
    else
      rootFrequencies = new FullFrequenciesSet(alphabet, "RootFreq");
  }
  else if(freqName == "Init")
  {
    vector<double> rootFreq(alphabet->getSize());
    if(AlphabetTools::isNucleicAlphabet(alphabet))
    {
      rootFreq[ 0] = ApplicationTools::getDoubleParameter("ancA", args, 0.25);
      rootFreq[ 1] = ApplicationTools::getDoubleParameter("ancC", args, 0.25);
      rootFreq[ 2] = ApplicationTools::getDoubleParameter("ancG", args, 0.25);
      rootFreq[ 3] = ApplicationTools::getDoubleParameter("ancT", args, 0.25);
      double theta  = (rootFreq[1] + rootFreq[2]) / (rootFreq[0] + rootFreq[1] + rootFreq[2] + rootFreq[3]);
      double theta1 = rootFreq[0] / (rootFreq[0] + rootFreq[3]);
      double theta2 = rootFreq[2] / (rootFreq[1] + rootFreq[2]);
      rootFrequencies = new FullNAFrequenciesSet(dynamic_cast<const NucleicAlphabet *>(alphabet), theta, theta1, theta2, "RootFreq");
    }
    else if(AlphabetTools::isProteicAlphabet(alphabet))
    {
      rootFreq[ 0] = ApplicationTools::getDoubleParameter("ancA", args, 0.05);
      rootFreq[ 1] = ApplicationTools::getDoubleParameter("ancR", args, 0.05);
      rootFreq[ 2] = ApplicationTools::getDoubleParameter("ancN", args, 0.05);
      rootFreq[ 3] = ApplicationTools::getDoubleParameter("ancD", args, 0.05);
      rootFreq[ 4] = ApplicationTools::getDoubleParameter("ancC", args, 0.05);
      rootFreq[ 5] = ApplicationTools::getDoubleParameter("ancQ", args, 0.05);
      rootFreq[ 6] = ApplicationTools::getDoubleParameter("ancE", args, 0.05);
      rootFreq[ 7] = ApplicationTools::getDoubleParameter("ancG", args, 0.05);
      rootFreq[ 8] = ApplicationTools::getDoubleParameter("ancH", args, 0.05);
      rootFreq[ 9] = ApplicationTools::getDoubleParameter("ancI", args, 0.05);
      rootFreq[10] = ApplicationTools::getDoubleParameter("ancL", args, 0.05);
      rootFreq[11] = ApplicationTools::getDoubleParameter("ancK", args, 0.05);
      rootFreq[12] = ApplicationTools::getDoubleParameter("ancM", args, 0.05);
      rootFreq[13] = ApplicationTools::getDoubleParameter("ancF", args, 0.05);
      rootFreq[14] = ApplicationTools::getDoubleParameter("ancP", args, 0.05);
      rootFreq[15] = ApplicationTools::getDoubleParameter("ancS", args, 0.05);
      rootFreq[16] = ApplicationTools::getDoubleParameter("ancT", args, 0.05);
      rootFreq[17] = ApplicationTools::getDoubleParameter("ancW", args, 0.05);
      rootFreq[18] = ApplicationTools::getDoubleParameter("ancY", args, 0.05);
      rootFreq[19] = ApplicationTools::getDoubleParameter("ancV", args, 0.05);
      rootFrequencies = new FullProteinFrequenciesSet(dynamic_cast<const ProteicAlphabet *>(alphabet), rootFreq, "RootFreq");
    }
    else throw Exception("Init root frequencies is only available with Nucleic and Proteic alphabet.");
  }
  else if(freqName == "ObservedGC" && data)
  {
    if(!AlphabetTools::isNucleicAlphabet(alphabet)) throw Exception("Error, unvalid option " + freqName + " with non-nucleic alphabet.");
    map<int, double> freqs = SequenceContainerTools::getFrequencies(*data);
    double theta  = (freqs[1] + freqs[2]) / (freqs[0] + freqs[1] + freqs[2] + freqs[3]);
    if(verbose) ApplicationTools::displayResult("Ancestral theta", TextTools::toString(theta));
    rootFrequencies = new GCFrequenciesSet(dynamic_cast<const NucleicAlphabet *>(alphabet), theta, "RootFreq");
  }
  else if(freqName == "BalancedGC")
  {
    if(!AlphabetTools::isNucleicAlphabet(alphabet)) throw Exception("Error, unvalid option " + freqName + " with non-nucleic alphabet.");
    rootFrequencies = new GCFrequenciesSet(dynamic_cast<const NucleicAlphabet *>(alphabet), "RootFreq");
  }
  else if(freqName == "InitGC")
  {
    if(!AlphabetTools::isNucleicAlphabet(alphabet)) throw Exception("Error, unvalid option " + freqName + " with non-nucleic alphabet.");
    double theta = ApplicationTools::getDoubleParameter("ancTheta", args, 0.5);
    rootFrequencies = new GCFrequenciesSet(dynamic_cast<const NucleicAlphabet *>(alphabet), theta, "RootFreq");
  }
  else throw Exception("Unvalid method for root frequencies: " + freqName);
    
  if(rateFreqs.size() > 0)
  {
    rootFrequencies = new MarkovModulatedFrequenciesSet(rootFrequencies, rateFreqs);
  }
  return rootFrequencies;
}

/******************************************************************************/

SubstitutionModelSet * PhylogeneticsApplicationTools::getSubstitutionModelSet(
  const Alphabet * alphabet,
  const SiteContainer * data,
  map<string, string> & params,
  const string & suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  if(!ApplicationTools::parameterExists("nonhomogeneous.number_of_models", params))
    throw Exception("You must specify this parameter: nonhomogeneous.number_of_models .");
  unsigned int nbModels = ApplicationTools::getParameter<unsigned int>("nonhomogeneous.number_of_models", params, 1, suffix, suffixIsOptional, false);
  if(nbModels == 0)
    throw Exception("The number of models can't be 0 !");
  if(verbose) ApplicationTools::displayResult("Number of distinct models", TextTools::toString(nbModels));
  
  //Deal with root frequencies, and build a new model set object:
  vector<double> rateFreqs;
  string tmpDesc = ApplicationTools::getStringParameter("model1", params, "JC69()", suffix, suffixIsOptional, verbose);
  map<string, string> unparsedParameterValues;
  SubstitutionModel* tmp = getSubstitutionModelDefaultInstance(alphabet, tmpDesc, unparsedParameterValues, true, true, verbose);
  if(tmp->getNumberOfStates() != alphabet->getSize())
  {
    //Markov-Modulated Markov Model...
    unsigned int n =(unsigned int)(tmp->getNumberOfStates() / alphabet->getSize());
    rateFreqs = vector<double>(n, 1./(double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
  }
  delete tmp;
  
  FrequenciesSet* rootFrequencies = getFrequenciesSet(alphabet, data, params, rateFreqs, suffix, suffixIsOptional, verbose);

  SubstitutionModelSet* modelSet = new SubstitutionModelSet(alphabet, rootFrequencies);
  
  // Now parse all models:
  map<string, double> existingParameters;
  for(unsigned int i = 0; i < nbModels; i++)
  {
    string prefix = "model" + TextTools::toString(i+1);
    string modelDesc = ApplicationTools::getStringParameter(prefix, params, "JC69()", suffix, suffixIsOptional, verbose);
    map<string, string> unparsedParameterValues;
    SubstitutionModel* model = getSubstitutionModelDefaultInstance(alphabet, modelDesc, unparsedParameterValues, true, true, verbose);
    prefix += ".";
    
    vector<string> specificParameters, sharedParameters;
    setSubstitutionModelParametersInitialValues(model,
        unparsedParameterValues, prefix, data,
        existingParameters, specificParameters, sharedParameters,
        verbose);
    vector<int> nodesId = ApplicationTools::getVectorParameter<int>(prefix + "nodes_id", params, ',', ':', TextTools::toString(i), suffix, suffixIsOptional, true);
    if(verbose) ApplicationTools::displayResult("Model" + TextTools::toString(i+1) + " is associated to", TextTools::toString(nodesId.size()) + " node(s).");
    //Add model and specific parameters:
    modelSet->addModel(model, nodesId, specificParameters);
    //Now set shared parameters:
    for(unsigned int j = 0; j < sharedParameters.size(); j++)
    {
      string pName = sharedParameters[j];
      string::size_type index = pName.find(".");
      if(index == string::npos) throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelSet. Bad parameter name: " + pName);
      string name = pName.substr(index + 1) + "_" + pName.substr(5, index - 5);
      modelSet->setParameterToModel(modelSet->getParameterIndex(name), modelSet->getNumberOfModels() - 1);
    }
  }

  return modelSet;
}

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

  if(distName == "Invariant")
  {
    //We have to parse the nested distribution first:
    string nestedDistDescription = args["dist"];
    if(TextTools::isEmpty(nestedDistDescription))
      throw Exception("PhylogeneticsApplicationTools::getRateDistributionDefaultInstance. Missing argument 'dist' for distribution 'Invariant'.");
    if(verbose)
      ApplicationTools::displayResult("Invariant Mixed distribution" , distName);
    map<string, string> unparsedParameterValuesNested;
    DiscreteDistribution* nestedDistribution = getRateDistributionDefaultInstance(nestedDistDescription, unparsedParameterValuesNested, constDistAllowed, verbose);
    
    //Now we create the Invariant rate distribution:
    rDist = new InvariantMixedDiscreteDistribution(nestedDistribution, 0., 0.000001, "Invariant");

    //Then we update the parameter set:
    for(map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      unparsedParameterValues["Invariant." + it->first] = it->second;
    if(args.find("p") != args.end())
      unparsedParameterValues["Invariant.p"] = args["p"];
  }
  else if(distName == "Uniform")
  {
    if(!constDistAllowed) throw Exception("You can't use a constant distribution here!");
    rDist = new ConstantDistribution(1.);
  }
  else if(distName == "Gamma")
  {
    if(args.find("n") == args.end())
      throw Exception("Missing argument 'n' (number of classes) in Gamma distribution");
    unsigned int nbClasses = TextTools::to<unsigned int>(args["n"]);
    rDist = new GammaDiscreteDistribution(nbClasses, 1., 1., "Gamma.");
    rDist->aliasParameters("alpha", "beta");
    if(args.find("alpha") != args.end())
      unparsedParameterValues["Gamma.alpha"] = args["alpha"];
  }
  else
  {
    throw Exception("Unknown distribution: " + distName + ".");
  }
  if(verbose)
  {
    ApplicationTools::displayResult("Rate distribution", distName);
    ApplicationTools::displayResult("Number of classes", TextTools::toString(rDist->getNumberOfCategories()));
  }

  return rDist;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::setRateDistributionParametersInitialValues(
  DiscreteDistribution* rDist,
  map<string, string>& unparsedParameterValues,
  bool verbose) throw (Exception)
{
  ParameterList pl = rDist->getIndependentParameters();
  for(unsigned int i = 0; i < pl.size(); i++)
  {
		Parameter * p = pl[i];
		AutoParameter * ap = new AutoParameter(* p);
		ap->setMessageHandler(ApplicationTools::warning);
		pl[i] = ap;
		delete p;
	}

  for(unsigned int i = 0; i < pl.size(); i++)
  {
    const string pName = pl[i]->getName();
    double value = ApplicationTools::getDoubleParameter(pName, unparsedParameterValues, pl[i]->getValue()); 
    pl[i]->setValue(value);
    if(verbose)
      ApplicationTools::displayResult("Parameter found", pName + "=" + TextTools::toString(pl[i]->getValue()));
  }
  rDist->matchParametersValues(pl);
  if(verbose)
  {
    for(unsigned int c = 0; c < rDist -> getNumberOfCategories(); c++)
    {
      ApplicationTools::displayResult("- Category " + TextTools::toString(c)
          + " (Pr = " + TextTools::toString(rDist->getProbability(c)) +") rate", TextTools::toString(rDist->getCategory(c)));
    }
  }
}
 

/******************************************************************************/

DiscreteDistribution* PhylogeneticsApplicationTools::getRateDistribution(
  map<string, string> & params,
  const string & suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  string distDescription = ApplicationTools::getStringParameter("rate_distribution", params, "Uniform()", suffix, suffixIsOptional);
  map<string, string> unparsedParameterValues;
  DiscreteDistribution * rDist = getRateDistributionDefaultInstance(distDescription, unparsedParameterValues, verbose);
  setRateDistributionParametersInitialValues(rDist, unparsedParameterValues, verbose);
  return rDist;
}

/******************************************************************************/

TreeLikelihood* PhylogeneticsApplicationTools::optimizeParameters(
  TreeLikelihood* tl,
  const ParameterList& parameters,
  map<string, string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose)
  throw (Exception)
{
  bool optimize = ApplicationTools::getBooleanParameter("optimization", params, true, suffix, suffixIsOptional, false);
  if(!optimize) return tl;
  
  unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2, suffix, suffixIsOptional);
  
  string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional);
  ostream * messageHandler = 
    (mhPath == "none") ? NULL :
      (mhPath == "std") ? &cout :
        new ofstream(mhPath.c_str(), ios::out);
  if(verbose) ApplicationTools::displayResult("Message handler", mhPath);

  string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional);
  ostream * profiler = 
    (prPath == "none") ? NULL :
      (prPath == "std") ? &cout :
        new ofstream(prPath.c_str(), ios::out);
  if(profiler != NULL) (*profiler) << setprecision(20);
  if(verbose) ApplicationTools::displayResult("Profiler", prPath);

  bool scaleFirst = ApplicationTools::getBooleanParameter("optimization.scale_first", params, false, suffix, suffixIsOptional, false);
  if(scaleFirst)
  {
    // We scale the tree before optimizing each branch length separately:
    if(verbose) ApplicationTools::displayMessage("Scaling the tree before optimizing each branch length separately.");
    double tolerance = ApplicationTools::getDoubleParameter("optimization.scale_first.tolerance", params, .0001, suffix, suffixIsOptional, true);
    if(verbose) ApplicationTools::displayResult("Scaling tolerance", TextTools::toString(tolerance));
    int nbEvalMax = ApplicationTools::getIntParameter("optimization.scale_first.max_number_f_eval", params, 1000000, suffix, suffixIsOptional, true);
    if(verbose) ApplicationTools::displayResult("Scaling max # f eval", TextTools::toString(nbEvalMax));
    int n = OptimizationTools::optimizeTreeScale(
      tl,
      tolerance,
      nbEvalMax,
      messageHandler,
      profiler);
    if(verbose) ApplicationTools::displayMessage("Performed " + TextTools::toString(n) + " function evaluations.");
  }

  // Should I ignore some parameters?
  ParameterList parametersToEstimate = parameters;
  string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameter", params, "", suffix, suffixIsOptional, false);
  StringTokenizer st(paramListDesc, ",");
  while(st.hasMoreToken())
  {
    try
    {
      string param = st.nextToken();
      if(param == "BrLen")
      {
        vector<string> vs = tl->getBranchLengthsParameters().getParameterNames();
        parametersToEstimate.deleteParameters(vs);
        if (verbose)
          ApplicationTools::displayResult("Parameter ignored", string("Branch lengths"));
      }
      else if(param == "Ancient")
      {
        NonHomogeneousTreeLikelihood *nhtl = dynamic_cast<NonHomogeneousTreeLikelihood *>(tl);
        if(!nhtl) ApplicationTools::displayWarning("The 'Ancient' parameters do not exist in homogeneous models, and will be ignored.");
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
    catch(ParameterNotFoundException & pnfe)
    {
      ApplicationTools::displayWarning("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
    }
  }
  
  unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000, suffix, suffixIsOptional);
  if(verbose) ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));
  
  double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", params, .000001, suffix, suffixIsOptional);
  if(verbose) ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));
  
  bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", params, false, suffix, suffixIsOptional, false);
  if(verbose) ApplicationTools::displayResult("Optimize topology", optimizeTopo ? "yes" : "no");
  string nniMethod = ApplicationTools::getStringParameter("optimization.topology.algorithm_nni.method", params, "phyml", suffix, suffixIsOptional, false);
  string nniAlgo;
  if(nniMethod == "fast")
  {
    nniAlgo = NNITopologySearch::FAST;
  }
  else if(nniMethod == "better")
  {
    nniAlgo = NNITopologySearch::BETTER;
  }
  else if(nniMethod == "phyml")
  {
    nniAlgo = NNITopologySearch::PHYML;
  }
  else throw Exception("Unknown NNI algorithm: '" + nniMethod + "'.");


  string method = ApplicationTools::getStringParameter("optimization.method", params, "DB", suffix, suffixIsOptional, false);
  string order  = ApplicationTools::getStringParameter("optimization.method.derivatives", params, "newton", suffix, suffixIsOptional, false);
  string optMethod;
  if(order == "gradient")
  {
    optMethod = OptimizationTools::OPTIMIZATION_GRADIENT;
  }
  else if(order == "newton")
  {
    optMethod = OptimizationTools::OPTIMIZATION_NEWTON;
  }
  else throw Exception("Unknown derivatives algorithm: '" + order + "'.");
  if(verbose) ApplicationTools::displayResult("Optimization method", method);
  if(verbose) ApplicationTools::displayResult("Algorithm used for derivable parameters", order);
  
  unsigned int n = 0;
  if(method == "DB")
  {
    //Uses Newton-Brent method:
    
    unsigned int nstep = ApplicationTools::getParameter<unsigned int>("optimization.method_DB.nstep", params, 1, suffix, suffixIsOptional, false);
    if(optimizeTopo)
    {
      bool optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, false);
      unsigned int n   = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, "", true, false);
      double tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional);
      double tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, suffixIsOptional);
      tl = OptimizationTools::optimizeTreeNNI(
          dynamic_cast<NNIHomogeneousTreeLikelihood *>(tl), parametersToEstimate,
          optNumFirst, tolBefore, tolDuring, nbEvalMax, n, messageHandler, profiler, optVerbose, optMethod, nstep, nniAlgo);
    }

    if(verbose && nstep > 1) ApplicationTools::displayResult("# of precision steps", TextTools::toString(nstep));
    n = OptimizationTools::optimizeNumericalParameters(
      dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl), parametersToEstimate,
      NULL, nstep, tolerance, nbEvalMax, messageHandler, profiler, optVerbose, optMethod);    
  }
  else if(method == "fullD")
  {
    //Uses Newton-raphson alogrithm with numerical derivatives when required.
    
    if(optimizeTopo)
    {
      bool optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, false);
      unsigned int n   = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, "", true, false);
      double tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional);
      double tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, suffixIsOptional);
      tl = OptimizationTools::optimizeTreeNNI2(
          dynamic_cast<NNIHomogeneousTreeLikelihood *>(tl), parametersToEstimate,
          optNumFirst, tolBefore, tolDuring, nbEvalMax, n, messageHandler, profiler, optVerbose, optMethod, nniAlgo);
    }

    n = OptimizationTools::optimizeNumericalParameters2(
      dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl), parametersToEstimate,
      NULL, tolerance, nbEvalMax, messageHandler, profiler, optVerbose, optMethod);       
  }
  else throw Exception("Unknown optimization method: " + method);

  string finalMethod = ApplicationTools::getStringParameter("optimization.final", params, "none", suffix, suffixIsOptional, true);
  Optimizer * finalOptimizer  = NULL;
  if(finalMethod == "none") {}
  else if(finalMethod == "simplex")
  {
    finalOptimizer = new DownhillSimplexMethod(tl);
  }
  else if(finalMethod == "powell")
  {
    finalOptimizer = new PowellMultiDimensions(tl);
  }
  else throw Exception("Unknown final optimization method: " + finalMethod);

  if(finalOptimizer)
  {
    parametersToEstimate.matchParametersValues(tl->getParameters());
    if(verbose) ApplicationTools::displayResult("Final optimization step", finalMethod);
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
  
  if(verbose) ApplicationTools::displayResult("Performed", TextTools::toString(n) + " function evaluations.");
  return tl;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::optimizeParameters(
  DiscreteRatesAcrossSitesClockTreeLikelihood * tl,
  const ParameterList& parameters,
  map<string, string> & params,
  const string & suffix,
  bool suffixIsOptional,
  bool verbose)
  throw (Exception)
{
  bool optimize = ApplicationTools::getBooleanParameter("optimization", params, true, suffix, suffixIsOptional, false);
  if(!optimize) return;
  
  unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2, suffix, suffixIsOptional);
  
  string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional);
  ostream * messageHandler = 
    (mhPath == "none") ? NULL :
      (mhPath == "std") ? &cout :
        new ofstream(mhPath.c_str(), ios::out);
  if(verbose) ApplicationTools::displayResult("Message handler", mhPath);

  string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional);
  ostream * profiler = 
    (prPath == "none") ? NULL :
      (prPath == "std") ? &cout :
        new ofstream(prPath.c_str(), ios::out);
  if(profiler != NULL) (*profiler) << setprecision(20);
  if(verbose) ApplicationTools::displayResult("Profiler", prPath);

  ParameterList parametersToEstimate = parameters;

  // Should I ignore some parameters?
  string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameter", params, "", suffix, suffixIsOptional, false);
  StringTokenizer st(paramListDesc, ",");
  while(st.hasMoreToken())
  {
    try
    {
      string param = st.nextToken();
      if(param == "BrLen")
      {
        vector<string> vs = tl->getBranchLengthsParameters().getParameterNames();
        parametersToEstimate.deleteParameters(vs);
        if (verbose)
          ApplicationTools::displayResult("Parameter ignored", string("Branch lengths"));
      }
      else if(param == "Ancient")
      {
        NonHomogeneousTreeLikelihood *nhtl = dynamic_cast<NonHomogeneousTreeLikelihood *>(tl);
        if(!nhtl) ApplicationTools::displayWarning("The 'Ancient' parameters do not exist in homogeneous models, and will be ignored.");
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
    catch(ParameterNotFoundException & pnfe)
    {
      ApplicationTools::displayError("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
    }
  }
  
  unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000, suffix, suffixIsOptional);
  if(verbose) ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));
  
  double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", params, .000001, suffix, suffixIsOptional);
  if(verbose) ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));
  
  string method = ApplicationTools::getStringParameter("optimization.method", params, "DB", suffix, suffixIsOptional, false);
  string order  = ApplicationTools::getStringParameter("optimization.method.derivatives", params, "gradient", suffix, suffixIsOptional, false);
  string optMethod, derMethod;
  if(order == "gradient")
  {
    optMethod = OptimizationTools::OPTIMIZATION_GRADIENT;
  }
  else if(order == "newton")
  {
    optMethod = OptimizationTools::OPTIMIZATION_NEWTON;
  }
  else throw Exception("Option '" + order + "' is not known for 'optimization.method.derivatives'.");
  if(verbose) ApplicationTools::displayResult("Optimization method", method);
  if(verbose) ApplicationTools::displayResult("Algorithm used for derivable parameters", order);
  
  unsigned int n = 0;
  if(method == "DB")
  {
    //Uses Newton-Brent method:
    unsigned int nstep = ApplicationTools::getParameter<unsigned int>("optimization.method_DB.nstep", params, 1, suffix, suffixIsOptional, false);
    if(verbose && nstep > 1) ApplicationTools::displayResult("# of precision steps", TextTools::toString(nstep));
    n = OptimizationTools::optimizeNumericalParametersWithGlobalClock(
      tl,
      parametersToEstimate,
      NULL,
      nstep,
      tolerance,
      nbEvalMax,
      messageHandler,
      profiler,
      optVerbose,
      optMethod);    
  }
  else if(method == "fullD")
  {
    //Uses Newton-raphson alogrithm with numerical derivatives when required.
    n = OptimizationTools::optimizeNumericalParametersWithGlobalClock2(
      tl,
      parametersToEstimate,
      NULL,
      tolerance,
      nbEvalMax,
      messageHandler,
      profiler,
      optVerbose,
      optMethod);       
  }
  else throw Exception("Unknown optimization method: " + method);

  string finalMethod = ApplicationTools::getStringParameter("optimization.final", params, "none", suffix, suffixIsOptional, false);
  Optimizer * finalOptimizer  = NULL;
  if(finalMethod == "none") {}
  else if(finalMethod == "simplex")
  {
    finalOptimizer = new DownhillSimplexMethod(tl);
  }
  else if(finalMethod == "powell")
  {
    finalOptimizer = new PowellMultiDimensions(tl);
  }
  else throw Exception("Unknown final optimization method: " + finalMethod);

  if(finalOptimizer)
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
  
  if(verbose) ApplicationTools::displayResult("Performed", TextTools::toString(n) + " function evaluations.");
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printOptimizationHelp(bool topo, bool clock)
{
  if(!ApplicationTools::message) return;
  *ApplicationTools::message << "optimization                   | [yes/no] optimize parameters?" << endl;
  *ApplicationTools::message << "optimization.method            | [DB/fullD] method to use" << endl;
  *ApplicationTools::message << "optimization.method.derivatives| [gradient/newton] use Conjugate gradient" << endl;
  *ApplicationTools::message << "                               | or Newton-Raphson" << endl;
  *ApplicationTools::message << "optimization.method_DB.step    | Number of progressive step to perform." << endl;
  *ApplicationTools::message << "optimization.final             | [none|simplex|powell] final step." << endl;
  *ApplicationTools::message << "optimization.verbose           | [0,1,2] level of verbose" << endl;
  *ApplicationTools::message << "optimization.message_handler   | [none, std ot file path] where to dislay" << endl;
  *ApplicationTools::message << "                               | optimization messages" << endl;
  *ApplicationTools::message << "                               | (if std, uses 'cout' to display messages)." << endl;
  *ApplicationTools::message << "optimization.profiler          | [none, std ot file path] where to display" << endl;
  *ApplicationTools::message << "                               | optimization steps (if std, uses 'cout'" << endl;
  *ApplicationTools::message << "                               | to display optimization steps)." << endl;
  *ApplicationTools::message << "optimization.tolerance         | [double] tolerance parameter for stopping" << endl;
  *ApplicationTools::message << "                               | the estimation." << endl;
  *ApplicationTools::message << "optimization.max_number_f_eval | [int] max. # of likelihood computations." << endl;
  *ApplicationTools::message << "optimization.ignore_parameter  | [list] parameters to ignore during optimization." << endl;
  if(!clock)
  {
  *ApplicationTools::message << "optimization.scale_first       | [yes, no] tell if a global scale" << endl;
  *ApplicationTools::message << "                               | optimization must be done prior to" << endl;
  *ApplicationTools::message << "                               | separate estimation of branch lengths." << endl;
  *ApplicationTools::message << "optimization.scale_first       | " << endl;
  *ApplicationTools::message << "                     .tolerance| [double] tolerance parameter for global" << endl;
  *ApplicationTools::message << "                               | scale optimization." << endl;
  *ApplicationTools::message << "             .max_number_f_eval| [int] maximum number of computation for" << endl;
  *ApplicationTools::message << "                               | global scale optimization." << endl;
  *ApplicationTools::message << "_______________________________|__________________________________________" << endl;
  }
  if(topo && !clock)
  {
  *ApplicationTools::message << "optimization.topology          | [yes/no] Optimize tree topology?" << endl;
  *ApplicationTools::message << "optimization.topology.algorithm| [nni] Topology movements to use." << endl;
  *ApplicationTools::message << "optimization.topology.nstep    | estimate numerical parameters every 'n'" << endl;
  *ApplicationTools::message << "                               | topology movement rounds." << endl;
  *ApplicationTools::message << "optimization.topology.numfirst | [yes/no] Optimize num. parameters first?" << endl;
  *ApplicationTools::message << "optimization.topology.tolerance| " << endl;
  *ApplicationTools::message << "                        .before| Tolerance for prior estimation." << endl;
  *ApplicationTools::message << "                        .during| Tolerance during estimation." << endl;
  }
  *ApplicationTools::message << "_______________________________|__________________________________________" << endl;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::writeTree(
  const TreeTemplate<Node> & tree,
  map<string, string> & params,
  const string & suffix,
  bool verbose) throw (Exception)
{
  string file = ApplicationTools::getAFilePath("output.tree.file", params, true, false, suffix, false);
  Newick newick;
  newick.write(tree, file, true);
  if(verbose) ApplicationTools::displayResult("Wrote tree to file ", file);
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printOutputTreeHelp()
{
  if(!ApplicationTools::message) return;
  *ApplicationTools::message << "Output tree parameters:" << endl;
  *ApplicationTools::message << "output.tree.file              | file where to write the tree" << endl;
  *ApplicationTools::message << "______________________________|___________________________________________" << endl;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printParameters(const SubstitutionModel* model, ostream& out)
{
  const UserProteinSubstitutionModel * trial1 = dynamic_cast<const UserProteinSubstitutionModel *>(model);
  if(trial1)
  {
    out << "model = Empirical(file=" << trial1->getPath() << ")" << endl;
  }
  else
  {
    const UserProteinSubstitutionModelF * trial2 = dynamic_cast<const UserProteinSubstitutionModelF *>(model);
    if(trial2)
    {
      out << "model = Empirical+F(file=" << trial2->getPath() << ")" << endl;
    }
    else
    {
      out << "model = " << model->getName() << "(";
      ParameterList pl = model->getParameters();
      for(unsigned int i = 0; i < pl.size(); i++)
      {
        if (i > 0) out << ", ";
        out << model->getParameterNameWithoutNamespace(pl[i]->getName()) << "=" << pl[i]->getValue();
      }
      out << ")" << endl;
    }
  }
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printParameters(const SubstitutionModelSet* modelSet, ostream& out)
{
  out << "nonhomogeneous = general" << endl;
  out << "nonhomogeneous.number_of_models = " << modelSet->getNumberOfModels() << endl;

  //Get the parameter links:
  map< unsigned int, vector<string> > modelLinks; // for each model index, stores the list of global parameters.
  map< string, vector<unsigned int> > parameterLinks; // for each parameter name, stores the list of model indices.
  ParameterList pl = modelSet->getParameters();
  ParameterList plroot = modelSet->getRootFrequenciesParameters();
  for(unsigned int i = 0; i < pl.size(); i++)
  {
    if(!plroot.hasParameter(pl[i]->getName()))
    {
      string name = pl[i]->getName();
      vector<unsigned int> models = modelSet->getModelsWithParameter(name);
      for(unsigned int j = 0; j < models.size(); j++)
      {
        modelLinks[models[j]].push_back(name);
        parameterLinks[name].push_back(models[j]);
      }
    }
  }

  //Loop over all models:
  for(unsigned int i = 0; i < modelSet->getNumberOfModels(); i++)
  {
    const SubstitutionModel* model = modelSet->getModel(i);
    out << endl;
    const UserProteinSubstitutionModel * trial1 = dynamic_cast<const UserProteinSubstitutionModel *>(model);
    if(trial1)
    {
      out << "model" << (i+1) << " = Empirical(file=" << trial1->getPath() << ")" << endl;
    }
    else
    {
      const UserProteinSubstitutionModelF * trial2 = dynamic_cast<const UserProteinSubstitutionModelF *>(model);
      if(trial2)
      {
        out << "model" << (i+1) << " = Empirical+F(file=" << trial2->getPath() << ")" << endl;
      }
      else
      {
        out << "model" << (i+1) << " = " << model->getName() << "(";
        vector<string> names = modelLinks[i];
        for(unsigned int j = 0; j < names.size(); j++)
        {
          if (j > 0) out << ", ";
          const string name = names[j];
          if(parameterLinks[name].size() == 1)
          {
            out << modelSet->getParameterModelName(name) << "=" << modelSet->getParameterValue(name);
          }
          else
          {
            //size must be > 1 !
            if(parameterLinks[name][0] == i)
            {
              out << modelSet->getParameterNameWithoutNamespace(modelSet->getParameterModelName(name)) << "=" << modelSet->getParameterValue(name);
            }
            else
            {
              out << modelSet->getParameterNameWithoutNamespace(modelSet->getParameterModelName(name)) << "=model" << (parameterLinks[name][0] + 1) << "." << modelSet->getParameterModelName(name);
            }
          }
        }
        out << ")" << endl;
      }
    }
    vector<int> ids = modelSet->getNodesWithModel(i);
    out << "model" << (i+1) << ".nodes_id = " << ids[0];
    for(unsigned int j = 1; j < ids.size(); j++)
      out << "," << ids[j];
    out << endl;
  }
 
  //Root frequencies:
  out << endl;
  out << "# Root frequencies:" << endl;
  if(plroot.size() == 1 && plroot[0]->getName() == "RootFreqtheta")
  {
    out << "nonhomogeneous.root_freq = InitGC(ancTheta=" << plroot[0]->getValue() << ")" << endl;
  }
  else
  {
    out << "nonhomogeneous.root_freq = Init(";
    vector<double> rootFreqs;
    try
    {
      const MarkovModulatedFrequenciesSet* mmFreqSet = dynamic_cast<const MarkovModulatedFrequenciesSet *>(modelSet->getRootFrequenciesSet());
      if(!mmFreqSet) throw Exception("");
      rootFreqs = mmFreqSet->getStatesFrequenciesSet()->getFrequencies();
    }
    catch(exception& e)
    {
      rootFreqs = modelSet->getRootFrequencies();
    }
    for(unsigned int i = 0; i < rootFreqs.size(); i++)
    {
      if (i > 0) out << ", ";
      out << "anc" << modelSet->getAlphabet()->intToChar((int)i) << "=" << rootFreqs[i];
    }
    out << ")" << endl;
  }
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printParameters(const DiscreteDistribution* rDist, ostream& out)
{
  out << "rate_distribution = ";
  const InvariantMixedDiscreteDistribution * invar = dynamic_cast<const InvariantMixedDiscreteDistribution *>(rDist);
  const DiscreteDistribution * dist = rDist;
  if (invar)
  {
    dist = invar->getVariableSubDistribution();
    cout << "Invariant(dist=";
  }
  const DiscreteDistribution* test;
  test = dynamic_cast<const ConstantDistribution *>(dist);
  if(test) out << "Uniform()";
  else
  {
    test = dynamic_cast<const GammaDiscreteDistribution *>(dist);
    if(test) out << "Gamma(n=" << rDist->getNumberOfCategories() << ", alpha=" << rDist->getParameter("alpha").getValue() << ")";
    else throw Exception("PhylogeneticsApplicationTools::printParameters(DiscreteDistribution). Unsupported distribution.");
  }
  if (invar)
    out << ")";
  out << endl;
}

/******************************************************************************/

