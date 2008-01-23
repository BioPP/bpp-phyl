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
  *ApplicationTools::message << "tree.file                     | file where to write the tree" << endl;
  *ApplicationTools::message << "______________________________|___________________________________________" << endl;
}

/******************************************************************************/

SubstitutionModel * PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance(
  const Alphabet * alphabet,
  const SiteContainer * data,
  map<string, string> & params,
  const string & prefix,
  const string & suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  SubstitutionModel * model = NULL;
  string covarionName = "none";

  if(AlphabetTools::isNucleicAlphabet(alphabet))
  {
    string modelName = ApplicationTools::getStringParameter(prefix + "name", params, "JCnuc", suffix, suffixIsOptional);
    if(modelName.find("+") != string::npos)
    {
      StringTokenizer st(modelName, "+");
      modelName = st.nextToken();
      covarionName = st.nextToken();
    }
    const NucleicAlphabet * alpha = dynamic_cast<const NucleicAlphabet *>(alphabet);
    bool useObsFreq = ApplicationTools::getBooleanParameter(prefix + "use_observed_freq", params, false, suffix, suffixIsOptional, false);
    
    if(verbose) ApplicationTools::displayResult("Substitution model" + suffix, modelName);

    if(modelName == "GTR")
    {
      model = new GTR(alpha);
      if(useObsFreq && data != NULL) model->setFreqFromData(*data);
    }
    else if(modelName == "TN93")
    {
      model = new TN93(alpha);
      if(useObsFreq && data != NULL) model->setFreqFromData(*data);
    }
    else if(modelName == "HKY85")
    {
      model = new HKY85(alpha);
      if(useObsFreq && data != NULL) model->setFreqFromData(*data);
    }
    else if(modelName == "F84")
    {
      model = new F84(alpha);
      if(useObsFreq && data != NULL) model->setFreqFromData(*data);
    }
    else if(modelName == "T92")
    {
      model = new T92(alpha);
      if(useObsFreq && data != NULL) model->setFreqFromData(*data);
    }
    else if(modelName == "K80")
    {
      model = new K80(alpha);
    }
    else if(modelName == "JCnuc")
    {
      model = new JCnuc(alpha);
    }
    else
    {
      throw Exception("Model '" + modelName + "' unknown.");
    }
  }
  else
  { 
    // Alphabet supposed to be proteic!
    string modelName = ApplicationTools::getStringParameter(prefix + "name", params, "JCprot", suffix, suffixIsOptional);
    if(modelName.find("+") != string::npos)
    {
      StringTokenizer st(modelName, "+");
      modelName = st.nextToken();
      covarionName = st.nextToken();
    }
    const ProteicAlphabet * alpha = dynamic_cast<const ProteicAlphabet *>(alphabet);
    bool useObsFreq = ApplicationTools::getBooleanParameter(prefix + "use_observed_freq", params, false, suffix, suffixIsOptional);
    
    if(modelName == "JCprot")
    {
      model = new JCprot(alpha);
    }
    else if(modelName == "DSO78")
    {
      model = new DSO78(alpha);
    }
    else if(modelName == "JTT92")
    {
      model = new JTT92(alpha);
    }
    else if(modelName == "empirical")
    {
      string file = ApplicationTools::getAFilePath("model_empirical.file", params, true, true, suffix, true);
      model = new UserProteinSubstitutionModel(alpha, file);
    }
    else
    {
      throw Exception("Model '" + modelName + "' unknown.");
    }
    if(useObsFreq && data != NULL) model->setFreqFromData(*data);
    if(verbose)
    {
      ApplicationTools::displayResult("Substitution model", modelName + (useObsFreq && (model != NULL) ? "-F" : ""));
    }
  }

  if(covarionName == "none") {}
  else if(covarionName == "G2001")
  {
    if(verbose)
    {
      ApplicationTools::displayResult("Covarion model" , covarionName);
    }
    DiscreteDistribution * rDist = getRateDistributionDefaultInstance(params, false, prefix, suffix, suffixIsOptional, verbose);
    ReversibleSubstitutionModel * tmp = dynamic_cast<ReversibleSubstitutionModel *>(model);
    model = new G2001(tmp, rDist); //The instance will delete the rDist object.
  }
  else if(covarionName == "TS98")
  {
    if(verbose)
    {
      ApplicationTools::displayResult("Covarion model" , covarionName);
    }
    ReversibleSubstitutionModel * tmp = dynamic_cast<ReversibleSubstitutionModel *>(model);
    model = new TS98(tmp);
  }
  else
  {
    throw Exception("Process unknown: " + covarionName + ".");
  }
  
  return model;
}

void PhylogeneticsApplicationTools::setSubstitutionModelParametersInitialValues(
  SubstitutionModel * model,
  map<string, string> & params,
  const string & prefix,
  const string & suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  ParameterList pl = model->getParameters();
  bool useObsFreq = ApplicationTools::getBooleanParameter(prefix + "use_observed_freq", params, false, suffix, suffixIsOptional, false);
  for(unsigned int i = 0; i < pl.size(); i++)
  {
    const string pName = pl[i]->getName();
    if(!useObsFreq || (pName.substr(0,5) != "theta"))
    {
      double value = ApplicationTools::getDoubleParameter(prefix + pName, params, pl[i]->getValue(), suffix, suffixIsOptional); 
      pl[i]->setValue(value);
    }
    if(verbose) ApplicationTools::displayResult(prefix + pName, TextTools::toString(pl[i]->getValue()));
  }
  model->matchParametersValues(pl);
}
 
SubstitutionModel * PhylogeneticsApplicationTools::getSubstitutionModel(
  const Alphabet * alphabet,
  const SiteContainer * data,
  map<string, string> & params,
  const string & suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  SubstitutionModel * model = getSubstitutionModelDefaultInstance(alphabet, data, params, "model.", suffix, suffixIsOptional, verbose);
  setSubstitutionModelParametersInitialValues(model, params, "model.", suffix, suffixIsOptional, verbose);
  return model;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printSubstitutionModelHelp()
{
  if(!ApplicationTools::message) return;
  *ApplicationTools::message << "Substitution Model:" << endl;
  *ApplicationTools::message << "model                   | Nucleotides (N): [JCnuc, K80, T92, F84, HKY85, TN93," << endl;
  *ApplicationTools::message << "                        | GTR]" << endl;
  *ApplicationTools::message << "                        | Proteins (P): [JCprot, DSO78, JTT92, empirical]" << endl;
  *ApplicationTools::message << "model.kappa             | kappa(N)  parameter in Q matrix" << endl;
  *ApplicationTools::message << "model.kappa1            | kappa1(N) parameter in Q matrix" << endl;
  *ApplicationTools::message << "model.kappa2            | kappa2(N) parameter in Q matrix" << endl;
  *ApplicationTools::message << "model.a,b,c,d,e,f       | GTR rates parameter in Q matrix" << endl;
  *ApplicationTools::message << "model.theta             | piG + piC" << endl;
  *ApplicationTools::message << "model.theta1            | piA / (piA + piT)" << endl;
  *ApplicationTools::message << "model.theta2            | piG / (piC + piG)" << endl;
  *ApplicationTools::message << "model.use_observed_freq | (N,P) Tell if the observed frequencies must be used." << endl; 
  *ApplicationTools::message << "model_empirical.file    | (P) The path toward data file to use (PAML format)." << endl; 
  *ApplicationTools::message << "________________________|_____________________________________________________" << endl;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::setSubstitutionModelParametersInitialValues(
    SubstitutionModel * model,
    map<string, double> & existingParams,
    vector<string> & specificParams,
    vector<string> & sharedParams,
    map<string, string> & params,
    const string & prefix,
    const string & suffix,
    bool suffixIsOptional,
    bool verbose) throw (Exception)
{
  ParameterList pl = model->getParameters();
  bool useObsFreq = ApplicationTools::getBooleanParameter(prefix + "use_observed_freq", params, false, suffix, suffixIsOptional);
  for(unsigned int i = 0; i < pl.size(); i++)
  {
    const string pName = pl[i]->getName();
    if(useObsFreq && (pName == "piA" || pName == "piC" || pName == "piG" || pName == "piT")) continue;
    string value = ApplicationTools::getStringParameter(prefix + pName, params, TextTools::toString(pl[i]->getValue()), suffix, suffixIsOptional);
    if(value.size() > 5 && value.substr(0, 5) == "model")
    {
      if(existingParams.find(value) != existingParams.end())
      {
        pl[i]->setValue(existingParams[value]);
        sharedParams.push_back(value);
      }
      else
      {
        throw Exception("Error, unknown parameter" + prefix + pName);
      }
    }
    else
    {
      double value2 = TextTools::toDouble(value);
      existingParams[prefix + pName] = value2;
      specificParams.push_back(pName);
      pl[i]->setValue(value2);
    }
    if(verbose) ApplicationTools::displayResult(prefix + pName, TextTools::toString(pl[i]->getValue()));
  }
  model->matchParametersValues(pl);
}
 
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
  string freqOpt = ApplicationTools::getStringParameter("nonhomogeneous.root_freq", params, "observed", suffix, suffixIsOptional);
  if(verbose) ApplicationTools::displayResult("Ancestral frequences method", freqOpt);
  if(freqOpt == "observed" && data)
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
    else
      rootFrequencies = new FullFrequenciesSet(alphabet, rootFreq, "RootFreq");
  }
  else if(freqOpt == "balanced")
  {
    if(AlphabetTools::isNucleicAlphabet(alphabet))
      rootFrequencies = new FullNAFrequenciesSet(dynamic_cast<const NucleicAlphabet*>(alphabet), "RootFreq");
    else
      rootFrequencies = new FullFrequenciesSet(alphabet, "RootFreq");
  }
  else if(freqOpt == "init")
  {
    vector<double> rootFreq(alphabet->getSize());
    if(AlphabetTools::isNucleicAlphabet(alphabet))
    {
      rootFreq[ 0] = ApplicationTools::getDoubleParameter("model.ancA", params, 0.25);
      rootFreq[ 1] = ApplicationTools::getDoubleParameter("model.ancC", params, 0.25);
      rootFreq[ 2] = ApplicationTools::getDoubleParameter("model.ancG", params, 0.25);
      rootFreq[ 3] = ApplicationTools::getDoubleParameter("model.ancT", params, 0.25);
      double theta  = (rootFreq[1] + rootFreq[2]) / (rootFreq[0] + rootFreq[1] + rootFreq[2] + rootFreq[3]);
      double theta1 = rootFreq[0] / (rootFreq[0] + rootFreq[3]);
      double theta2 = rootFreq[2] / (rootFreq[1] + rootFreq[2]);
      rootFrequencies = new FullNAFrequenciesSet(dynamic_cast<const NucleicAlphabet *>(alphabet), theta, theta1, theta2, "RootFreq");
    }
    else if(AlphabetTools::isProteicAlphabet(alphabet))
    {
      rootFreq[ 0] = ApplicationTools::getDoubleParameter("model.ancA", params, 0.05);
      rootFreq[ 1] = ApplicationTools::getDoubleParameter("model.ancR", params, 0.05);
      rootFreq[ 2] = ApplicationTools::getDoubleParameter("model.ancN", params, 0.05);
      rootFreq[ 3] = ApplicationTools::getDoubleParameter("model.ancD", params, 0.05);
      rootFreq[ 4] = ApplicationTools::getDoubleParameter("model.ancC", params, 0.05);
      rootFreq[ 5] = ApplicationTools::getDoubleParameter("model.ancQ", params, 0.05);
      rootFreq[ 6] = ApplicationTools::getDoubleParameter("model.ancE", params, 0.05);
      rootFreq[ 7] = ApplicationTools::getDoubleParameter("model.ancG", params, 0.05);
      rootFreq[ 8] = ApplicationTools::getDoubleParameter("model.ancH", params, 0.05);
      rootFreq[ 9] = ApplicationTools::getDoubleParameter("model.ancI", params, 0.05);
      rootFreq[10] = ApplicationTools::getDoubleParameter("model.ancL", params, 0.05);
      rootFreq[11] = ApplicationTools::getDoubleParameter("model.ancK", params, 0.05);
      rootFreq[12] = ApplicationTools::getDoubleParameter("model.ancM", params, 0.05);
      rootFreq[13] = ApplicationTools::getDoubleParameter("model.ancF", params, 0.05);
      rootFreq[14] = ApplicationTools::getDoubleParameter("model.ancP", params, 0.05);
      rootFreq[15] = ApplicationTools::getDoubleParameter("model.ancS", params, 0.05);
      rootFreq[16] = ApplicationTools::getDoubleParameter("model.ancT", params, 0.05);
      rootFreq[17] = ApplicationTools::getDoubleParameter("model.ancW", params, 0.05);
      rootFreq[18] = ApplicationTools::getDoubleParameter("model.ancY", params, 0.05);
      rootFreq[19] = ApplicationTools::getDoubleParameter("model.ancV", params, 0.05);
      rootFrequencies = new FullProteinFrequenciesSet(dynamic_cast<const ProteicAlphabet *>(alphabet), rootFreq, "RootFreq");
    }
    else throw Exception("Init root frequencies is only available with Nucleic and Proteic alphabet.");
  }
  else if(freqOpt == "observedGC" && data)
  {
    if(!AlphabetTools::isNucleicAlphabet(alphabet)) throw Exception("Error, unvalid option " + freqOpt + " with non-nucleic alphabet.");
    map<int, double> freqs = SequenceContainerTools::getFrequencies(*data);
    double theta  = (freqs[1] + freqs[2]) / (freqs[0] + freqs[1] + freqs[2] + freqs[3]);
    if(verbose) ApplicationTools::displayResult("Ancestral theta", TextTools::toString(theta));
    rootFrequencies = new GCFrequenciesSet(dynamic_cast<const NucleicAlphabet *>(alphabet), theta, "RootFreq");
  }
  else if(freqOpt == "balancedGC")
  {
    if(!AlphabetTools::isNucleicAlphabet(alphabet)) throw Exception("Error, unvalid option " + freqOpt + " with non-nucleic alphabet.");
    rootFrequencies = new GCFrequenciesSet(dynamic_cast<const NucleicAlphabet *>(alphabet), "RootFreq");
  }
  else if(freqOpt == "initGC")
  {
    if(!AlphabetTools::isNucleicAlphabet(alphabet)) throw Exception("Error, unvalid option " + freqOpt + " with non-nucleic alphabet.");
    double theta = ApplicationTools::getDoubleParameter("model.ancTheta", params, 0.5);
    rootFrequencies = new GCFrequenciesSet(dynamic_cast<const NucleicAlphabet *>(alphabet), theta, "RootFreq");
  }
  else throw Exception("Unvalid method for root frequencies: " + freqOpt);
    
  if(rateFreqs.size() > 0)
    rootFrequencies = new MarkovModulatedFrequenciesSet(rootFrequencies, rateFreqs);
  return rootFrequencies;
}

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
  SubstitutionModel * tmp = getSubstitutionModelDefaultInstance(alphabet, data, params, string("model1."), suffix, suffixIsOptional, verbose);
  if(tmp->getNumberOfStates() != alphabet->getSize())
  {
    //Markov-Modulated Markov Model...
    unsigned int n =(unsigned int)(tmp->getNumberOfStates() / alphabet->getSize());
    vector<double> ratesFreq(n, 1./(double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                                               // we should assume a rate distribution for the root also!!!  
  }
  delete tmp;
  
  FrequenciesSet * rootFrequencies = getFrequenciesSet(alphabet, data, params, rateFreqs, suffix, suffixIsOptional, verbose);

  SubstitutionModelSet * modelSet = new SubstitutionModelSet(alphabet, rootFrequencies);
  
  // Now parse all models:
  map<string, double> existingParameters;
  for(unsigned int i = 0; i < nbModels; i++)
  {
    string prefix = "model" + TextTools::toString(i+1) + ".";
    SubstitutionModel * model = getSubstitutionModelDefaultInstance(alphabet, data, params, prefix, suffix, suffixIsOptional, verbose);
    vector<string> specificParameters, sharedParameters;
    setSubstitutionModelParametersInitialValues(model, existingParameters, specificParameters, sharedParameters, params, prefix, suffix, suffixIsOptional, verbose);
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

void PhylogeneticsApplicationTools::printCovarionModelHelp()
{
  if(!ApplicationTools::message) return;
  *ApplicationTools::message << "Covarion Model:" << endl;
  *ApplicationTools::message << "model                   | [+G2001|+TS98]" << endl;
  *ApplicationTools::message << "model.nu                | covarion rate parameter in G2001 model." << endl;
  *ApplicationTools::message << "model.rate_distribution | [gamma] covarion rate distribution in G2001 model." << endl;
  *ApplicationTools::message << "model.rate_distribution |" << endl;
  *ApplicationTools::message << "                 .alpha | alpha parameter of the gamma law." << endl;
  *ApplicationTools::message << "model.rate_distribution |" << endl;
  *ApplicationTools::message << "         .classes_number| number of classes for rate distribution." << endl;
  *ApplicationTools::message << "model.nu                | covarion rate parameter in G2001 model." << endl;
  *ApplicationTools::message << "model.nu                | covarion rate parameter in G2001 model." << endl;
  *ApplicationTools::message << "model.s1                | covarion rate parameter in TS98 model." << endl;
  *ApplicationTools::message << "model.s2                | covarion rate parameter in TS98 model." << endl;
  *ApplicationTools::message << "________________________|________________________________________________" << endl;
}

/******************************************************************************/

DiscreteDistribution * PhylogeneticsApplicationTools::getRateDistributionDefaultInstance(
  map<string, string> & params,
  bool constDistAllowed,
  const string & prefix,
  const string & suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  string defaultDist = constDistAllowed ? "constant" : "gamma";
  string distributionType = ApplicationTools::getStringParameter(prefix + "rate_distribution", params, defaultDist, suffix, suffixIsOptional);
  string mixedType = "";
  DiscreteDistribution * rDist;
  string::size_type x = distributionType.find("+");
  if(x != string::npos)
  {
    mixedType = distributionType.substr(x+1);
    distributionType = distributionType.substr(0, x);
  }
  if(distributionType == "constant")
  {
    if(!constDistAllowed) throw Exception("You can't use a constant distribution here!");
    rDist = new ConstantDistribution(1.);
  }
  else if(distributionType == "gamma")
  {
    int nbClasses = ApplicationTools::getIntParameter(prefix + "rate_distribution.classes_number", params, 4, suffix, suffixIsOptional);
    rDist = new GammaDiscreteDistribution(nbClasses);
  }
  else
  {
    throw Exception("Distribution unknown: " + distributionType + ".");
  }
  if(mixedType == "invariant")
  {
    rDist = new InvariantMixedDiscreteDistribution(rDist, 0., 0.000001, true);
  }
  if(verbose)
  {
    ApplicationTools::displayResult("Rate distribution", distributionType + (mixedType != "" ? "+" + mixedType : ""));
    ApplicationTools::displayResult("Number of classes", TextTools::toString(rDist->getNumberOfCategories()));
  }

  return rDist;
}

void PhylogeneticsApplicationTools::setRateDistributionParametersInitialValues(
  DiscreteDistribution * rDist,
  map<string, string> & params,
  const string & prefix,
  const string & suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  ParameterList pl = rDist->getParameters();
  for(unsigned int i = 0; i < pl.size(); i++)
  {
    const string pName = pl[i]->getName();
    double value = ApplicationTools::getDoubleParameter(prefix + "rate_distribution." + pName, params, pl[i]->getValue(), suffix, suffixIsOptional); 
    if(pName == "alpha" && value < 0)
    {
      throw Exception("Alpha parameter in gamma distribution of rates must be > 0, found " + TextTools::toString(value) + ".");
    }
    pl[i]->setValue(value);
    if(verbose) ApplicationTools::displayResult(prefix + "rate_distribution." + pName, TextTools::toString(pl[i]->getValue()));
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
 

DiscreteDistribution * PhylogeneticsApplicationTools::getRateDistribution(
  map<string, string> & params,
  const string & suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  DiscreteDistribution * rDist = getRateDistributionDefaultInstance(params, true, "", suffix, suffixIsOptional, verbose);
  setRateDistributionParametersInitialValues(rDist, params, "", suffix, suffixIsOptional, verbose);
  return rDist;
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printRateDistributionHelp()
{
  if(!ApplicationTools::message) return;
  *ApplicationTools::message << "Rate distribution parameters:" << endl;
  *ApplicationTools::message << "rate_distribution               | constant or gamma, add '+invariant' for invariants." << endl;
  *ApplicationTools::message << "rate_distribution.classes_number| discrete approximation: number of" << endl;
  *ApplicationTools::message << "                                | categories (default to 4)." << endl;
  *ApplicationTools::message << "rate_distribution.alpha         | the gamma law's alpha parameter." << endl;
  *ApplicationTools::message << "rate_distribution.p             | the proportion of invariant (default to 0)." << endl;
  *ApplicationTools::message << "________________________________|_________________________________________" << endl;
}

/******************************************************************************/

TreeLikelihood * PhylogeneticsApplicationTools::optimizeParameters(
  TreeLikelihood * tl,
  map<string, string> & params,
  const string & suffix,
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
        for(unsigned int i = 0; i < vs.size(); i++)
        {
          dynamic_cast<AbstractTreeLikelihood *>(tl)->ignoreParameter(vs[i]);
        }
      }
      else if(param == "Ancient")
      {
        vector<string> vs = dynamic_cast<NonHomogeneousTreeLikelihood *>(tl)->getRootFrequenciesParameters().getParameterNames();
        for(unsigned int i = 0; i < vs.size(); i++)
        {
          dynamic_cast<AbstractTreeLikelihood *>(tl)->ignoreParameter(vs[i]);
        }
      }
      else dynamic_cast<AbstractTreeLikelihood *>(tl)->ignoreParameter(param);
    } 
    catch(ParameterNotFoundException pnfe)
    {
      ApplicationTools::displayError("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
    }
    catch(exception e)
    {
      ApplicationTools::displayError("DEBUB: ERROR!!! This functionality can only be used with HomogeneousTreeLikelihood for now.");
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
          dynamic_cast<NNIHomogeneousTreeLikelihood *>(tl),
          optNumFirst, tolBefore, tolDuring, nbEvalMax, n, messageHandler, profiler, optVerbose, optMethod, nstep, nniAlgo);
    }

    if(verbose && nstep > 1) ApplicationTools::displayResult("# of precision steps", TextTools::toString(nstep));
    n = OptimizationTools::optimizeNumericalParameters(
      dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl),
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
          dynamic_cast<NNIHomogeneousTreeLikelihood *>(tl),
          optNumFirst, tolBefore, tolDuring, nbEvalMax, n, messageHandler, profiler, optVerbose, optMethod, nniAlgo);
    }

    n = OptimizationTools::optimizeNumericalParameters2(
      dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(tl), NULL, tolerance, nbEvalMax, messageHandler, profiler, optVerbose, optMethod);       
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
    if(verbose) ApplicationTools::displayResult("Final optimization step", finalMethod);
    finalOptimizer->setProfiler(profiler);
    finalOptimizer->setMessageHandler(messageHandler);
    finalOptimizer->setMaximumNumberOfEvaluations(nbEvalMax);
    finalOptimizer->getStopCondition()->setTolerance(tolerance);
    finalOptimizer->setVerbose(verbose);
    finalOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    finalOptimizer->init(tl->getParameters());
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
        for(unsigned int i = 0; i < vs.size(); i++)
        {
          dynamic_cast<AbstractTreeLikelihood *>(tl)->ignoreParameter(vs[i]);
        }
      }
      else if(param == "Ancient")
      {
        vector<string> vs = dynamic_cast<NonHomogeneousTreeLikelihood *>(tl)->getRootFrequenciesParameters().getParameterNames();
        for(unsigned int i = 0; i < vs.size(); i++)
        {
          dynamic_cast<AbstractTreeLikelihood *>(tl)->ignoreParameter(vs[i]);
        }
      }
      else dynamic_cast<AbstractTreeLikelihood *>(tl)->ignoreParameter(param);
    } 
    catch(ParameterNotFoundException pnfe)
    {
      ApplicationTools::displayError("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
    }
    catch(exception e)
    {
      ApplicationTools::displayError("DEBUB: ERROR!!! This functionality can only be used with HomogeneousTreeLikelihood for now.");
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
    ApplicationTools::displayResult("Final optimization step", finalMethod);
    finalOptimizer->setProfiler(profiler);
    finalOptimizer->setMessageHandler(messageHandler);
    finalOptimizer->setMaximumNumberOfEvaluations(nbEvalMax);
    finalOptimizer->getStopCondition()->setTolerance(tolerance);
    finalOptimizer->setVerbose(verbose);
    finalOptimizer->setConstraintPolicy(AutoParameter::CONSTRAINTS_AUTO);
    finalOptimizer->init(tl->getParameters());
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

