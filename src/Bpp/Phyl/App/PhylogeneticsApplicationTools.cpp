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
#include "../Model.all"
#include "../Likelihood.all"
#include "../OptimizationTools.h"
#include "../Tree.h"
#include "../Io/Newick.h"
#include "../Io/NexusIoTree.h"
#include "../Io/Nhx.h"

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

using namespace std;

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
  else throw Exception("Unknow format for tree reading: " + format);
  Tree* tree = treeReader->read(treeFilePath);
  delete treeReader;

  if (verbose) ApplicationTools::displayResult("Tree file", treeFilePath);
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
  else throw Exception("Unknow format for tree reading: " + format);
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

   bool word = ((modelName == "Word") || (modelName == "Triplet") || (modelName == "CodonNeutral")
                || (modelName == "CodonAsynonymous"));

   bool wordfreq = ((modelName == "CodonAsynonymousFrequencies")
                    || (modelName == "CodonNeutralFrequencies"));

  // //////////////////////////////////
  // / MIXED MODELS
  // ////////////////////////////////

  if (modelName == "MixedModel")
  {
    map<string, string> unparsedParameterValuesNested;
    if (args.find("model") == args.end())
      throw Exception("The argument 'model' is missing from MixedSubstitutionModel description");
    string nestedModelDescription = args["model"];
    SubstitutionModel* pSM = getSubstitutionModelDefaultInstance(alphabet,
                                                                 nestedModelDescription,
                                                                 unparsedParameterValuesNested,
                                                                 allowCovarions,
                                                                 allowGaps,
                                                                 verbose);
    map<string,DiscreteDistribution*> mdist;
    map<string, string> unparsedParameterValuesNested2, unparsedParameterValuesNested3;

    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin();
         it != unparsedParameterValuesNested.end(); it++)
    {
      if (it->second.find("(") != string::npos)
      {
        unparsedParameterValuesNested3.clear();

        mdist[pSM->getParameterNameWithoutNamespace(it->first)] = getDistributionDefaultInstance(it->second, unparsedParameterValuesNested3,true);
        for (map<string, string>::iterator it2 = unparsedParameterValuesNested3.begin();
             it2 != unparsedParameterValuesNested3.end(); it2++)
        {
          unparsedParameterValuesNested2[it->first + "_" + it2->first] = it2->second;
        }
      }
      else
        unparsedParameterValuesNested2[it->first] = it->second;
    }

    for (map<string, string>::iterator it = unparsedParameterValuesNested2.begin();
         it != unparsedParameterValuesNested2.end(); it++)
    {
      unparsedParameterValues[it->first] = it->second;
    }

    model = new MixtureOfASubstitutionModel(alphabet,pSM,mdist);

    vector<string> v = model->getParameters().getParameterNames();

    for (map<string,DiscreteDistribution*>::iterator it = mdist.begin();
         it != mdist.end(); it++)
    {
      delete it->second;
    }

    if (verbose)
      ApplicationTools::displayResult("Mixture Of A Substitution Model", nestedModelDescription );
  }
  else if (modelName == "Mixture")
    {
      vector<string> v_nestedModelDescription;
      vector<SubstitutionModel*> v_pSM;
      map<string, string> unparsedParameterValuesNested;

      if (args.find("model1") == args.end()) {
        throw Exception("Missing argument 'model1' for model " + modelName + ".");
      }
      unsigned int nbmodels = 0;
      
      while (args.find("model" + TextTools::toString(nbmodels+1)) != args.end()){
        v_nestedModelDescription.push_back(args["model" + TextTools::toString(++nbmodels)]);
      }
    
      if (nbmodels < 2)
        throw Exception("Missing nested models for model " + modelName + ".");

      for (unsigned i = 0; i < v_nestedModelDescription.size(); i++)
        {
          unparsedParameterValuesNested.clear();
          model = getSubstitutionModelDefaultInstance(alphabet, v_nestedModelDescription[i], unparsedParameterValuesNested, false, false, false);
          for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++){
            unparsedParameterValues[modelName + "." + TextTools::toString(i+1) + "_" + it->first] = it->second;
          }
          v_pSM.push_back(model);
        }

      model = new MixtureOfSubstitutionModels(alphabet, v_pSM);
      if (verbose)
        ApplicationTools::displayResult("Mixture Of Substitution Models", modelName );
    }
  
  // /////////////////////////////////
  // / WORDS and CODONS defined by models
  // ///////////////////////////////

  else if (word)
  {
    vector<string> v_nestedModelDescription;
    vector<SubstitutionModel*> v_pSM;
    const WordAlphabet* pWA;

    string s, nestedModelDescription;
    unsigned int nbmodels;

    if ((modelName == "Word" && !AlphabetTools::isWordAlphabet(alphabet)) ||
        (modelName != "Word" && !AlphabetTools::isCodonAlphabet(alphabet)))
      throw Exception("Bad alphabet type "
                      + alphabet->getAlphabetType() + " for  model " + modelName + ".");

    pWA = dynamic_cast<const WordAlphabet*>(alphabet);

    if (args.find("model") != args.end())
    {
      nestedModelDescription = args["model"];
      if (modelName == "Word")
      {
        v_nestedModelDescription.push_back(nestedModelDescription);
        nbmodels = pWA->getLength();
      }
      else
      {
        v_nestedModelDescription.push_back(nestedModelDescription);
        nbmodels = 3;
      }
    }
    else
    {
      if (args.find("model1") == args.end())
      {
        throw Exception("Missing argument 'model' or 'model1' for model " + modelName + ".");
      }
      nbmodels = 0;

      while (args.find("model" + TextTools::toString(nbmodels+1)) != args.end())
      {
        v_nestedModelDescription.push_back(args["model" + TextTools::toString(++nbmodels)]);
      }
    }

    if (nbmodels < 2)
      throw Exception("Missing nested models for model " + modelName + ".");

    if (pWA->getLength() != nbmodels)
      throw Exception("Bad alphabet type "
                      + alphabet->getAlphabetType() + " for  model " + modelName + ".");

    map<string, string> unparsedParameterValuesNested;

    if (v_nestedModelDescription.size() != nbmodels)
    {
      model = getSubstitutionModelDefaultInstance(pWA->getNAlphabet(0), v_nestedModelDescription[0], unparsedParameterValuesNested, false, false, false);
      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedParameterValues[modelName + "._" + it->first] = it->second;
      }
      v_pSM.push_back(model);
    }
    else
    {
      for (unsigned i = 0; i < v_nestedModelDescription.size(); i++)
      {
        unparsedParameterValuesNested.clear();
        model = getSubstitutionModelDefaultInstance(pWA->getNAlphabet(i), v_nestedModelDescription[i], unparsedParameterValuesNested, false, false, false);
        for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
        {
          unparsedParameterValues[modelName + "." + TextTools::toString(i+1) + "_" + it->first] = it->second;
        }
        v_pSM.push_back(model);
      }
    }

    // /////////////////////////////////
    // / WORD
    // ///////////////////////////////

    if (modelName == "Word")
    {
      model = (v_nestedModelDescription.size() != nbmodels)
              ? new WordReversibleSubstitutionModel(v_pSM[0], nbmodels)
              : new WordReversibleSubstitutionModel(v_pSM);
    }

    // /////////////////////////////////
    // / TRIPLET
    // ///////////////////////////////

    else if (modelName == "Triplet")
    {
      if (dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0])==NULL)
        throw Exception("Non simple NucleotideSubstitutionModel imbedded in " + modelName + " model.");
      
      if (v_nestedModelDescription.size() != 3)
        model = new TripletReversibleSubstitutionModel(
                                                       dynamic_cast<const CodonAlphabet*>(pWA),
                                                       dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]));
      else {
        if (dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1])==NULL || dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2])==NULL)
        throw Exception("Non simple NucleotideSubstitutionModel imbedded in " + modelName + " model.");

        model = new TripletReversibleSubstitutionModel(
                                                       dynamic_cast<const CodonAlphabet*>(pWA),
                                                       dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                       dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                                                       dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]));
      }
    }

    // /////////////////////////////////
    // / CODON NEUTRAL
    // ///////////////////////////////

    else if (modelName == "CodonNeutral")
    {
      if (dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0])==NULL)
        throw Exception("Non simple NucleotideSubstitutionModel imbedded in " + modelName + " model.");
      
      if (v_nestedModelDescription.size() != 3)
        model = new CodonNeutralReversibleSubstitutionModel(
                                                       dynamic_cast<const CodonAlphabet*>(pWA),
                                                       dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]));
      else {
        if (dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1])==NULL || dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2])==NULL)
        throw Exception("Non simple NucleotideSubstitutionModel imbedded in " + modelName + " model.");

        model = new CodonNeutralReversibleSubstitutionModel(
                                                       dynamic_cast<const CodonAlphabet*>(pWA),
                                                       dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                       dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                                                       dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]));
      }
    }

    // /////////////////////////////////
    // / CODON ASYNONYMOUS
    // ///////////////////////////////

    else if (modelName == "CodonAsynonymous")
    {
      if (args.find("genetic_code") == args.end())
        args["genetic_code"]=pWA->getAlphabetType();
      
      GeneticCode* pgc = SequenceApplicationTools::getGeneticCode(dynamic_cast<const NucleicAlphabet*>(pWA->getNAlphabet(0)),args["genetic_code"]);
      if (pgc->getSourceAlphabet()->getAlphabetType() != pWA->getAlphabetType())
        throw Exception("Mismatch between genetic code and codon alphabet");

      AlphabetIndex2<double>* pai2;

      if (args.find("aadistance") == args.end())
        pai2 = 0;
      else
        pai2 = SequenceApplicationTools::getAADistance(args["aadistance"]);

      if (dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0])==NULL)
        throw Exception("Non simple NucleotideSubstitutionModel imbedded in " + modelName + " model.");
      
      if (v_nestedModelDescription.size() != 3)
        model = new CodonAsynonymousReversibleSubstitutionModel(pgc,
                                                                dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]), pai2);
      else {
        if (dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1])==NULL || dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2])==NULL)
        throw Exception("Non simple NucleotideSubstitutionModel imbedded in " + modelName + " model.");

        model = new CodonAsynonymousReversibleSubstitutionModel(
                                                                pgc,
                                                                dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                                dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                                                                dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]), pai2);
      }
    }
  }

  // /////////////////////////////////
  // / CODON MODELS with FREQUENCIES
  // ///////////////////////////////

  else if (wordfreq)
  {
    if (!AlphabetTools::isCodonAlphabet(alphabet))
      throw Exception("Alphabet should be Codon Alphabet.");

    const CodonAlphabet* pCA = dynamic_cast<const CodonAlphabet*>(alphabet);

    if (args.find("frequencies") == args.end())
      throw Exception("Missing equilibrium frequencies.");

    map<string, string> unparsedParameterValuesNested;

    FrequenciesSet* pFS = getFrequenciesSetDefaultInstance(pCA, args["frequencies"], unparsedParameterValuesNested);

    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedParameterValues[modelName + "." + it->first] = it->second;
    }

    // /////////////////////////////////
    // / CODONNEUTRALFREQUENCIES
    // ///////////////////////////////

    if (modelName == "CodonNeutralFrequencies")
    {
      model = new CodonNeutralFrequenciesReversibleSubstitutionModel(pCA, pFS);

      // for description
      modelName += args["frequencies"];
    }

    // /////////////////////////////////
    // / CODONASYNONYMOUSFREQUENCIES
    // ///////////////////////////////

    else if (modelName == "CodonAsynonymousFrequencies")
    {
      if (args.find("genetic_code") == args.end())
        args["genetic_code"] = pCA->getAlphabetType();

      GeneticCode* pgc = SequenceApplicationTools::getGeneticCode(dynamic_cast<const NucleicAlphabet*>(pCA->getNAlphabet(0)),args["genetic_code"]);
      if (pgc->getSourceAlphabet()->getAlphabetType() != pCA->getAlphabetType())
        throw Exception("Mismatch between genetic code and codon alphabet");

      AlphabetIndex2<double>* pai2;

      if (args.find("aadistance") == args.end())
        pai2 = 0;
      else
        pai2  = SequenceApplicationTools::getAADistance(args["aadistance"]);

      model = new CodonAsynonymousFrequenciesReversibleSubstitutionModel(pgc, pFS, pai2);
    }

  }

  // //////////////////////////////////////
  // predefined codon models
  // //////////////////////////////////////

  else if ((modelName == "MG94") || (modelName == "YN98") ||
           (modelName == "GY94") || (modelName.substr(0,5)=="YNGKP"))
  {
    if (!AlphabetTools::isCodonAlphabet(alphabet))
      throw Exception("Alphabet should be Codon Alphabet.");

    const CodonAlphabet* pCA = (const CodonAlphabet*)(alphabet);

    if (args.find("genetic_code") == args.end())
      args["genetic_code"]=pCA->getAlphabetType();

    GeneticCode* pgc = SequenceApplicationTools::getGeneticCode(dynamic_cast<const NucleicAlphabet*>(pCA->getNAlphabet(0)),args["genetic_code"]);
    if (pgc->getSourceAlphabet()->getAlphabetType() != pCA->getAlphabetType())
      throw Exception("Mismatch between genetic code and codon alphabet");


    FrequenciesSet* codonFreqs;
    
    if (args.find("frequencies") != args.end()){

      map<string, string> unparsedParameterValuesNested;
    
      codonFreqs = getFrequenciesSetDefaultInstance(pCA, args["frequencies"], unparsedParameterValuesNested);
      

    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
        unparsedParameterValues[modelName + "." + it->first] = it->second;
    }
    else {
      string freqOpt = ApplicationTools::getStringParameter("codon_freqs", args, "F0");
      short opt = 0;
      if (freqOpt == "F0")
        opt = FrequenciesSet::F0;
      else if (freqOpt == "F1X4")
        opt = FrequenciesSet::F1X4;
      else if (freqOpt == "F3X4")
        opt = FrequenciesSet::F3X4;
      else if (freqOpt == "F61")
        opt = FrequenciesSet::F61;
      else
        throw Exception("Unvalid codon frequency option. Should be one of F0, F1X4, F3X4 or F61");
      
      codonFreqs = FrequenciesSet::getFrequenciesSetForCodons(opt, *pgc);
    }
    
    if (modelName == "MG94")
      model = new MG94(pgc, codonFreqs);
    else if (modelName == "GY94")
      model = new GY94(pgc, codonFreqs);
    else if ((modelName == "YN98") || (modelName == "YNGKP_M0"))
      model = new YN98(pgc, codonFreqs);
    else if (modelName == "YNGKP_M1")
      model = new YNGKP_M1(pgc, codonFreqs);
     else if (modelName == "YNGKP_M2")
      model = new YNGKP_M2(pgc, codonFreqs);
     else if (modelName == "YNGKP_M3")
       if (args.find("n") == args.end())
         model = new YNGKP_M3(pgc, codonFreqs);
       else
         model = new YNGKP_M3(pgc, codonFreqs,TextTools::to<unsigned int>(args["n"]));
     else if ((modelName == "YNGKP_M7") || modelName == "YNGKP_M8"){
       if (args.find("n") == args.end())
         throw Exception("Missing argument 'n' (number of classes) in " + modelName + " distribution");
       unsigned int nbClasses = TextTools::to<unsigned int>(args["n"]);

       if (modelName == "YNGKP_M7")
         model = new YNGKP_M7(pgc, codonFreqs,nbClasses);
       else if (modelName == "YNGKP_M8")
         model = new YNGKP_M8(pgc, codonFreqs,nbClasses);
     }
    else
      throw Exception("Unknown Codon model: " + modelName);

  }

  // /////////////////////////////////
  // / RE08
  // ///////////////////////////////

  else if (modelName == "RE08")
  {
    if (!allowGaps)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance. No Gap model allowed here.");

    // We have to parse the nested model first:
    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance. Missing argument 'model' for model 'RE08'.");
    if (verbose)
      ApplicationTools::displayResult("Gap model", modelName);
    map<string, string> unparsedParameterValuesNested;
    SubstitutionModel* nestedModel = getSubstitutionModelDefaultInstance(alphabet, nestedModelDescription, unparsedParameterValuesNested, allowCovarions, false, verbose);

    // Now we create the RE08 substitution model:
    ReversibleSubstitutionModel* tmp = dynamic_cast<ReversibleSubstitutionModel*>(nestedModel);
    model = new RE08(tmp);

    // Then we update the parameter set:
    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedParameterValues["RE08.model_" + it->first] = it->second;
    }
  }

  // /////////////////////////////////
  // / TS98
  // ///////////////////////////////

  else if (modelName == "TS98")
  {
    if (!allowCovarions)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance. No Covarion model allowed here.");

    // We have to parse the nested model first:
    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance. Missing argument 'model' for model 'TS98'.");
    if (verbose)
      ApplicationTools::displayResult("Covarion model", modelName);
    map<string, string> unparsedParameterValuesNested;
    SubstitutionModel* nestedModel = getSubstitutionModelDefaultInstance(alphabet, nestedModelDescription, unparsedParameterValuesNested, false, allowGaps, verbose);

    // Now we create the TS98 substitution model:
    ReversibleSubstitutionModel* tmp = dynamic_cast<ReversibleSubstitutionModel*>(nestedModel);
    model = new TS98(tmp);

    // Then we update the parameter set:
    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedParameterValues["TS98.model_" + it->first] = it->second;
    }
  }

  // /////////////////////////////////
  // / G01
  // ///////////////////////////////

  else if (modelName == "G01")
  {
    if (!allowCovarions)
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance. No Covarion model allowed here.");

    // We have to parse the nested model first:
    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance. Missing argument 'model' for model 'G01'.");
    string nestedRateDistDescription = args["rdist"];
    if (TextTools::isEmpty(nestedRateDistDescription))
      throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelDefaultInstance. Missing argument 'rdist' for model 'G01'.");
    if (verbose)
      ApplicationTools::displayResult("Covarion model", modelName);

    map<string, string> unparsedParameterValuesNestedModel;
    SubstitutionModel* nestedModel = getSubstitutionModelDefaultInstance(alphabet, nestedModelDescription, unparsedParameterValuesNestedModel, false, allowGaps, verbose);
    map<string, string> unparsedParameterValuesNestedDist;
    DiscreteDistribution* nestedRDist = getRateDistributionDefaultInstance(nestedRateDistDescription, unparsedParameterValuesNestedDist, false, verbose);

    // Now we create the G01 substitution model:
    ReversibleSubstitutionModel* tmp = dynamic_cast<ReversibleSubstitutionModel*>(nestedModel);
    model = new G2001(tmp, nestedRDist);

    // Then we update the parameter set:
    for (map<string, string>::iterator it = unparsedParameterValuesNestedModel.begin(); it != unparsedParameterValuesNestedModel.end(); it++)
    {
      unparsedParameterValues["G01.model_" + it->first] = it->second;
    }
    for (map<string, string>::iterator it = unparsedParameterValuesNestedDist.begin(); it != unparsedParameterValuesNestedDist.end(); it++)
    {
      unparsedParameterValues["G01.rdist_" + it->first] = it->second;
    }
  }
  else
  {
    // This is a 'simple' model...
    if (AlphabetTools::isNucleicAlphabet(alphabet))
    {
      const NucleicAlphabet* alpha = dynamic_cast<const NucleicAlphabet*>(alphabet);

      // /////////////////////////////////
      // / GTR
      // ///////////////////////////////

      if (modelName == "GTR")
      {
        model = new GTR(alpha);
      }


      // /////////////////////////////////
      // / L95
      // ///////////////////////////////

      else if (modelName == "SSR")
        {
          model = new SSR(alpha);
        }

      // /////////////////////////////////
      // / L95
      // ///////////////////////////////

      else if (modelName == "L95")
      {
        model = new L95(alpha);
      }

      // /////////////////////////////////
      // / TN93
      // //////////////////////////////

      else if (modelName == "TN93")
      {
        model = new TN93(alpha);
      }

      // /////////////////////////////////
      // / HKY85
      // ///////////////////////////////

      else if (modelName == "HKY85")
      {
        model = new HKY85(alpha);
      }

      // /////////////////////////////////
      // / F84
      // ///////////////////////////////

      else if (modelName == "F84")
      {
        model = new F84(alpha);
      }

      // /////////////////////////////////
      // / T92
      // ///////////////////////////////

      else if (modelName == "T92")
      {
        model = new T92(alpha);
      }

      // /////////////////////////////////
      // / K80
      // ///////////////////////////////

      else if (modelName == "K80")
      {
        model = new K80(alpha);
      }


      // /////////////////////////////////
      // / JC69
      // ///////////////////////////////

      else if (modelName == "JC69")
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
      const ProteicAlphabet* alpha = dynamic_cast<const ProteicAlphabet*>(alphabet);

      if (modelName == "JC69+F")
        model = new JCprot(alpha, new FullProteinFrequenciesSet(alpha), true);
      else if (modelName == "DSO78+F")
        model = new DSO78(alpha, new FullProteinFrequenciesSet(alpha), true);
      else if (modelName == "JTT92+F")
        model = new JTT92(alpha, new FullProteinFrequenciesSet(alpha), true);
      else if (modelName == "LG08+F")
        model = new LG08(alpha, new FullProteinFrequenciesSet(alpha), true);
      else if (modelName == "WAG01+F")
        model = new WAG01(alpha, new FullProteinFrequenciesSet(alpha), true);
      else if (modelName == "Empirical+F")
      {
        string prefix = args["name"];
        if (TextTools::isEmpty(prefix))
          throw Exception("'name' argument missing for user-defined substitution model.");
        model = new UserProteinSubstitutionModel(alpha, args["file"], new FullProteinFrequenciesSet(alpha), prefix + "+F.", true);
      }
      else if (modelName == "JC69")
        model = new JCprot(alpha);
      else if (modelName == "DSO78")
        model = new DSO78(alpha);
      else if (modelName == "JTT92")
        model = new JTT92(alpha);
      else if (modelName == "LG08")
        model = new LG08(alpha);
      else if (modelName == "WAG01")
        model = new WAG01(alpha);
      else if (modelName == "LLG08_EHO")
        model = new LLG08_EHO(alpha);
      else if (modelName == "LLG08_EX2")
        model = new LLG08_EX2(alpha);
      else if (modelName == "LLG08_EX3")
        model = new LLG08_EX3(alpha);
      else if (modelName == "LLG08_UL2")
        model = new LLG08_UL2(alpha);
      else if (modelName == "LLG08_UL3")
        model = new LLG08_UL3(alpha);
      else if (modelName == "Empirical")
      {
        string prefix = args["name"];
        if (TextTools::isEmpty(prefix))
          throw Exception("'name' argument missing for user-defined substitution model.");
        model = new UserProteinSubstitutionModel(alpha, args["file"], prefix);
      }
      else if (modelName == "JC69+COA")
        model = new COA(alpha, "JC69");
      else if (modelName == "DSO78+COA")
        model = new COA(alpha, "DSO78");
      else if (modelName == "JTT92+COA")
        model = new COA(alpha, "JTT92");
      else if (modelName == "LG08+COA")
        model = new COA(alpha, "LG08");
      else if (modelName == "WAG01+COA")
        model = new COA(alpha, "WAG01");
      else if (modelName == "Empirical+COA")
        {
        string prefix = args["name"];
        if (TextTools::isEmpty(prefix))
          throw Exception("'name' argument missing for user-defined substitution model.");
        model = new COA(alpha, args["file"]);
        }
      else
        throw Exception("Model '" + modelName + "' unknown.");
    }
    if (verbose)
      ApplicationTools::displayResult("Substitution model", modelName);
  }

  //Update parameter args:
  vector<string> pnames = model->getParameters().getParameterNames();

  for (unsigned int i = 0; i < pnames.size(); i++)
  {
    string name = model->getParameterNameWithoutNamespace(pnames[i]);
    if (args.find(name) != args.end()) {
      unparsedParameterValues[modelName + "." + name] = args[name];
    }
  }

  // Now look if some parameters are aliased:
  ParameterList pl = model->getIndependentParameters();
  string pname, pval, pname2;
  for (unsigned int i = 0; i < pl.size(); i++)
  {
    pname = model->getParameterNameWithoutNamespace(pl[i].getName());

    if (args.find(pname) == args.end()) continue;
    pval = args[pname];

    if ((pval.length() >= 5 && pval.substr(0, 5) == "model") ||
        (pval.find("(") != string::npos))
      continue;
    bool found = false;
    for (unsigned int j = 0; j < pl.size() && !found; j++)
    {
      pname2 = model->getParameterNameWithoutNamespace(pl[j].getName());

      //if (j == i || args.find(pname2) == args.end()) continue; Julien 03/03/2010: This extra condition prevents complicated (nested) models to work properly...
      if (j == i) continue;
      if (pval == pname2)
      {
        // This is an alias...
        // NB: this may throw an exception if uncorrect! We leave it as is for now :s
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
    modelDescription = ApplicationTools::getStringParameter("model", params, "CodonNeutral(model=JC69)", suffix, suffixIsOptional, verbose);
  else if (AlphabetTools::isWordAlphabet(alphabet))
    modelDescription = ApplicationTools::getStringParameter("model", params, "Word(model=JC69)", suffix, suffixIsOptional, verbose);
  else
    modelDescription = ApplicationTools::getStringParameter("model", params, "JC69", suffix, suffixIsOptional, verbose);
  
  map<string, string> unparsedParameterValues;
  SubstitutionModel* model = getSubstitutionModelDefaultInstance(alphabet, modelDescription, unparsedParameterValues, true, true, verbose);
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
  
  bool useObsFreq = ApplicationTools::getBooleanParameter(model->getNamespace() + "useObservedFreqs", unparsedParameterValues, false, "", true, false);
  if (verbose) ApplicationTools::displayResult("Use observed frequencies for model", useObsFreq ? "yes" : "no");
  if (useObsFreq && data != 0)
  {
   unsigned int psi = ApplicationTools::getParameter<unsigned int>(model->getNamespace() + "useObservedFreqs.pseudoCount", unparsedParameterValues, 0);
   model->setFreqFromData(*data, psi);
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
   posp=pName.rfind(".");
   if (!useObsFreq || (pName.substr(posp+1,5) != "theta"))
     {
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
  bool useObsFreq = ApplicationTools::getBooleanParameter(model->getNamespace() + "useObservedFreqs", unparsedParameterValues, false, "","", false);
  if (verbose) ApplicationTools::displayResult("Use observed frequencies for model", useObsFreq ? "yes" : "no");
  if (useObsFreq && data != 0)
  {
   unsigned int psi = ApplicationTools::getParameter<unsigned int>(model->getNamespace() + "useObservedFreqs.pseudoCount", unparsedParameterValues, 0);
   model->setFreqFromData(*data, psi);
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
   string value;
    if (!useObsFreq || (model->getParameterNameWithoutNamespace(pName).substr(0,5) != "theta"))
    {
      value = ApplicationTools::getStringParameter(pName, unparsedParameterValues, TextTools::toString(pl[i].getValue()));
      if (value.size() > 5 && value.substr(0, 5) == "model")
      {
        if (existingParams.find(value) != existingParams.end())
        {
          pl[i].setValue(existingParams[value]);
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
      unsigned int psc=0;
      if (unparsedParameterValues.find("pseudoCount") != unparsedParameterValues.end())
        psc=TextTools::toInt(unparsedParameterValues["pseudoCount"]);
      
      map<int, double> freqs;
      SequenceContainerTools::getFrequencies(*data, freqs, psc);

      pFS->setFrequenciesFromMap(freqs);
    }
    else if (init == "balanced")
    {
      // Nothing to do here, this is the default instanciation.
    }
    else
      throw Exception("Unknown init argument");
  }
  else if (unparsedParameterValues.find("values") != unparsedParameterValues.end())
  {
    // Initialization using the "values" argument
    vector<double> frequencies;
    string rf = unparsedParameterValues["values"];
    StringTokenizer strtok(rf.substr(1, rf.length() - 2), ",");
    while (strtok.hasMoreToken())
      frequencies.push_back(TextTools::toDouble(strtok.nextToken()));
    pFS->setFrequencies(frequencies);
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
   string freqName;
   map<string, string> args;
   KeyvalTools::parseProcedure(freqDescription, freqName, args);
   FrequenciesSet* pFS;

  if (freqName.substr(0, 4) == "Full")
  {
    if (AlphabetTools::isNucleicAlphabet(alphabet))
    {
      pFS = new FullNucleotideFrequenciesSet(dynamic_cast<const NucleicAlphabet*>(alphabet));
    }
    else if (AlphabetTools::isProteicAlphabet(alphabet))
    {
      pFS = new FullProteinFrequenciesSet(dynamic_cast<const ProteicAlphabet*>(alphabet));
    }
    else if (AlphabetTools::isCodonAlphabet(alphabet))
    {
      pFS = new FullCodonFrequenciesSet(dynamic_cast<const CodonAlphabet*>(alphabet));
    }
    else
    {
      pFS = new FullFrequenciesSet(alphabet);
    }
  
    // Update parameter values:
    if (args.find("theta") != args.end())
      unparsedParameterValues["Full.theta"] = args["theta"];
    for (unsigned int i = 1; i < alphabet->getSize() ; i++)
    {
      if (args.find("theta" + TextTools::toString(i) ) != args.end())
        {
        unparsedParameterValues["Full.theta" + TextTools::toString(i) ] = args["theta" + TextTools::toString(i) ];
        }
    }
  }
  else if (freqName == "GC")
  {
    if (!AlphabetTools::isNucleicAlphabet(alphabet))
      throw Exception("Error, unvalid frequencies " + freqName + " with non-nucleic alphabet.");

    pFS = new GCFrequenciesSet(dynamic_cast<const NucleicAlphabet*>(alphabet));

    if (args.find("theta") != args.end())
      unparsedParameterValues["GC.theta"] = args["theta"];
  }

  // INDEPENDENTWORD
  else if (freqName == "Word")
  {
    if (!AlphabetTools::isWordAlphabet(alphabet))
      throw Exception("PhylogeneticsApplicationTools::getFrequenciesSetDefaultInstance.\n\t Bad alphabet type "
                      + alphabet->getAlphabetType() + " for frequenciesset " + freqName + ".");

    const WordAlphabet* pWA = dynamic_cast<const WordAlphabet*>(alphabet);

    if (args.find("frequency") != args.end())
    {
      string sAFS = args["frequency"];

      unsigned int nbfreq = pWA->getLength();
      FrequenciesSet* pFS2;
      string st = "";
      for (unsigned i = 0; i < nbfreq; i++)
      {
        st += TextTools::toString(i+1);
      }

      map<string, string> unparsedParameterValuesNested;
      pFS2 = getFrequenciesSetDefaultInstance(pWA->getNAlphabet(0), sAFS, unparsedParameterValuesNested);
      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedParameterValues["Word." + st + "_" + it->first] = it->second;
      }
      pFS = new WordFromUniqueFrequenciesSet(pWA,pFS2);
    }
    else
    {
      if (args.find("frequency1") == args.end())
        throw Exception("PhylogeneticsApplicationTools::getFrequenciesSetDefaultInstance. Missing argument 'frequency' or 'frequency1' for frequencies set 'Word'.");
      vector<string> v_sAFS;
      vector<FrequenciesSet*> v_AFS;
      unsigned int nbfreq = 1;

      while (args.find("frequency" + TextTools::toString(nbfreq)) != args.end())
      {
        v_sAFS.push_back(args["frequency" + TextTools::toString(nbfreq++)]);
      }

      if (v_sAFS.size() != pWA->getLength())
        throw Exception("PhylogeneticsApplicationTools::getFrequenciesSetDefaultInstance. Number of frequencies (" + TextTools::toString(v_sAFS.size()) + ") does not match length of the words (" + TextTools::toString(pWA->getLength()) + ")");

      map<string, string> unparsedParameterValuesNested;
      for (unsigned i = 0; i < v_sAFS.size(); i++)
      {
        unparsedParameterValuesNested.clear();
        pFS = getFrequenciesSetDefaultInstance(pWA->getNAlphabet(i), v_sAFS[i], unparsedParameterValuesNested);
        for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
        {
          unparsedParameterValues["Word." + TextTools::toString(i+1) + "_" + it->first] = it->second;
        }
        v_AFS.push_back(pFS);
      }

      pFS = new WordFromIndependentFrequenciesSet(pWA, v_AFS);
    }
  }
  else
    throw Exception("Unknown frequency option: " + freqName);

  // Forward arguments:
  if (args.find("init") != args.end())
    unparsedParameterValues["init"] = args["init"];
  if (args.find("pseudoCount") != args.end())
    unparsedParameterValues["pseudoCount"] = args["pseudoCount"];
  if (args.find("values") != args.end())
    unparsedParameterValues["values"] = args["values"];

  return pFS;
}


/******************************************************************************/

SubstitutionModelSet* PhylogeneticsApplicationTools::getSubstitutionModelSet(
  const Alphabet* alphabet,
  const SiteContainer* data,
  map<string, string>& params,
  const string& suffix,
  bool suffixIsOptional,
  bool verbose) throw (Exception)
{
  if (!ApplicationTools::parameterExists("nonhomogeneous.number_of_models", params))
    throw Exception("You must specify this parameter: nonhomogeneous.number_of_models .");
  unsigned int nbModels = ApplicationTools::getParameter<unsigned int>("nonhomogeneous.number_of_models", params, 1, suffix, suffixIsOptional, false);
  if (nbModels == 0)
    throw Exception("The number of models can't be 0 !");

  if (verbose) ApplicationTools::displayResult("Number of distinct models", TextTools::toString(nbModels));


  // ///////////////////////////////////////////
  // Build a new model set object:

  vector<double> rateFreqs;
  string tmpDesc;
  if (AlphabetTools::isCodonAlphabet(alphabet))
    tmpDesc = ApplicationTools::getStringParameter("model1", params, "CodonNeutral(model=JC69)", suffix, suffixIsOptional, verbose);
  else if (AlphabetTools::isWordAlphabet(alphabet))
    tmpDesc = ApplicationTools::getStringParameter("model1", params, "Word(model=JC69)", suffix, suffixIsOptional, verbose);
  else
    tmpDesc = ApplicationTools::getStringParameter("model1", params, "JC69", suffix, suffixIsOptional, verbose);
  

  map<string, string> tmpUnparsedParameterValues;
  auto_ptr<SubstitutionModel> tmp(getSubstitutionModelDefaultInstance(alphabet, tmpDesc, tmpUnparsedParameterValues, true, true, 0));
  if (tmp->getNumberOfStates() != alphabet->getSize())
  {
    // Markov-Modulated Markov Model...
    unsigned int n = (unsigned int)(tmp->getNumberOfStates() / alphabet->getSize());
    rateFreqs = vector<double>(n, 1. / (double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
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
  
  SubstitutionModelSet* modelSet = stationarity ?
    new SubstitutionModelSet(alphabet, true) : //Stationarity assumed.
    new SubstitutionModelSet(alphabet, rootFrequencies);

  // //////////////////////////////////////
  // Now parse all models:

  map<string, double> existingParameters;

  for (unsigned int i = 0; i < nbModels; i++)
  {
    string prefix = "model" + TextTools::toString(i + 1);
    string modelDesc;
    if (AlphabetTools::isCodonAlphabet(alphabet))
      modelDesc = ApplicationTools::getStringParameter(prefix, params, "CodonNeutral(model=JC69)", suffix, suffixIsOptional, verbose);
    else
      if (AlphabetTools::isWordAlphabet(alphabet))
        modelDesc = ApplicationTools::getStringParameter(prefix, params, "Word(model=JC69)", suffix, suffixIsOptional, verbose);
      else
        modelDesc = ApplicationTools::getStringParameter(prefix, params, "JC69", suffix, suffixIsOptional, verbose);
    

    map<string, string> unparsedParameterValues;
    SubstitutionModel* model = getSubstitutionModelDefaultInstance(alphabet, modelDesc, unparsedParameterValues, true, true, verbose);
    prefix += ".";

    vector<string> specificParameters, sharedParameters;
    setSubstitutionModelParametersInitialValues(model,
                                                unparsedParameterValues, prefix, data,
                                                existingParameters, specificParameters, sharedParameters,
                                                verbose);
    vector<int> nodesId = ApplicationTools::getVectorParameter<int>(prefix + "nodes_id", params, ',', ':', TextTools::toString(i), suffix, suffixIsOptional, true);
    if (verbose) ApplicationTools::displayResult("Model" + TextTools::toString(i + 1) + " is associated to", TextTools::toString(nodesId.size()) + " node(s).");
    // Add model and specific parameters:
    modelSet->addModel(model, nodesId, specificParameters);
    // Now set shared parameters:
    for (unsigned int j = 0; j < sharedParameters.size(); j++)
    {
      string pName = sharedParameters[j];
      string::size_type index = pName.find(".");
      if (index == string::npos) throw Exception("PhylogeneticsApplicationTools::getSubstitutionModelSet. Bad parameter name: " + pName);
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
    if (!constDistAllowed) throw Exception("You can't use a constant distribution here!");
    rDist = new ConstantDistribution(1.,true);
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
    vector<string> v=rDist->getParameters().getParameterNames();

    for (unsigned int i = 0; i < v.size(); i++)
      unparsedParameterValues[v[i]] = TextTools::toString(rDist->getParameterValue(rDist->getParameterNameWithoutNamespace(v[i])));
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

    if (v_nestedDistrDescr.size()!=probas.size())
      throw Exception("Number of distributions do not fit the number of probabilities");
    map<string, string> unparsedParameterValuesNested;

    for (unsigned i = 0; i < v_nestedDistrDescr.size(); i++)
    {
      unparsedParameterValuesNested.clear();
      pdd = getDistributionDefaultInstance(v_nestedDistrDescr[i], unparsedParameterValuesNested, false);
          
      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedParameterValues[distName + "." + TextTools::toString(i+1) + "_" + it->first] = it->second;
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
      if (args.find("alpha") == args.end())
        throw Exception("Missing argument 'alpha' (shape) in Gamma distribution");
      if (args.find("beta") == args.end())
        throw Exception("Missing argument 'beta' (scale) in Gamma distribution");
      rDist = new GammaDiscreteDistribution(nbClasses, TextTools::to<double>(args["alpha"]),
                                            TextTools::to<double>(args["beta"]));
      unparsedParameterValues["Gamma.alpha"] = args["alpha"];
      unparsedParameterValues["Gamma.beta"] = args["beta"];
    }
    else if (distName == "Gaussian")
      {
        if (args.find("mu") == args.end())
          throw Exception("Missing argument 'mu' (mean) in Gaussian distribution");
        if (args.find("sigma") == args.end())
          throw Exception("Missing argument 'sigma' (standard deviation) in Gaussian distribution");
        rDist = new GaussianDiscreteDistribution(nbClasses, TextTools::to<double>(args["mu"]),
                                              TextTools::to<double>(args["sigma"]));
        unparsedParameterValues["Gaussian.mu"] = args["mu"];
        unparsedParameterValues["Gaussian.sigma"] = args["sigma"];
      }
    else if (distName == "Beta")
      {
        if (args.find("alpha") == args.end())
          throw Exception("Missing argument 'alpha' in Beta distribution");
        if (args.find("beta") == args.end())
          throw Exception("Missing argument 'beta' in Beta distribution");
        rDist = new BetaDiscreteDistribution(nbClasses, TextTools::to<double>(args["alpha"]),
                                                 TextTools::to<double>(args["beta"]));
        unparsedParameterValues["Beta.alpha"] = args["alpha"];
        unparsedParameterValues["Beta.beta"] = args["beta"];
      }
    else if (distName == "Exponential")
    {
      if (args.find("lambda") == args.end())
        throw Exception("Missing argument 'lambda' in Exponential distribution");
      if (args.find("median") == args.end())
        rDist = new ExponentialDiscreteDistribution(nbClasses,
                                                    TextTools::to<double>(args["lambda"]));
      else
        rDist = new ExponentialDiscreteDistribution(nbClasses,
                                                    TextTools::to<double>(args["lambda"]), true);

      unparsedParameterValues["Exponential.lambda"] = args["lambda"];
    }
    else if (distName == "TruncExponential")
    {
      if (args.find("lambda") == args.end())
        throw Exception("Missing argument 'lambda' in Truncated Exponential distribution");
      if (args.find("tp") == args.end())
        throw Exception("Missing argument 'tp' (truncation point) in Truncated Exponential distribution");
      if (args.find("median") == args.end())
        rDist = new TruncatedExponentialDiscreteDistribution(nbClasses,
                                                             TextTools::to<double>(args["lambda"]),
                                                             TextTools::to<double>(args["tp"]));
      else
        rDist = new TruncatedExponentialDiscreteDistribution(nbClasses,
                                                             TextTools::to<double>(args["lambda"]),
                                                             TextTools::to<double>(args["tp"]),
                                                             true);

      unparsedParameterValues["TruncExponential.lambda"] = args["lambda"];
      unparsedParameterValues["TruncExponential.tp"] = args["tp"];
    }
    else if (distName == "Uniform")
        {
          if (args.find("begin") == args.end())
            throw Exception("Missing argument 'begin' (mean) in Uniform distribution");
          if (args.find("end") == args.end())
            throw Exception("Missing argument 'end' (standard deviation) in Uniform distribution");
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
      vector<string> v=pMDD->getParameters().getParameterNames();

      for (unsigned int i = 0; i < v.size(); i++)
        unparsedParameterValues[v[i]] = TextTools::toString(pMDD->getParameterValue(pMDD->getParameterNameWithoutNamespace(v[i])));
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
  if (optimization == "None") return tl;
  string optName;
  map<string, string> optArgs;
  KeyvalTools::parseProcedure(optimization, optName, optArgs);

  unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2, suffix, suffixIsOptional);

  string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional);
  OutputStream* messageHandler =
    (mhPath == "none") ? 0 :
    (mhPath == "std") ? ApplicationTools::message :
    new StlOutputStream(auto_ptr<ostream>(new ofstream(mhPath.c_str(), ios::out)));
  if (verbose) ApplicationTools::displayResult("Message handler", mhPath);

  string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional);
  OutputStream* profiler =
    (prPath == "none") ? 0 :
    (prPath == "std") ? ApplicationTools::message :
    new StlOutputStream(auto_ptr<ostream>(new ofstream(prPath.c_str(), ios::out)));
  if (profiler) profiler->setPrecision(20);
  if (verbose) ApplicationTools::displayResult("Profiler", prPath);

  bool scaleFirst = ApplicationTools::getBooleanParameter("optimization.scale_first", params, false, suffix, suffixIsOptional, false);
  if (scaleFirst)
  {
    // We scale the tree before optimizing each branch length separately:
    if (verbose) ApplicationTools::displayMessage("Scaling the tree before optimizing each branch length separately.");
    double tolerance = ApplicationTools::getDoubleParameter("optimization.scale_first.tolerance", params, .0001, suffix, suffixIsOptional, true);
    if (verbose) ApplicationTools::displayResult("Scaling tolerance", TextTools::toString(tolerance));
    int nbEvalMax = ApplicationTools::getIntParameter("optimization.scale_first.max_number_f_eval", params, 1000000, suffix, suffixIsOptional, true);
    if (verbose) ApplicationTools::displayResult("Scaling max # f eval", TextTools::toString(nbEvalMax));
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
        if (!nhtl) ApplicationTools::displayWarning("The 'Ancient' parameters do not exist in homogeneous models, and will be ignored.");
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
      ApplicationTools::displayWarning("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
    }
  }

  unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000, suffix, suffixIsOptional);
  if (verbose) ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));

  double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", params, .000001, suffix, suffixIsOptional);
  if (verbose) ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));

  bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", params, false, suffix, suffixIsOptional, false);
  if (verbose) ApplicationTools::displayResult("Optimize topology", optimizeTopo ? "yes" : "no");
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
  else throw Exception("Unknown NNI algorithm: '" + nniMethod + "'.");


  string order  = ApplicationTools::getStringParameter("derivatives", optArgs, "Newton", "", true, false);
  string optMethod;
  if (order == "Gradient")
  {
    optMethod = OptimizationTools::OPTIMIZATION_GRADIENT;
  }
  else if (order == "Newton")
  {
    optMethod = OptimizationTools::OPTIMIZATION_NEWTON;
  }
  else throw Exception("Unknown derivatives algorithm: '" + order + "'.");
  if (verbose) ApplicationTools::displayResult("Optimization method", optName);
  if (verbose) ApplicationTools::displayResult("Algorithm used for derivable parameters", order);

  //See if we should reparametrize:
  bool reparam = ApplicationTools::getBooleanParameter("optimization.reparametrization", params, false);
  if (verbose) ApplicationTools::displayResult("Reparametrization", (reparam ? "yes" : "no"));

  unsigned int n = 0;
  if (optName == "DB")
  {
    // Uses Newton-Brent method:

    unsigned int nstep = ApplicationTools::getParameter<unsigned int>("nstep", optArgs, 1, "", true, false);
    if (optimizeTopo)
    {
      bool        optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, false);
      unsigned int topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, "", true, false);
      double        tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional);
      double        tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, suffixIsOptional);
      tl = OptimizationTools::optimizeTreeNNI(
        dynamic_cast<NNIHomogeneousTreeLikelihood*>(tl), parametersToEstimate,
        optNumFirst, tolBefore, tolDuring, nbEvalMax, topoNbStep, messageHandler, profiler,
        reparam, optVerbose, optMethod, nstep, nniAlgo);
    }

    if (verbose && nstep > 1) ApplicationTools::displayResult("# of precision steps", TextTools::toString(nstep));
    parametersToEstimate.matchParametersValues(tl->getParameters());
    n = OptimizationTools::optimizeNumericalParameters(
      dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(tl), parametersToEstimate,
      0, nstep, tolerance, nbEvalMax, messageHandler, profiler, reparam, optVerbose, optMethod);
  }
  else if (optName == "FullD")
  {
    // Uses Newton-raphson alogrithm with numerical derivatives when required.

    if (optimizeTopo)
    {
      bool        optNumFirst = ApplicationTools::getBooleanParameter("optimization.topology.numfirst", params, true, suffix, suffixIsOptional, false);
      unsigned int topoNbStep = ApplicationTools::getParameter<unsigned int>("optimization.topology.nstep", params, 1, "", true, false);
      double        tolBefore = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.before", params, 100, suffix, suffixIsOptional);
      double        tolDuring = ApplicationTools::getDoubleParameter("optimization.topology.tolerance.during", params, 100, suffix, suffixIsOptional);
      tl = OptimizationTools::optimizeTreeNNI2(
        dynamic_cast<NNIHomogeneousTreeLikelihood*>(tl), parametersToEstimate,
        optNumFirst, tolBefore, tolDuring, nbEvalMax, topoNbStep, messageHandler, profiler,
        reparam, optVerbose, optMethod, nniAlgo);
    }

    parametersToEstimate.matchParametersValues(tl->getParameters());
    n = OptimizationTools::optimizeNumericalParameters2(
      dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(tl), parametersToEstimate,
      0, tolerance, nbEvalMax, messageHandler, profiler, reparam, optVerbose, optMethod);
  }
  else throw Exception("Unknown optimization method: " + optName);

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
  else throw Exception("Unknown final optimization method: " + finalMethod);

  if (finalOptimizer)
  {
    parametersToEstimate.matchParametersValues(tl->getParameters());
    if (verbose) ApplicationTools::displayResult("Final optimization step", finalMethod);
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

  if (verbose) ApplicationTools::displayResult("Performed", TextTools::toString(n) + " function evaluations.");
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
  if (optimization == "None") return;
  string optName;
  map<string, string> optArgs;
  KeyvalTools::parseProcedure(optimization, optName, optArgs);

  unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2, suffix, suffixIsOptional);

  string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false, suffix, suffixIsOptional);
  OutputStream* messageHandler =
    (mhPath == "none") ? 0 :
    (mhPath == "std") ? ApplicationTools::message :
    new StlOutputStream(auto_ptr<ostream>(new ofstream(mhPath.c_str(), ios::out)));
  if (verbose) ApplicationTools::displayResult("Message handler", mhPath);

  string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false, suffix, suffixIsOptional);
  OutputStream* profiler =
    (prPath == "none") ? 0 :
    (prPath == "std") ? ApplicationTools::message :
    new StlOutputStream(auto_ptr<ostream>(new ofstream(prPath.c_str(), ios::out)));
  if (profiler) profiler->setPrecision(20);
  if (verbose) ApplicationTools::displayResult("Profiler", prPath);

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
        if (!nhtl) ApplicationTools::displayWarning("The 'Ancient' parameters do not exist in homogeneous models, and will be ignored.");
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
  if (verbose) ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));

  double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", params, .000001, suffix, suffixIsOptional);
  if (verbose) ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));

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
  else throw Exception("Option '" + order + "' is not known for 'optimization.method.derivatives'.");
  if (verbose) ApplicationTools::displayResult("Optimization method", optName);
  if (verbose) ApplicationTools::displayResult("Algorithm used for derivable parameters", order);

  unsigned int n = 0;
  if (optName == "DB")
  {
    // Uses Newton-Brent method:
    unsigned int nstep = ApplicationTools::getParameter<unsigned int>("nstep", optArgs, 1, "", true, false);
    if (verbose && nstep > 1) ApplicationTools::displayResult("# of precision steps", TextTools::toString(nstep));
    n = OptimizationTools::optimizeNumericalParametersWithGlobalClock(
      tl,
      parametersToEstimate,
      0,
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
      0,
      tolerance,
      nbEvalMax,
      messageHandler,
      profiler,
      optVerbose,
      optMethod);
  }
  else throw Exception("Unknown optimization method: " + optName);

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
  else throw Exception("Unknown final optimization method: " + finalMethod);

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

  if (verbose) ApplicationTools::displayResult("Performed", TextTools::toString(n) + " function evaluations.");
}

/******************************************************************************/

void PhylogeneticsApplicationTools::checkEstimatedParameters(const ParameterList& pl)
{
  for (unsigned int i = 0; i < pl.size(); ++i) {
    const Constraint* constraint = pl[i].getConstraint();
    if (constraint) {
      double value = pl[i].getValue();
      if (!constraint->isCorrect(value - 1e-6) || !constraint->isCorrect(value + 1e-6)) {
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
    treeWriter = new Nhx();
  else throw Exception("Unknow format for tree writing: " + format);
  if (!checkOnly)
    treeWriter->write(tree, file, true);
  delete treeWriter;
  if (verbose) ApplicationTools::displayResult("Wrote tree to file ", file);
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
  else throw Exception("Unknow format for tree writing: " + format);
  if (!checkOnly)
    treeWriter->write(trees, file, true);
  delete treeWriter;
  if (verbose) ApplicationTools::displayResult("Wrote trees to file ", file);
}

/******************************************************************************/

void PhylogeneticsApplicationTools::describeParameters_(const ParameterAliasable* parametrizable, OutputStream& out, map<string, string>& globalAliases, const vector<string>& names, bool printLocalAliases)
{
  ParameterList pl = parametrizable->getIndependentParameters().subList(names);
  unsigned int p = out.getPrecision();
  out.setPrecision(12);
  for (unsigned int i = 0; i < pl.size(); i++)
  {
    if (i > 0) out << ", ";
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
  }
  out.setPrecision(p);
}

/******************************************************************************/

void PhylogeneticsApplicationTools::describeSubstitutionModel_(const SubstitutionModel* model, OutputStream& out, map<string, string>& globalAliases)
{
  const UserProteinSubstitutionModel* trial1 = dynamic_cast<const UserProteinSubstitutionModel*>(model);
  if (trial1)
  {
    out << "Empirical";
    vector<string> pnames = trial1->getParameters().getParameterNames();
    if (pnames.size() > 0)
      out << "+F";
    out << "(file=" << trial1->getPath();
    for (unsigned int i = 0; i < pnames.size(); i++)
      out << ", " << pnames[i] << "=" << trial1->getParameterValue(pnames[i]);
    out << ")";
    out.endLine();
  }
  else
  {
    const MarkovModulatedSubstitutionModel* trial3 = dynamic_cast<const MarkovModulatedSubstitutionModel*>(model);
    if (trial3)
    {
      out << trial3->getName() << "(model=";
      const SubstitutionModel* nestedModel = trial3->getNestedModel();
      describeSubstitutionModel_(nestedModel, out, globalAliases);
      out << ", ";
      vector<string> names;
      const G2001* trial4 = dynamic_cast<const G2001*>(model);
      if (trial4)
      {
        // Also print distribution here:
        out << "rdist=";
        const DiscreteDistribution* nestedDist = trial4->getRateDistribution();
        describeDiscreteDistribution_(nestedDist, out, globalAliases);
        out << ", ";
        names.push_back(trial4->getParameter("nu").getName());
      }
      const TS98* trial5 = dynamic_cast<const TS98*>(model);
      if (trial5)
      {
        names.push_back(trial5->getParameter("s1").getName());
        names.push_back(trial5->getParameter("s2").getName());
      }
      describeParameters_(trial3, out, globalAliases, names);
      out << ")";
    }
    else
    {
      const RE08* trial4 = dynamic_cast<const RE08*>(model);
      if (trial4)
      {
        out << trial4->getName() << "(model=";
        const SubstitutionModel* nestedModel = trial4->getNestedModel();
        describeSubstitutionModel_(nestedModel, out, globalAliases);
        out << ", ";
        vector<string> names;
        names.push_back(trial4->getParameter("lambda").getName());
        names.push_back(trial4->getParameter("mu").getName());
        describeParameters_(trial4, out, globalAliases, names);
        out << ")";
      }
      else
      {
        out << model->getName() << "(";
        describeParameters_(model, out, globalAliases, model->getIndependentParameters().getParameterNames());
        out << ")";
      }
    }
  }
}

/******************************************************************************/


void PhylogeneticsApplicationTools::describeFrequenciesSet_(const FrequenciesSet* pfreqset, OutputStream& out)
{
  if (!pfreqset)
  {
    out << "None";
    return;
  }
  out << pfreqset->getName() << "(";
  ParameterList pl = pfreqset->getParameters();
  unsigned int p = out.getPrecision();
  out.setPrecision(12);
  for (unsigned int i = 0; i < pl.size(); i++)
  {
    if (i > 0) out << ", ";
    string pname = pfreqset->getParameterNameWithoutNamespace(pl[i].getName());
    (out << pname << "=").enableScientificNotation(false) << pl[i].getValue();
  }
  out << ")";
  out.setPrecision(p);
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printParameters(const SubstitutionModel* model, OutputStream& out)
{
  out << "model = ";
  map<string, string> globalAliases;
  describeSubstitutionModel_(model, out, globalAliases);
  out.endLine();
}

/******************************************************************************/

void PhylogeneticsApplicationTools::printParameters(const SubstitutionModelSet* modelSet, OutputStream& out)
{
  (out << "nonhomogeneous = general").endLine();
  (out << "nonhomogeneous.number_of_models = " << modelSet->getNumberOfModels()).endLine();

  // Get the parameter links:
  map< unsigned int, vector<string> > modelLinks; // for each model index, stores the list of global parameters.
  map< string, vector<unsigned int> > parameterLinks; // for each parameter name, stores the list of model indices.
  ParameterList pl = modelSet->getParameters();
  ParameterList plroot = modelSet->getRootFrequenciesParameters();
  for (unsigned int i = 0; i < pl.size(); i++)
  {
    if (!plroot.hasParameter(pl[i].getName()))
    {
      string name = pl[i].getName();
      vector<unsigned int> models = modelSet->getModelsWithParameter(name);
      for (unsigned int j = 0; j < models.size(); j++)
      {
        modelLinks[models[j]].push_back(name);
        parameterLinks[name].push_back(models[j]);
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
        if (parameterLinks[name][0] != i) // Otherwise, this is the 'reference' value
        {
          globalAliases[modelSet->getParameterModelName(name)] = "model" + TextTools::toString(parameterLinks[name][0] + 1) + "." + modelSet->getParameterModelName(name);
        }
      }
    }

    // Now print it:
    out.endLine() << "model" << (i + 1) << " = ";
    describeSubstitutionModel_(model, out, globalAliases);
    out.endLine();
    vector<int> ids = modelSet->getNodesWithModel(i);
    out << "model" << (i + 1) << ".nodes_id = " << ids[0];
    for (unsigned int j = 1; j < ids.size(); j++)
    {
      out << "," << ids[j];
    }
    out.endLine();
  }

  // Root frequencies:
  out.endLine();
  (out << "# Root frequencies:").endLine();
  out << "nonhomogeneous.root_freq = ";
  describeFrequenciesSet_(modelSet->getRootFrequenciesSet(), out);
}

/******************************************************************************/

void PhylogeneticsApplicationTools::describeDiscreteDistribution_(const DiscreteDistribution* rDist, OutputStream& out, map<string, string>& globalAliases)
{
  const InvariantMixedDiscreteDistribution* invar = dynamic_cast<const InvariantMixedDiscreteDistribution*>(rDist);
  const DiscreteDistribution* test = rDist;
  if (invar)
  {
    test = invar->getVariableSubDistribution();
    out << "Invariant(dist=";
    describeDiscreteDistribution_(test, out, globalAliases);
    out << ", ";
    vector<string> names;
    names.push_back(invar->getParameter("p").getName());
    describeParameters_(invar, out, globalAliases, names);
    out << ")";
  }
  else
  {
    test = dynamic_cast<const ConstantDistribution*>(rDist);
    if (test) out << "Uniform()";
    else
    {
      test = dynamic_cast<const GammaDiscreteDistribution*>(rDist);
      if (test)
      {
        out << "Gamma(n=" << rDist->getNumberOfCategories() << ", ";
        describeParameters_(rDist, out, globalAliases, rDist->getIndependentParameters().getParameterNames(), false);
        out << ")";
      }
      else throw Exception("PhylogeneticsApplicationTools::printParameters(DiscreteDistribution). Unsupported distribution.");
    }
  }
}

void PhylogeneticsApplicationTools::printParameters(const DiscreteDistribution* rDist, OutputStream& out)
{
  out << "rate_distribution = ";
  map<string, string> globalAliases;
  describeDiscreteDistribution_(rDist, out, globalAliases);
  out.endLine();
}

/******************************************************************************/

