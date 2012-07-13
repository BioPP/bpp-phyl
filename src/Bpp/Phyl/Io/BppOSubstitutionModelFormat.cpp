//
// File: BppOSubstitutionModelFormatFormat.cpp
// Created by: Laurent Guéguen
// Created on: mercredi 4 juillet 2012, à 13h 58
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

#include "BppOSubstitutionModelFormat.h"

#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/KeyvalTools.h>

#include "../Model.all"
#include "../App/PhylogeneticsApplicationTools.h"

#include "BppOFrequenciesSetFormat.h"

#include <Bpp/Seq/App/SequenceApplicationTools.h>

#include <Bpp/Numeric/Prob.all>

using namespace bpp;

// From the STL:
#include <iomanip>

using namespace std;

SubstitutionModel* BppOSubstitutionModelFormat::read(const Alphabet* alphabet,
                                                     const std::string& modelDescription,
                                                     std::map<std::string, std::string>& unparsedParameterValues,
                                                     bool allowCovarions,
                                                     bool allowMixed,
                                                     bool allowGaps,
                                                     bool verbose)
{
  SubstitutionModel* model = 0;
  string modelName = "";
  map<string, string> args;
  BppOFrequenciesSetFormat* bIOFreq=new BppOFrequenciesSetFormat();

  KeyvalTools::parseProcedure(modelDescription, modelName, args);

  // //////////////////////////////////
  // / MIXED MODELS
  // ////////////////////////////////

  if ((modelName == "MixedModel" || (modelName == "Mixture")) && allowMixed)
    model=readMixed(alphabet, modelDescription, unparsedParameterValues, allowCovarions, allowGaps, verbose);
  
  
  // /////////////////////////////////
  // / WORDS and CODONS
  // ///////////////////////////////

  else if ((modelName == "Word") || (modelName == "Triplet") || (modelName.substr(0,5) == "Codon"))
    model=readWord(alphabet, modelDescription, unparsedParameterValues, allowCovarions, allowMixed, allowGaps, verbose);


// //////////////////////////////////////
// PREDEFINED CODON MODELS
// //////////////////////////////////////

  else if ((modelName == "MG94") || (modelName == "YN98") ||
           (modelName == "GY94") || (modelName.substr(0, 5) == "YNGKP"))
    {
      if (!AlphabetTools::isCodonAlphabet(alphabet))
        throw Exception("Alphabet should be Codon Alphabet.");

      const CodonAlphabet* pCA = (const CodonAlphabet*)(alphabet);

      if (args.find("genetic_code") == args.end())
        args["genetic_code"] = pCA->getAlphabetType();

      GeneticCode* pgc = SequenceApplicationTools::getGeneticCode(dynamic_cast<const NucleicAlphabet*>(pCA->getNAlphabet(0)), args["genetic_code"]);
      if (pgc->getSourceAlphabet()->getAlphabetType() != pCA->getAlphabetType())
        throw Exception("Mismatch between genetic code and codon alphabet");

      string freqOpt = ApplicationTools::getStringParameter("frequencies", args, "F0", "", true, verbose);
      map<string, string> unparsedParameterValuesNested;
      
      FrequenciesSet* codonFreqs = bIOFreq->read(pCA, freqOpt, unparsedParameterValuesNested);

      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
        {
          unparsedParameterValues[modelName + "." + it->first] = it->second;
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
          model = new YNGKP_M3(pgc, codonFreqs, TextTools::to<unsigned int>(args["n"]));
      else if ((modelName == "YNGKP_M7") || modelName == "YNGKP_M8")
        {
          if (args.find("n") == args.end())
            throw Exception("Missing argument 'n' (number of classes) in " + modelName + " distribution");
          unsigned int nbClasses = TextTools::to<unsigned int>(args["n"]);

          if (modelName == "YNGKP_M7")
            model = new YNGKP_M7(pgc, codonFreqs, nbClasses);
          else if (modelName == "YNGKP_M8")
            model = new YNGKP_M8(pgc, codonFreqs, nbClasses);
        }
      else
        throw Exception("Unknown Codon model: " + modelName);
    }

  
  ////////////////////////////////////                                                                                                                                             
  // gBGC                                                                                                                                                                           
  ////////////////////////////////////
  
  else if (modelName == "gBGC")
    {
      // We have to parse the nested model first:
      string nestedModelDescription = args["model"];
      if (TextTools::isEmpty(nestedModelDescription))
        throw Exception("PhylogeneticsApplicationTools::read. Missing argument 'model' for model 'gBGC'.");
      if (verbose)
        ApplicationTools::displayResult("Biased gene conversion", modelName);
      map<string, string> unparsedParameterValuesNested;
      SubstitutionModel* nestedModel = read(alphabet, nestedModelDescription, unparsedParameterValuesNested, allowCovarions, allowMixed, false, verbose);

      // Now we create the RE08 substitution model:
      model = new gBGC(dynamic_cast<const NucleicAlphabet*>(alphabet), dynamic_cast<NucleotideSubstitutionModel*>(nestedModel));

      // Then we update the parameter set:
      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
        {
          unparsedParameterValues["gBGC." + it->first] = it->second;
        }
    }
  
  ////////////////////////////////////
  // YpR
  ////////////////////////////////////
  
  else if (modelName == "YpR_Sym")
    {
      if (alphabet->getAlphabetType()!="RNY alphabet")
        throw Exception("Mismatch alphabet: " + alphabet->getAlphabetType() + " for model: " + modelName);
      const RNY* prny=dynamic_cast<const RNY*>(alphabet);

      string nestedModelDescription = args["model"];
      if (TextTools::isEmpty(nestedModelDescription))
        throw Exception("PhylogeneticsApplicationTools::read. Missing argument 'model' for model 'YpR_sym'.");
      if (verbose)
        ApplicationTools::displayResult("Symetric YpR model" , modelName);
      map<string, string> unparsedParameterValuesNested;
      SubstitutionModel* nestedModel = read(&prny->getLetterAlphabet(), nestedModelDescription, unparsedParameterValuesNested, false, false, false, verbose);

      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
        unparsedParameterValues["YpR_Sym." + it->first] = it->second;

      model=new YpR_Sym(prny, nestedModel);
    }
  else if (modelName=="YpR_Gen"){
    if (alphabet->getAlphabetType()!="RNY alphabet")
      throw Exception("Mismatch alphabet: " + alphabet->getAlphabetType() + " for model: " + modelName);
    const RNY* prny=dynamic_cast<const RNY*>(alphabet);

    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("PhylogeneticsApplicationTools::read. Missing argument 'model' for model 'YpR_gen'.");
    if (verbose)
      ApplicationTools::displayResult("General YpR model" , modelName);
    map<string, string> unparsedParameterValuesNested;
    SubstitutionModel* nestedModel = read(&prny->getLetterAlphabet(), nestedModelDescription, unparsedParameterValuesNested, false, false, false,verbose);

    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      unparsedParameterValues["YpR_Gen." + it->first] = it->second;

    model=new YpR_Gen(prny, nestedModel);
  }


  // /////////////////////////////////
  // / RE08
  // ///////////////////////////////

  else if (modelName == "RE08")
    {
      if (!allowGaps)
        throw Exception("PhylogeneticsApplicationTools::read. No Gap model allowed here.");

      // We have to parse the nested model first:
      string nestedModelDescription = args["model"];
      if (TextTools::isEmpty(nestedModelDescription))
        throw Exception("PhylogeneticsApplicationTools::read. Missing argument 'model' for model 'RE08'.");
      if (verbose)
        ApplicationTools::displayResult("Gap model", modelName);
      map<string, string> unparsedParameterValuesNested;
      SubstitutionModel* nestedModel = read(alphabet, nestedModelDescription, unparsedParameterValuesNested, allowCovarions, allowMixed, false, verbose);

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
        throw Exception("PhylogeneticsApplicationTools::read. No Covarion model allowed here.");

      // We have to parse the nested model first:
      string nestedModelDescription = args["model"];
      if (TextTools::isEmpty(nestedModelDescription))
        throw Exception("PhylogeneticsApplicationTools::read. Missing argument 'model' for model 'TS98'.");
      if (verbose)
        ApplicationTools::displayResult("Covarion model", modelName);
      map<string, string> unparsedParameterValuesNested;
      SubstitutionModel* nestedModel = read(alphabet, nestedModelDescription, unparsedParameterValuesNested, false, allowMixed, allowGaps, verbose);

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
        throw Exception("PhylogeneticsApplicationTools::read. No Covarion model allowed here.");

      // We have to parse the nested model first:
      string nestedModelDescription = args["model"];
      if (TextTools::isEmpty(nestedModelDescription))
        throw Exception("PhylogeneticsApplicationTools::read. Missing argument 'model' for model 'G01'.");
      string nestedRateDistDescription = args["rdist"];
      if (TextTools::isEmpty(nestedRateDistDescription))
        throw Exception("PhylogeneticsApplicationTools::read. Missing argument 'rdist' for model 'G01'.");
      if (verbose)
        ApplicationTools::displayResult("Covarion model", modelName);

      map<string, string> unparsedParameterValuesNestedModel;
      SubstitutionModel* nestedModel = read(alphabet, nestedModelDescription, unparsedParameterValuesNestedModel, false, allowMixed, allowGaps, verbose);
      map<string, string> unparsedParameterValuesNestedDist;
      DiscreteDistribution* nestedRDist = PhylogeneticsApplicationTools::getRateDistributionDefaultInstance(nestedRateDistDescription, unparsedParameterValuesNestedDist, false, verbose);

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
          // / SSR
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
          // / RN95
          // ///////////////////////////////

          else if (modelName == "RN95")
            {
              model = new RN95(alpha);
            }

          // /////////////////////////////////
          // / RN95s
          // ///////////////////////////////

          else if (modelName == "RN95s")
            {
              model = new RN95s(alpha);
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
          else if (AlphabetTools::isBinaryAlphabet(alphabet))
            {
              const BinaryAlphabet* balpha = dynamic_cast<const BinaryAlphabet*>(alphabet);

              if (modelName == "Binary")
                model = new BinarySubstitutionModel(balpha);
            }
          else
            throw Exception("Model '" + modelName + "' unknown.");
        }
      if (verbose)
        ApplicationTools::displayResult("Substitution model", modelName);
    }

  // Update parameter args:
  vector<string> pnames = model->getParameters().getParameterNames();

  string pref=model->getNamespace();

  for (unsigned int i = 0; i < pnames.size(); i++)
    {
      string name = model->getParameterNameWithoutNamespace(pnames[i]);
      if (args.find(name) != args.end())
        unparsedParameterValues[pref + name] = args[name];
    }

  // Now look if some parameters are aliased:
  ParameterList pl = model->getIndependentParameters();
  string pname, pval, pname2;
  for (unsigned int i = 0; i < pl.size(); i++)
    {
      pname = model->getParameterNameWithoutNamespace(pl[i].getName());

      if (args.find(pname) == args.end())
        continue;
      pval = args[pname];

      if ((pval.length() >= 5 && pval.substr(0, 5) == "model") ||
          (pval.find("(") != string::npos))
        continue;
      bool found = false;
      for (unsigned int j = 0; j < pl.size() && !found; j++)
        {
          pname2 = model->getParameterNameWithoutNamespace(pl[j].getName());

          // if (j == i || args.find(pname2) == args.end()) continue; Julien 03/03/2010: This extra condition prevents complicated (nested) models to work properly...
          if (j == i)
            continue;
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

  // 2 following tests be removed in a later version
  if (args.find("useObservedFreqs") != args.end())
    throw Exception("useObservedFreqs argument is obsolete. Please use 'initFreqs=observed' instead.");
  if (args.find("useObservedFreqs.pseudoCount") != args.end())
    throw Exception("useObservedFreqs.pseudoCount argument is obsolete. Please use 'initFreqs.observedPseudoCount' instead.");

  
  if (args.find("initFreqs") != args.end())
    unparsedParameterValues[pref + "initFreqs"] = args["initFreqs"];
  if (args.find("initFreqs.observedPseudoCount") != args.end())
    unparsedParameterValues[pref + "initFreqs.observedPseudoCount"] = args["initFreqs.observedPseudoCount"];

  delete bIOFreq;
  return model;

}


SubstitutionModel* BppOSubstitutionModelFormat::readWord(const Alphabet* alphabet,
                                                              const std::string& modelDescription,
                                                              std::map<std::string, std::string>& unparsedParameterValues,
                                                              bool allowCovarions,
                                                              bool allowMixed,
                                                              bool allowGaps,
                                                              bool verbose)
{
  BppOFrequenciesSetFormat* bIOFreq=new BppOFrequenciesSetFormat();
  SubstitutionModel* model = 0;
  string modelName = "";
  map<string, string> args;

  KeyvalTools::parseProcedure(modelDescription, modelName, args);

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



  if (args.find("model") != args.end()){
      v_nestedModelDescription.push_back(args["model"]);
      nbmodels = (modelName == "Word")?pWA->getLength():3;
  }
  else
    {
      if (args.find("model1") == args.end())
          throw Exception("Missing argument 'model' or 'model1' for model " + modelName + ".");

      nbmodels = 0;

      while (args.find("model" + TextTools::toString(nbmodels + 1)) != args.end())
          v_nestedModelDescription.push_back(args["model" + TextTools::toString(++nbmodels)]);
    }

  if (nbmodels < 2)
    throw Exception("Missing nested models for model " + modelName + ".");

  if (pWA->getLength() != nbmodels)
    throw Exception("Bad alphabet type "
                    + alphabet->getAlphabetType() + " for  model " + modelName + ".");

  
  map<string, string> unparsedParameterValuesNested;

  if (v_nestedModelDescription.size() != nbmodels)
    {
      model = read(pWA->getNAlphabet(0), v_nestedModelDescription[0], unparsedParameterValuesNested, false, true, false, false);
      string pref = "";
      for (unsigned int i = 0; i < nbmodels; i++)
          pref += TextTools::toString(i + 1);

      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
          unparsedParameterValues[modelName + "." + pref + "_" + it->first] = it->second;

      v_pSM.push_back(model);
    }
  else
    {
      for (unsigned i = 0; i < v_nestedModelDescription.size(); i++)
        {
          unparsedParameterValuesNested.clear();
          model = read(pWA->getNAlphabet(i), v_nestedModelDescription[i], unparsedParameterValuesNested, false, true, false, false);
          for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
              unparsedParameterValues[modelName + "." + TextTools::toString(i + 1) + "_" + it->first] = it->second;

          v_pSM.push_back(model);
        }
    }

  // /////////////////////////////////
  // / WORD
  // ///////////////////////////////

  if (modelName == "Word")
    {
      model = (v_nestedModelDescription.size() != nbmodels)
        ? new WordSubstitutionModel(v_pSM[0], nbmodels)
        : new WordSubstitutionModel(v_pSM);
    }

  // /////////////////////////////////
  // / CODON
  // ///////////////////////////////

  else {

    const CodonAlphabet* pCA=dynamic_cast<const CodonAlphabet*>(pWA);
    if (pCA==0)
      throw Exception("Non codon Alphabet fo model" + modelName + " model.");
      
    AlphabetIndex2<double>* pai2(0);
    GeneticCode* pgc(0);
    FrequenciesSet* pFS(0);
    
    if ((dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]) == 0) ||
        ((v_nestedModelDescription.size() == 3) && 
         (dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]) == 0 || dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]) == 0)))
      throw Exception("Non simple NucleotideSubstitutionModel imbedded in " + modelName + " model.");


    if (modelName.find("Dist")!=string::npos)
      {
        if (args.find("genetic_code") == args.end())
          args["genetic_code"] = pCA->getAlphabetType();

        pgc = SequenceApplicationTools::getGeneticCode(dynamic_cast<const NucleicAlphabet*>(pCA->getNAlphabet(0)), args["genetic_code"]);
        if (pgc->getSourceAlphabet()->getAlphabetType() != pCA->getAlphabetType())
          throw Exception("Mismatch between genetic code and codon alphabet");

        pai2 = (args.find("aadistance") == args.end())?0:SequenceApplicationTools::getAADistance(args["aadistance"]);
      }
    

    if (modelName.find("Freq")!=string::npos){
      if (args.find("frequencies") == args.end())
        throw Exception("Missing equilibrium frequencies.");

      unparsedParameterValuesNested.clear();
      pFS = bIOFreq->read(pCA, args["frequencies"], unparsedParameterValuesNested);

      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
        unparsedParameterValues[modelName + "." + it->first] = it->second;
    }
      

    ////
    
    if (modelName == "Triplet")
      model = (v_nestedModelDescription.size() != 3)
        ?new TripletSubstitutionModel(
                                      pCA, 
                                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]))
        :new TripletSubstitutionModel(
                                      pCA, 
                                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]));
    
    else if (modelName == "CodonRate")
      model = (v_nestedModelDescription.size() != 3)
        ?new CodonRateSubstitutionModel(
                                        pCA, 
                                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]))
        :new CodonRateSubstitutionModel(
                                        pCA, 
                                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]));

    
    else if (modelName == "CodonDistance")
      {
        if (v_nestedModelDescription.size() != 3)
          model = new CodonDistanceSubstitutionModel(pgc,
                                                     dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]), pai2);
        else
          model = new CodonDistanceSubstitutionModel(
                                                     pgc,
                                                     dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                     dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                                                     dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]), pai2);
      }

    else if (modelName == "CodonRateFreq")
      {
        if (v_nestedModelDescription.size() != 3)
          model = new CodonRateFrequenciesSubstitutionModel(pCA,
                                                            dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]), pFS);
        else
          model = new CodonRateFrequenciesSubstitutionModel(pCA,
                                                            dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                            dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                                                            dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]),
                                                            pFS);
      }
    
    else if (modelName == "CodonDistFreq")
      {
        if (v_nestedModelDescription.size() != 3)
          model = new CodonDistanceFrequenciesSubstitutionModel(pgc,
                                                                dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                                pFS,
                                                                pai2);
        else
          model = new CodonDistanceFrequenciesSubstitutionModel(
                                                                pgc,
                                                                dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                                dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                                                                dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]),
                                                                pFS,
                                                                pai2);
      }

    else if (modelName == "CodonDistPhasFreq")
      {
        if (v_nestedModelDescription.size() != 3)
          model = new CodonDistancePhaseFrequenciesSubstitutionModel(pgc,
                                                                     dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                                     pFS,
                                                                     pai2);
        else
          model = new CodonDistancePhaseFrequenciesSubstitutionModel(
                                                                     pgc,
                                                                     dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                                     dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                                                                     dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]),
                                                                     pFS,
                                                                     pai2);
      }
  }
  delete bIOFreq;
  return model;
}


MixedSubstitutionModel* BppOSubstitutionModelFormat::readMixed(const Alphabet* alphabet,
                                                               const std::string& modelDescription,
                                                               std::map<std::string, std::string>& unparsedParameterValues,
                                                               bool allowCovarions,
                                                               bool allowGaps,
                                                               bool verbose)
{
  MixedSubstitutionModel* model(0);

  string modelName = "";
  map<string, string> args;
  KeyvalTools::parseProcedure(modelDescription, modelName, args);

  map<string, string> unparsedParameterValuesNested;
  SubstitutionModel* pSM;

  if (modelName == "MixedModel"){

    if (args.find("model") == args.end())
      throw Exception("The argument 'model' is missing from MixedSubstitutionModel description");
    string nestedModelDescription = args["model"];
    pSM = read(alphabet,
               nestedModelDescription,
               unparsedParameterValuesNested,
               allowCovarions,
               true,
               allowGaps,
               verbose);
  
    map<string, DiscreteDistribution*> mdist;
    map<string, string> unparsedParameterValuesNested2, unparsedParameterValuesNested3;

    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin();
         it != unparsedParameterValuesNested.end();
         it++){
      if (it->second.find("(") != string::npos){
        unparsedParameterValuesNested3.clear();
        mdist[pSM->getParameterNameWithoutNamespace(it->first)] = PhylogeneticsApplicationTools::getDistributionDefaultInstance(it->second, unparsedParameterValuesNested3, true);
        for (map<string, string>::iterator it2 = unparsedParameterValuesNested3.begin();
             it2 != unparsedParameterValuesNested3.end();
             it2++)
          unparsedParameterValuesNested2[it->first + "_" + it2->first] = it2->second;

      }
      else
        unparsedParameterValuesNested2[it->first] = it->second;
    }

    for (map<string, string>::iterator it = unparsedParameterValuesNested2.begin();
         it != unparsedParameterValuesNested2.end();
         it++)
      unparsedParameterValues[it->first] = it->second;

    
    int fi(-1), ti(-1);

    if (args.find("from") != args.end())
      fi = alphabet->charToInt(args["from"]);
    if (args.find("to") != args.end())
      ti = alphabet->charToInt(args["to"]);

    model = new MixtureOfASubstitutionModel(alphabet, pSM, mdist, fi, ti);

    vector<string> v = model->getParameters().getParameterNames();

    for (map<string, DiscreteDistribution*>::iterator it = mdist.begin();
         it != mdist.end(); it++)
      delete it->second;
    
    if (verbose)
      ApplicationTools::displayResult("Mixture Of A Substitution Model", nestedModelDescription );
  }

  
  else if (modelName == "Mixture") {
    vector<string> v_nestedModelDescription;
    vector<SubstitutionModel*> v_pSM;
      
    if (args.find("model1") == args.end())
      {
        throw Exception("Missing argument 'model1' for model " + modelName + ".");
      }
    unsigned int nbmodels = 0;
      
    while (args.find("model" + TextTools::toString(nbmodels + 1)) != args.end())
      {
        v_nestedModelDescription.push_back(args["model" + TextTools::toString(++nbmodels)]);
      }
      
    if (nbmodels < 2)
      throw Exception("Missing nested models for model " + modelName + ".");

    for (unsigned i = 0; i < v_nestedModelDescription.size(); i++)
      {
        unparsedParameterValuesNested.clear();
        pSM = read(alphabet, v_nestedModelDescription[i], unparsedParameterValuesNested, false, true, false, false);
        for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
          {
            unparsedParameterValues[modelName + "." + TextTools::toString(i + 1) + "_" + it->first] = it->second;
          }
        v_pSM.push_back(pSM);
      }
      
    model = new MixtureOfSubstitutionModels(alphabet, v_pSM);
    if (verbose)
      ApplicationTools::displayResult("Mixture Of Substitution Models", modelName );
  }

  else
    throw Exception("Unknown model name for mixture " + modelName);

  return model;
}


void BppOSubstitutionModelFormat::write(const SubstitutionModel& model,
                                        OutputStream& out,
                                        std::map<std::string, std::string>& globalAliases,
                                        std::vector<std::string>& writtenNames) const
{

  BppOFrequenciesSetFormat* bIOFreq=new BppOFrequenciesSetFormat();
  // Mixed Model

  if (dynamic_cast<const MixedSubstitutionModel*>(&model)!=NULL){
    write(*dynamic_cast<const MixedSubstitutionModel*>(&model), out, globalAliases, writtenNames);
    return;
  }
  
  // Is it a protein user defined model?
  const UserProteinSubstitutionModel* userModel = dynamic_cast<const UserProteinSubstitutionModel*>(&model);
  if (userModel)
    {
      out << "Empirical";
      vector<string> pnames = userModel->getParameters().getParameterNames();
      if (pnames.size() > 0)
        out << "+F";
      out << "(file=" << userModel->getPath();
      for (unsigned int i = 0; i < pnames.size(); i++)
        {
          if (find(writtenNames.begin(),writtenNames.end(),pnames[i])==writtenNames.end()){
            out << ", " << pnames[i] << "=" << userModel->getParameterValue(pnames[i]);
            writtenNames.push_back(pnames[i]);
          }
        }
      out << ")";
      out.endLine();
      return;
    }

  // Is it a markov-modulated model?
  const MarkovModulatedSubstitutionModel* mmModel = dynamic_cast<const MarkovModulatedSubstitutionModel*>(&model);
  if (mmModel)
    {
      out << mmModel->getName() << "(model=";
      const SubstitutionModel* nestedModel = mmModel->getNestedModel();
      write(*nestedModel, out, globalAliases, writtenNames);
      out << ", ";

      const G2001* gModel = dynamic_cast<const G2001*>(&model);
      if (gModel)
        {
          // Also print distribution here:
          out << "rdist=";
          const DiscreteDistribution* nestedDist = gModel->getRateDistribution();
          PhylogeneticsApplicationTools::describeDiscreteDistribution_(nestedDist, out, globalAliases, writtenNames);
          out << ", ";
        }
      PhylogeneticsApplicationTools::describeParameters_(mmModel, out, globalAliases, model.getIndependentParameters().getParameterNames(), writtenNames);
      out << ")";
      return;
    }

  // Is it a model with gaps?
  const RE08* reModel = dynamic_cast<const RE08*>(&model);
  if (reModel)
    {
      out << reModel->getName() << "(model=";
      const SubstitutionModel* nestedModel = reModel->getNestedModel();
      write(*nestedModel, out, globalAliases, writtenNames);
      out << ", ";
      PhylogeneticsApplicationTools::describeParameters_(reModel, out, globalAliases, model.getIndependentParameters().getParameterNames(), writtenNames);
      out << ")";
      return;
    }

  // Is it a YpR model?
  const YpR* yprModel= dynamic_cast<const YpR*>(&model);
  if (yprModel)
    {
      out << yprModel->getName() << "(model=";
      const SubstitutionModel* nestedModel = yprModel->getNestedModel();
      write(*nestedModel, out, globalAliases, writtenNames);
      out << ", ";
      PhylogeneticsApplicationTools::describeParameters_(yprModel, out, globalAliases, model.getIndependentParameters().getParameterNames(), writtenNames);
      out << ")";
      return;
    }
  // Regular model

  out << model.getName() << "(";

  const FrequenciesSet* pfs=model.getFrequenciesSet();
  if (pfs){ 
    out << "frequencies=";
    bIOFreq->write(pfs, out, writtenNames);
  }
  PhylogeneticsApplicationTools::describeParameters_(&model, out, globalAliases, model.getIndependentParameters().getParameterNames(), writtenNames, true, pfs!=NULL);
  
  out << ")";
  delete bIOFreq;
}



void BppOSubstitutionModelFormat::write(const MixedSubstitutionModel& model,
                                        OutputStream& out,
                                        std::map<std::string, std::string>& globalAliases,
                                        std::vector<std::string>& writtenNames) const
{
  bool flag(false);

  
  if (dynamic_cast<const MixtureOfSubstitutionModels*>(&model)!=NULL){
    const MixtureOfSubstitutionModels* pMS=dynamic_cast<const MixtureOfSubstitutionModels*>(&model);

    for (unsigned int i=0; i< pMS->getNumberOfModels(); i++){
      const SubstitutionModel* eM=pMS->getNModel(i);

      vector<string> vpl=eM->getIndependentParameters().getParameterNames();
      for (unsigned j=0;j<vpl.size();j++)
        if (eM->getParameterNameWithoutNamespace(vpl[j])=="rate")
          writtenNames.push_back(vpl[j]);
    }

    out << "Mixture(";
    for (unsigned int i=0; i< pMS->getNumberOfModels(); i++){
      if (i!=0)
        out << ", ";
      out << "model" + TextTools::toString(i+1) + "=";
      write(*pMS->getNModel(i),out,globalAliases, writtenNames);
    }
  }
  else {
    const MixtureOfASubstitutionModel* pMS=dynamic_cast<const MixtureOfASubstitutionModel*>(&model);
    out << "MixedModel(model= ";
    const SubstitutionModel* eM=pMS->getNModel(0);

    ParameterList pl=eM->getIndependentParameters();
    vector<string> vpl=pl.getParameterNames();

    
    out << eM->getName() << "(";
    for (unsigned j=0;j<vpl.size();j++){
      if (find(writtenNames.begin(),writtenNames.end(),vpl[j])==writtenNames.end()){
        if (eM->getParameterNameWithoutNamespace(vpl[j])!="rate"){
          if (flag)
            out << ",";
          else
            flag=true;
          
          out << eM->getParameterNameWithoutNamespace(vpl[j]) << "=" ;
          const DiscreteDistribution* pDD=pMS->getDistribution(vpl[j]);
          if (pDD && dynamic_cast<const ConstantDistribution*>(pDD)==NULL)
            PhylogeneticsApplicationTools::describeDiscreteDistribution_(pDD, out, globalAliases, writtenNames);
          else
            out << pl[j].getValue();
        }
        writtenNames.push_back(vpl[j]);
      }
    }
    out << ")";
    if (pMS->from()!=-1)
      out << ", from=" << model.getAlphabet()->intToChar(pMS->from()) << ", to=" << model.getAlphabet()->intToChar(pMS->to());
  }

  PhylogeneticsApplicationTools::describeParameters_(&model, out, globalAliases, model.getIndependentParameters().getParameterNames(), writtenNames, true, true);

  out << ")";
}
  
 
