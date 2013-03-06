//
// File: BppOSubstitutionModelFormat.cpp
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
#include "BppORateDistributionFormat.h"

#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/KeyvalTools.h>

#include "../Model/Codon/MG94.h"
#include "../Model/Codon/GY94.h"
#include "../Model/Codon/YNGKP_M1.h"
#include "../Model/Codon/YNGKP_M2.h"
#include "../Model/Codon/YNGKP_M3.h"
#include "../Model/Codon/YNGKP_M7.h"
#include "../Model/Codon/YNGKP_M8.h"
#include "../Model/Codon/YN98.h"
#include "../Model/Codon/TripletSubstitutionModel.h"
#include "../Model/Codon/CodonRateSubstitutionModel.h"
#include "../Model/Codon/CodonDistanceSubstitutionModel.h"
#include "../Model/Codon/CodonRateFrequenciesSubstitutionModel.h"
#include "../Model/Codon/CodonDistanceFrequenciesSubstitutionModel.h"
#include "../Model/Codon/CodonDistancePhaseFrequenciesSubstitutionModel.h"
#include "../Model/Codon/CodonDistanceFitnessPhaseFrequenciesSubstitutionModel.h"
#include "../Model/RE08.h"
#include "../Model/TS98.h"
#include "../Model/G2001.h"
#include "../Model/Nucleotide/F84.h"
#include "../Model/Nucleotide/NucleotideSubstitutionModel.h"
#include "../Model/Nucleotide/gBGC.h"
#include "../Model/Nucleotide/RN95.h"
#include "../Model/Nucleotide/GTR.h"
#include "../Model/Nucleotide/RN95s.h"
#include "../Model/Nucleotide/HKY85.h"
#include "../Model/Nucleotide/SSR.h"
#include "../Model/Nucleotide/JCnuc.h"
#include "../Model/Nucleotide/T92.h"
#include "../Model/Nucleotide/K80.h"
#include "../Model/Nucleotide/TN93.h"
#include "../Model/Nucleotide/L95.h"
#include "../Model/Nucleotide/YpR.h"
#include "../Model/Protein/CoalaCore.h"
#include "../Model/Protein/LLG08_EX2.h"
#include "../Model/Protein/Coala.h"
#include "../Model/Protein/LLG08_EX3.h"
#include "../Model/Protein/DSO78.h"
#include "../Model/Protein/LLG08_UL2.h"
#include "../Model/Protein/JCprot.h"
#include "../Model/Protein/LLG08_UL3.h"
#include "../Model/Protein/JTT92.h"
#include "../Model/Protein/ProteinSubstitutionModel.h"
#include "../Model/Protein/LG08.h" 
#include "../Model/Protein/UserProteinSubstitutionModel.h"
#include "../Model/Protein/LGL08_CAT.h"
#include "../Model/Protein/WAG01.h"
#include "../Model/Protein/LLG08_EHO.h"
#include "../Model/BinarySubstitutionModel.h"

#include "../App/PhylogeneticsApplicationTools.h"

#include "BppOFrequenciesSetFormat.h"

#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>

#include <Bpp/Io/OutputStream.h>
#include <Bpp/Io/BppOParametrizableFormat.h>
#include <Bpp/Io/BppODiscreteDistributionFormat.h>

//From Numeric

#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/AutoParameter.h>

using namespace bpp;

// From the STL:
#include <iomanip>

using namespace std;

unsigned char BppOSubstitutionModelFormat::DNA = 1;
unsigned char BppOSubstitutionModelFormat::RNA = 2;
unsigned char BppOSubstitutionModelFormat::NUCLEOTIDE = 1 | 2;
unsigned char BppOSubstitutionModelFormat::PROTEIN = 4;
unsigned char BppOSubstitutionModelFormat::CODON = 8;
unsigned char BppOSubstitutionModelFormat::WORD = 16;
unsigned char BppOSubstitutionModelFormat::BINARY = 32;
unsigned char BppOSubstitutionModelFormat::ALL = 1 | 2 | 4 | 8 | 16 | 32;


SubstitutionModel* BppOSubstitutionModelFormat::read(const Alphabet* alphabet,
                                                     const std::string& modelDescription,
                                                     const SiteContainer* data,
                                                     bool parseArguments)
{
  unparsedArguments_.clear();
  auto_ptr<SubstitutionModel> model;
  string modelName = "";
  map<string, string> args;
  KeyvalTools::parseProcedure(modelDescription, modelName, args);

  // //////////////////////////////////
  // / MIXED MODELS
  // ////////////////////////////////

  if ((modelName == "MixedModel" || (modelName == "Mixture")) && allowMixed_)
    model.reset(readMixed_(alphabet, modelDescription, data));


  // /////////////////////////////////
  // / WORDS and CODONS
  // ///////////////////////////////

  else if ((modelName == "Word") || (modelName == "Triplet") || (modelName.substr(0, 5) == "Codon"))
    model.reset(readWord_(alphabet, modelDescription, data));


  // //////////////////////////////////////
  // PREDEFINED CODON MODELS
  // //////////////////////////////////////

  else if (((modelName == "MG94") || (modelName == "YN98") ||
            (modelName == "GY94") || (modelName.substr(0, 5) == "YNGKP")) && (alphabetCode_ & CODON))
  {
    if (!(alphabetCode_ & CODON))
      throw Exception("BppOSubstitutionModelFormat::read. Codon alphabet not supported.");
    if (!AlphabetTools::isCodonAlphabet(alphabet))
      throw Exception("Alphabet should be Codon Alphabet.");

    const CodonAlphabet* pCA = dynamic_cast<const CodonAlphabet*>(alphabet);

    if (args.find("genetic_code") == args.end())
      args["genetic_code"] = pCA->getAlphabetType();

    auto_ptr<GeneticCode> pgc(SequenceApplicationTools::getGeneticCode(dynamic_cast<const NucleicAlphabet*>(pCA->getNAlphabet(0)), args["genetic_code"]));
    if (pgc->getSourceAlphabet()->getAlphabetType() != pCA->getAlphabetType())
      throw Exception("Mismatch between genetic code and codon alphabet");

    string freqOpt = ApplicationTools::getStringParameter("frequencies", args, "F0", "", true, verbose_);
    BppOFrequenciesSetFormat freqReader(BppOFrequenciesSetFormat::ALL, verbose_);
    auto_ptr<FrequenciesSet> codonFreqs(freqReader.read(pCA, freqOpt, data, false));
    map<string, string> unparsedParameterValuesNested(freqReader.getUnparsedArguments());

    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedArguments_[modelName + "." + it->first] = it->second;
    }

    if (modelName == "MG94")
      model.reset(new MG94(pgc.release(), codonFreqs.release()));
    else if (modelName == "GY94")
      model.reset(new GY94(pgc.release(), codonFreqs.release()));
    else if ((modelName == "YN98") || (modelName == "YNGKP_M0"))
      model.reset(new YN98(pgc.release(), codonFreqs.release()));
    else if (modelName == "YNGKP_M1")
      model.reset(new YNGKP_M1(pgc.release(), codonFreqs.release()));
    else if (modelName == "YNGKP_M2")
      model.reset(new YNGKP_M2(pgc.release(), codonFreqs.release()));
    else if (modelName == "YNGKP_M3")
      if (args.find("n") == args.end())
        model.reset(new YNGKP_M3(pgc.release(), codonFreqs.release()));
      else
        model.reset(new YNGKP_M3(pgc.release(), codonFreqs.release(), TextTools::to<unsigned int>(args["n"])));
    else if ((modelName == "YNGKP_M7") || modelName == "YNGKP_M8")
    {
      if (args.find("n") == args.end())
        throw Exception("Missing argument 'n' (number of classes) in " + modelName + " distribution");
      unsigned int nbClasses = TextTools::to<unsigned int>(args["n"]);
      if (verbose_)
        ApplicationTools::displayResult("Number of classes in model", nbClasses);

      if (modelName == "YNGKP_M7")
        model.reset(new YNGKP_M7(pgc.release(), codonFreqs.release(), nbClasses));
      else if (modelName == "YNGKP_M8")
        model.reset(new YNGKP_M8(pgc.release(), codonFreqs.release(), nbClasses));
    }
    else
      throw Exception("Unknown Codon model: " + modelName);
  }


  // //////////////////////////////////
  // gBGC
  // //////////////////////////////////

  else if (modelName == "gBGC")
  {
    if (!(alphabetCode_ & NUCLEOTIDE))
      throw Exception("BppOSubstitutionModelFormat::read. Nucleotide alphabet not supported.");
    // We have to parse the nested model first:
    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'model' for model 'gBGC'.");
    if (verbose_)
      ApplicationTools::displayResult("Biased gene conversion", modelName);
    BppOSubstitutionModelFormat nestedReader(NUCLEOTIDE, true, true, false, verbose_);
    auto_ptr<NucleotideSubstitutionModel> nestedModel(dynamic_cast<NucleotideSubstitutionModel*>(nestedReader.read(alphabet, nestedModelDescription, data, false)));
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    // Now we create the gBGC substitution model:
    model.reset(new gBGC(dynamic_cast<const NucleicAlphabet*>(alphabet), nestedModel.release()));

    // Then we update the parameter set:
    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedArguments_["gBGC." + it->first] = it->second;
    }
  }

  // //////////////////////////////////
  // YpR
  // //////////////////////////////////

  else if (modelName == "YpR_Sym")
  {
    if (!(alphabetCode_ & NUCLEOTIDE))
      throw Exception("BppOSubstitutionModelFormat::read. Nucleotide alphabet not supported.");
    if (alphabet->getAlphabetType() != "RNY alphabet")
      throw Exception("Mismatch alphabet: " + alphabet->getAlphabetType() + " for model: " + modelName);
    const RNY* prny = dynamic_cast<const RNY*>(alphabet);

    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'model' for model 'YpR_sym'.");
    if (verbose_)
      ApplicationTools::displayResult("Symetric YpR model", modelName);
    BppOSubstitutionModelFormat nestedReader(NUCLEOTIDE, false, false, false, verbose_);
    auto_ptr<NucleotideSubstitutionModel> nestedModel(dynamic_cast<NucleotideSubstitutionModel*>(nestedReader.read(&prny->getLetterAlphabet(), nestedModelDescription, data, false)));
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedArguments_["YpR_Sym." + it->first] = it->second;
    }

    model.reset(new YpR_Sym(prny, nestedModel.release()));
  }
  else if (modelName == "YpR_Gen")
  {
    if (!(alphabetCode_ & NUCLEOTIDE))
      throw Exception("BppOSubstitutionModelFormat::read. Nucleotide alphabet not supported.");
    if (alphabet->getAlphabetType() != "RNY alphabet")
      throw Exception("Mismatch alphabet: " + alphabet->getAlphabetType() + " for model: " + modelName);
    const RNY* prny = dynamic_cast<const RNY*>(alphabet);

    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'model' for model 'YpR_gen'.");
    if (verbose_)
      ApplicationTools::displayResult("General YpR model", modelName);
    BppOSubstitutionModelFormat nestedReader(NUCLEOTIDE, false, false, false, verbose_);
    auto_ptr<NucleotideSubstitutionModel> nestedModel(dynamic_cast<NucleotideSubstitutionModel*>(nestedReader.read(&prny->getLetterAlphabet(), nestedModelDescription, data, false)));
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedArguments_["YpR_Gen." + it->first] = it->second;
    }

    model.reset(new YpR_Gen(prny, nestedModel.release()));
  }


  // /////////////////////////////////
  // / RE08
  // ///////////////////////////////

  else if (modelName == "RE08")
  {
    if (!allowGaps_)
      throw Exception("BppOSubstitutionModelFormat::read. No Gap model allowed here.");

    // We have to parse the nested model first:
    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'model' for model 'RE08'.");
    if (verbose_)
      ApplicationTools::displayResult("Gap model", modelName);
    BppOSubstitutionModelFormat nestedReader(ALL, allowCovarions_, false, false, verbose_);
    auto_ptr<ReversibleSubstitutionModel> nestedModel(dynamic_cast<ReversibleSubstitutionModel*>(nestedReader.read(alphabet, nestedModelDescription, data, false)));
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    // Now we create the RE08 substitution model:
    model.reset(new RE08(nestedModel.release()));

    // Then we update the parameter set:
    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedArguments_["RE08.model_" + it->first] = it->second;
    }
  }

  // /////////////////////////////////
  // / TS98
  // ///////////////////////////////

  else if (modelName == "TS98")
  {
    if (!allowCovarions_)
      throw Exception("BppOSubstitutionModelFormat::read. No Covarion model allowed here.");

    // We have to parse the nested model first:
    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'model' for model 'TS98'.");
    if (verbose_)
      ApplicationTools::displayResult("Covarion model", modelName);
    BppOSubstitutionModelFormat nestedReader(ALL, false, allowMixed_, allowGaps_, false);
    auto_ptr<ReversibleSubstitutionModel> nestedModel(dynamic_cast<ReversibleSubstitutionModel*>(nestedReader.read(alphabet, nestedModelDescription, data, false)));
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    // Now we create the TS98 substitution model:
    model.reset(new TS98(nestedModel.release()));

    // Then we update the parameter set:
    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedArguments_["TS98.model_" + it->first] = it->second;
    }
  }

  // /////////////////////////////////
  // / G01
  // ///////////////////////////////

  else if (modelName == "G01")
  {
    if (!allowCovarions_)
      throw Exception("BppOSubstitutionModelFormat::read. No Covarion model allowed here.");

    // We have to parse the nested model first:
    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'model' for model 'G01'.");
    string nestedRateDistDescription = args["rdist"];
    if (TextTools::isEmpty(nestedRateDistDescription))
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'rdist' for model 'G01'.");
    if (verbose_)
      ApplicationTools::displayResult("Covarion model", modelName);
    BppOSubstitutionModelFormat nestedReader(ALL, false, allowMixed_, allowGaps_, verbose_);
    auto_ptr<ReversibleSubstitutionModel> nestedModel(dynamic_cast<ReversibleSubstitutionModel*>(nestedReader.read(alphabet, nestedModelDescription, data, false)));
    map<string, string> unparsedParameterValuesNestedModel(nestedReader.getUnparsedArguments());
    BppORateDistributionFormat rateReader(false);
    auto_ptr<DiscreteDistribution> nestedRDist(rateReader.read(nestedRateDistDescription, false));
    map<string, string> unparsedParameterValuesNestedDist(rateReader.getUnparsedArguments());

    // Now we create the TS98 substitution model:
    model.reset(new G2001(nestedModel.release(), nestedRDist.release()));

    // Then we update the parameter set:
    for (map<string, string>::iterator it = unparsedParameterValuesNestedModel.begin(); it != unparsedParameterValuesNestedModel.end(); it++)
    {
      unparsedArguments_["G01.model_" + it->first] = it->second;
    }
    for (map<string, string>::iterator it = unparsedParameterValuesNestedDist.begin(); it != unparsedParameterValuesNestedDist.end(); it++)
    {
      unparsedArguments_["G01.rdist_" + it->first] = it->second;
    }
  }
  else
  {
    // This is a 'simple' model...
    if (AlphabetTools::isNucleicAlphabet(alphabet))
    {
      if (!(alphabetCode_ & NUCLEOTIDE))
        throw Exception("BppOSubstitutionModelFormat::read. Nucleotide alphabet not supported.");
      const NucleicAlphabet* alpha = dynamic_cast<const NucleicAlphabet*>(alphabet);

      // /////////////////////////////////
      // / GTR
      // ///////////////////////////////

      if (modelName == "GTR")
      {
        model.reset(new GTR(alpha));
      }


      // /////////////////////////////////
      // / SSR
      // ///////////////////////////////

      else if (modelName == "SSR")
      {
        model.reset(new SSR(alpha));
      }

      // /////////////////////////////////
      // / L95
      // ///////////////////////////////

      else if (modelName == "L95")
      {
        model.reset(new L95(alpha));
      }

      // /////////////////////////////////
      // / RN95
      // ///////////////////////////////

      else if (modelName == "RN95")
      {
        model.reset(new RN95(alpha));
      }

      // /////////////////////////////////
      // / RN95s
      // ///////////////////////////////

      else if (modelName == "RN95s")
      {
        model.reset(new RN95s(alpha));
      }

      // /////////////////////////////////
      // / TN93
      // //////////////////////////////

      else if (modelName == "TN93")
      {
        model.reset(new TN93(alpha));
      }

      // /////////////////////////////////
      // / HKY85
      // ///////////////////////////////

      else if (modelName == "HKY85")
      {
        model.reset(new HKY85(alpha));
      }

      // /////////////////////////////////
      // / F84
      // ///////////////////////////////

      else if (modelName == "F84")
      {
        model.reset(new F84(alpha));
      }

      // /////////////////////////////////
      // / T92
      // ///////////////////////////////

      else if (modelName == "T92")
      {
        model.reset(new T92(alpha));
      }

      // /////////////////////////////////
      // / K80
      // ///////////////////////////////

      else if (modelName == "K80")
      {
        model.reset(new K80(alpha));
      }


      // /////////////////////////////////
      // / JC69
      // ///////////////////////////////

      else if (modelName == "JC69")
      {
        model.reset(new JCnuc(alpha));
      }
      else
      {
        throw Exception("Model '" + modelName + "' unknown.");
      }
    }
    else
    {
      if (!(alphabetCode_ & PROTEIN))
        throw Exception("BppOSubstitutionModelFormat::read. Protein alphabet not supported.");
      const ProteicAlphabet* alpha = dynamic_cast<const ProteicAlphabet*>(alphabet);

      if (modelName.find("+F")!=string::npos){
        string freqOpt = ApplicationTools::getStringParameter("frequencies", args, "Full", "", true, verbose_);
        BppOFrequenciesSetFormat freqReader(BppOFrequenciesSetFormat::ALL, false);
        auto_ptr<FrequenciesSet> protFreq(freqReader.read(alpha, freqOpt, data, true));
        map<string, string> unparsedParameterValuesNested(freqReader.getUnparsedArguments());

        for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
          {
            unparsedArguments_[modelName + "." + it->first] = it->second;
          }
        
        if (modelName == "JC69+F")
          model.reset(new JCprot(alpha, dynamic_cast<ProteinFrequenciesSet*>(protFreq.release()), true));
        else if (modelName == "DSO78+F")
          model.reset(new DSO78(alpha, dynamic_cast<ProteinFrequenciesSet*>(protFreq.release()), true));
        else if (modelName == "JTT92+F")
          model.reset(new JTT92(alpha, dynamic_cast<ProteinFrequenciesSet*>(protFreq.release()), true));
        else if (modelName == "LG08+F")
          model.reset(new LG08(alpha, dynamic_cast<ProteinFrequenciesSet*>(protFreq.release()), true));
        else if (modelName == "WAG01+F")
          model.reset(new WAG01(alpha, dynamic_cast<ProteinFrequenciesSet*>(protFreq.release()), true));
        else if (modelName == "Empirical+F")
          {
            string prefix = args["name"];
            if (TextTools::isEmpty(prefix))
              throw Exception("'name' argument missing for user-defined substitution model.");
            model.reset(new UserProteinSubstitutionModel(alpha, args["file"], dynamic_cast<ProteinFrequenciesSet*>(protFreq.release()), prefix + "+F.", true));
          }
      }
      else if (modelName == "JC69")
        model.reset(new JCprot(alpha));
      else if (modelName == "DSO78")
        model.reset(new DSO78(alpha));
      else if (modelName == "JTT92")
        model.reset(new JTT92(alpha));
      else if (modelName == "LG08")
        model.reset(new LG08(alpha));
      else if (modelName == "WAG01")
        model.reset(new WAG01(alpha));
      else if (modelName == "LLG08_EHO")
        model.reset(new LLG08_EHO(alpha));
      else if (modelName == "LLG08_EX2")
        model.reset(new LLG08_EX2(alpha));
      else if (modelName == "LLG08_EX3")
        model.reset(new LLG08_EX3(alpha));
      else if (modelName == "LLG08_UL2")
        model.reset(new LLG08_UL2(alpha));
      else if (modelName == "LLG08_UL3")
        model.reset(new LLG08_UL3(alpha));
      else if (modelName == "LGL08_CAT")
      {
        unsigned int nbCat = TextTools::toInt(args["nbCat"]);
        model.reset(new LGL08_CAT(alpha, nbCat));
      }
      else if (modelName == "Empirical")
      {
        string prefix = args["name"];
        if (TextTools::isEmpty(prefix))
          throw Exception("'name' argument missing for user-defined substitution model.");
        model.reset(new UserProteinSubstitutionModel(alpha, args["file"], prefix));
      }
      else if (modelName == "COaLA")
      {
        string nbrOfParametersPerBranch = args["nbrAxes"];
        if (TextTools::isEmpty(nbrOfParametersPerBranch))
          throw Exception("'nbrAxes' argument missing to define the number of axis of the Correspondence Analysis.");
        string exchangeability = args["exch"];
        if (TextTools::isEmpty(exchangeability))
          throw Exception("'exch' argument missing to define exchangeability.");
        string file = args["file"];
        if (exchangeability == "Empirical" && TextTools::isEmpty(file))
          throw Exception("'file' argument missing to specify the exchangeabilities of the user-defined empirical model.");
        model.reset(new Coala(alpha, TextTools::toInt(nbrOfParametersPerBranch), exchangeability, file));
        dynamic_cast<Coala*>(model.get())->setFreqFromData(*data);
      }

      else if (AlphabetTools::isBinaryAlphabet(alphabet))
      {
        if (!(alphabetCode_ & BINARY))
          throw Exception("BppOSubstitutionModelFormat::read. Binary alphabet not supported.");
        const BinaryAlphabet* balpha = dynamic_cast<const BinaryAlphabet*>(alphabet);

        if (modelName == "Binary")
          model.reset(new BinarySubstitutionModel(balpha));
      }
      else
        throw Exception("Model '" + modelName + "' unknown.");
    }
    if (verbose_)
      ApplicationTools::displayResult("Substitution model", modelName);
  }

  // Update parameter args:
  vector<string> pnames = model->getParameters().getParameterNames();

  string pref = model->getNamespace();

  for (size_t i = 0; i < pnames.size(); i++)
  {
    string name = model->getParameterNameWithoutNamespace(pnames[i]);
    if (args.find(name) != args.end())
      unparsedArguments_[pref + name] = args[name];
  }

  // Now look if some parameters are aliased:
  ParameterList pl = model->getIndependentParameters();
  string pname, pval, pname2;
  for (size_t i = 0; i < pl.size(); i++)
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
        if (verbose_)
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
    unparsedArguments_[pref + "initFreqs"] = args["initFreqs"];
  if (args.find("initFreqs.observedPseudoCount") != args.end())
    unparsedArguments_[pref + "initFreqs.observedPseudoCount"] = args["initFreqs.observedPseudoCount"];

  if (parseArguments)
    initialize_(*model, data);

  return model.release();
}


SubstitutionModel* BppOSubstitutionModelFormat::readWord_(const Alphabet* alphabet, const std::string& modelDescription, const SiteContainer* data)
{
  auto_ptr<SubstitutionModel> model;
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

  if (args.find("model") != args.end())
  {
    v_nestedModelDescription.push_back(args["model"]);
    nbmodels = (modelName == "Word") ? pWA->getLength() : 3;
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


  if (v_nestedModelDescription.size() != nbmodels)
  {
    BppOSubstitutionModelFormat nestedReader(alphabetCode_, false, true, false, false);
    model.reset(nestedReader.read(pWA->getNAlphabet(0), v_nestedModelDescription[0], data, false));
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
    string pref = "";
    for (unsigned int i = 0; i < nbmodels; i++)
    {
      pref += TextTools::toString(i + 1);
    }

    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedArguments_[modelName + "." + pref + "_" + it->first] = it->second;
    }

    v_pSM.push_back(model.release());
  }
  else
  {
    for (unsigned i = 0; i < v_nestedModelDescription.size(); i++)
    {
      BppOSubstitutionModelFormat nestedReader(alphabetCode_, false, true, false, false);
      model.reset(nestedReader.read(pWA->getNAlphabet(i), v_nestedModelDescription[i], data, false));
      map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_[modelName + "." + TextTools::toString(i + 1) + "_" + it->first] = it->second;
      }

      v_pSM.push_back(model.release());
    }
  }

  // /////////////////////////////////
  // / WORD
  // ///////////////////////////////

  if (modelName == "Word")
  {
    model.reset((v_nestedModelDescription.size() != nbmodels)
                ? new WordSubstitutionModel(v_pSM[0], nbmodels)
                : new WordSubstitutionModel(v_pSM));
  }

  // /////////////////////////////////
  // / CODON
  // ///////////////////////////////

  else
  {
    const CodonAlphabet* pCA = dynamic_cast<const CodonAlphabet*>(pWA);
    if (pCA == 0)
      throw Exception("Non codon Alphabet fo model" + modelName + " model.");

    auto_ptr< AlphabetIndex2 > pai2;
    auto_ptr<GeneticCode> pgc;
    auto_ptr<FrequenciesSet> pFS;

    if ((dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]) == 0) ||
        ((v_nestedModelDescription.size() == 3) &&
         (dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]) == 0 || dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]) == 0)))
      throw Exception("Non simple NucleotideSubstitutionModel imbedded in " + modelName + " model.");


    if (modelName.find("Dist") != string::npos)
    {
      if (args.find("genetic_code") == args.end())
        args["genetic_code"] = pCA->getAlphabetType();

      pgc.reset(SequenceApplicationTools::getGeneticCode(dynamic_cast<const NucleicAlphabet*>(pCA->getNAlphabet(0)), args["genetic_code"]));
      if (pgc->getSourceAlphabet()->getAlphabetType() != pCA->getAlphabetType())
        throw Exception("Mismatch between genetic code and codon alphabet");

      pai2.reset((args.find("aadistance") == args.end()) ? 0 : SequenceApplicationTools::getAlphabetIndex2(&AlphabetTools::PROTEIN_ALPHABET, args["aadistance"]));
    }


    if (modelName.find("Freq") != string::npos)
    {
      if (args.find("frequencies") == args.end())
        throw Exception("Missing equilibrium frequencies.");

      BppOFrequenciesSetFormat bIOFreq(alphabetCode_, verbose_);
      pFS.reset(bIOFreq.read(pCA, args["frequencies"], data, false));
      map<string, string> unparsedParameterValuesNested(bIOFreq.getUnparsedArguments());

      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_[modelName + "." + it->first] = it->second;
      }
    }


    // //

    if (modelName == "Triplet")
      model.reset((v_nestedModelDescription.size() != 3)
                  ? new TripletSubstitutionModel(
                    pCA,
                    dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]))
                  : new TripletSubstitutionModel(
                    pCA,
                    dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                    dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                    dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2])));

    else if (modelName == "CodonRate")
      model.reset((v_nestedModelDescription.size() != 3)
                  ? new CodonRateSubstitutionModel(
                    pCA,
                    dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]))
                  : new CodonRateSubstitutionModel(
                    pCA,
                    dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                    dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                    dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2])));


    else if (modelName == "CodonDist")
    {
      if (v_nestedModelDescription.size() != 3)
        model.reset(new CodonDistanceSubstitutionModel(pgc.release(), dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]), pai2.release()));
      else
        model.reset(new CodonDistanceSubstitutionModel(
                      pgc.release(),
                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]), pai2.release()));
    }

    else if (modelName == "CodonRateFreq")
    {
      if (v_nestedModelDescription.size() != 3)
        model.reset(new CodonRateFrequenciesSubstitutionModel(pCA, dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]), pFS.release()));
      else
        model.reset(new CodonRateFrequenciesSubstitutionModel(pCA,
                                                              dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                              dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                                                              dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]),
                                                              pFS.release()));
    }

    else if (modelName == "CodonDistFreq")
    {
      if (v_nestedModelDescription.size() != 3)
        model.reset(new CodonDistanceFrequenciesSubstitutionModel(pgc.release(),
                                                                  dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                                  pFS.release(),
                                                                  pai2.release()));
      else
        model.reset(new CodonDistanceFrequenciesSubstitutionModel(
                      pgc.release(),
                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]),
                      pFS.release(),
                      pai2.release()));
    }

    else if (modelName == "CodonDistPhasFreq")
    {
      if (v_nestedModelDescription.size() != 3)
        model.reset(new CodonDistancePhaseFrequenciesSubstitutionModel(pgc.release(),
                                                                       dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                                       pFS.release(),
                                                                       pai2.release()));
      else
        model.reset(new CodonDistancePhaseFrequenciesSubstitutionModel(
                      pgc.release(),
                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]),
                      pFS.release(),
                      pai2.release()));
    }
    else if (modelName == "CodonDistFitPhasFreq")
    {
      if (args.find("fitness") == args.end())
        throw Exception("Missing fitness in model " + modelName + ".");

      BppOFrequenciesSetFormat bIOFreq(alphabetCode_, verbose_);
      auto_ptr<FrequenciesSet> pFit(bIOFreq.read(pCA, args["fitness"], data, false));
      map<string, string> unparsedParameterValuesNested(bIOFreq.getUnparsedArguments());

      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_[modelName + ".fit_" + it->first] = it->second;
      }

      if (v_nestedModelDescription.size() != 3)
      {
        model.reset(new CodonDistanceFitnessPhaseFrequenciesSubstitutionModel(pgc.release(),
                                                                              dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                                              pFit.release(),
                                                                              pFS.release(),
                                                                              pai2.release()));
      }
      else
        model.reset(new CodonDistanceFitnessPhaseFrequenciesSubstitutionModel(
                      pgc.release(),
                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]),
                      pFit.release(),
                      pFS.release(),
                      pai2.release()));
    }
  }
  return model.release();
}


MixedSubstitutionModel* BppOSubstitutionModelFormat::readMixed_(const Alphabet* alphabet, const std::string& modelDescription, const SiteContainer* data)
{
  auto_ptr<MixedSubstitutionModel> model;

  string modelName = "";
  map<string, string> args;
  KeyvalTools::parseProcedure(modelDescription, modelName, args);
  auto_ptr<SubstitutionModel> pSM;

  if (modelName == "MixedModel")
  {
    if (args.find("model") == args.end())
      throw Exception("The argument 'model' is missing from MixedSubstitutionModel description");
    string nestedModelDescription = args["model"];
    BppOSubstitutionModelFormat nestedReader(alphabetCode_, allowCovarions_, true, allowGaps_, false);
    pSM.reset(nestedReader.read(alphabet, nestedModelDescription, data, false));
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    map<string, DiscreteDistribution*> mdist;
    map<string, string> unparsedParameterValuesNested2;

    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin();
         it != unparsedParameterValuesNested.end();
         it++)
    {
      if (it->second.find("(") != string::npos)
      {
        BppODiscreteDistributionFormat bIO(false);
        mdist[pSM->getParameterNameWithoutNamespace(it->first)] = bIO.read(it->second, false);
        map<string, string> unparsedParameterValuesNested3(bIO.getUnparsedArguments());
        for (map<string, string>::iterator it2 = unparsedParameterValuesNested3.begin();
             it2 != unparsedParameterValuesNested3.end();
             it2++)
        {
          unparsedParameterValuesNested2[it->first + "_" + it2->first] = it2->second;
        }
      }
      else
        unparsedParameterValuesNested2[it->first] = it->second;
    }

    for (map<string, string>::iterator it = unparsedParameterValuesNested2.begin();
         it != unparsedParameterValuesNested2.end();
         it++)
    {
      unparsedArguments_[it->first] = it->second;
    }


    int fi(-1), ti(-1);

    if (args.find("from") != args.end())
      fi = alphabet->charToInt(args["from"]);
    if (args.find("to") != args.end())
      ti = alphabet->charToInt(args["to"]);

    string sModN=pSM->getName();
    model.reset(new MixtureOfASubstitutionModel(alphabet, pSM.release(), mdist, fi, ti));

    vector<string> v = model->getParameters().getParameterNames();

    for (map<string, DiscreteDistribution*>::iterator it = mdist.begin();
         it != mdist.end(); it++)
    {
      delete it->second;
    }

    if (verbose_)
      {
        ApplicationTools::displayResult("Mixture Of A Substitution Model", sModN);
        ApplicationTools::displayResult("Number of classes", model->getNumberOfModels());
      }
  }


  else if (modelName == "Mixture")
  {
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
      BppOSubstitutionModelFormat nestedReader(alphabetCode_, false, true, false, false);
      pSM.reset(nestedReader.read(alphabet, v_nestedModelDescription[i], data, false));
      map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_[modelName + "." + TextTools::toString(i + 1) + "_" + it->first] = it->second;
      }
      v_pSM.push_back(pSM.release());
    }

    model.reset(new MixtureOfSubstitutionModels(alphabet, v_pSM));
    if (verbose_)
        ApplicationTools::displayResult("Mixture Of Substitution Models", modelName );
  }
  else
    throw Exception("Unknown model name for mixture " + modelName);

  return model.release();
}


void BppOSubstitutionModelFormat::write(const SubstitutionModel& model,
                                        OutputStream& out,
                                        std::map<std::string, std::string>& globalAliases,
                                        std::vector<std::string>& writtenNames) const
{
  bool comma = false;

  //  Mixed Model that are defined as "Mixture" and "Mixed"

  if ((dynamic_cast<const MixedSubstitutionModel*>(&model) != NULL) && (dynamic_cast<const AbstractBiblioMixedSubstitutionModel*>(&model) == NULL))
  {
    writeMixed_(*dynamic_cast<const MixedSubstitutionModel*>(&model), out, globalAliases, writtenNames);
    return;
  }

  out << model.getName() + "(";

  // Is it a protein user defined model?
  const UserProteinSubstitutionModel* userModel = dynamic_cast<const UserProteinSubstitutionModel*>(&model);
  if (userModel)
  {
    out << "file=" << userModel->getPath();
    comma = true;
  }

  // Is it a markov-modulated model?
  const MarkovModulatedSubstitutionModel* mmModel = dynamic_cast<const MarkovModulatedSubstitutionModel*>(&model);
  if (mmModel)
  {
    out << "model=";
    const SubstitutionModel* nestedModel = mmModel->getNestedModel();
    write(*nestedModel, out, globalAliases, writtenNames);

    const G2001* gModel = dynamic_cast<const G2001*>(&model);
    if (gModel)
    {
      // Also print distribution here:
      out << ",rdist=";
      const DiscreteDistribution* nestedDist = gModel->getRateDistribution();
      const BppODiscreteDistributionFormat* bIO = new BppODiscreteDistributionFormat();

      bIO->write(*nestedDist, out, globalAliases, writtenNames);
      delete bIO;
    }
    comma = true;
  }

  // Is it a model with gaps?
  const RE08* reModel = dynamic_cast<const RE08*>(&model);
  if (reModel)
  {
    out << "model=";
    const SubstitutionModel* nestedModel = reModel->getNestedModel();
    write(*nestedModel, out, globalAliases, writtenNames);
    comma = true;
  }

  // Is it a YpR model?
  const YpR* yprModel = dynamic_cast<const YpR*>(&model);
  if (yprModel)
  {
    out << "model=";
    const SubstitutionModel* nestedModel = yprModel->getNestedModel();
    write(*nestedModel, out, globalAliases, writtenNames);
    comma = true;
  }

  // Is it a word model?

  const AbstractWordSubstitutionModel* wM = dynamic_cast<const AbstractWordSubstitutionModel*>(&model);
  if (wM)
  {
    size_t nmod = wM->getNumberOfModels();
    const SubstitutionModel* mod0 = wM->getNModel(0);
    if (nmod == 1)
    {
      out << "model=";
      write(*mod0, out, globalAliases, writtenNames);
    }
    else
    {
      const SubstitutionModel* mod1 = wM->getNModel(1);
      if (mod1 == mod0)
      {
        out << "model=";
        write(*mod0, out, globalAliases, writtenNames);
      }
      else
      {
        out << "model1=";
        write(*mod0, out, globalAliases, writtenNames);
        for (unsigned int i = 1; i < nmod; i++)
        {
          out << ",model" + TextTools::toString(i + 1) + "=";
          write(*wM->getNModel(i), out, globalAliases, writtenNames);
        }
      }
    }
    comma = true;
  }

  // Is it a COaLA model ?
  const Coala* coalaModel = dynamic_cast<const Coala*>(&model);
  if (coalaModel)
  {
    out << "exch=" << coalaModel->getExch() << ",nbrAxes=" << coalaModel->getNbrOfAxes();
    comma = true;
  }

  // Regular model
  const FrequenciesSet* pfs = model.getFrequenciesSet();
  if (pfs)
  {
    if (comma)
      out << ",";
    out << "frequencies=";

    BppOFrequenciesSetFormat bIOFreq(alphabetCode_, false);
    bIOFreq.write(pfs, out, writtenNames);
    comma = true;
  }

  // Specific case of CodonFitnessSubstitutionModel

  const CodonDistanceFitnessPhaseFrequenciesSubstitutionModel* pCF = dynamic_cast<const CodonDistanceFitnessPhaseFrequenciesSubstitutionModel*>(&model);
  if (pCF)
  {
    if (comma)
      out << ",";
    out << "fitness=";

    BppOFrequenciesSetFormat bIOFreq(alphabetCode_, false);
    bIOFreq.write(pCF->getFitness(), out, writtenNames);
    comma = true;
  }

  // and the other parameters

  BppOParametrizableFormat bIO;
  bIO.write(&model, out, globalAliases, model.getIndependentParameters().getParameterNames(), writtenNames, true, comma);
  out << ")";
}


void BppOSubstitutionModelFormat::writeMixed_(const MixedSubstitutionModel& model,
                                              OutputStream& out,
                                              std::map<std::string, std::string>& globalAliases,
                                              std::vector<std::string>& writtenNames) const
{
  if (dynamic_cast<const MixtureOfSubstitutionModels*>(&model) != NULL)
  {
    const MixtureOfSubstitutionModels* pMS = dynamic_cast<const MixtureOfSubstitutionModels*>(&model);

    for (unsigned int i = 0; i < pMS->getNumberOfModels(); i++)
    {
      const SubstitutionModel* eM = pMS->getNModel(i);

      vector<string> vpl = eM->getIndependentParameters().getParameterNames();
      for (unsigned j = 0; j < vpl.size(); j++)
      {
        if (eM->getParameterNameWithoutNamespace(vpl[j]) == "rate")
          writtenNames.push_back(vpl[j]);
      }
    }

    out << "Mixture(";
    for (unsigned int i = 0; i < pMS->getNumberOfModels(); i++)
    {
      if (i != 0)
        out << ", ";
      out << "model" + TextTools::toString(i + 1) + "=";
      write(*pMS->getNModel(i), out, globalAliases, writtenNames);
    }
  }
  else
  {
    const MixtureOfASubstitutionModel* pMS = dynamic_cast<const MixtureOfASubstitutionModel*>(&model);
    out << "MixedModel(model=";
    const SubstitutionModel* eM = pMS->getNModel(0);

    ParameterList pl = eM->getIndependentParameters();
    vector<string> vpl = pl.getParameterNames();

    for (unsigned j = 0; j < vpl.size(); j++)
    {
      if (find(writtenNames.begin(), writtenNames.end(), vpl[j]) == writtenNames.end())
      {
        if (eM->getParameterNameWithoutNamespace(vpl[j]) == "rate")
          writtenNames.push_back(vpl[j]);
        else
        {
          const DiscreteDistribution* pDD = pMS->getDistribution(vpl[j]);
          if (pDD && dynamic_cast<const ConstantDistribution*>(pDD) == NULL)
          {
            const BppODiscreteDistributionFormat* bIO = new BppODiscreteDistributionFormat();
            StdStr sout;
            bIO->write(*pDD, sout, globalAliases, writtenNames);
            globalAliases[vpl[j]] = sout.str();
            delete bIO;
          }
        }
      }
    }

    write(*eM, out, globalAliases, writtenNames);

    if (pMS->from() != -1)
      out << ",from=" << model.getAlphabet()->intToChar(pMS->from()) << ",to=" << model.getAlphabet()->intToChar(pMS->to());
  }

  const BppOParametrizableFormat* bIO = new BppOParametrizableFormat();
  bIO->write(&model, out, globalAliases, model.getIndependentParameters().getParameterNames(), writtenNames, true, true);
  delete bIO;

  out << ")";
}

void BppOSubstitutionModelFormat::initialize_(
  SubstitutionModel& model,
  const SiteContainer* data) throw (Exception)
{
  string initFreqs = ApplicationTools::getStringParameter(model.getNamespace() + "initFreqs", unparsedArguments_, "", "", true, false);
  if (verbose_)
    ApplicationTools::displayResult("External frequencies initialization for model", (initFreqs == "") ? "None" : initFreqs);

  if (initFreqs != "")
  {
    if (initFreqs == "observed")
    {
      if (!data)
        throw Exception("BppOSubstitutionModelFormat::initialize_(). Missing data for observed frequencies");
      unsigned int psi = ApplicationTools::getParameter<unsigned int>(model.getNamespace() + "initFreqs.observedPseudoCount", unparsedArguments_, 0);
      model.setFreqFromData(*data, psi);
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
      model.setFreq(frequencies);
    }
    else
      throw Exception("Unknown initFreqs argument");
  }

  ParameterList pl = model.getIndependentParameters();
  for (size_t i = 0; i < pl.size(); i++)
  {
    AutoParameter ap(pl[i]);
    ap.setMessageHandler(ApplicationTools::warning);
    pl.setParameter(i, ap);
  }
  size_t posp;
  for (size_t i = 0; i < pl.size(); i++)
  {
    const string pName = pl[i].getName();
    posp = pName.rfind(".");
    bool test1 = (initFreqs == "");
    bool test2 = (model.getParameterNameWithoutNamespace(pName).size() < posp + 6) || (model.getParameterNameWithoutNamespace(pName).substr(posp + 1, 5) != "theta");
    bool test3 = (unparsedArguments_.find(pName) != unparsedArguments_.end());
    if (test1 || test2 || test3)
    {
      if (!test1 && !test2 && test3)
        ApplicationTools::displayWarning("Warning, initFreqs argument is set and a value is set for parameter " + pName);
      double value = ApplicationTools::getDoubleParameter(pName, unparsedArguments_, pl[i].getValue());
      pl[i].setValue(value);
    }
    if (verbose_)
      ApplicationTools::displayResult("Parameter found", pName + "=" + TextTools::toString(pl[i].getValue()));
  }

  model.matchParametersValues(pl);
}


