// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Io/BppODiscreteDistributionFormat.h>
#include <Bpp/Io/BppOParametrizableFormat.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Io/OutputStream.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/TextTools.h>

#include "../App/PhylogeneticsApplicationTools.h"
#include "../Model/AbstractBiblioMixedTransitionModel.h"
#include "../Model/BinarySubstitutionModel.h"
#include "../Model/D1WalkSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonAAFitnessSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonAARateSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonBGCSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonClusterAASubstitutionModel.h"
#include "../Model/Codon/AbstractCodonCpGSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonDistanceSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonDistanceSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonFitnessSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonFrequenciesSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonPhaseFrequenciesSubstitutionModel.h"
#include "../Model/Codon/CodonAdHocSubstitutionModel.h"
#include "../Model/Codon/CodonSameAARateSubstitutionModel.h"
#include "../Model/Codon/DFP07.h"
#include "../Model/Codon/GY94.h"
#include "../Model/Codon/KCM.h"
#include "../Model/Codon/KroneckerCodonDistanceFrequenciesSubstitutionModel.h"
#include "../Model/Codon/KroneckerCodonDistanceSubstitutionModel.h"
#include "../Model/Codon/MG94.h"
#include "../Model/Codon/SENCA.h"
#include "../Model/Codon/TripletSubstitutionModel.h"
#include "../Model/Codon/YN98.h"
#include "../Model/Codon/YNGP_M10.h"
#include "../Model/Codon/YNGP_M7.h"
#include "../Model/Codon/YNGP_M8.h"
#include "../Model/Codon/YNGP_M9.h"
#include "../Model/EquiprobableSubstitutionModel.h"
#include "../Model/FromMixtureSubstitutionModel.h"
#include "../Model/G2001.h"
#include "../Model/InMixedSubstitutionModel.h"
#include "../Model/IntegrationOfSubstitutionModel.h"
#include "../Model/KroneckerWordSubstitutionModel.h"
#include "../Model/MixtureOfATransitionModel.h"
#include "../Model/MixtureOfTransitionModels.h"
#include "../Model/Nucleotide/F81.h"
#include "../Model/Nucleotide/F84.h"
#include "../Model/Nucleotide/GTR.h"
#include "../Model/Nucleotide/HKY85.h"
#include "../Model/Nucleotide/JCnuc.h"
#include "../Model/Nucleotide/K80.h"
#include "../Model/Nucleotide/L95.h"
#include "../Model/Nucleotide/NucleotideSubstitutionModel.h"
#include "../Model/Nucleotide/RN95.h"
#include "../Model/Nucleotide/RN95s.h"
#include "../Model/Nucleotide/SSR.h"
#include "../Model/Nucleotide/T92.h"
#include "../Model/Nucleotide/TN93.h"
#include "../Model/Nucleotide/YpR.h"
#include "../Model/Nucleotide/gBGC.h"
#include "../Model/OneChangeRegisterTransitionModel.h"
#include "../Model/OneChangeTransitionModel.h"
#include "../Model/POMO.h"
#include "../Model/Protein/Coala.h"
#include "../Model/Protein/CoalaCore.h"
#include "../Model/Protein/DSO78.h"
#include "../Model/Protein/JCprot.h"
#include "../Model/Protein/JTT92.h"
#include "../Model/Protein/LG08.h"
#include "../Model/Protein/LGL08_CAT.h"
#include "../Model/Protein/ProteinSubstitutionModel.h"
#include "../Model/Protein/UserProteinSubstitutionModel.h"
#include "../Model/Protein/WAG01.h"
#include "../Model/RE08.h"
#include "../Model/RegisterRatesSubstitutionModel.h"
#include "../Model/TS98.h"
#include "../Model/TwoParameterBinarySubstitutionModel.h"
#include "../Model/WordSubstitutionModel.h"
#include "BppOFrequencySetFormat.h"
#include "BppORateDistributionFormat.h"
#include "BppOSubstitutionModelFormat.h"
#include "BppOTransitionModelFormat.h"

// From Numeric

#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/AutoParameter.h>


#include <Bpp/Text/StringTokenizer.h>


using namespace bpp;

// From the STL:
#include <iomanip>
#include <set>
#include <memory>

using namespace std;

unsigned char BppOSubstitutionModelFormat::DNA = 1;
unsigned char BppOSubstitutionModelFormat::RNA = 2;
unsigned char BppOSubstitutionModelFormat::NUCLEOTIDE = 1 | 2;
unsigned char BppOSubstitutionModelFormat::PROTEIN = 4;
unsigned char BppOSubstitutionModelFormat::CODON = 8;
unsigned char BppOSubstitutionModelFormat::WORD = 16;
unsigned char BppOSubstitutionModelFormat::BINARY = 32;
unsigned char BppOSubstitutionModelFormat::INTEGER = 64;
unsigned char BppOSubstitutionModelFormat::ALL = 1 | 2 | 4 | 8 | 16 | 32 | 64;


unique_ptr<SubstitutionModelInterface> BppOSubstitutionModelFormat::readSubstitutionModel(
    shared_ptr<const Alphabet> alphabet,
    const string& modelDescription,
    const std::map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mData,
    size_t nData,
    bool parseArguments)
{
  unparsedArguments_.clear();
  unique_ptr<SubstitutionModelInterface> model;
  string modelName = "";
  map<string, string> args;
  KeyvalTools::parseProcedure(modelDescription, modelName, args);


  // ///////////////////////////////
  // LINKED WITH MIXTURES
  // //////////////////////////////

  if (modelName == "InMixed")
  {
    if (args.find("model") == args.end())
      throw Exception("'model' argument missing to define the InMixed model.");

    string nameMod = "";
    size_t numMod = 0;

    if (args.find("numMod") == args.end())
    {
      if (args.find("nameMod") == args.end())
        throw Exception("'numMod' and 'nameMod' arguments missing to define the InMixed submodel.");
      else
        nameMod = args["nameMod"];
    }
    else
      numMod = (size_t)TextTools::toInt(args["numMod"]);

    string modelDesc2 = args["model"];
    BppOTransitionModelFormat nestedReader(alphabetCode_, false, true, false, false, warningLevel_);
    nestedReader.setGeneticCode(geneticCode_); // This uses the same instance as the o

    unique_ptr<TransitionModelInterface> nestedModelTmp = nestedReader.readTransitionModel(alphabet, modelDesc2, mData, nData, false);
    unique_ptr<MixedTransitionModelInterface> nestedModel(dynamic_cast<MixedTransitionModelInterface*>(nestedModelTmp.release()));

    // Check that nestedModel is fine and has subModel of given name
    // or number

    if (!nestedModel)
      throw Exception("Unknown mixed model " + modelDesc2 + ".");

    if (nameMod != "")
    {
      try
      {
        nestedModel->model(nameMod);
      }
      catch (NullPointerException&)
      {
        throw Exception("BppOSubstitutionModelFormat::read. " + nestedModel->getName() + "argument for model 'InMixed' has no submodel with name " + nameMod + ".");
      }
      model = make_unique<InMixedSubstitutionModel>(std::move(nestedModel), nameMod, modelDesc2);
    }
    else
    {
      if (numMod == 0 || (nestedModel->getNumberOfModels() < numMod))
        throw Exception("BppOSubstitutionModelFormat::read. " + nestedModel->getName() + "argument for model 'InMixed' has no submodel with number " + TextTools::toString(numMod) + ".");
      model = make_unique<InMixedSubstitutionModel>(std::move(nestedModel), numMod - 1, modelDesc2);
    }

    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
    for (auto& it :unparsedParameterValuesNested)
    {
      unparsedArguments_[it.first] = it.second;
    }

    if (verbose_)
      ApplicationTools::displayResult("Number of submodel", TextTools::toString(numMod));
  }


  else if (modelName == "FromRegister")
  {
    // We have to parse the nested model first:
    if (args.find("model") == args.end())
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'model' for model 'FromRegister'.");

    string nestedModelDescription = args["model"];
    BppOSubstitutionModelFormat nestedReader(ALL, false, allowMixed_, allowGaps_, verbose_, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_);

    auto nestedModel = nestedReader.readSubstitutionModel(alphabet, nestedModelDescription, mData, nData, false);
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    // We look for the register:
    if (args.find("register") == args.end())
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'register' for model 'FromRegister'.");

    string registerDescription = args["register"];
    shared_ptr<AlphabetIndex2> weights = nullptr;
    shared_ptr<AlphabetIndex2> distances = nullptr;

    auto reg = PhylogeneticsApplicationTools::getSubstitutionRegister(registerDescription, nestedModel->getStateMap(), geneticCode_, weights, distances, verbose_);

    // is it normalized (default : false)
    bool isNorm = false;

    if (args.find("isNormalized") != args.end() && args["isNormalized"] == "true")
      isNorm = true;

    model = make_unique<RegisterRatesSubstitutionModel>(std::move(nestedModel), *reg, isNorm);

    if (TextTools::isEmpty(registerDescription))
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'register' for model 'FromRegister'.");

    // Then we update the parameter set:
    for (auto& it : unparsedParameterValuesNested)
    {
      unparsedArguments_["FromRegister." + it.first] = it.second;
    }
  }

  else if (modelName == "POMO")
  {
    auto allelic = dynamic_pointer_cast<const AllelicAlphabet>(alphabet);
    if (!allelic)
      throw Exception("BppOSubstitutionModelFormat;;read. POMO model with no allelic alphabet.");

    // We have to parse the nested model first:
    if (args.find("model") == args.end())
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'model' for model 'POMO'.");

    map<string, string> unparsedParameterValuesNested;
    // fitness
    unique_ptr<FrequencySetInterface> nestedFreq(nullptr);

    if (args.find("fitness") != args.end())
    {
      string nestedFreqDescription = args["fitness"];
      BppOFrequencySetFormat nestedFreqReader(ALL, false, warningLevel_);

      nestedFreq = nestedFreqReader.readFrequencySet(allelic->getStateAlphabet(), nestedFreqDescription, mData, nData, false);
      unparsedParameterValuesNested = nestedFreqReader.getUnparsedArguments();

      for (auto& it : unparsedParameterValuesNested)
      {
        unparsedParameterValuesNested["fit_" + it.first] = it.second;
      }
    }

    // model

    string nestedModelDescription = args["model"];
    BppOSubstitutionModelFormat nestedReader(ALL, false, allowMixed_, allowGaps_, verbose_, warningLevel_);
    nestedReader.setGeneticCode(geneticCode_);

    auto nestedModel = nestedReader.readSubstitutionModel(allelic->getStateAlphabet(), nestedModelDescription, mData, nData, false);
    unparsedParameterValuesNested.insert(nestedReader.getUnparsedArguments().begin(),
        nestedReader.getUnparsedArguments().end());


    model = make_unique<POMO>(allelic, std::move(nestedModel), std::move(nestedFreq));

    // Then we update the parameter set:
    for (auto& it : unparsedParameterValuesNested)
    {
      unparsedArguments_["POMO." + it.first] = it.second;
    }
  }


  // /////////////////////////////////
  // / WORDS and CODONS
  // ///////////////////////////////

  else if ((modelName == "Word") || (modelName.substr(0, 4) == "Kron") || (modelName == "Triplet") || (modelName.substr(0, 5) == "Codon") || (modelName == "SENCA") )
    model = readWord_(alphabet, modelDescription, mData, nData);


  // //////////////////
  // PREDEFINED CODON MODELS

  else if (((modelName == "MG94") || (modelName == "YN98") || (modelName == "YNGP_M0") ||
      (modelName == "GY94") ||  (modelName.substr(0, 3) == "KCM") || (modelName == "SameAARate"))
      && (alphabetCode_ & CODON))
  {
    if (!(alphabetCode_ & CODON))
      throw Exception("BppOSubstitutionModelFormat::read. Codon alphabet not supported.");
    if (!geneticCode_)
      throw Exception("BppOSubstitutionModelFormat::readSubstitionModel(). No genetic code specified! Consider using 'setGeneticCode'.");

    if (!AlphabetTools::isCodonAlphabet(*alphabet))
      throw Exception("Alphabet should be Codon Alphabet.");

    auto pCA = dynamic_pointer_cast<const CodonAlphabet>(alphabet);

    if (args.find("genetic_code") != args.end())
    {
      ApplicationTools::displayWarning("'genetic_code' argument is no longer supported inside model description, and has been supersided by a global 'genetic_code' option.");
      throw Exception("BppOSubstitutionModelFormat::read. Deprecated 'genetic_code' argument.");
    }

    if (geneticCode_->getSourceAlphabet()->getAlphabetType() != pCA->getAlphabetType())
      throw Exception("Mismatch between genetic code and codon alphabet");

    unique_ptr<CodonFrequencySetInterface> codonFreqs(nullptr);

    if (args.find("frequencies") != args.end())
    {
      string freqOpt = ApplicationTools::getStringParameter("frequencies", args, "F0", "", true, warningLevel_);
      BppOFrequencySetFormat freqReader(BppOFrequencySetFormat::ALL, false, warningLevel_);
      freqReader.setGeneticCode(geneticCode_); // This uses the same instance as the one that will be used by the model.

      auto tmpFreqs = freqReader.readFrequencySet(pCA, freqOpt, mData, nData);
      codonFreqs = unique_ptr<CodonFrequencySetInterface>(dynamic_cast<CodonFrequencySetInterface*>(tmpFreqs.release()));
      for (const auto& it : freqReader.getUnparsedArguments())
      {
        unparsedArguments_[modelName + "." + it.first] = it.second;
      }
    }
    else
    // codonFreqs compulsory for all models but SameAARate
    if (modelName != "SameAARate")
      throw Exception("Missing 'frequencies' for model " + modelName);

    if (modelName == "MG94")
      model = make_unique<MG94>(geneticCode_, std::move(codonFreqs));
    else if (modelName == "GY94")
      model = make_unique<GY94>(geneticCode_, std::move(codonFreqs));
    else if ((modelName == "YN98") || (modelName == "YNGP_M0"))
      model = make_unique<YN98>(geneticCode_, std::move(codonFreqs));
    else if (modelName == "KCM7")
      model = make_unique<KCM>(geneticCode_, true);
    else if (modelName == "KCM19")
      model = make_unique<KCM>(geneticCode_, false);
    else if (modelName == "SameAARate")
    {
      if (args.find("protmodel") == args.end())
        throw Exception("Missing 'protmodel in model " + modelName + ".");

      BppOSubstitutionModelFormat nestedProtReader(PROTEIN, false, allowMixed_, allowGaps_, verbose_, warningLevel_);
      auto tmpProtModel = nestedProtReader.readSubstitutionModel(geneticCode_->getTargetAlphabet(), args["protmodel"], mData, nData, false);
      auto nestedProtModel = unique_ptr<ProteinSubstitutionModelInterface>(dynamic_cast<ProteinSubstitutionModelInterface*>(tmpProtModel.release()));

      auto unparsedParameterValuesNested  = nestedProtReader.getUnparsedArguments();
      unparsedArguments_.insert(unparsedParameterValuesNested.begin(), unparsedParameterValuesNested.end());

      if (args.find("codonmodel") == args.end())
        throw Exception("Missing 'codonmodel in model " + modelName + ".");

      BppOSubstitutionModelFormat nestedCodonReader(CODON, false, allowMixed_, allowGaps_, verbose_, warningLevel_);
      nestedCodonReader.setGeneticCode(geneticCode_); // This uses the same instance as the o

      auto tmpSubstitutionModel = nestedCodonReader.readSubstitutionModel(alphabet, args["codonmodel"], mData, nData, false);
      auto nestedCodonModel = unique_ptr<CodonSubstitutionModelInterface>(dynamic_cast<CodonSubstitutionModelInterface*>(tmpSubstitutionModel.release()));

      unparsedParameterValuesNested = nestedCodonReader.getUnparsedArguments();
      for (const auto& it:unparsedParameterValuesNested)
      {
        unparsedArguments_["SameAARate." + it.first] = it.second;
      }

      model = make_unique<CodonSameAARateSubstitutionModel>(std::move(nestedProtModel), std::move(nestedCodonModel), std::move(codonFreqs), geneticCode_);
    }
    else
      throw Exception("Unknown Codon model: " + modelName);
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
    auto prny = dynamic_pointer_cast<const RNY>(alphabet);

    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'model' for model 'YpR_sym'.");
    if (verbose_)
      ApplicationTools::displayResult("Symmetric YpR model", modelName);
    BppOSubstitutionModelFormat nestedReader(NUCLEOTIDE, false, false, false, verbose_, warningLevel_);
    auto tmpModel = nestedReader.readSubstitutionModel(prny->getLetterAlphabet(), nestedModelDescription, mData, nData, false);
    unique_ptr<NucleotideSubstitutionModelInterface> nestedModel(dynamic_cast<NucleotideSubstitutionModelInterface*>(tmpModel.release()));
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
    for (auto& it : unparsedParameterValuesNested)
    {
      unparsedArguments_["YpR_Sym." + it.first] = it.second;
    }

    model = make_unique<YpR_Sym>(prny, std::move(nestedModel));
  }
  else if (modelName == "YpR_Gen")
  {
    if (!(alphabetCode_ & NUCLEOTIDE))
      throw Exception("BppOSubstitutionModelFormat::read. Nucleotide alphabet not supported.");
    if (alphabet->getAlphabetType() != "RNY alphabet")
      throw Exception("Mismatch alphabet: " + alphabet->getAlphabetType() + " for model: " + modelName);
    auto prny = dynamic_pointer_cast<const RNY>(alphabet);

    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'model' for model 'YpR_gen'.");
    if (verbose_)
      ApplicationTools::displayResult("General YpR model", modelName);
    BppOSubstitutionModelFormat nestedReader(NUCLEOTIDE, false, false, false, verbose_, warningLevel_);
    auto tmpModel = nestedReader.readSubstitutionModel(prny->getLetterAlphabet(), nestedModelDescription, mData, nData, false);
    unique_ptr<NucleotideSubstitutionModelInterface> nestedModel(dynamic_cast<NucleotideSubstitutionModelInterface*>(tmpModel.release()));
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    for (auto& it : unparsedParameterValuesNested)
    {
      unparsedArguments_["YpR_Gen." + it.first] = it.second;
    }

    model = make_unique<YpR_Gen>(prny, std::move(nestedModel));
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
    BppOSubstitutionModelFormat nestedReader(ALL, allowCovarions_, false, false, verbose_, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_);
    auto tmpModel = nestedReader.readSubstitutionModel(alphabet, nestedModelDescription, mData, nData, false);

    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    // Now we create the RE08 substitution model:
    if (AlphabetTools::isNucleicAlphabet(*alphabet))
    {
      if (!(alphabetCode_ & NUCLEOTIDE))
        throw Exception("BppOSubstitutionModelFormat::read. Nucleic alphabet not supported.");
      unique_ptr<NucleotideReversibleSubstitutionModelInterface> nestedModel(dynamic_cast<NucleotideReversibleSubstitutionModelInterface*>(tmpModel.release()));
      model = make_unique<RE08Nucleotide>(std::move(nestedModel));
      if (!model)
        throw Exception("BppOSubstitutionModelFormat::readSubstitionModel(). Invalid submodel, must be 'reversible' and 'nucleotide'.");
    }
    else if (AlphabetTools::isProteicAlphabet(*alphabet))
    {
      if (!(alphabetCode_ & PROTEIN))
        throw Exception("BppOSubstitutionModelFormat::read. Protein alphabet not supported.");
      unique_ptr<ProteinReversibleSubstitutionModelInterface> nestedModel(dynamic_cast<ProteinReversibleSubstitutionModelInterface*>(tmpModel.release()));
      model = make_unique<RE08Protein>(std::move(nestedModel));
      if (!model)
        throw Exception("BppOSubstitutionModelFormat::readSubstitionModel(). Invalid submodel, must be 'reversible' and 'protein'.");
    }
    else if (AlphabetTools::isCodonAlphabet(*alphabet))
    {
      if (!(alphabetCode_ & CODON))
        throw Exception("BppOSubstitutionModelFormat::read. Codon alphabet not supported.");
      unique_ptr<CodonReversibleSubstitutionModelInterface> nestedModel(dynamic_cast<CodonReversibleSubstitutionModelInterface*>(tmpModel.release()));
      model = make_unique<RE08Codon>(std::move(nestedModel));
      if (!model)
        throw Exception("BppOSubstitutionModelFormat::readSubstitionModel(). Invalid submodel, must be 'reversible' and 'codon'.");
    }
    else
    {
      unique_ptr<ReversibleSubstitutionModelInterface> nestedModel(dynamic_cast<ReversibleSubstitutionModelInterface*>(tmpModel.release()));
      model = make_unique<RE08>(std::move(nestedModel));
    }

    // Then we update the parameter set:
    for (auto& it : unparsedParameterValuesNested)
    {
      unparsedArguments_["RE08.model_" + it.first] = it.second;
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
    BppOSubstitutionModelFormat nestedReader(ALL, false, allowMixed_, allowGaps_, false, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_);
    auto tmpModel = nestedReader.readSubstitutionModel(alphabet, nestedModelDescription, mData, nData, false);
    unique_ptr<ReversibleSubstitutionModelInterface> nestedModel(dynamic_cast<ReversibleSubstitutionModelInterface*>(tmpModel.release()));
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    // Now we create the TS98 substitution model:
    model = make_unique<TS98>(std::move(nestedModel));

    // Then we update the parameter set:
    for (auto& it : unparsedParameterValuesNested)
    {
      unparsedArguments_["TS98.model_" + it.first] = it.second;
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
    BppOSubstitutionModelFormat nestedReader(ALL, false, allowMixed_, allowGaps_, verbose_, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_);
    auto tmpModel = nestedReader.readSubstitutionModel(alphabet, nestedModelDescription, mData, nData, false);
    unique_ptr<ReversibleSubstitutionModelInterface> nestedModel(dynamic_cast<ReversibleSubstitutionModelInterface*>(tmpModel.release()));
    map<string, string> unparsedParameterValuesNestedModel(nestedReader.getUnparsedArguments());

    BppORateDistributionFormat rateReader(false);
    auto nestedRDist = rateReader.readDiscreteDistribution(nestedRateDistDescription, false);
    map<string, string> unparsedParameterValuesNestedDist(rateReader.getUnparsedArguments());

    // Now we create the G01 substitution model:
    model = make_unique<G2001>(std::move(nestedModel), std::move(nestedRDist));

    // Then we update the parameter set:
    for (auto& it : unparsedParameterValuesNestedModel)
    {
      unparsedArguments_["G01.model_" + it.first] = it.second;
    }
    for (auto& it : unparsedParameterValuesNestedDist)
    {
      unparsedArguments_["G01.rdist_" + it.first] = it.second;
    }
  }

  ///////////////////////////////////////////////////////
  ///
  ///
  ///  SIMPLE MODELS
  ///
  ///
  //////////////////////////////////////////////////////////

  else
  {
    // This is a 'simple' model...

    ////////////////////////////
    //// Equiprobable (most simple)
    if (modelName == "Equi")
      model.reset(new EquiprobableSubstitutionModel(alphabet));
    ///////////////////////////////////////////////////
    // nucleotidic model
    else if (AlphabetTools::isNucleicAlphabet(*alphabet))
    {
      if (!(alphabetCode_ & NUCLEOTIDE))
        throw Exception("BppOSubstitutionModelFormat::read. Nucleotide alphabet not supported.");
      auto alpha = dynamic_pointer_cast<const NucleicAlphabet>(alphabet);

      // //////////////////////////////////
      // gBGC
      // //////////////////////////////////
      if (modelName.find("+gBGC") != string::npos)
      {
        string subModName = modelName.substr(0, modelName.find("+gBGC"));

        if (verbose_)
          ApplicationTools::displayResult("Biased gene conversion for", subModName);

        string nestedModelDescription = subModName;

        if (modelDescription.find_first_of("(") != string::npos)
        {
          string::size_type begin = modelDescription.find_first_of("(");
          string::size_type end = modelDescription.find_last_of(")");

          nestedModelDescription += modelDescription.substr(begin, end - begin + 1);
        }

        BppOSubstitutionModelFormat nestedReader(NUCLEOTIDE, true, true, false, verbose_, warningLevel_);
        auto tmpModel = nestedReader.readSubstitutionModel(alphabet, nestedModelDescription, mData, nData, false);
        unique_ptr<NucleotideSubstitutionModelInterface> nestedModel(dynamic_cast<NucleotideSubstitutionModelInterface*>(tmpModel.release()));
        map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

        // Now we create the gBGC substitution model:
        model = make_unique<gBGC>(alpha, std::move(nestedModel));

        // Then we update the parameter set:
        for (auto& it : unparsedParameterValuesNested)
        {
          unparsedArguments_["gBGC." + it.first] = it.second;
        }
      }


      // /////////////////////////////////
      // / GTR
      // ///////////////////////////////

      else if (modelName == "GTR")
      {
        model = make_unique<GTR>(alpha);
      }


      // /////////////////////////////////
      // / SSR
      // ///////////////////////////////

      else if (modelName == "SSR")
      {
        model = make_unique<SSR>(alpha);
      }

      // /////////////////////////////////
      // / L95
      // ///////////////////////////////

      else if (modelName == "L95")
      {
        model = make_unique<L95>(alpha);
      }

      // /////////////////////////////////
      // / RN95
      // ///////////////////////////////

      else if (modelName == "RN95")
      {
        model = make_unique<RN95>(alpha);
      }

      // /////////////////////////////////
      // / RN95s
      // ///////////////////////////////

      else if (modelName == "RN95s")
      {
        model = make_unique<RN95s>(alpha);
      }

      // /////////////////////////////////
      // / TN93
      // //////////////////////////////

      else if (modelName == "TN93")
      {
        model = make_unique<TN93>(alpha);
      }

      // /////////////////////////////////
      // / HKY85
      // ///////////////////////////////

      else if (modelName == "HKY85")
      {
        model = make_unique<HKY85>(alpha);
      }

      // /////////////////////////////////
      // / F81
      // ///////////////////////////////

      else if (modelName == "F81")
      {
        model = make_unique<F81>(alpha);
      }
      // /////////////////////////////////
      // / F84
      // ///////////////////////////////

      else if (modelName == "F84")
      {
        model = make_unique<F84>(alpha);
      }

      // /////////////////////////////////
      // / T92
      // ///////////////////////////////

      else if (modelName == "T92")
      {
        model = make_unique<T92>(alpha);
      }

      // /////////////////////////////////
      // / K80
      // ///////////////////////////////

      else if (modelName == "K80")
      {
        model = make_unique<K80>(alpha);
      }


      // /////////////////////////////////
      // / JC69
      // ///////////////////////////////

      else if (modelName == "JC69")
      {
        model = make_unique<JCnuc>(alpha);
      }
      else
        throw Exception("Model '" + modelName + "' unknown, or does not fit nucleic alphabet.");
    }
    else if (AlphabetTools::isProteicAlphabet(*alphabet))
    {
      ////////////////////////////////////
      //// Proteic model

      if (!(alphabetCode_ & PROTEIN))
        throw Exception("BppOSubstitutionModelFormat::read. Protein alphabet not supported.");
      auto alpha = dynamic_pointer_cast<const ProteicAlphabet>(alphabet);

      if (modelName.find("+F") != string::npos)
      {
        string freqOpt = ApplicationTools::getStringParameter("frequencies", args, "Full", "", true, warningLevel_);
        BppOFrequencySetFormat freqReader(BppOFrequencySetFormat::ALL, false, warningLevel_);
        auto tmpFreq = freqReader.readFrequencySet(alpha, freqOpt, mData, nData, true);
        unique_ptr<ProteinFrequencySetInterface> protFreq(dynamic_cast<ProteinFrequencySetInterface*>(tmpFreq.release()));

        map<string, string> unparsedParameterValuesNested(freqReader.getUnparsedArguments());

        for (auto& it : unparsedParameterValuesNested)
        {
          unparsedArguments_[modelName + "." + it.first] = it.second;
        }

        if (modelName == "JC69+F")
          model = make_unique<JCprot>(alpha, std::move(protFreq));
        else if (modelName == "DSO78+F")
          model = make_unique<DSO78>(alpha, std::move(protFreq));
        else if (modelName == "JTT92+F")
          model = make_unique<JTT92>(alpha, std::move(protFreq));
        else if (modelName == "LG08+F")
          model = make_unique<LG08>(alpha, std::move(protFreq));
        else if (modelName == "WAG01+F")
          model = make_unique<WAG01>(alpha, std::move(protFreq));
        else if (modelName == "Empirical+F")
        {
          string prefix = args["name"];
          if (TextTools::isEmpty(prefix))
            throw Exception("'name' argument missing for user-defined substitution model.");
          string fname = args["file"];
          if (TextTools::isEmpty(fname))
            throw Exception("'file' argument missing for user-defined substitution model.");
          model = make_unique<UserProteinSubstitutionModel>(alpha, args["file"], std::move(protFreq), prefix + "+F.");
        }
      }
      else if (modelName == "JC69")
        model = make_unique<JCprot>(alpha);
      else if (modelName == "DSO78")
        model = make_unique<DSO78>(alpha);
      else if (modelName == "JTT92")
        model = make_unique<JTT92>(alpha);
      else if (modelName == "LG08")
        model = make_unique<LG08>(alpha);
      else if (modelName == "WAG01")
        model = make_unique<WAG01>(alpha);
      // submodels of mixture models
      else if (modelName.substr(0, 9) == "LGL08_CAT")
      {
        string subModelName = modelName.substr(10);

        size_t posp = modelDescription.find("(");

        string modelDesc2 = modelName.substr(0, 9) + modelDescription.substr(posp);

        BppOTransitionModelFormat nestedReader(PROTEIN, true, true, false, verbose_, warningLevel_);

        auto tmpModel = nestedReader.readTransitionModel(alphabet, modelDesc2, mData, nData, false);
        auto nestedModel = unique_ptr<MixedTransitionModelInterface>(dynamic_cast<MixedTransitionModelInterface*>(tmpModel.release()));

        // Check that nestedModel is fine and has subModel of given name
        if (!nestedModel)
          throw Exception("Unknown model " + modelName + ".");

        try
        {
          nestedModel->model(subModelName);
        }
        catch (NullPointerException&)
        {
          throw Exception("BppOSubstitutionModelFormat::read. " + nestedModel->getName() + "argument for model 'FromMixture' has no submodel with name " + subModelName + ".");
        }

        // Now we create the FromMixture substitution model:
        model = make_unique<FromMixtureSubstitutionModel>(*nestedModel, subModelName, modelDesc2);
      }
      else if (modelName == "Empirical")
      {
        string prefix = args["name"];
        if (TextTools::isEmpty(prefix))
          throw Exception("'name' argument missing for user-defined substitution model.");
        model = make_unique<UserProteinSubstitutionModel>(alpha, args["file"], prefix);
      }
      else if (modelName == "Coala")
      {
        string exchangeability = args["exch"];
        if (TextTools::isEmpty(exchangeability))
          throw Exception("BppOSubstitutionModelFormat::read. missing argument 'exch' for model 'Coala'.");
        string prefix = args["name"];
        if (exchangeability == "Empirical" && TextTools::isEmpty(prefix))
          throw Exception("'name' argument missing to specify the exchangeabilities of the user-defined empirical model.");
        BppOSubstitutionModelFormat nestedReader(PROTEIN, false, allowMixed_, allowGaps_, verbose_, warningLevel_);
        auto tmpModel = nestedReader.readSubstitutionModel(alphabet, exchangeability, mData, nData, false);
        unique_ptr<ProteinSubstitutionModelInterface> nestedModel(dynamic_cast<ProteinSubstitutionModelInterface*>(tmpModel.release()));
        string nbrOfParametersPerBranch = args["nbrAxes"];
        if (TextTools::isEmpty(nbrOfParametersPerBranch))
          throw Exception("'nbrAxes' argument missing to define the number of axis of the Correspondence Analysis.");
        // Now we create the Coala model
        model = make_unique<Coala>(alpha, *nestedModel, TextTools::to<unsigned int>(nbrOfParametersPerBranch));
        if (nData)
          model->setFreqFromData(*mData.at(nData));
      }
      else
        throw Exception("Model '" + modelName + "' is unknown, or does not fit proteic alphabet.");
    }
    else if (AlphabetTools::isBinaryAlphabet(*alphabet))
    {
      if (!(alphabetCode_ & BINARY))
        throw Exception("BppOSubstitutionModelFormat::read. Binary alphabet not supported.");

      auto balpha = dynamic_pointer_cast<const BinaryAlphabet>(alphabet);

      if (modelName == "Binary")
        model = make_unique<BinarySubstitutionModel>(balpha);
      else if (modelName == "TwoParameterBinary")
        model = make_unique<TwoParameterBinarySubstitutionModel>(balpha);
      else
        throw Exception("Model '" + modelName + "' unknown, or does not fit binary alphabet.");
    }
    else if (AlphabetTools::isIntegerAlphabet(*alphabet))
    {
      if (!(alphabetCode_ & INTEGER))
        throw Exception("BppOSubstitutionModelFormat::read. Integer alphabet not supported.");

      auto balpha = dynamic_pointer_cast<const IntegerAlphabet>(alphabet);

      if (modelName == "D1Walk")
        model.reset(new D1WalkSubstitutionModel(balpha));
      else
        throw Exception("Model '" + modelName + "' unknown, or does not fit integer alphabet.");
    }
    else
      throw Exception("Model '" + modelName + "' unknown, or does not fit " + alphabet->getAlphabetType() + " alphabet.");
  }

  if (verbose_)
    ApplicationTools::displayResult("Substitution model", modelName);

  updateParameters_(*model, args);


  if (parseArguments)
  {
    if (nData)
      initialize_(*model, mData.at(nData));
    else
      initialize_(*model, 0);
  }
  
  return model;
}


void BppOSubstitutionModelFormat::updateParameters_(
    BranchModelInterface& model,
    map<std::string, std::string>& args)
{
  // Update parameter args:
  vector<string> pnames = model.getParameters().getParameterNames();

  string pref = model.getNamespace();

  for (auto pname : pnames)
  {
    string name = model.getParameterNameWithoutNamespace(pname);
    if (args.find(name) != args.end())
      unparsedArguments_[pref + name] = args[name];
  }

  // Now look if some parameters are aliased:
  ParameterList pl = model.getIndependentParameters();
  string pname, pval, pname2;
  for (size_t i = 0; i < pl.size(); ++i)
  {
    pname = model.getParameterNameWithoutNamespace(pl[i].getName());

    if (args.find(pname) == args.end())
      continue;
    pval = args[pname];

    if (((pval.rfind("_") != string::npos) && (TextTools::isDecimalInteger(pval.substr(pval.rfind("_") + 1, string::npos)))) ||
        (pval.find("(") != string::npos))
      continue;

    bool found = false;
    for (size_t j = 0; j < pl.size() && !found; ++j)
    {
      pname2 = model.getParameterNameWithoutNamespace(pl[j].getName());

      // if (j == i || args.find(pname2) == args.end()) continue; Julien 03/03/2010: This extra condition prevents complicated (nested) models to work properly...
      if (j == i)
        continue;
      if (pval == pname2)
      {
        // This is an alias...
        // NB: this may throw an exception if incorrect! We leave it as is for now :s
        model.aliasParameters(pname2, pname);
        if (verbose_)
          ApplicationTools::displayResult("Parameter alias found", pname + "->" + pname2);
        found = true;
      }
    }
    // if (!TextTools::isDecimalNumber(pval) && !found)
    //   throw Exception("Incorrect parameter syntax: parameter " + pval + " was not found and can't be used as a value for parameter " + pname + ".");
  }

  // 2 following tests be removed in a later version
  if (args.find("useObservedFreqs") != args.end())
    throw Exception("useObservedFreqs argument is obsolete. Please use 'initFreqs=observed' instead.");
  if (args.find("useObservedFreqs.pseudoCount") != args.end())
    throw Exception("useObservedFreqs.pseudoCount argument is obsolete. Please use 'initFreqs.observedPseudoCount' instead.");

  if (args.find("initFreqs") != args.end())
    // Specific case since already used
    if (pref!="Coala.")
      unparsedArguments_[pref + "initFreqs"] = args["initFreqs"];
  if (args.find("initFreqs.observedPseudoCount") != args.end())
    unparsedArguments_[pref + "initFreqs.observedPseudoCount"] = args["initFreqs.observedPseudoCount"];

  if (args.find("rate") != args.end())
  {
    model.addRateParameter();
    unparsedArguments_[pref + "rate"] = args["rate"];
  }
}


unique_ptr<SubstitutionModelInterface> BppOSubstitutionModelFormat::readWord_(
    std::shared_ptr<const Alphabet> alphabet,
    const std::string& modelDescription,
    const std::map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mData,
    size_t nData)
{
  unique_ptr<SubstitutionModelInterface> model;
  string modelName = "";
  map<string, string> args;
  KeyvalTools::parseProcedure(modelDescription, modelName, args);

  vector<string> v_nestedModelDescription;
  vector< unique_ptr<SubstitutionModelInterface>> v_pSM;
  shared_ptr<const CoreWordAlphabet> pWA;

  string s, nestedModelDescription;
  unsigned int nbmodels;

  if (((modelName == "Word" || modelName == "Kron") && !AlphabetTools::isWordAlphabet(*alphabet)) ||
      ((!(modelName == "Word" || modelName == "Kron")) && !AlphabetTools::isCodonAlphabet(*alphabet)))
    throw Exception("Bad alphabet type "
          + alphabet->getAlphabetType() + " for  model " + modelName + ".");

  pWA = dynamic_pointer_cast<const CoreWordAlphabet>(alphabet);


  ////////////////////////////////////
  /// Reading the submodels

  if (args.find("model") != args.end())
  {
    v_nestedModelDescription.push_back(args["model"]);
    nbmodels = (modelName == "Word" || modelName == "Kron") ? pWA->getLength() : 3;
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
    BppOSubstitutionModelFormat nestedReader(NUCLEOTIDE, false, true, false, false, warningLevel_);
    model = nestedReader.readSubstitutionModel(pWA->getNAlphabet(0), v_nestedModelDescription[0], mData, nData, false);
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
    string pref = "";
    for (unsigned int i = 0; i < nbmodels; ++i)
    {
      pref += TextTools::toString(i + 1);
    }

    for (auto& it : unparsedParameterValuesNested)
    {
      unparsedArguments_[modelName + "." + pref + "_" + it.first] = it.second;
    }

    v_pSM.push_back(std::move(model));
  }
  else
  {
    for (unsigned i = 0; i < v_nestedModelDescription.size(); ++i)
    {
      BppOSubstitutionModelFormat nestedReader(alphabetCode_, false, true, false, false, warningLevel_);
      model = nestedReader.readSubstitutionModel(pWA->getNAlphabet(i), v_nestedModelDescription[i], mData, nData, false);
      map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
      for (auto& it : unparsedParameterValuesNested)
      {
        unparsedArguments_[modelName + "." + TextTools::toString(i + 1) + "_" + it.first] = it.second;
      }

      v_pSM.push_back(std::move(model));
    }
  }

  // //////////////////////////////
  // In case of a Kron... model, mgmt of the positions

  // check the vector of simultaneous changing positions

  vector<set<size_t>> vposKron;
  if (modelName.substr(0, 4) == "Kron")
  {
    if (args.find("positions") != args.end())
    {
      StringTokenizer nst(args["positions"], "+");

      while (nst.hasMoreToken())
      {
        string spos = nst.nextToken();
        StringTokenizer nst2(spos, "*");

        set<size_t> ss;

        while (nst2.hasMoreToken())
        {
          ss.insert((size_t)TextTools::toInt(nst2.nextToken()));
        }

        vposKron.push_back(ss);
      }
    }
  }

  // /////////////////////////////////
  // / WORD
  // ///////////////////////////////

  if (modelName == "Word")
  {
    if (v_nestedModelDescription.size() != nbmodels)
    {
      model = make_unique<WordSubstitutionModel>(std::move(v_pSM[0]), nbmodels);
    }
    else
    {
      ModelList ml(v_pSM);
      model = make_unique<WordSubstitutionModel>(ml);
    }
  }

  // /////////////////////////////////
  // / KRON
  // ///////////////////////////////

  else if (modelName == "Kron")
  {
    if (vposKron.size() == 0)
    {
      if (v_nestedModelDescription.size() != nbmodels)
      {
        model = make_unique<KroneckerWordSubstitutionModel>(std::move(v_pSM[0]), nbmodels);
      }
      else
      {
        ModelList ml(v_pSM);
        model = make_unique<KroneckerWordSubstitutionModel>(ml);
      }
    }
    else
    {
      if (v_nestedModelDescription.size() != nbmodels)
      {
        model = make_unique<KroneckerWordSubstitutionModel>(std::move(v_pSM[0]), nbmodels, vposKron);
      }
      else
      {
        ModelList ml(v_pSM);
        model = make_unique<KroneckerWordSubstitutionModel>(ml, vposKron);
      }
    }
  }

  // /////////////////////////////////
  // / CODON
  // ///////////////////////////////

  else
  {
    auto pCA = dynamic_pointer_cast<const CodonAlphabet>(pWA);
    if (pCA == 0)
      throw Exception("Non codon Alphabet for model" + modelName + ".");

    unique_ptr<AlphabetIndex2> pai2;
    unique_ptr<CodonFrequencySetInterface> pFS;

    if ((!dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[0].get())) ||
        ((v_nestedModelDescription.size() == 3) &&
        (!dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[1].get()) || !dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[2].get()))))
      throw Exception("Non simple NucleotideSubstitutionModel embedded in " + modelName + " model.");

    if (args.find("genetic_code") != args.end())
    {
      ApplicationTools::displayWarning("'genetic_code' argument is no longer supported inside model description, and has been supersided by a global 'genetic_code' option.");
      throw Exception("BppOSubstitutionModelFormat::read. Deprecated 'genetic_code' argument.");
    }

    if (!geneticCode_)
      throw Exception("BppOSubstitutionModelFormat::readWord_(). No genetic code specified! Consider using 'setGeneticCode'.");


    ///////////////////////////////////
    /// Dist

    if ((modelName.find("Dist") != string::npos) || (modelName == "SENCA"))
      pai2 = (args.find("aadistance") == args.end()) ? nullptr : SequenceApplicationTools::getAlphabetIndex2(AlphabetTools::PROTEIN_ALPHABET, args["aadistance"]);

    //////////////////////////////////
    /// Freq

    if (modelName.find("Freq") != string::npos)
    {
      if (args.find("frequencies") == args.end())
        throw Exception("Missing equilibrium frequencies.");

      BppOFrequencySetFormat bIOFreq(alphabetCode_, verbose_, warningLevel_);
      bIOFreq.setGeneticCode(geneticCode_); // This uses the same instance as the one that will be used by the model
      auto tmp = bIOFreq.readFrequencySet(pCA, args["frequencies"], mData, nData);
      pFS = unique_ptr<CodonFrequencySetInterface>(dynamic_cast<CodonFrequencySetInterface*>(tmp.release()));
      if (!pFS)
        throw IOException("BppOSubstitutionModelFormat::readWord_(). Incorrect codon frequencies.");
      map<string, string> unparsedParameterValuesNested(bIOFreq.getUnparsedArguments());

      for (auto& it : unparsedParameterValuesNested)
      {
        unparsedArguments_[modelName + "." + it.first] = it.second;
      }
    }

    // //////////////
    // // Triplet

    if (modelName == "Triplet")
      model = (v_nestedModelDescription.size() != 3)
                  ? make_unique<TripletSubstitutionModel>(
            pCA,
            unique_ptr<NucleotideSubstitutionModelInterface>(
            dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[0].release())
            )
            )
                  : make_unique<TripletSubstitutionModel>(
            pCA,
            unique_ptr<NucleotideSubstitutionModelInterface>(
            dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[0].release())
            ),
            unique_ptr<NucleotideSubstitutionModelInterface>(
            dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[1].release())
            ),
            unique_ptr<NucleotideSubstitutionModelInterface>(
            dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[2].release())
            )
            );

    else if (modelName.find("Codon") != string::npos)
    {
      vector< unique_ptr<CoreCodonSubstitutionModelInterface>> vCSM;
      string name = "Codon";
      map<string, string> unparsedParameterValuesNested;

      if (modelName.find("Dist") != string::npos)
      {
        name += "Dist";

        vCSM.push_back(make_unique<AbstractCodonDistanceSubstitutionModel>(std::move(pai2), geneticCode_, ""));
      }
      else if (modelName.find("BGC") != string::npos)
      {
        name += "BGC";

        vCSM.push_back(make_unique<AbstractCodonBGCSubstitutionModel>(geneticCode_, ""));
      }
      else if (modelName.find("Prot") != string::npos)
      {
        name += "Prot";

        if (args.find("protmodel") == args.end())
          throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'protmodel' for codon model argument 'Prot'.");

        nestedModelDescription = args["protmodel"];
        BppOSubstitutionModelFormat nestedReader(PROTEIN, false, false, allowGaps_, verbose_, warningLevel_);

        shared_ptr<SubstitutionModelInterface> tmpModel = nestedReader.readSubstitutionModel(geneticCode_->getTargetAlphabet(), nestedModelDescription, mData, nData, false);
        auto nestedModel = dynamic_pointer_cast<ProteinSubstitutionModelInterface>(tmpModel);

        unparsedParameterValuesNested.insert(nestedReader.getUnparsedArguments().begin(), nestedReader.getUnparsedArguments().end());

        vCSM.push_back(make_unique<AbstractCodonAARateSubstitutionModel>(nestedModel, geneticCode_, ""));
      }
      if (modelName.find("AAClust") != string::npos)
      {
        name += "AAClust";

        // Initialization using the "assign" argument
        vector<uint> partition;
        if (args.find("partition") != args.end())
        {
          string rf = args["partition"];

          StringTokenizer strtok(rf.substr(1, rf.length() - 2), ",");
          while (strtok.hasMoreToken())
            partition.push_back(TextTools::to<uint>(strtok.nextToken()));
        }

        unique_ptr<AbstractCodonClusterAASubstitutionModel> aca = partition.size() != 0 ?
            make_unique<AbstractCodonClusterAASubstitutionModel>(geneticCode_, "", partition) :
            make_unique<AbstractCodonClusterAASubstitutionModel>(geneticCode_, "");

        vCSM.push_back(std::move(aca));
      }

      /// Default name in none used before
      if (vCSM.size() == 0)
        name += "Rate";

      if (modelName.find("CpG") != string::npos)
      {
        name += "CpG";
        vCSM.push_back(make_unique<AbstractCodonCpGSubstitutionModel>(pCA, ""));
      }

      if (modelName.find("AAFit") != string::npos)
      {
        name += "AAFit";

        if (args.find("fitness") == args.end())
          throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'fitness' for codon model argument 'AAFit'.");

        string nestedFreqDescription = args["fitness"];
        BppOFrequencySetFormat nestedReader(PROTEIN, verbose_, warningLevel_);

        auto nestedFreq = nestedReader.readFrequencySet(geneticCode_->getTargetAlphabet(), nestedFreqDescription, mData, nData, false);

        for (auto it : nestedReader.getUnparsedArguments())
        {
          unparsedParameterValuesNested["fit_" + it.first] = it.second;
        }

        auto aca = make_unique<AbstractCodonAAFitnessSubstitutionModel>(std::move(nestedFreq), geneticCode_, "");

        if (args.find("Ns") != args.end())
          aca->addNsParameter();

        vCSM.push_back(std::move(aca));
      }
      else if (modelName.find("Fit") != string::npos)
      {
        if (args.find("fitness") == args.end())
          throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'fitness' for codon model argument 'Fit'.");
        string nestedFreqDescription = args["fitness"];

        BppOFrequencySetFormat nestedReader(alphabetCode_, verbose_, warningLevel_);
        nestedReader.setGeneticCode(geneticCode_);

        auto nestedFreq = nestedReader.readFrequencySet(alphabet, nestedFreqDescription, mData, nData, false);

        for (auto it : nestedReader.getUnparsedArguments())
        {
          unparsedParameterValuesNested["fit_" + it.first] = it.second;
        }

        vCSM.push_back(make_unique<AbstractCodonFitnessSubstitutionModel>(std::move(nestedFreq), geneticCode_, ""));
      }

      if (modelName.find("PhasFreq") != string::npos)
      {
        name += "PhasFreq";
        vCSM.push_back(make_unique<AbstractCodonPhaseFrequenciesSubstitutionModel>(std::move(pFS), ""));
      }
      else if (modelName.find("Freq") != string::npos)
      {
        name += "Freq";
        vCSM.push_back(make_unique<AbstractCodonFrequenciesSubstitutionModel>(std::move(pFS), ""));
      }

      // Then we update the parameter set:
      for (auto it : unparsedParameterValuesNested)
      {
        unparsedArguments_[name + "." + it.first] = it.second;
      }

      model = (v_nestedModelDescription.size() != 3)
                  ? make_unique<CodonAdHocSubstitutionModel>(
            geneticCode_,
            unique_ptr<NucleotideSubstitutionModelInterface>(
            dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[0].release())
            ),
            vCSM,
            name)
                  : make_unique<CodonAdHocSubstitutionModel>(
            geneticCode_,
            unique_ptr<NucleotideSubstitutionModelInterface>(
            dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[0].release())
            ),
            unique_ptr<NucleotideSubstitutionModelInterface>(
            dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[1].release())
            ),
            unique_ptr<NucleotideSubstitutionModelInterface>(
            dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[2].release())
            ),
            vCSM,
            name);
    }
    else if (modelName == "KronDistFreq")
    {
      if (v_nestedModelDescription.size() != 3)
      {
        if (vposKron.size() == 0)
          model = make_unique<KroneckerCodonDistanceFrequenciesSubstitutionModel>(geneticCode_,
                unique_ptr<NucleotideSubstitutionModelInterface>(
                dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[0].release())),
                std::move(pFS), std::move(pai2));
        else
          model = make_unique<KroneckerCodonDistanceFrequenciesSubstitutionModel>(geneticCode_,
                unique_ptr<NucleotideSubstitutionModelInterface>(
                dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[0].release())),
                vposKron, std::move(pFS), std::move(pai2));
      }
      else
      {
        if (vposKron.size() != 0)
          model = make_unique<KroneckerCodonDistanceFrequenciesSubstitutionModel>(
                geneticCode_,
                unique_ptr<NucleotideSubstitutionModelInterface>(
                dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[0].release())
                ),
                unique_ptr<NucleotideSubstitutionModelInterface>(
                dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[1].release())
                ),
                unique_ptr<NucleotideSubstitutionModelInterface>(
                dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[2].release())
                ),
                std::move(pFS),
                std::move(pai2));
        else
          model = make_unique<KroneckerCodonDistanceFrequenciesSubstitutionModel>(
                geneticCode_,
                unique_ptr<NucleotideSubstitutionModelInterface>(
                dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[0].release())
                ),
                unique_ptr<NucleotideSubstitutionModelInterface>(
                dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[1].release())
                ),
                unique_ptr<NucleotideSubstitutionModelInterface>(
                dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[2].release())
                ),
                vposKron,
                std::move(pFS),
                std::move(pai2));
      }
    }
    else if (modelName == "KronDist")
    {
      if (v_nestedModelDescription.size() != 3)
      {
        if (vposKron.size() == 0)
          model = make_unique<KroneckerCodonDistanceSubstitutionModel>(geneticCode_,
                unique_ptr<NucleotideSubstitutionModelInterface>(dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[0].release())),
                std::move(pai2));
        else
          model = make_unique<KroneckerCodonDistanceSubstitutionModel>(geneticCode_,
                unique_ptr<NucleotideSubstitutionModelInterface>(dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[0].release())),
                vposKron, std::move(pai2));
      }
      else
      {
        if (vposKron.size() != 0)
          model = make_unique<KroneckerCodonDistanceSubstitutionModel>(
                geneticCode_,
                unique_ptr<NucleotideSubstitutionModelInterface>(
                dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[0].release())
                ),
                unique_ptr<NucleotideSubstitutionModelInterface>(
                dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[1].release())
                ),
                unique_ptr<NucleotideSubstitutionModelInterface>(
                dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[2].release())
                ),
                std::move(pai2));
        else
          model = make_unique<KroneckerCodonDistanceSubstitutionModel>(
                geneticCode_,
                unique_ptr<NucleotideSubstitutionModelInterface>(
                dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[0].release())
                ),
                unique_ptr<NucleotideSubstitutionModelInterface>(
                dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[1].release())
                ),
                unique_ptr<NucleotideSubstitutionModelInterface>(
                dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[2].release())
                ),
                vposKron,
                std::move(pai2));
      }
    }

    else if (modelName == "SENCA")
    {
      if (args.find("fitness") == args.end())
        throw Exception("Missing fitness in model " + modelName + ".");

      BppOFrequencySetFormat bIOFreq(alphabetCode_, verbose_, warningLevel_);
      bIOFreq.setGeneticCode(geneticCode_);

      auto pFit(bIOFreq.readFrequencySet(pCA, args["fitness"], mData, nData, false));
      map<string, string> unparsedParameterValuesNested(bIOFreq.getUnparsedArguments());

      for (auto& it : unparsedParameterValuesNested)
      {
        unparsedArguments_[modelName + ".fit_" + it.first] = it.second;
      }

      if (v_nestedModelDescription.size() != 3)
      {
        model = make_unique<SENCA>(geneticCode_,
              unique_ptr<NucleotideSubstitutionModelInterface>(
              dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[0].release())
              ),
              std::move(pFit), std::move(pai2));
      }
      else
        model = make_unique<SENCA>(
              geneticCode_,
              unique_ptr<NucleotideSubstitutionModelInterface>(
              dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[0].release())
              ),
              unique_ptr<NucleotideSubstitutionModelInterface>(
              dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[1].release())
              ),
              unique_ptr<NucleotideSubstitutionModelInterface>(
              dynamic_cast<NucleotideSubstitutionModelInterface*>(v_pSM[2].release())
              ),
              std::move(pFit),
              std::move(pai2));
    }
  }

  return model;
}

void BppOSubstitutionModelFormat::write(const BranchModelInterface& model,
    OutputStream& out,
    std::map<std::string, std::string>& globalAliases,
    std::vector<std::string>& writtenNames) const
{
  bool comma = false;

  //  Mixed Model that are defined as "Mixture" and "Mixed"

  if ((dynamic_cast<const MixedTransitionModelInterface*>(&model) != nullptr) && (dynamic_cast<const AbstractBiblioMixedTransitionModel*>(&model) == nullptr))
  {
    writeMixed_(dynamic_cast<const MixedTransitionModelInterface&>(model), out, globalAliases, writtenNames);
    return;
  }

  auto name = model.getName();
  auto parend = (name.rfind(")") == name.size() - 1);

  out << (parend ? name.substr(0, name.size() - 1) : name + "(");

  // Is it a protein user defined model?
  try
  {
    const UserProteinSubstitutionModel& userModel = dynamic_cast<const UserProteinSubstitutionModel&>(model);
    out << "file=" << userModel.getPath();
    string ns = userModel.getNamespace();

    if (TextTools::hasSubstring(ns, "+F.") )
      ns = ns.substr(0, ns.size() - 3);
    else
      ns = ns.substr(0, ns.size() - 1);

    out << ", name=" << ns;
    comma = true;
  }
  catch (bad_cast&)
  {}

  // Is it a markov-modulated model?
  try
  {
    const MarkovModulatedSubstitutionModel& mmModel = dynamic_cast<const MarkovModulatedSubstitutionModel&>(model);
    out << "model=";
    write(mmModel.nestedModel(), out, globalAliases, writtenNames);

    try
    {
      const G2001& gModel = dynamic_cast<const G2001&>(model);
      // Also print distribution here:
      out << ", rdist=";
      const DiscreteDistributionInterface& nestedDist = gModel.rateDistribution();
      const BppODiscreteDistributionFormat bIO;
      bIO.writeDiscreteDistribution(nestedDist, out, globalAliases, writtenNames);
    }
    catch (bad_cast&)
    {}
    comma = true;
  }
  catch (bad_cast&)
  {}

  // Is it a model with gaps?
  try
  {
    const RE08& reModel = dynamic_cast<const RE08&>(model);
    out << "model=";
    write(reModel.nestedModel(), out, globalAliases, writtenNames);
    comma = true;
  }
  catch (bad_cast&)
  {}

  // Is it a YpR model?
  try
  {
    const YpR& yprModel = dynamic_cast<const YpR&>(model);
    out << "model=";
    write(yprModel.nestedModel(), out, globalAliases, writtenNames);
    comma = true;
  }
  catch (bad_cast&)
  {}

  // Is it a OneChange model?
  try
  {
    const OneChangeTransitionModel& oneChangeTransitionModel = dynamic_cast<const OneChangeTransitionModel&>(model);
    out << "model=";
    write(oneChangeTransitionModel.model(), out, globalAliases, writtenNames);
    comma = true;
  }
  catch (bad_cast&)
  {
    // Is it a model with register?
    try
    {
      const OneChangeRegisterTransitionModel& oneChangeRegisterTransitionModel = dynamic_cast<const OneChangeRegisterTransitionModel&>(model);
      out << "model=";
      write(oneChangeRegisterTransitionModel.model(), out, globalAliases, writtenNames);
      comma = true;
      out << ", register=";
      out << oneChangeRegisterTransitionModel.getRegisterName();
      out << ", numReg=" << VectorTools::paste(oneChangeRegisterTransitionModel.getRegisterNumbers(), "+");
    }
    catch (bad_cast&)
    {}

    try
    {
      const RegisterRatesSubstitutionModel& registerRatesSubstitutionModel = dynamic_cast<const RegisterRatesSubstitutionModel&>(model);
      out << "model=";
      write(registerRatesSubstitutionModel.model(), out, globalAliases, writtenNames);
      comma = true;
      out << ", register=";
      out << registerRatesSubstitutionModel.getRegisterName();
    }
    catch (bad_cast&)
    {}
  }
  // Is it an Integration model
  try {
    const IntegrationOfSubstitutionModel& integr = dynamic_cast<const IntegrationOfSubstitutionModel&>(model);
    out << "model=";
    write(integr.model(), out, globalAliases, writtenNames);

    out << ", k=";
    out << integr.k();
    out << ",zetas=(";
    const auto& zetas=integr.getZetas();
    comma = false;    
    for (const auto& zeta:zetas)
    {
      if (!comma)
        comma=true;
      else
        out << ",";
      out << zeta;
    }
    out << ")";
    for (size_t i=1;i<zetas.size();i++)
      writtenNames.push_back("Integrate.theta"+TextTools::toString(i));
  } catch(bad_cast&) {}
  
  // Is it a model with register?
  try {
    const OneChangeRegisterTransitionModel& oneChangeRegisterTransitionModel = dynamic_cast<const OneChangeRegisterTransitionModel&>(model);
    out << "model=";
    write(oneChangeRegisterTransitionModel.model(), out, globalAliases, writtenNames);
    comma = true;
    out << ", register=";
    out << oneChangeRegisterTransitionModel.getRegisterName();
    out << ", numReg=" << VectorTools::paste(oneChangeRegisterTransitionModel.getRegisterNumbers(), "+");
  } catch(bad_cast&) {}

  try {
    const RegisterRatesSubstitutionModel& registerRatesSubstitutionModel = dynamic_cast<const RegisterRatesSubstitutionModel&>(model);
    out << "model=";
    write(registerRatesSubstitutionModel.model(), out, globalAliases, writtenNames);
    comma = true;
    out << ", register=";
    out << registerRatesSubstitutionModel.getRegisterName();
  } catch(bad_cast&) {}

  // Is it a gBGC model?
  try
  {
    const gBGC& gbgcModel = dynamic_cast<const gBGC&>(model);
    StdStr sout;
    write(gbgcModel.nestedModel(), sout, globalAliases, writtenNames);

    string ss = sout.str();
    auto begin = ss.find_first_of("(");
    auto end = ss.find_last_of(")");

    out << ss.substr(begin + 1, end - begin - 1);

    comma = true;
  }
  catch (bad_cast&)
  {}

  // Is it a word model?
  try
  {
    const AbstractWordSubstitutionModel& wM = dynamic_cast<const AbstractWordSubstitutionModel&>(model);
    size_t nmod = wM.getNumberOfModels();
    const TransitionModelInterface& mod0 = wM.nModel(0);
    if (nmod == 1)
    {
      out << "model=";
      write(mod0, out, globalAliases, writtenNames);
    }
    else
    {
      const TransitionModelInterface& mod1 = wM.nModel(1);
      if (&mod1 == &mod0)
      {
        out << "model=";
        write(mod0, out, globalAliases, writtenNames);
      }
      else
      {
        out << "model1=";
        write(mod0, out, globalAliases, writtenNames);
        for (unsigned int i = 1; i < nmod; ++i)
        {
          out << ", model" + TextTools::toString(i + 1) + "=";
          write(wM.nModel(i), out, globalAliases, writtenNames);
        }
      }
    }
    comma = true;
  }
  catch (bad_cast&)
  {}

  // Is it a COaLA model ?
  try
  {
    const Coala& coalaModel = dynamic_cast<const Coala&>(model);
    out << "exch=" << coalaModel.getExch() << ", nbrAxes=" << coalaModel.getNbrOfAxes();
    comma = true;
  }
  catch (bad_cast&)
  {}

  // Is it a InMixed model?
  try
  {
    const InMixedSubstitutionModel& inModel = dynamic_cast<const InMixedSubstitutionModel&>(model);
    out << "model=";
    write(inModel.mixedModel(), out, globalAliases, writtenNames);
    out << ", numMod=" << TextTools::toString(inModel.getSubModelNumber());
    comma = true;
  }
  catch (bad_cast&)
  {}

  // Is it a POMO model?
  const POMO* pomo = dynamic_cast<const POMO*>(&model);
  if (pomo)
  {
    out << "model=";
    write(pomo->mutationModel(), out, globalAliases, writtenNames);

    if (pomo->hasFitness())
    {
      out << ", fitness=";

      BppOFrequencySetFormat bIOFreq(alphabetCode_, false, warningLevel_);
      bIOFreq.writeFrequencySet(pomo->fitness(), out, globalAliases, writtenNames);
    }
    comma = true;
  }

  // Is it a model with FrequencySet?

  try
  {
    auto& pfs = model.frequencySet();

    StdStr freqdesc;
    BppOFrequencySetFormat bIOFreq(alphabetCode_, false, warningLevel_);
    bIOFreq.writeFrequencySet(pfs, freqdesc, globalAliases, writtenNames);
    string fd = freqdesc.str();
    if (fd.size() != 0)
    {
      if (comma)
        out << ", ";
      out << "frequencies=" + fd;
    }
    comma = true;
  }
  catch (exception&)
  {}

  // Is it a codon model with Protein Model or partition in it?
  try
  {
    auto& casm = dynamic_cast<const CodonAdHocSubstitutionModel&>(model);
    for (size_t i = 0; i < casm.getNumberOfModels(); ++i)
    {
      try
      {
        auto& acr = dynamic_cast<const AbstractCodonAARateSubstitutionModel&>(casm.layerModel(i));
        if (comma)
          out << ", ";
        out << "protmodel=";

        write(acr.aaModel(), out, globalAliases, writtenNames);
        comma = true;
      }
      catch (bad_cast&)
      {}

      try
      {
        auto& acf = dynamic_cast<const AbstractCodonAAFitnessSubstitutionModel&>(casm.layerModel(i));
        if (comma)
          out << ", ";
        out << "fitness=";

        BppOFrequencySetFormat bIOFreq(PROTEIN, false, warningLevel_);
        bIOFreq.writeFrequencySet(acf.aaFitness(), out, globalAliases, writtenNames);
        comma = true;
      }
      catch (bad_cast&)
      {}

      try
      {
        auto& acc = dynamic_cast<const AbstractCodonClusterAASubstitutionModel&>(casm.layerModel(i));
        if (comma)
          out << ", ";
        out << "partition=(";
        const vector<uint>& vass = acc.getAssign();

        for (size_t j = 0; j < vass.size(); ++j)
        {
          if (j != 0)
            out << ", ";
          out << vass[j];
        }
        out << ")";
        comma = true;
      }
      catch (bad_cast&)
      {}
    }
  }
  catch (bad_cast&)
  {}

  // Specific case of SENCA

  try
  {
    auto& pCF = dynamic_cast<const SENCA&>(model);
    if (comma)
      out << ", ";
    out << "fitness=";

    BppOFrequencySetFormat bIOFreq(alphabetCode_, false, warningLevel_);
    bIOFreq.writeFrequencySet(pCF.fitness(), out, globalAliases, writtenNames);
    comma = true;
  }
  catch (bad_cast&)
  {}

  // and bibliomixed models

  try
  {
    auto& pM7 = dynamic_cast<const YNGP_M7&>(model);
    if (comma)
      out << ", ";
    out << "n=" << pM7.getNumberOfModels();
    comma = true;
  }
  catch (bad_cast&)
  {}

  try
  {
    auto& pM8 = dynamic_cast<const YNGP_M8&>(model);
    if (comma)
      out << ", ";
    out << "n=" << pM8.getNumberOfModels() - 1;
    comma = true;
  }
  catch (bad_cast&)
  {}

  try
  {
    auto& pM9 = dynamic_cast<const YNGP_M9&>(model);
    if (comma)
      out << ", ";
    out << "nbeta=" << pM9.getNBeta() << ", ngamma=" << pM9.getNGamma();

    comma = true;
  }
  catch (bad_cast&)
  {}

  try
  {
    auto& pM10 = dynamic_cast<const YNGP_M10&>(model);
    if (comma)
      out << ", ";
    out << "nbeta=" << pM10.getNBeta() << ", ngamma=" << pM10.getNGamma();

    comma = true;
  }
  catch (bad_cast&)
  {}

  try
  {
    auto& pLGL = dynamic_cast<const LGL08_CAT&>(model);
    if (comma)
      out << ", ";
    out << "nbCat=" << pLGL.getNumberOfCategories();

    comma = true;
  }
  catch (bad_cast&)
  {}

  try
  {
    auto& pDFP = dynamic_cast<const DFP07&>(model);
    if (comma)
      out << ", ";

    out << "protmodel=";
    write(pDFP.proteinModel(), out, globalAliases, writtenNames);

    comma = true;
  }
  catch (bad_cast&)
  {}

  try
  {
    auto& pSameAA = dynamic_cast<const CodonSameAARateSubstitutionModel&>(model);
    if (comma)
      out << ", ";

    out << "codonmodel=";
    write(pSameAA.codonModel(), out, globalAliases, writtenNames);

    out << ", protmodel=";
    write(pSameAA.proteinModel(), out, globalAliases, writtenNames);

    comma = true;
  }
  catch (bad_cast&)
  {}


  // case of Biblio models, update writtenNames

  try
  {
    auto& absm = dynamic_cast<const AbstractBiblioTransitionModel&>(model);
    size_t wNs = writtenNames.size();

    for (size_t i = 0; i < wNs; ++i)
    {
      try
      {
        writtenNames.push_back(absm.getNamespace() + absm.getParNameFromPmodel(writtenNames[i]));
      }
      catch (Exception&)
      {}
    }
  }
  catch (bad_cast&)
  {}

  BppOParametrizableFormat bIO;
  bIO.write(model, out, globalAliases, model.getIndependentParameters().getParameterNames(), writtenNames, true, comma);
  out << ")";
}


void BppOSubstitutionModelFormat::writeMixed_(const MixedTransitionModelInterface& model,
    OutputStream& out,
    std::map<std::string, std::string>& globalAliases,
    std::vector<std::string>& writtenNames) const
{
  try
  {
    auto& pMS = dynamic_cast<const MixtureOfTransitionModels&>(model);

    out << "Mixture(";
    for (unsigned int i = 0; i < pMS.getNumberOfModels(); ++i)
    {
      if (i != 0)
        out << ", ";
      out << "model" + TextTools::toString(i + 1) + "=";
      write(pMS.nModel(i), out, globalAliases, writtenNames);
    }
  }
  catch (bad_cast&)
  {
    try
    {
      auto& pMS = dynamic_cast<const MixtureOfATransitionModel&>(model);
      out << "MixedModel(model=";
      auto& eM = pMS.model(0);

      ParameterList pl = eM.getIndependentParameters();
      vector<string> vpl = pl.getParameterNames();

      for (auto& pn : vpl)
      {
        if (find(writtenNames.begin(), writtenNames.end(), pn) == writtenNames.end())
        {
          if (pMS.hasDistribution(pn))
          {
            auto& pDD = pMS.distribution(pn);
            if (!dynamic_cast<const ConstantDistribution*>(&pDD))
            {
              const BppODiscreteDistributionFormat bIO;
              StdStr sout;
              bIO.writeDiscreteDistribution(pDD, sout, globalAliases, writtenNames);
              globalAliases[pn] = sout.str();
            }
          }
        }
      }

      write(eM, out, globalAliases, writtenNames);

      if (pMS.from() != -1)
        out << ", from=" << model.getAlphabet()->intToChar(pMS.from()) << ", to=" << model.getAlphabet()->intToChar(pMS.to());
    }
    catch (bad_cast&)
    {}
  }

  const BppOParametrizableFormat bIO;
  bIO.write(model, out, globalAliases, model.getIndependentParameters().getParameterNames(), writtenNames, true, true);
  out << ")";
}

void BppOSubstitutionModelFormat::initialize_(
    BranchModelInterface& model,
    std::shared_ptr<const AlignmentDataInterface> data)
{
  string initFreqs = ApplicationTools::getStringParameter(model.getNamespace() + "initFreqs", unparsedArguments_, "", "", true, warningLevel_);
  if (verbose_)
    ApplicationTools::displayResult("External model frequencies init", (initFreqs == "") ? "None" : initFreqs);

  if (initFreqs != "")
  {
    try
    {
      auto& tmodel = dynamic_cast<TransitionModelInterface&>(model);
      if (initFreqs == "observed")
      {
        if (!data)
          throw Exception("BppOSubstitutionModelFormat::initialize_(). Missing data for observed frequencies");
        unsigned int psi = ApplicationTools::getParameter<unsigned int>(model.getNamespace() + "initFreqs.observedPseudoCount", unparsedArguments_, 0, "", true, warningLevel_);
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
    catch (bad_cast&)
    {
      ApplicationTools::displayMessage("Frequencies initialization not possible for model " + model.getName());
    }
    unparsedArguments_.erase(unparsedArguments_.find(model.getNamespace() + "initFreqs"));
    if (unparsedArguments_.find(model.getNamespace() + "initFreqs.observedPseudoCount") != unparsedArguments_.end())
      unparsedArguments_.erase(unparsedArguments_.find(model.getNamespace() + "initFreqs.observedPseudoCount"));
  }

  ParameterList pl = model.getIndependentParameters();
  for (size_t i = 0; i < pl.size(); ++i)
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
    bool test2 = (pName.size() < posp + 6) || (pName.substr(posp + 1, 5) != "theta");
    bool test3 = (unparsedArguments_.find(pName) != unparsedArguments_.end());
    try
    {
      if (test1 || test2 || test3)
      {
        if (!test1 && !test2 && test3)
          ApplicationTools::displayWarning("Warning, initFreqs argument is set and a value is set for parameter " + pName);

        double value = ApplicationTools::getDoubleParameter(pName, unparsedArguments_, pl[i].getValue(), "", true, warningLevel_);
        pl[i].setValue(value);

        if (unparsedArguments_.find(pName) != unparsedArguments_.end())
          unparsedArguments_.erase(unparsedArguments_.find(pName));
      }
      if (verbose_)
        ApplicationTools::displayResult("Parameter found", pName + "=" + TextTools::toString(pl[i].getValue()));
    }

    catch (Exception& e)
    {}
  }

  model.matchParametersValues(pl);
}
