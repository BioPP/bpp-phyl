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
#include "../Model/Codon/DFP07.h"
#include "../Model/Codon/RELAX.h"
#include "../Model/Codon/YNGP_M1.h"
#include "../Model/Codon/YNGP_M10.h"
#include "../Model/Codon/YNGP_M2.h"
#include "../Model/Codon/YNGP_M3.h"
#include "../Model/Codon/YNGP_M7.h"
#include "../Model/Codon/YNGP_M8.h"
#include "../Model/Codon/YNGP_M9.h"
#include "../Model/MixedTransitionModel.h"
#include "../Model/MixtureOfATransitionModel.h"
#include "../Model/MixtureOfTransitionModels.h"
#include "../Model/OneChangeRegisterTransitionModel.h"
#include "../Model/OneChangeTransitionModel.h"
#include "../Model/Protein/LG10_EX_EHO.h"
#include "../Model/Protein/LGL08_CAT.h"
#include "../Model/Protein/LLG08_EHO.h"
#include "../Model/Protein/LLG08_EX2.h"
#include "../Model/Protein/LLG08_EX3.h"
#include "../Model/Protein/LLG08_UL2.h"
#include "../Model/Protein/LLG08_UL3.h"
#include "BppOFrequencySetFormat.h"
#include "BppORateDistributionFormat.h"
#include "BppOTransitionModelFormat.h"

// From Numeric

#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/AutoParameter.h>


#include <Bpp/Text/StringTokenizer.h>

using namespace bpp;

// From the STL:
#include <iomanip>
#include <set>

using namespace std;

unique_ptr<TransitionModelInterface> BppOTransitionModelFormat::readTransitionModel(
    std::shared_ptr<const Alphabet> alphabet,
    const std::string& modelDescription,
    const std::map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mData,
    size_t nData,
    bool parseArguments)
{
  unparsedArguments_.clear();
  unique_ptr<TransitionModelInterface> model;
  string modelName = "";
  map<string, string> args;
  KeyvalTools::parseProcedure(modelDescription, modelName, args);

  // get data number

  if (args.find("data") != args.end())
    nData = TextTools::to<size_t>(args["data"]);

  // //////////////////////////////////
  // / MIXED MODELS
  // ////////////////////////////////

  if ((modelName == "MixedModel" || (modelName == "Mixture")) && allowMixed_)
    model = readMixed_(alphabet, modelDescription, mData, nData);
  else if (modelName == "OneChange")
  {
    // We have to parse the nested model first:
    if (args.find("model") == args.end())
      throw Exception("BppOTransitionModelFormat::read. Missing argument 'model' for model 'OneChange'.");
    string nestedModelDescription = args["model"];
    BppOSubstitutionModelFormat nestedReader(ALL, false, allowMixed_, allowGaps_, verbose_, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_);

    auto nestedModel = nestedReader.readSubstitutionModel(alphabet, nestedModelDescription, mData, nData, false);
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    // We look for the register:
    if (args.find("register") == args.end())
      model = make_unique<OneChangeTransitionModel>(std::move(nestedModel));
    else
    {
      shared_ptr<AlphabetIndex2> weights;
      shared_ptr<AlphabetIndex2> distances;
      string registerDescription = args["register"];
      auto reg = PhylogeneticsApplicationTools::getSubstitutionRegister(registerDescription, nestedModel->getStateMap(), geneticCode_, weights, distances);

      if (args.find("numReg") == args.end())
        throw Exception("Missing argument 'numReg' (number of event for register in model " + modelName);

      vector<size_t> vNumRegs;

      StringTokenizer nst(args["numReg"], "+");

      bool out = true;

      while (nst.hasMoreToken())
      {
        size_t n = TextTools::to<size_t>(nst.nextToken());
        vNumRegs.push_back(n);
        if (verbose_)
        {
          ApplicationTools::displayResult(out ? "Register types" : "", reg->getTypeName(n));
          out = false;
        }
      }

      model = make_unique<OneChangeRegisterTransitionModel>(std::move(nestedModel), *reg, vNumRegs);
    }

    // Then we update the parameter set:
    for (auto& it : unparsedParameterValuesNested)
    {
      unparsedArguments_["OneChange." + it.first] = it.second;
    }
  }
  // //////////////////
  // PREDEFINED CODON MODELS
  else if (((modelName.substr(0, 4) == "YNGP") || (modelName == "DFP07") || (modelName == "RELAX")) && (alphabetCode_ & CODON))
  {
    if (!(alphabetCode_ & CODON))
      throw Exception("BppOTransitionModelFormat::read. Codon alphabet not supported.");
    if (!geneticCode_)
      throw Exception("BppOTransitionModelFormat::readTransitionModel(). No genetic code specified! Consider using 'setGeneticCode'.");

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


    unique_ptr<CodonFrequencySetInterface> codonFreqs;

    if (args.find("frequencies") != args.end())
    {
      string freqOpt = ApplicationTools::getStringParameter("frequencies", args, "F0", "", true, warningLevel_);
      BppOFrequencySetFormat freqReader(BppOFrequencySetFormat::ALL, verbose_, warningLevel_);
      freqReader.setGeneticCode(geneticCode_); // This uses the same instance as the one that will be used by the model.

      auto tmpP = freqReader.readFrequencySet(pCA, freqOpt, mData, nData, false);
      codonFreqs = unique_ptr<CodonFrequencySetInterface>(dynamic_cast<CodonFrequencySetInterface*>(tmpP.release()));
      for (const auto& it : freqReader.getUnparsedArguments())
      {
        unparsedArguments_[modelName + "." + it.first] = it.second;
      }
    }
    else
      throw Exception("Missing 'frequencies' for model " + modelName);

    if (modelName == "YNGP_M1")
      model = make_unique<YNGP_M1>(geneticCode_, std::move(codonFreqs));
    else if (modelName == "YNGP_M2")
      model = make_unique<YNGP_M2>(geneticCode_, std::move(codonFreqs));
    else if (modelName == "RELAX")
      model = make_unique<RELAX>(geneticCode_, std::move(codonFreqs));
    else if (modelName == "YNGP_M3")
      if (args.find("n") == args.end())
        model = make_unique<YNGP_M3>(geneticCode_, std::move(codonFreqs));
      else
        model = make_unique<YNGP_M3>(geneticCode_, std::move(codonFreqs), TextTools::to<unsigned int>(args["n"]));
    else if ((modelName == "YNGP_M7") || modelName == "YNGP_M8" || modelName == "YNGP_M8a")
    {
      if (args.find("n") == args.end())
        throw Exception("Missing argument 'n' (number of classes) in " + modelName + " distribution");
      unsigned int nbClasses = TextTools::to<unsigned int>(args["n"]);
      if (verbose_)
        ApplicationTools::displayResult("Number of classes in model", nbClasses);

      if (modelName == "YNGP_M7")
        model = make_unique<YNGP_M7>(geneticCode_, std::move(codonFreqs), nbClasses);
      else if (modelName == "YNGP_M8")
        model = make_unique<YNGP_M8>(geneticCode_, std::move(codonFreqs), nbClasses);
      else if (modelName == "YNGP_M8a")
        model = make_unique<YNGP_M8>(geneticCode_, std::move(codonFreqs), nbClasses, true);
    }
    else if (modelName == "YNGP_M9" || modelName == "YNGP_M10")
    {
      if (args.find("nbeta") == args.end())
        throw Exception("Missing argument 'nbeta' (number of classes of beta distribution) in " + modelName + " distribution");
      unsigned int nbBeta = TextTools::to<unsigned int>(args["nbeta"]);
      if (args.find("ngamma") == args.end())
        throw Exception("Missing argument 'ngamma' (number of classes of gamma distribution) in " + modelName + " distribution");
      unsigned int nbGamma = TextTools::to<unsigned int>(args["ngamma"]);
      if (verbose_)
        ApplicationTools::displayResult("Number of classes in model", nbBeta + nbGamma);

      if (modelName == "YNGP_M9")
        model = make_unique<YNGP_M9>(geneticCode_, std::move(codonFreqs), nbBeta, nbGamma);
      else
        model = make_unique<YNGP_M10>(geneticCode_, std::move(codonFreqs), nbBeta, nbGamma);
    }
    else if (modelName == "DFP07")
    {
      if (args.find("protmodel") == args.end())
        throw Exception("Missing 'protmodel in model " + modelName + ".");

      BppOSubstitutionModelFormat nestedProtReader(PROTEIN, false, allowMixed_, allowGaps_, verbose_, warningLevel_);
      auto tmpP = nestedProtReader.readSubstitutionModel(
            geneticCode_->getTargetAlphabet(),
            args["protmodel"], mData, nData, false);
      auto nestedProtModel = unique_ptr<ProteinSubstitutionModelInterface>(
            dynamic_cast<ProteinSubstitutionModelInterface*>(tmpP.release())
            );

      auto unparsedParameterValuesNested  = nestedProtReader.getUnparsedArguments();
      unparsedArguments_.insert(unparsedParameterValuesNested.begin(), unparsedParameterValuesNested.end());

      model = make_unique<DFP07>(geneticCode_, std::move(nestedProtModel), std::move(codonFreqs));
    }
  }
  else if (AlphabetTools::isProteicAlphabet(*alphabet))
  {
    if (!(alphabetCode_ & PROTEIN))
      throw Exception("BppOTransitionModelFormat::read. Protein alphabet not supported.");
    auto alpha = dynamic_pointer_cast<const ProteicAlphabet>(alphabet);

    if (modelName == "LLG08_EHO")
      model = make_unique<LLG08_EHO>(alpha);
    else if (modelName == "LLG08_EX2")
      model = make_unique<LLG08_EX2>(alpha);
    else if (modelName == "LLG08_EX3")
      model = make_unique<LLG08_EX3>(alpha);
    else if (modelName == "LLG08_UL2")
      model = make_unique<LLG08_UL2>(alpha);
    else if (modelName == "LLG08_UL3")
      model = make_unique<LLG08_UL3>(alpha);
    else if (modelName == "LG10_EX_EHO")
      model = make_unique<LG10_EX_EHO>(alpha);
    else if (modelName == "LGL08_CAT")
    {
      if (args.find("nbCat") == args.end())
        throw Exception("'nbCat' argument is compulsory for model 'LGL08_CAT'");

      unsigned int nbCat = TextTools::to<unsigned int>(args["nbCat"]);
      if (nbCat == 0)
        throw Exception("nbCat argument has to be 10, 20, 30, 40, 50 or 60.");
      model = make_unique<LGL08_CAT>(alpha, nbCat);
    }
  }

  if (!model)
    model = readSubstitutionModel(alphabet, modelDescription, mData, nData, parseArguments);
  else
  {
    if (verbose_)
      ApplicationTools::displayResult("Transition model", modelName);

    updateParameters_(*model, args);

    if (parseArguments)
    {
      if (nData)
        initialize_(*model, mData.at(nData));
      else
        initialize_(*model, 0);
    }
  }

  return model;
}

unique_ptr<MixedTransitionModelInterface> BppOTransitionModelFormat::readMixed_(
    std::shared_ptr<const Alphabet> alphabet,
    const std::string& modelDescription,
    const std::map<size_t, std::shared_ptr<const AlignmentDataInterface>>& mData,
    size_t nData)
{
  unique_ptr<MixedTransitionModelInterface> model;

  string modelName = "";
  map<string, string> args;
  KeyvalTools::parseProcedure(modelDescription, modelName, args);
  unique_ptr<TransitionModelInterface> pSM;

  if (modelName == "MixedModel")
  {
    if (args.find("model") == args.end())
      throw Exception("The argument 'model' is missing from MixedModel description");
    string nestedModelDescription = args["model"];
    BppOTransitionModelFormat nestedReader(alphabetCode_, allowCovarions_, true, allowGaps_, false, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_); // This uses the same
    // instance as the one
    // that will be used
    // by the model.
    pSM = nestedReader.readTransitionModel(alphabet, nestedModelDescription, mData, nData, false);

    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    map<string, unique_ptr<DiscreteDistributionInterface>> mdist;
    map<string, string> unparsedParameterValuesNested2;

    for (auto& it : unparsedParameterValuesNested)
    {
      if (it.second.find("(") != string::npos)
      {
        BppODiscreteDistributionFormat bIO(false);
        mdist[pSM->getParameterNameWithoutNamespace(it.first)] = bIO.readDiscreteDistribution(it.second, false);
        map<string, string> unparsedParameterValuesNested3(bIO.getUnparsedArguments());
        for (auto& it2 : unparsedParameterValuesNested3)
        {
          unparsedParameterValuesNested2[it.first + "_" + it2.first] = it2.second;
        }
      }
      else
        unparsedParameterValuesNested2[it.first] = it.second;
    }

    for (auto& it : unparsedParameterValuesNested2)
    {
      unparsedArguments_[it.first] = it.second;
    }


    int fi(-1), ti(-1);

    if (args.find("from") != args.end())
      fi = alphabet->charToInt(args["from"]);
    if (args.find("to") != args.end())
      ti = alphabet->charToInt(args["to"]);

    string sModN = pSM->getName();
    model = make_unique<MixtureOfATransitionModel>(alphabet, std::move(pSM), mdist, fi, ti);

    vector<string> v = model->getParameters().getParameterNames();

    if (verbose_)
    {
      ApplicationTools::displayResult("Mixture Of A TransitionModel Model", sModN);
      ApplicationTools::displayResult("Number of classes", model->getNumberOfModels());
    }
  }
  else if (modelName == "Mixture")
  {
    vector<string> v_nestedModelDescription;
    vector< std::unique_ptr<TransitionModelInterface>> v_pSM;

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

    for (unsigned i = 0; i < v_nestedModelDescription.size(); ++i)
    {
      BppOTransitionModelFormat nestedReader(alphabetCode_, false, true, false, false, warningLevel_);
      if (geneticCode_)
        nestedReader.setGeneticCode(geneticCode_); // This uses the same instance as the one that will be used by the model.

      pSM = nestedReader.readTransitionModel(alphabet, v_nestedModelDescription[i], mData, nData, false);

      map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
      for (auto& it : unparsedParameterValuesNested)
      {
        unparsedArguments_[modelName + "." + TextTools::toString(i + 1) + "_" + it.first] = it.second;
      }

      v_pSM.push_back(std::move(pSM));
    }

    model = make_unique<MixtureOfTransitionModels>(alphabet, v_pSM);
    if (verbose_)
      ApplicationTools::displayResult("Mixture Of TransitionModel Models", modelName );
  }
  else
    throw Exception("Unknown model name for mixture " + modelName);

  return model;
}
