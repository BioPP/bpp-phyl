//
// File: BppOTransitionModelFormat.cpp
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

#include "BppOTransitionModelFormat.h"
#include "BppORateDistributionFormat.h"

#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/KeyvalTools.h>

#include "../Model/OneChangeTransitionModel.h"
#include "../Model/OneChangeRegisterTransitionModel.h"

#include "../App/PhylogeneticsApplicationTools.h"

#include "BppOFrequencySetFormat.h"

#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>

#include <Bpp/Io/OutputStream.h>
#include <Bpp/Io/BppOParametrizableFormat.h>
#include <Bpp/Io/BppODiscreteDistributionFormat.h>

#include "../Model/MixedTransitionModel.h"
#include "../Model/MixtureOfATransitionModel.h"
#include "../Model/MixtureOfTransitionModels.h"

#include "../Model/Codon/DFP07.h"
#include "../Model/Codon/YNGP_M1.h"
#include "../Model/Codon/YNGP_M2.h"
#include "../Model/Codon/YNGP_M3.h"
#include "../Model/Codon/YNGP_M7.h"
#include "../Model/Codon/YNGP_M8.h"
#include "../Model/Codon/YNGP_M9.h"
#include "../Model/Codon/YNGP_M10.h"
#include "../Model/Codon/RELAX.h"
#include "../Model/Protein/LLG08_EX2.h"
#include "../Model/Protein/LLG08_EX3.h"
#include "../Model/Protein/LLG08_UL2.h"
#include "../Model/Protein/LLG08_UL3.h"
#include "../Model/Protein/LGL08_CAT.h"
#include "../Model/Protein/LLG08_EHO.h"
#include "../Model/Protein/LG10_EX_EHO.h"

//From Numeric

#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/AutoParameter.h>


#include <Bpp/Text/StringTokenizer.h>

using namespace bpp;

// From the STL:
#include <iomanip>
#include <set>

using namespace std;

TransitionModel* BppOTransitionModelFormat::readTransitionModel(
  const Alphabet* alphabet,
    const std::string& modelDescription,
    const AlignedValuesContainer* data,
    bool parseArguments)
{
  unparsedArguments_.clear();
  unique_ptr<TransitionModel> model;
  string modelName = "";
  map<string, string> args;
  KeyvalTools::parseProcedure(modelDescription, modelName, args);

  // //////////////////////////////////
  // / MIXED MODELS
  // ////////////////////////////////

  if ((modelName == "MixedModel" || (modelName == "Mixture")) && allowMixed_)
    model.reset(readMixed_(alphabet, modelDescription, data));
  else if (modelName == "OneChange")
  {
    // We have to parse the nested model first:
    if (args.find("model")==args.end())
      throw Exception("BppOTransitionModelFormat::read. Missing argument 'model' for model 'OneChange'.");
    string nestedModelDescription = args["model"];
    BppOSubstitutionModelFormat nestedReader(ALL, false, allowMixed_, allowGaps_, verbose_, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_);
    
    SubstitutionModel* nestedModel=nestedReader.readSubstitutionModel(alphabet, nestedModelDescription, data, false);
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    // We look for the register:
    if (args.find("register")==args.end())
      model.reset(new OneChangeTransitionModel(*nestedModel));
    else
    {
      AlphabetIndex2* weights=0;
      AlphabetIndex2* distances=0;
      
      string registerDescription = args["register"];
      unique_ptr<SubstitutionRegister> reg(PhylogeneticsApplicationTools::getSubstitutionRegister(registerDescription, nestedModel->getStateMap(), geneticCode_, weights, distances));

      if (args.find("numReg") == args.end())
        throw Exception("Missing argument 'numReg' (number of event for register in model " + modelName);

      vector<size_t> vNumRegs;
      
      StringTokenizer nst(args["numReg"], "+");

      bool out=true;
      
      while (nst.hasMoreToken())
      {
        size_t n=TextTools::to<size_t>(nst.nextToken());
        vNumRegs.push_back(n);
        if (verbose_)
        {
          ApplicationTools::displayResult(out?"Register types":"", reg->getTypeName(n));
          out=false;
        }
      }

      model.reset(new OneChangeRegisterTransitionModel(*nestedModel, *reg, vNumRegs));
    }

    // Then we update the parameter set:
    for (auto&  it:unparsedParameterValuesNested)
    {
      unparsedArguments_["OneChange." + it.first] = it.second;
    }
    delete nestedModel;
  }
  // //////////////////
  // PREDEFINED CODON MODELS
  else if (((modelName.substr(0, 4) == "YNGP") || (modelName == "DFP07") || (modelName == "RELAX")) && (alphabetCode_ & CODON))
  {
    if (!(alphabetCode_ & CODON))
      throw Exception("BppOTransitionModelFormat::read. Codon alphabet not supported.");
    if (!geneticCode_)
      throw Exception("BppOTransitionModelFormat::readTransitionModel(). No genetic code specified! Consider using 'setGeneticCode'.");
  
    if (!AlphabetTools::isCodonAlphabet(alphabet))
      throw Exception("Alphabet should be Codon Alphabet.");
    const CodonAlphabet* pCA = dynamic_cast<const CodonAlphabet*>(alphabet);

    if (args.find("genetic_code") != args.end()) {
      ApplicationTools::displayWarning("'genetic_code' argument is no longer supported inside model description, and has been supersided by a global 'genetic_code' option.");
      throw Exception("BppOSubstitutionModelFormat::read. Deprecated 'genetic_code' argument.");
    }

    if (geneticCode_->getSourceAlphabet()->getAlphabetType() != pCA->getAlphabetType())
      throw Exception("Mismatch between genetic code and codon alphabet");


    shared_ptr<CodonFrequencySet> codonFreqs(0);

    if (args.find("frequencies")!=args.end())
    {
      string freqOpt = ApplicationTools::getStringParameter("frequencies", args, "F0", "", true, warningLevel_);
      BppOFrequencySetFormat freqReader(BppOFrequencySetFormat::ALL, verbose_, warningLevel_);
      freqReader.setGeneticCode(geneticCode_); //This uses the same instance as the one that will be used by the model.

      codonFreqs = std::dynamic_pointer_cast<CodonFrequencySet>(freqReader.readFrequencySet(pCA, freqOpt, data, false));
      for (const auto& it:freqReader.getUnparsedArguments())
        unparsedArguments_[modelName + "." + it.first] = it.second;
    }
    else 
      // codonFreqs compulsory for all models but SameAARate
      if (modelName!="DFP07")
        throw Exception("Missing 'frequencies' for model " + modelName);

    if (modelName == "YNGP_M1")
      model.reset(new YNGP_M1(geneticCode_, codonFreqs));
    else if (modelName == "YNGP_M2")
      model.reset(new YNGP_M2(geneticCode_, codonFreqs));
    else if (modelName == "RELAX")
      model.reset(new RELAX(geneticCode_, codonFreqs));
    else if (modelName == "YNGP_M3")
      if (args.find("n") == args.end())
        model.reset(new YNGP_M3(geneticCode_, codonFreqs));
      else
        model.reset(new YNGP_M3(geneticCode_, codonFreqs, TextTools::to<unsigned int>(args["n"])));
    else if ((modelName == "YNGP_M7") || modelName == "YNGP_M8")
    {
      if (args.find("n") == args.end())
        throw Exception("Missing argument 'n' (number of classes) in " + modelName + " distribution");
      unsigned int nbClasses = TextTools::to<unsigned int>(args["n"]);
      if (verbose_)
        ApplicationTools::displayResult("Number of classes in model", nbClasses);

      if (modelName == "YNGP_M7")
        model.reset(new YNGP_M7(geneticCode_, codonFreqs, nbClasses));
      else if (modelName == "YNGP_M8")
        model.reset(new YNGP_M8(geneticCode_, codonFreqs, nbClasses));
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
        model.reset(new YNGP_M9(geneticCode_, codonFreqs, nbBeta, nbGamma));
      else
        model.reset(new YNGP_M10(geneticCode_, codonFreqs, nbBeta, nbGamma));
    }
    else if (modelName == "DFP07")
    {
      if (args.find("protmodel") == args.end())
        throw Exception("Missing 'protmodel in model " + modelName + ".");

      BppOSubstitutionModelFormat nestedProtReader(PROTEIN, false, allowMixed_, allowGaps_, verbose_, warningLevel_);
      auto  nestedProtModel = std::shared_ptr<ProteinSubstitutionModel>(dynamic_cast<ProteinSubstitutionModel*>(nestedProtReader.readSubstitutionModel(geneticCode_->getTargetAlphabet(), args["protmodel"], data, false)));

      auto unparsedParameterValuesNested  = nestedProtReader.getUnparsedArguments();
      unparsedArguments_.insert(unparsedParameterValuesNested.begin(), unparsedParameterValuesNested.end());

      model.reset(new DFP07(geneticCode_, nestedProtModel, codonFreqs));
    }
  }
  else if (AlphabetTools::isProteicAlphabet(alphabet))
  {
    if (!(alphabetCode_ & PROTEIN))
      throw Exception("BppOTransitionModelFormat::read. Protein alphabet not supported.");
    const ProteicAlphabet* alpha = dynamic_cast<const ProteicAlphabet*>(alphabet);

    if (modelName == "LLG08_EHO")
      model.reset(new LLG08_EHO(alpha));
    else if (modelName == "LLG08_EX2")
      model.reset(new LLG08_EX2(alpha));
    else if (modelName == "LLG08_EX3")
      model.reset(new LLG08_EX3(alpha));
    else if (modelName == "LLG08_UL2")
      model.reset(new LLG08_UL2(alpha));
    else if (modelName == "LLG08_UL3")
      model.reset(new LLG08_UL3(alpha));
    else if (modelName == "LG10_EX_EHO")
      model.reset(new LG10_EX_EHO(alpha));	
    else if (modelName == "LGL08_CAT")
    {
      if (args.find("nbCat")==args.end())
        throw Exception("'nbCat' argument is compulsory for model 'LGL08_CAT'");
          
      unsigned int nbCat = TextTools::to<unsigned int>(args["nbCat"]);
      model.reset(new LGL08_CAT(alpha, nbCat));
    }
  }

  if (!model)
    model.reset(readSubstitutionModel(alphabet, modelDescription, data, parseArguments));
  else
  {
    if (verbose_)
      ApplicationTools::displayResult("Transition model", modelName);
  
    updateParameters_(model.get(), args);
  
    if (parseArguments)
      initialize_(*model, data);
  }
  
  return model.release();
}

MixedTransitionModel* BppOTransitionModelFormat::readMixed_(const Alphabet* alphabet, const std::string& modelDescription, const AlignedValuesContainer* data)
{
  unique_ptr<MixedTransitionModel> model;

  string modelName = "";
  map<string, string> args;
  KeyvalTools::parseProcedure(modelDescription, modelName, args);
  unique_ptr<TransitionModel> pSM;

  if (modelName == "MixedModel")
  {
    if (args.find("model") == args.end())
      throw Exception("The argument 'model' is missing from MixedModel description");
    string nestedModelDescription = args["model"];
    BppOTransitionModelFormat nestedReader(alphabetCode_, allowCovarions_, true, allowGaps_, false, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_); //This uses the same
    //instance as the one
    //that will be used
    //by the model.
    pSM.reset(nestedReader.readTransitionModel(alphabet, nestedModelDescription, data, false));

    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    map<string, DiscreteDistribution*> mdist;
    map<string, string> unparsedParameterValuesNested2;

    for (auto& it:unparsedParameterValuesNested)
    {
      if (it.second.find("(") != string::npos)
      {
        BppODiscreteDistributionFormat bIO(false);
        mdist[pSM->getParameterNameWithoutNamespace(it.first)] = bIO.readDiscreteDistribution(it.second, false);
        map<string, string> unparsedParameterValuesNested3(bIO.getUnparsedArguments());
        for (auto& it2:unparsedParameterValuesNested3)
          unparsedParameterValuesNested2[it.first + "_" + it2.first] = it2.second;
      }
      else
        unparsedParameterValuesNested2[it.first] = it.second;
    }

    for (auto& it:unparsedParameterValuesNested2)
      unparsedArguments_[it.first] = it.second;


    int fi(-1), ti(-1);

    if (args.find("from") != args.end())
      fi = alphabet->charToInt(args["from"]);
    if (args.find("to") != args.end())
      ti = alphabet->charToInt(args["to"]);

    string sModN=pSM->getName();
    model.reset(new MixtureOfATransitionModel(alphabet, pSM.release(), mdist, fi, ti));

    vector<string> v = model->getParameters().getParameterNames();

    for (auto&  it:mdist)
      delete it.second;

    if (verbose_)
    {
      ApplicationTools::displayResult("Mixture Of A TransitionModel Model", sModN);
      ApplicationTools::displayResult("Number of classes", model->getNumberOfModels());
    }
  }
  else if (modelName == "Mixture")
  {
    vector<string> v_nestedModelDescription;
    vector<TransitionModel*> v_pSM;

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
      BppOTransitionModelFormat nestedReader(alphabetCode_, false, true, false, false, warningLevel_);
      if (geneticCode_)
        nestedReader.setGeneticCode(geneticCode_); //This uses the same instance as the one that will be used by the model.

      pSM.reset(nestedReader.readTransitionModel(alphabet, v_nestedModelDescription[i], data, false));

      map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
      for (auto& it:unparsedParameterValuesNested)
        unparsedArguments_[modelName + "." + TextTools::toString(i + 1) + "_" + it.first] = it.second;

      v_pSM.push_back(pSM.release());
    }

    model.reset(new MixtureOfTransitionModels(alphabet, v_pSM));
    if (verbose_)
      ApplicationTools::displayResult("Mixture Of TransitionModel Models", modelName );
  }
  else
    throw Exception("Unknown model name for mixture " + modelName);

  return model.release();
}

