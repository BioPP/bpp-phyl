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

#include "../Model/WordSubstitutionModel.h"
#include "../Model/KroneckerWordSubstitutionModel.h"
#include "../Model/Codon/MG94.h"
#include "../Model/Codon/GY94.h"
#include "../Model/Codon/YN98.h"
#include "../Model/Codon/TripletSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonDistanceSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonAARateSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonAAFitnessSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonClusterAASubstitutionModel.h"
#include "../Model/Codon/AbstractCodonCpGSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonBGCSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonDistanceSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonFitnessSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonFrequenciesSubstitutionModel.h"
#include "../Model/Codon/AbstractCodonPhaseFrequenciesSubstitutionModel.h"
#include "../Model/Codon/CodonAdHocSubstitutionModel.h"
#include "../Model/Codon/KroneckerCodonDistanceFrequenciesSubstitutionModel.h"
#include "../Model/Codon/KroneckerCodonDistanceSubstitutionModel.h"
#include "../Model/Codon/KCM.h"
#include "../Model/Codon/SENCA.h"
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
#include "../Model/Protein/Coala.h"
#include "../Model/Protein/DSO78.h"
#include "../Model/Protein/JCprot.h"
#include "../Model/Protein/JTT92.h"
#include "../Model/Protein/ProteinSubstitutionModel.h"
#include "../Model/Protein/LG08.h" 
#include "../Model/Protein/UserProteinSubstitutionModel.h"
#include "../Model/Protein/WAG01.h"
#include "../Model/BinarySubstitutionModel.h"
#include "../Model/FromMixtureSubstitutionModel.h"
#include "../Model/AbstractBiblioMixedTransitionModel.h"
#include "../Model/InMixedSubstitutionModel.h"
#include "../Model/MixtureOfTransitionModels.h"
#include "../Model/MixtureOfATransitionModel.h"
#include "../Model/OneChangeTransitionModel.h"
#include "../Model/OneChangeRegisterTransitionModel.h"
#include "../Model/RegisterRatesSubstitutionModel.h"
#include "../Model/Codon/YNGP_M7.h"
#include "../Model/Codon/YNGP_M8.h"
#include "../Model/Codon/YNGP_M9.h"
#include "../Model/Codon/YNGP_M10.h"

#include "../App/PhylogeneticsApplicationTools.h"

#include "BppOFrequenciesSetFormat.h"
#include "BppOTransitionModelFormat.h"

#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>

#include <Bpp/Io/OutputStream.h>
#include <Bpp/Io/BppOParametrizableFormat.h>
#include <Bpp/Io/BppODiscreteDistributionFormat.h>

//From Numeric

#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/AutoParameter.h>


#include <Bpp/Text/StringTokenizer.h>


using namespace bpp;

// From the STL:
#include <iomanip>
#include <set>

using namespace std;

unsigned char BppOSubstitutionModelFormat::DNA = 1;
unsigned char BppOSubstitutionModelFormat::RNA = 2;
unsigned char BppOSubstitutionModelFormat::NUCLEOTIDE = 1 | 2;
unsigned char BppOSubstitutionModelFormat::PROTEIN = 4;
unsigned char BppOSubstitutionModelFormat::CODON = 8;
unsigned char BppOSubstitutionModelFormat::WORD = 16;
unsigned char BppOSubstitutionModelFormat::BINARY = 32;
unsigned char BppOSubstitutionModelFormat::ALL = 1 | 2 | 4 | 8 | 16 | 32;


SubstitutionModel* BppOSubstitutionModelFormat::readSubstitionModel(
  const Alphabet* alphabet,
  const std::string& modelDescription,
  const SiteContainer* data,
  bool parseArguments)
{
  unparsedArguments_.clear();
  unique_ptr<SubstitutionModel> model;
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

    string nameMod="";
    size_t numMod=0;
    
    if (args.find("numMod") == args.end())
    {
      if (args.find("nameMod") == args.end())
        throw Exception("'numMod' and 'nameMod' arguments missing to define the InMixed submodel.");
      else
        nameMod=args["nameMod"];
    }
    else
      numMod=(size_t)TextTools::toInt(args["numMod"]);
    
    string modelDesc2=args["model"];
    BppOTransitionModelFormat nestedReader(alphabetCode_, false, true, false, false, warningLevel_);
    nestedReader.setGeneticCode(geneticCode_); //This uses the same instance as the o

    MixedTransitionModel* nestedModel=dynamic_cast<MixedTransitionModel*>(nestedReader.readTransitionModel(alphabet, modelDesc2, data, false));

    // Check that nestedModel is fine and has subModel of given name
    // or number
    
    if (nestedModel==NULL)
      throw Exception("Unknown mixed model " + modelDesc2 + ".");
    
    if (nameMod!="")
    {
      if (nestedModel->getModel(nameMod)==NULL)
        throw Exception("BppOSubstitutionModelFormat::read. " + nestedModel->getName() + "argument for model 'InMixed' has no submodel with name " + nameMod + ".");
      model.reset(new InMixedSubstitutionModel(*nestedModel, nameMod, modelDesc2));
    }
    else
    {
      if (numMod==0 || (nestedModel->getNumberOfModels()<numMod))
        throw Exception("BppOSubstitutionModelFormat::read. " + nestedModel->getName() + "argument for model 'InMixed' has no submodel with number " + TextTools::toString(numMod) + ".");
      model.reset(new InMixedSubstitutionModel(*nestedModel, numMod-1, modelDesc2));      
    }

    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedArguments_[it->first] = it->second;
    }

    if (verbose_)
      ApplicationTools::displayResult("Number of submodel", TextTools::toString(numMod));

    delete nestedModel;
  }


  else if (modelName == "FromRegister")
  {
    // We have to parse the nested model first:
    if (args.find("model")==args.end())
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'model' for model 'FromRegister'.");

    string nestedModelDescription = args["model"];
    BppOSubstitutionModelFormat nestedReader(ALL, false, allowMixed_, allowGaps_, verbose_, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_);
    
    SubstitutionModel* nestedModel=nestedReader.readSubstitionModel(alphabet, nestedModelDescription, data, false);
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    // We look for the register:
    if (args.find("register")==args.end())
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'register' for model 'FromRegister'.");
    
    string registerDescription = args["register"];
    SubstitutionRegister* reg=PhylogeneticsApplicationTools::getSubstitutionRegister(registerDescription, nestedModel);

    // is it normalized (default : false)
    bool isNorm=false;
    
    if (args.find("isNormalized")!=args.end() && args["isNormalized"]=="true")
      isNorm=true;

    model.reset(new RegisterRatesSubstitutionModel(*nestedModel, *reg, isNorm));
    
    if (TextTools::isEmpty(registerDescription))
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'register' for model 'FromRegister'.");

    // Then we update the parameter set:
    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedArguments_["FromRegister."+it->first] = it->second;
    }
    
    delete nestedModel;
  }
  
  
  // /////////////////////////////////
  // / WORDS and CODONS
  // ///////////////////////////////

  else if ((modelName == "Word") || (modelName.substr(0,4) == "Kron") || (modelName == "Triplet") || (modelName.substr(0, 5) == "Codon") || (modelName == "SENCA") )
    model.reset(readWord_(alphabet, modelDescription, data));

  
  // //////////////////
  // PREDEFINED CODON MODELS
  
  else if (((modelName == "MG94") || (modelName == "YN98") ||
            (modelName == "GY94") ||  (modelName.substr(0,3) == "KCM"))
           && (alphabetCode_ & CODON))
  {
    if (!(alphabetCode_ & CODON))
      throw Exception("BppOSubstitutionModelFormat::read. Codon alphabet not supported.");
    if (!geneticCode_)
      throw Exception("BppOSubstitutionModelFormat::readSubstitionModel(). No genetic code specified! Consider using 'setGeneticCode'.");
    
    if (!AlphabetTools::isCodonAlphabet(alphabet))
      throw Exception("Alphabet should be Codon Alphabet.");

    const CodonAlphabet* pCA = dynamic_cast<const CodonAlphabet*>(alphabet);

    if (args.find("genetic_code") != args.end()) {
      ApplicationTools::displayWarning("'genetic_code' argument is no longer supported inside model description, and has been supersided by a global 'genetic_code' option.");
      throw Exception("BppOSubstitutionModelFormat::read. Deprecated 'genetic_code' argument.");
    }

    if (geneticCode_->getSourceAlphabet()->getAlphabetType() != pCA->getAlphabetType())
      throw Exception("Mismatch between genetic code and codon alphabet");

    string freqOpt = ApplicationTools::getStringParameter("frequencies", args, "F0", "", true, warningLevel_);
    BppOFrequenciesSetFormat freqReader(BppOFrequenciesSetFormat::ALL, verbose_, warningLevel_);
    freqReader.setGeneticCode(geneticCode_); //This uses the same instance as the one that will be used by the model.
    unique_ptr<FrequenciesSet> codonFreqs(freqReader.read(pCA, freqOpt, data, false));
    map<string, string> unparsedParameterValuesNested(freqReader.getUnparsedArguments());

    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedArguments_[modelName + "." + it->first] = it->second;
    }

    if (modelName == "MG94")
      model.reset(new MG94(geneticCode_, codonFreqs.release()));
    else if (modelName == "GY94")
      model.reset(new GY94(geneticCode_, codonFreqs.release()));
    else if ((modelName == "YN98") || (modelName == "YNGP_M0"))
      model.reset(new YN98(geneticCode_, codonFreqs.release()));
    else if (modelName == "KCM7")
      model.reset(new KCM(geneticCode_, true));
    else if (modelName == "KCM19")
      model.reset(new KCM(geneticCode_, false));
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
    const RNY* prny = dynamic_cast<const RNY*>(alphabet);

    string nestedModelDescription = args["model"];
    if (TextTools::isEmpty(nestedModelDescription))
      throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'model' for model 'YpR_sym'.");
    if (verbose_)
      ApplicationTools::displayResult("Symetric YpR model", modelName);
    BppOSubstitutionModelFormat nestedReader(NUCLEOTIDE, false, false, false, verbose_, warningLevel_);
    unique_ptr<NucleotideSubstitutionModel> nestedModel(dynamic_cast<NucleotideSubstitutionModel*>(nestedReader.readSubstitionModel(&prny->getLetterAlphabet(), nestedModelDescription, data, false)));
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
    BppOSubstitutionModelFormat nestedReader(NUCLEOTIDE, false, false, false, verbose_, warningLevel_);
    unique_ptr<NucleotideSubstitutionModel> nestedModel(dynamic_cast<NucleotideSubstitutionModel*>(nestedReader.readSubstitionModel(&prny->getLetterAlphabet(), nestedModelDescription, data, false)));
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
    BppOSubstitutionModelFormat nestedReader(ALL, allowCovarions_, false, false, verbose_, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_);
    unique_ptr<ReversibleSubstitutionModel> nestedModel(dynamic_cast<ReversibleSubstitutionModel*>(nestedReader.readSubstitionModel(alphabet, nestedModelDescription, data, false)));
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    // Now we create the RE08 substitution model:
    if (AlphabetTools::isNucleicAlphabet(alphabet))
    {
      if (!(alphabetCode_ & NUCLEOTIDE))
        throw Exception("BppOSubstitutionModelFormat::read. Nucleic alphabet not supported.");
      model.reset(new RE08Nucleotide(dynamic_cast<NucleotideReversibleSubstitutionModel*>(nestedModel.release())));
      if (!model.get())
        throw Exception("BppOSubstitutionModelFormat::readSubstitionModel(). Invalid submodel, must be 'reversible' and 'nucleotide'.");
    }
    else if (AlphabetTools::isProteicAlphabet(alphabet))
    {
      if (!(alphabetCode_ & PROTEIN))
        throw Exception("BppOSubstitutionModelFormat::read. Protein alphabet not supported.");
      model.reset(new RE08Protein(dynamic_cast<ProteinReversibleSubstitutionModel*>(nestedModel.release())));
      if (!model.get())
        throw Exception("BppOSubstitutionModelFormat::readSubstitionModel(). Invalid submodel, must be 'reversible' and 'protein'.");
    }
    else if (AlphabetTools::isCodonAlphabet(alphabet))
    {
      if (!(alphabetCode_ & CODON))
        throw Exception("BppOSubstitutionModelFormat::read. Codon alphabet not supported.");
      model.reset(new RE08Codon(dynamic_cast<CodonReversibleSubstitutionModel*>(nestedModel.release())));
      if (!model.get())
        throw Exception("BppOSubstitutionModelFormat::readSubstitionModel(). Invalid submodel, must be 'reversible' and 'codon'.");
    }
    else
    {
      model.reset(new RE08(nestedModel.release()));
    }

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
    BppOSubstitutionModelFormat nestedReader(ALL, false, allowMixed_, allowGaps_, false, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_);
    unique_ptr<ReversibleSubstitutionModel> nestedModel(dynamic_cast<ReversibleSubstitutionModel*>(nestedReader.readSubstitionModel(alphabet, nestedModelDescription, data, false)));
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
    BppOSubstitutionModelFormat nestedReader(ALL, false, allowMixed_, allowGaps_, verbose_, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_);
    unique_ptr<ReversibleSubstitutionModel> nestedModel(dynamic_cast<ReversibleSubstitutionModel*>(nestedReader.readSubstitionModel(alphabet, nestedModelDescription, data, false)));
    map<string, string> unparsedParameterValuesNestedModel(nestedReader.getUnparsedArguments());
    BppORateDistributionFormat rateReader(false);
    unique_ptr<DiscreteDistribution> nestedRDist(rateReader.read(nestedRateDistDescription, false));
    map<string, string> unparsedParameterValuesNestedDist(rateReader.getUnparsedArguments());

    // Now we create the G01 substitution model:
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
    if (AlphabetTools::isNucleicAlphabet(alphabet))
    {
      if (!(alphabetCode_ & NUCLEOTIDE))
        throw Exception("BppOSubstitutionModelFormat::read. Nucleotide alphabet not supported.");
      const NucleicAlphabet* alpha = dynamic_cast<const NucleicAlphabet*>(alphabet);

      // //////////////////////////////////
      // gBGC
      // //////////////////////////////////
      if (modelName.find("+gBGC") != string::npos)
      {
        string subModName=modelName.substr(0,modelName.find("+gBGC"));

        if (verbose_)
          ApplicationTools::displayResult("Biased gene conversion for", subModName);

        string nestedModelDescription = subModName;
        
        if (modelDescription.find_first_of("(")!=string::npos)
        {
          string::size_type begin = modelDescription.find_first_of("(");
          string::size_type end = modelDescription.find_last_of(")");

          nestedModelDescription += modelDescription.substr(begin,end-begin+1);
        }
        
        BppOSubstitutionModelFormat nestedReader(NUCLEOTIDE, true, true, false, verbose_, warningLevel_);
        unique_ptr<NucleotideSubstitutionModel> nestedModel(dynamic_cast<NucleotideSubstitutionModel*>(nestedReader.readSubstitionModel(alphabet, nestedModelDescription, data, false)));
        map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

        // Now we create the gBGC substitution model:
        model.reset(new gBGC(alpha, nestedModel.release()));

        // Then we update the parameter set:
        for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
        {
          unparsedArguments_["gBGC." + it->first] = it->second;
        }
      }


      // /////////////////////////////////
      // / GTR
      // ///////////////////////////////

      else if (modelName == "GTR")
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
        throw Exception("Model '" + modelName + "' unknown, or does not fit nucleic alphabet.");
    }
    else if (AlphabetTools::isProteicAlphabet(alphabet))
    {
      if (!(alphabetCode_ & PROTEIN))
        throw Exception("BppOSubstitutionModelFormat::read. Protein alphabet not supported.");
      const ProteicAlphabet* alpha = dynamic_cast<const ProteicAlphabet*>(alphabet);
      
      if (modelName.find("+F") != string::npos) {
        string freqOpt = ApplicationTools::getStringParameter("frequencies", args, "Full", "", true, warningLevel_);
        BppOFrequenciesSetFormat freqReader(BppOFrequenciesSetFormat::ALL, false, warningLevel_);
        unique_ptr<FrequenciesSet> protFreq(freqReader.read(alpha, freqOpt, data, true));
  
        map<string, string> unparsedParameterValuesNested(freqReader.getUnparsedArguments());

        for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
        {
          unparsedArguments_[modelName + "." + it->first] = it->second;
        }

        if (modelName == "JC69+F")
          model.reset(new JCprot(alpha, dynamic_cast<ProteinFrequenciesSet*>(protFreq.release())));
        else if (modelName == "DSO78+F")
          model.reset(new DSO78(alpha, dynamic_cast<ProteinFrequenciesSet*>(protFreq.release())));
        else if (modelName == "JTT92+F")
          model.reset(new JTT92(alpha, dynamic_cast<ProteinFrequenciesSet*>(protFreq.release())));
        else if (modelName == "LG08+F")
          model.reset(new LG08(alpha, dynamic_cast<ProteinFrequenciesSet*>(protFreq.release())));
        else if (modelName == "WAG01+F")
          model.reset(new WAG01(alpha, dynamic_cast<ProteinFrequenciesSet*>(protFreq.release())));
        else if (modelName == "Empirical+F")
        {
          string prefix = args["name"];
          if (TextTools::isEmpty(prefix))
            throw Exception("'name' argument missing for user-defined substitution model.");
          string fname = args["file"];
          if (TextTools::isEmpty(fname))
            throw Exception("'file' argument missing for user-defined substitution model.");
          model.reset(new UserProteinSubstitutionModel(alpha, args["file"], dynamic_cast<ProteinFrequenciesSet*>(protFreq.release()), prefix + "+F."));
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
      // submodels of mixture models
      else if (modelName.substr(0,9) == "LGL08_CAT")
      {
        string subModelName = modelName.substr(10);

        size_t posp=modelDescription.find("(");

        string modelDesc2=modelName.substr(0,9)+modelDescription.substr(posp);
        
        BppOTransitionModelFormat nestedReader(PROTEIN, true, true, false, verbose_, warningLevel_);
        
        MixedTransitionModel* nestedModel=dynamic_cast<MixedTransitionModel*>(nestedReader.readTransitionModel(alphabet, modelDesc2, data, false));
      
        // Check that nestedModel is fine and has subModel of given name
        if (nestedModel==NULL)
          throw Exception("Unknown model " + modelName + ".");
        
        if (nestedModel->getModel(subModelName)==NULL)
          throw Exception("BppOSubstitutionModelFormat::read. " + nestedModel->getName() + "argument for model 'FromMixture' has no submodel with name " + subModelName + ".");

        // Now we create the FromMixture substitution model:
        model.reset(new FromMixtureSubstitutionModel(*nestedModel, subModelName, modelDesc2));
        
        delete nestedModel;
      }      
      else if (modelName == "Empirical")
      {
        string prefix = args["name"];
        if (TextTools::isEmpty(prefix))
          throw Exception("'name' argument missing for user-defined substitution model.");
        model.reset(new UserProteinSubstitutionModel(alpha, args["file"], prefix));
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
        unique_ptr<ProteinSubstitutionModel> nestedModel(dynamic_cast<ProteinSubstitutionModel*>(nestedReader.readSubstitionModel(alphabet, exchangeability, data, false)));
        string nbrOfParametersPerBranch = args["nbrAxes"];
        if (TextTools::isEmpty(nbrOfParametersPerBranch))
          throw Exception("'nbrAxes' argument missing to define the number of axis of the Correspondence Analysis.");
        //Now we create the Coala model:
        model.reset(new Coala(alpha, *nestedModel, TextTools::to<unsigned int>(nbrOfParametersPerBranch)));
        model->setFreqFromData(*data);
      }
      else
        throw Exception("Model '" + modelName + "' is unknown, or does not fit proteic alphabet.");      
    }
    else if (AlphabetTools::isBinaryAlphabet(alphabet))
    {
      if (!(alphabetCode_ & BINARY))
        throw Exception("BppOSubstitutionModelFormat::read. Binary alphabet not supported.");
      
      const BinaryAlphabet* balpha = dynamic_cast<const BinaryAlphabet*>(alphabet);

      if (modelName == "Binary")
        model.reset(new BinarySubstitutionModel(balpha));
      else
        throw Exception("Model '" + modelName + "' unknown, or does not fit binary alphabet.");
    }
    else
      throw Exception("Model '" + modelName + "' unknown, or does not fit " + alphabet->getAlphabetType() + " alphabet.");
    
  }
    
  if (verbose_)
    ApplicationTools::displayResult("Substitution model", modelName);
  
  updateParameters_(model.get(), args);
  
  if (parseArguments)
    initialize_(*model, data);
  
  return model.release();  
}

  
void BppOSubstitutionModelFormat::updateParameters_(TransitionModel* model, std::map<std::string, std::string>& args)
{
  
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

    if (((pval.rfind("_") != string::npos) && (TextTools::isDecimalInteger(pval.substr(pval.rfind("_")+1,string::npos)))) ||
        (pval.find("(") != string::npos))
      continue;
    bool found = false;
    for (size_t j = 0; j < pl.size() && !found; j++)
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
    // if (!TextTools::isDecimalNumber(pval) && !found)
    //   throw Exception("Incorrect parameter syntax: parameter " + pval + " was not found and can't be used as a value for parameter " + pname + ".");
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

  if (args.find("rate") != args.end())
  {
    model->addRateParameter();
    unparsedArguments_[pref+"rate"] = args["rate"];
  }  
}


SubstitutionModel* BppOSubstitutionModelFormat::readWord_(const Alphabet* alphabet, const std::string& modelDescription, const SiteContainer* data)
{
  unique_ptr<SubstitutionModel> model;
  string modelName = "";
  map<string, string> args;
  KeyvalTools::parseProcedure(modelDescription, modelName, args);

  vector<string> v_nestedModelDescription;
  vector<SubstitutionModel*> v_pSM;
  const CoreWordAlphabet* pWA;

  string s, nestedModelDescription;
  unsigned int nbmodels;

  if (((modelName == "Word" || modelName == "Kron") && !AlphabetTools::isWordAlphabet(alphabet)) ||
      ((!(modelName == "Word" || modelName == "Kron")) && !AlphabetTools::isCodonAlphabet(alphabet)))
    throw Exception("Bad alphabet type "
                    + alphabet->getAlphabetType() + " for  model " + modelName + ".");

  pWA = dynamic_cast<const CoreWordAlphabet*>(alphabet);

  
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
    BppOSubstitutionModelFormat nestedReader(alphabetCode_, false, true, false, false, warningLevel_);
    model.reset(nestedReader.readSubstitionModel(pWA->getNAlphabet(0), v_nestedModelDescription[0], data, false));
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
      BppOSubstitutionModelFormat nestedReader(alphabetCode_, false, true, false, false, warningLevel_);
      model.reset(nestedReader.readSubstitionModel(pWA->getNAlphabet(i), v_nestedModelDescription[i], data, false));
      map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_[modelName + "." + TextTools::toString(i + 1) + "_" + it->first] = it->second;
      }

      v_pSM.push_back(model.release());
    }
  }

  // //////////////////////////////
  // In case of a Kron... model, mgmt of the positions

  // check the vector of simultaneous changing positions
    
  vector<set<size_t> > vposKron;
  if (modelName.substr(0,4) == "Kron")
  {
    if (args.find("positions")!=args.end())
    {
      StringTokenizer nst(args["positions"], "+");
    
      while (nst.hasMoreToken())
      {
        string spos=nst.nextToken();
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
    vector<SubstitutionModel*> vsm;
    
    for (auto mod : v_pSM)
      if (dynamic_cast<SubstitutionModel*>(mod)==0)
        throw Exception("Need Markov SubstitutionModel for Words");
      else
        vsm.push_back(dynamic_cast<SubstitutionModel*>(mod));
    
    if (v_nestedModelDescription.size() != nbmodels) {
      model.reset(new WordSubstitutionModel(vsm[0], nbmodels));
    } else {
      ModelList ml(vsm);
      model.reset(new WordSubstitutionModel(ml));
    }
  }

  // /////////////////////////////////
  // / KRON
  // ///////////////////////////////

  else if (modelName == "Kron")
  {
    vector<SubstitutionModel*> vsm;
    
    for (auto mod : v_pSM)
      if (dynamic_cast<SubstitutionModel*>(mod)==0)
        throw Exception("Need Markov SubstitutionModel for Krons");
      else
        vsm.push_back(dynamic_cast<SubstitutionModel*>(mod));

    if (vposKron.size()==0)
    {
      if (v_nestedModelDescription.size() != nbmodels) {
        model.reset(new KroneckerWordSubstitutionModel(vsm[0], nbmodels));
      } else {
        ModelList ml(vsm);
        model.reset(new KroneckerWordSubstitutionModel(ml));
      }
    }
    else
    {
      if (v_nestedModelDescription.size() != nbmodels) {
        model.reset(new KroneckerWordSubstitutionModel(vsm[0], nbmodels, vposKron));
      } else {
        ModelList ml(vsm);
        model.reset(new KroneckerWordSubstitutionModel(ml, vposKron));
      }
    }
  }

  // /////////////////////////////////
  // / CODON
  // ///////////////////////////////

  else
  {
    const CodonAlphabet* pCA = dynamic_cast<const CodonAlphabet*>(pWA);
    if (pCA == 0)
      throw Exception("Non codon Alphabet for model" + modelName + ".");

    unique_ptr< AlphabetIndex2 > pai2;
    unique_ptr<FrequenciesSet> pFS;

    if ((dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]) == 0) ||
        ((v_nestedModelDescription.size() == 3) &&
         (dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]) == 0 || dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]) == 0)))
      throw Exception("Non simple NucleotideSubstitutionModel imbedded in " + modelName + " model.");

    if (args.find("genetic_code") != args.end()) {
      ApplicationTools::displayWarning("'genetic_code' argument is no longer supported inside model description, and has been supersided by a global 'genetic_code' option.");
      throw Exception("BppOSubstitutionModelFormat::read. Deprecated 'genetic_code' argument.");
    }

    if (!geneticCode_)
      throw Exception("BppOSubstitutionModelFormat::readWord_(). No genetic code specified! Consider using 'setGeneticCode'.");


    ///////////////////////////////////
    /// Dist
    
    if ((modelName.find("Dist") != string::npos) || (modelName=="SENCA"))
      pai2.reset((args.find("aadistance") == args.end()) ? 0 : SequenceApplicationTools::getAlphabetIndex2(&AlphabetTools::PROTEIN_ALPHABET, args["aadistance"]));
    
    //////////////////////////////////
    /// Freq
    
    if (modelName.find("Freq") != string::npos)
    {
      if (args.find("frequencies") == args.end())
        throw Exception("Missing equilibrium frequencies.");

      BppOFrequenciesSetFormat bIOFreq(alphabetCode_, verbose_, warningLevel_);
      bIOFreq.setGeneticCode(geneticCode_); //This uses the same instance as the one that will be used by the model
      pFS.reset(bIOFreq.read(pCA, args["frequencies"], data, false));
      map<string, string> unparsedParameterValuesNested(bIOFreq.getUnparsedArguments());

      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_[modelName + "." + it->first] = it->second;
      }
    }

    // //////////////
    // // Triplet

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

    else if (modelName.find("Codon")!=string::npos)
    {
      vector<CoreCodonSubstitutionModel*> vCSM;
      string name="Codon";
      map<string, string> unparsedParameterValuesNested;

      if (modelName.find("Dist")!=string::npos)
      {
        name +="Dist";
        
        vCSM.push_back(new AbstractCodonDistanceSubstitutionModel(pai2.release(), geneticCode_, ""));
      }
      else if (modelName.find("BGC")!=string::npos)
      {
        name +="BGC";
        
        vCSM.push_back(new AbstractCodonBGCSubstitutionModel(geneticCode_, ""));
      }
      else if (modelName.find("Prot")!=string::npos)
      {
        name+="Prot";

        if (args.find("protmodel")==args.end())
          throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'protmodel' for codon model argument 'Prot'.");

        nestedModelDescription = args["protmodel"];
        BppOSubstitutionModelFormat nestedReader(PROTEIN, false, false, allowGaps_, verbose_, warningLevel_);
        
        shared_ptr<ProteinSubstitutionModel> nestedModel(dynamic_cast<ProteinSubstitutionModel*>(nestedReader.readSubstitionModel(geneticCode_->getTargetAlphabet(), nestedModelDescription, data, false)));
        
        unparsedParameterValuesNested.insert(nestedReader.getUnparsedArguments().begin(),nestedReader.getUnparsedArguments().end());

        vCSM.push_back(new AbstractCodonAARateSubstitutionModel(nestedModel, geneticCode_, ""));
      }
      if (modelName.find("AAClust")!=string::npos)
      {
        name+="AAClust";

        // Initialization using the "assign" argument
        vector<uint> partition;
        if (args.find("partition")!=args.end())
        {
          string rf = args["partition"];

          StringTokenizer strtok(rf.substr(1, rf.length() - 2), ",");
          while (strtok.hasMoreToken())
            partition.push_back(TextTools::to<uint>(strtok.nextToken()));
        }

        AbstractCodonClusterAASubstitutionModel* aca=partition.size()!=0?
          new AbstractCodonClusterAASubstitutionModel(geneticCode_, "", partition):
          new AbstractCodonClusterAASubstitutionModel(geneticCode_, "");

        vCSM.push_back(aca);
      }

      /// Default name in none used before
      if (vCSM.size()==0)
        name+="Rate";
      
      if (modelName.find("CpG")!=string::npos)
      {
        name+="CpG";
        vCSM.push_back(new AbstractCodonCpGSubstitutionModel(*pCA,""));
      }
      
      if (modelName.find("AAFit")!=string::npos)
      {
        name+="AAFit";

        if (args.find("fitness")==args.end())
          throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'fitness' for codon model argument 'AAFit'.");

        string nestedFreqDescription = args["fitness"];
        BppOFrequenciesSetFormat nestedReader(PROTEIN, verbose_, warningLevel_);

        FrequenciesSet* nestedFreq=nestedReader.read(geneticCode_->getTargetAlphabet(), nestedFreqDescription, data, false);
        
        for (auto  it : nestedReader.getUnparsedArguments())
          unparsedParameterValuesNested["fit_" + it.first] = it.second;

        AbstractCodonAAFitnessSubstitutionModel* aca=new AbstractCodonAAFitnessSubstitutionModel(nestedFreq, geneticCode_, "");

        if (args.find("Ns")!=args.end())
          aca->addNsParameter();
        
        vCSM.push_back(aca);
      }
      else if (modelName.find("Fit")!=string::npos)
      {
        if (args.find("fitness") == args.end())
          throw Exception("BppOSubstitutionModelFormat::read. Missing argument 'fitness' for codon model argument 'Fit'.");
        string nestedFreqDescription = args["fitness"];
        
        BppOFrequenciesSetFormat nestedReader(alphabetCode_, verbose_, warningLevel_);
        nestedReader.setGeneticCode(geneticCode_);
      
        FrequenciesSet* nestedFreq=nestedReader.read(alphabet, nestedFreqDescription, data, false);
        
        for (auto  it : nestedReader.getUnparsedArguments())
          unparsedParameterValuesNested["fit_" + it.first] = it.second;
        
        vCSM.push_back(new AbstractCodonFitnessSubstitutionModel(nestedFreq, geneticCode_, ""));
      }

      if (modelName.find("PhasFreq")!=string::npos)
      {
        name+="PhasFreq";
        vCSM.push_back(new AbstractCodonPhaseFrequenciesSubstitutionModel(pFS.release(),""));
      }
      else if (modelName.find("Freq")!=string::npos)
      {
        name+="Freq";
        vCSM.push_back(new AbstractCodonFrequenciesSubstitutionModel(pFS.release(),""));
      }

      // Then we update the parameter set:
      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_[name+"."+it->first] = it->second;
      }
    
      model.reset((v_nestedModelDescription.size() != 3)
                  ? new CodonAdHocSubstitutionModel(
                    geneticCode_,
                    dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                    vCSM,
                    name)
                  : new CodonAdHocSubstitutionModel(
                    geneticCode_,
                    dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                    dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                    dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]),
                    vCSM,
                    name));
      
    }
    else if (modelName == "KronDistFreq")
    {
      if (v_nestedModelDescription.size() != 3){
        
        if (vposKron.size()==0)
          model.reset(new KroneckerCodonDistanceFrequenciesSubstitutionModel(geneticCode_,
                                                                             dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                                             pFS.release(),
                                                                             pai2.release()));
        else
          model.reset(new KroneckerCodonDistanceFrequenciesSubstitutionModel(geneticCode_,
                                                                             dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                                             vposKron,
                                                                             pFS.release(),
                                                                             pai2.release()));

      }
      else
      {
        if (vposKron.size()!=0)
          model.reset(new KroneckerCodonDistanceFrequenciesSubstitutionModel(
                        geneticCode_,
                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]),
                        pFS.release(),
                        pai2.release()));
        else
          model.reset(new KroneckerCodonDistanceFrequenciesSubstitutionModel(
                        geneticCode_,
                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]),
                        vposKron, 
                        pFS.release(),
                        pai2.release()));

      }
    }
    else if (modelName == "KronDist")
    {
      if (v_nestedModelDescription.size() != 3){
        
        if (vposKron.size()==0)
          model.reset(new KroneckerCodonDistanceSubstitutionModel(geneticCode_,
                                                                  dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                                  pai2.release()));
        else
          model.reset(new KroneckerCodonDistanceSubstitutionModel(geneticCode_,
                                                                  dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                                                                  vposKron,
                                                                  pai2.release()));

      }
      else
      {
        if (vposKron.size()!=0)
          model.reset(new KroneckerCodonDistanceSubstitutionModel(
                        geneticCode_,
                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]),
                        pai2.release()));
        else
          model.reset(new KroneckerCodonDistanceSubstitutionModel(
                        geneticCode_,
                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                        dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]),
                        vposKron, 
                        pai2.release()));

      }
    }
    
    else if (modelName == "SENCA")
    {
      if (args.find("fitness") == args.end())
        throw Exception("Missing fitness in model " + modelName + ".");

      BppOFrequenciesSetFormat bIOFreq(alphabetCode_, verbose_, warningLevel_);
      bIOFreq.setGeneticCode(geneticCode_);
      
      unique_ptr<FrequenciesSet> pFit(bIOFreq.read(pCA, args["fitness"], data, false));
      map<string, string> unparsedParameterValuesNested(bIOFreq.getUnparsedArguments());

      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_[modelName + ".fit_" + it->first] = it->second;
      }

      if (v_nestedModelDescription.size() != 3)
      {
        model.reset(new SENCA(geneticCode_,
                              dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                              pFit.release(),
                              pai2.release()));
      }
      else
        model.reset(new SENCA(
                      geneticCode_,
                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[0]),
                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[1]),
                      dynamic_cast<NucleotideSubstitutionModel*>(v_pSM[2]),
                      pFit.release(),
                      pai2.release()));
    }
  }

  return model.release();
}



void BppOSubstitutionModelFormat::write(const TransitionModel& model,
                                        OutputStream& out,
                                        std::map<std::string, std::string>& globalAliases,
                                        std::vector<std::string>& writtenNames) const
{
  bool comma = false;

  //  Mixed Model that are defined as "Mixture" and "Mixed"

  if ((dynamic_cast<const MixedTransitionModel*>(&model) != NULL) && (dynamic_cast<const AbstractBiblioMixedTransitionModel*>(&model) == NULL))
  {
    writeMixed_(*dynamic_cast<const MixedTransitionModel*>(&model), out, globalAliases, writtenNames);
    return;
  }

  out << model.getName() + "(";

  // Is it a protein user defined model?
  const UserProteinSubstitutionModel* userModel = dynamic_cast<const UserProteinSubstitutionModel*>(&model);
  if (userModel)
  {
    out << "file=" << userModel->getPath();
    string ns=userModel->getNamespace();
    
    if (TextTools::hasSubstring(ns, "+F.") )
      ns = ns.substr(0,ns.size()-3); 
    else  
      ns = ns.substr(0,ns.size()-1); 

    out << ",name=" << ns;
    comma = true;    
  }

  // Is it a markov-modulated model?
  const MarkovModulatedSubstitutionModel* mmModel = dynamic_cast<const MarkovModulatedSubstitutionModel*>(&model);
  if (mmModel)
  {
    out << "model=";
    const TransitionModel* nestedModel = mmModel->getNestedModel();
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
    const TransitionModel* nestedModel = reModel->getNestedModel();
    write(*nestedModel, out, globalAliases, writtenNames);
    comma = true;
  }

  // Is it a YpR model?
  const YpR* yprModel = dynamic_cast<const YpR*>(&model);
  if (yprModel)
  {
    out << "model=";
    const TransitionModel* nestedModel = yprModel->getNestedModel();
    write(*nestedModel, out, globalAliases, writtenNames);
    comma = true;
  }

  // Is it a OneChange model?
  const OneChangeTransitionModel* onechangetransitionmodel = dynamic_cast<const OneChangeTransitionModel*>(&model);
  if (onechangetransitionmodel)
  {
    out << "model=";
    const TransitionModel& nestedModel = onechangetransitionmodel->getModel();
    write(nestedModel, out, globalAliases, writtenNames);
    comma = true;
  }
  else
  {
    // Is it a model with register?
    const OneChangeRegisterTransitionModel* onechangeregistertransitionmodel = dynamic_cast<const OneChangeRegisterTransitionModel*>(&model);
    const RegisterRatesSubstitutionModel* regratessubsmodel = dynamic_cast<const RegisterRatesSubstitutionModel*>(&model);
    if (onechangeregistertransitionmodel || regratessubsmodel)
    {
      out << "model=";
      const TransitionModel& nestedModel = onechangeregistertransitionmodel?onechangeregistertransitionmodel->getModel():regratessubsmodel->getModel();
      
      write(nestedModel, out, globalAliases, writtenNames);
      comma = true;
      out << ", register=";
      if (onechangeregistertransitionmodel)
        out << onechangeregistertransitionmodel->getRegisterName();
      else
        out << regratessubsmodel->getRegisterName();
      
      if (onechangeregistertransitionmodel)
        out << ", numReg=" << VectorTools::paste(onechangeregistertransitionmodel->getRegisterNumbers(),"+");
    }
  }
  
  // Is it a gBGC model?
  const gBGC* gbgcModel = dynamic_cast<const gBGC*>(&model);
  if (gbgcModel)
  {
    StdStr sout;
    
    const TransitionModel* nestedModel = gbgcModel->getNestedModel();
    write(*nestedModel, sout, globalAliases, writtenNames);

    string ss=sout.str();

    string::size_type begin = ss.find_first_of("(");
    string::size_type end = ss.find_last_of(")");

    out << ss.substr(begin+1,end-begin-1);
    
    comma = true;
  }

  // Is it a word model?

  const AbstractWordSubstitutionModel* wM = dynamic_cast<const AbstractWordSubstitutionModel*>(&model);
  if (wM)
  {
    size_t nmod = wM->getNumberOfModels();
    const TransitionModel* mod0 = wM->getNModel(0);
    if (nmod == 1)
    {
      out << "model=";
      write(*mod0, out, globalAliases, writtenNames);
    }
    else
    {
      const TransitionModel* mod1 = wM->getNModel(1);
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

  // Is it a InMixed model?

  const InMixedSubstitutionModel* inModel = dynamic_cast<const InMixedSubstitutionModel*>(&model);
  if (inModel)
  {
    out << "model=";
    write(inModel->getMixedModel(), out, globalAliases, writtenNames);
    out << ", numMod=" << TextTools::toString(inModel->getSubModelNumber());
    comma = true;
  }
  
  // Is it a model with FrequenciesSet?
  
  const FrequenciesSet* pfs = model.getFrequenciesSet();
  auto paf=dynamic_cast<const AbstractWrappedModel*>(&model);
  
  if ((pfs!=0) && (!paf))
  {
    if (comma)
      out << ",";
    out << "frequencies=";

    BppOFrequenciesSetFormat bIOFreq(alphabetCode_, false, warningLevel_);
    bIOFreq.write(pfs, out, globalAliases, writtenNames);
    
    comma = true;
  }

  // Is it a codon model with Protein Model or partition in it? 
  const CodonAdHocSubstitutionModel* casm=dynamic_cast<const CodonAdHocSubstitutionModel*>(&model);
  if (casm)
  {
    for (size_t i=0; i<casm->getNumberOfModels(); i++)
    {
      const AbstractCodonAARateSubstitutionModel* acr=dynamic_cast<const AbstractCodonAARateSubstitutionModel*>(casm->getNModel(i).get());
      if (acr)
      {
        if (comma)
          out << ",";
        out << "protmodel=";
    
        write(*acr->getAAModel().get(), out, globalAliases, writtenNames);
        comma = true;
      }
      const AbstractCodonAAFitnessSubstitutionModel* acf=dynamic_cast<const AbstractCodonAAFitnessSubstitutionModel*>(casm->getNModel(i).get());
      if (acf)
      {
        if (comma)
          out << ",";
        out << "fitness=";
    
        BppOFrequenciesSetFormat bIOFreq(PROTEIN, false, warningLevel_);
        bIOFreq.write(&acf->getAAFitness(), out, globalAliases, writtenNames);
        comma = true;
      }
      const AbstractCodonClusterAASubstitutionModel* acc=dynamic_cast<const AbstractCodonClusterAASubstitutionModel*>(casm->getNModel(i).get());
      if (acc)
      {
        if (comma)
          out << ",";
        out << "partition=(";
        const vector<uint>& vass=acc->getAssign();
        
        for (size_t j=0; j<vass.size(); j++)
        {
          if (j != 0)
            out << ",";
          out << vass[j];
        }
        out << ")";
        comma = true;
      }
    }
  }

  // Specific case of SENCA

  const SENCA* pCF = dynamic_cast<const SENCA*>(&model);
  if (pCF)
  {
    if (comma)
      out << ",";
    out << "fitness=";

    BppOFrequenciesSetFormat bIOFreq(alphabetCode_, false, warningLevel_);
    bIOFreq.write(pCF->getFitness(), out, globalAliases, writtenNames);
    comma = true;
  }

  // and bibliomixed models

  const YNGP_M7* pM7 = dynamic_cast<const YNGP_M7*>(&model);
  if (pM7)
  {
    if (comma)
      out << ",";
    out << "n=" << pM7->getNumberOfModels();
    comma=true;
  }

  const YNGP_M8* pM8 = dynamic_cast<const YNGP_M8*>(&model);
  if (pM8)
  {
    if (comma)
      out << ",";
    out << "n=" << pM8->getNumberOfModels()-1;
    comma=true;
  }

  const YNGP_M9* pM9 = dynamic_cast<const YNGP_M9*>(&model);
  if (pM9)
  {
    if (comma)
      out << ",";
    out << "nbeta=" << pM9->getNBeta() << ",ngamma=" << pM9->getNGamma();
    
    comma=true;
  }

  const YNGP_M10* pM10 = dynamic_cast<const YNGP_M10*>(&model);
  if (pM10)
  {
    if (comma)
      out << ",";
    out << "nbeta=" << pM10->getNBeta() << ",ngamma=" << pM10->getNGamma();
    
    comma=true;
  }
  
  BppOParametrizableFormat bIO;

  // case of Biblio models, update writtenNames

  const AbstractBiblioSubstitutionModel* absm=dynamic_cast<const AbstractBiblioSubstitutionModel*>(&model);
  
  if (absm)
  {
    size_t wNs=writtenNames.size();
    
    for (size_t i=0;i<wNs; i++)
    {
      try
      {
        writtenNames.push_back(absm->getNamespace()+absm->getParNameFromPmodel(writtenNames[i]));
      }
      catch (Exception& e) {}
    }
  }
    
  bIO.write(&model, out, globalAliases, model.getIndependentParameters().getParameterNames(), writtenNames, true, comma);
  out << ")";
}


void BppOSubstitutionModelFormat::writeMixed_(const MixedTransitionModel& model,
                                              OutputStream& out,
                                              std::map<std::string, std::string>& globalAliases,
                                              std::vector<std::string>& writtenNames) const
{
  if (dynamic_cast<const MixtureOfTransitionModels*>(&model) != NULL)
  {
    const MixtureOfTransitionModels* pMS = dynamic_cast<const MixtureOfTransitionModels*>(&model);

    // for (unsigned int i = 0; i < pMS->getNumberOfModels(); i++)
    // {
    //   const SubstitutionModel* eM = pMS->getModel(i);

    //   vector<string> vpl = eM->getIndependentParameters().getParameterNames();
    //   for (unsigned j = 0; j < vpl.size(); j++)
    //   {
    //     if (eM->getParameterNameWithoutNamespace(vpl[j]) == "rate")
    //       writtenNames.push_back(vpl[j]);
    //   }
    // }

    out << "Mixture(";
    for (unsigned int i = 0; i < pMS->getNumberOfModels(); i++)
    {
      if (i != 0)
        out << ", ";
      out << "model" + TextTools::toString(i + 1) + "=";
      write(*pMS->AbstractMixedTransitionModel::getNModel(i), out, globalAliases, writtenNames);
    }
  }
  else
  {
    const MixtureOfATransitionModel* pMS = dynamic_cast<const MixtureOfATransitionModel*>(&model);
    out << "MixedModel(model=";
    const TransitionModel* eM = pMS->getModel(0);

    ParameterList pl = eM->getIndependentParameters();
    vector<string> vpl = pl.getParameterNames();

    for (unsigned j = 0; j < vpl.size(); j++)
    {
      if (find(writtenNames.begin(), writtenNames.end(), vpl[j]) == writtenNames.end())
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
  TransitionModel& model,
  const SiteContainer* data)
{
  string initFreqs = ApplicationTools::getStringParameter(model.getNamespace() + "initFreqs", unparsedArguments_, "", "", true, warningLevel_);
  if (verbose_)
    ApplicationTools::displayResult("External frequencies initialization for model", (initFreqs == "") ? "None" : initFreqs);

  if (initFreqs != "")
  {
    if (initFreqs == "observed")
    {
      if (!data)
        throw Exception("BppOSubstitutionModelFormat::initialize_(). Missing data for observed frequencies");
      unsigned int psi = ApplicationTools::getParameter<unsigned int>(model.getNamespace() + "initFreqs.observedPseudoCount", unparsedArguments_, 0, "", true, warningLevel_);
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

    unparsedArguments_.erase(unparsedArguments_.find(model.getNamespace() + "initFreqs"));
  }

  ParameterList pl = model.getIndependentParameters();
  for (size_t i = 0; i < pl.size(); i++)
  {
    AutoParameter ap(pl[i]);
    ap.setMessageHandler(ApplicationTools::warning.get());
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
    try {
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
  
    catch (Exception& e) {}
  }
  
  model.matchParametersValues(pl);
}


