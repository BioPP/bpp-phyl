//
// File: BppOFrequenciesSetFormatFormat.cpp
// Created by: Laurent Guéguen
// Created on: lundi 9 juillet 2012, à 12h 56
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

#include "BppOFrequenciesSetFormat.h"


#include "../Model/FrequenciesSet/NucleotideFrequenciesSet.h"
#include "../Model/FrequenciesSet/ProteinFrequenciesSet.h"
#include "../Model/FrequenciesSet/CodonFrequenciesSet.h"
#include "../Model/FrequenciesSet/MvaFrequenciesSet.h"

//From bpp-core:
#include <Bpp/Io/FileTools.h>
#include <Bpp/Io/BppOParametrizableFormat.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/AutoParameter.h>

//From bpp-seq:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

// From the STL:
#include <iomanip>

#include "BppOSubstitutionModelFormat.h"

using namespace std;

unsigned char BppOFrequenciesSetFormat::DNA = 1;
unsigned char BppOFrequenciesSetFormat::RNA = 2;
unsigned char BppOFrequenciesSetFormat::NUCLEOTIDE = 1 | 2;
unsigned char BppOFrequenciesSetFormat::PROTEIN = 4;
unsigned char BppOFrequenciesSetFormat::CODON = 8;
unsigned char BppOFrequenciesSetFormat::WORD = 16;
unsigned char BppOFrequenciesSetFormat::ALL = 1 | 2 | 4 | 8 | 16;

FrequenciesSet* BppOFrequenciesSetFormat::read(const Alphabet* alphabet, const std::string& freqDescription, const SiteContainer* data, bool parseArguments)
{
  unparsedArguments_.clear();
  string freqName;
  map<string, string> args;
  KeyvalTools::parseProcedure(freqDescription, freqName, args);
  auto_ptr<FrequenciesSet> pFS;

  if (freqName == "Fixed")
  {
    if (AlphabetTools::isNucleicAlphabet(alphabet))
    {
      if (alphabetCode_ & NUCLEOTIDE)
        pFS.reset(new FixedNucleotideFrequenciesSet(dynamic_cast<const NucleicAlphabet*>(alphabet)));
      else
        throw Exception("Nucleotide alphabet not supported.");
    }
    else if (AlphabetTools::isProteicAlphabet(alphabet))
    {
      if (alphabetCode_ & PROTEIN)
        pFS.reset(new FixedProteinFrequenciesSet(dynamic_cast<const ProteicAlphabet*>(alphabet)));
      else
        throw Exception("Protein alphabet not supported.");
    }
    else if (AlphabetTools::isCodonAlphabet(alphabet))
    {
      if (alphabetCode_ & CODON) {
        if (!geneticCode_)
          throw Exception("BppOFrequenciesSetFormat::read(). No genetic code specified! Consider using 'setGeneticCode'.");
        pFS.reset(new FixedCodonFrequenciesSet(geneticCode_));
      } else {
        throw Exception("Codon alphabet not supported.");
      }
    }
    else
    {
      pFS.reset(new FixedFrequenciesSet(new CanonicalStateMap(alphabet, false)));
    }
  }
  else if (freqName == "Full")
  {
    unsigned short method = 1;
    if (args.find("method") != args.end()){
      if (args["method"] == "local")
        method=2;
      else
        if (args["method"] == "binary")
          method=3;
    }
      
    if (AlphabetTools::isNucleicAlphabet(alphabet))
    {
      if (alphabetCode_ & NUCLEOTIDE)
        pFS.reset(new FullNucleotideFrequenciesSet(dynamic_cast<const NucleicAlphabet*>(alphabet)));
      else
        throw Exception("Nucleotide alphabet not supported.");
    }
    else if (AlphabetTools::isProteicAlphabet(alphabet))
    {
      if (alphabetCode_ & PROTEIN)
        pFS.reset(new FullProteinFrequenciesSet(dynamic_cast<const ProteicAlphabet*>(alphabet), false, method));
      else
        throw Exception("Protein alphabet not supported.");
    }
    else if (AlphabetTools::isCodonAlphabet(alphabet))
    {
      if (alphabetCode_ & CODON) {
        if (!geneticCode_)
          throw Exception("BppOFrequenciesSetFormat::read(). No genetic code specified! Consider using 'setGeneticCode'.");
        pFS.reset(new FullCodonFrequenciesSet(geneticCode_));
      } else {
        throw Exception("Codon alphabet not supported.");
      }
    }
    else
    {
      //NB: jdutheil 25/09/14 => gap models will not be supported before we add the appropriate option!
      pFS.reset(new FullFrequenciesSet(new CanonicalStateMap(alphabet, false), false, method));
    }
  }
  else if (freqName == "GC")
  {
    if (!AlphabetTools::isNucleicAlphabet(alphabet))
      throw Exception("Error, unvalid frequencies " + freqName + " with non-nucleic alphabet.");
    if (alphabetCode_ & NUCLEOTIDE)
      pFS.reset(new GCFrequenciesSet(dynamic_cast<const NucleicAlphabet*>(alphabet)));
    else
      throw Exception("Nucleotide alphabet not supported.");
  }

  // INDEPENDENTWORD
  else if (freqName == "Word")
  {
    if (!(alphabetCode_& WORD))
      throw Exception("Word alphabet not supported.");
    if (!AlphabetTools::isWordAlphabet(alphabet))
      throw Exception("BppOFrequenciesSetFormat::read(...).\n\t Bad alphabet type "
                      + alphabet->getAlphabetType() + " for frequencies set " + freqName + ".");

    const WordAlphabet* pWA = dynamic_cast<const WordAlphabet*>(alphabet);

    if (args.find("frequency") != args.end())
    {
      string sAFS = args["frequency"];

      unsigned int nbfreq = pWA->getLength();
      string st = "";
      for (unsigned i = 0; i < nbfreq; i++)
      {
        st += TextTools::toString(i + 1);
      }

      BppOFrequenciesSetFormat nestedReader(alphabetCode_, false, warningLevel_);
      auto_ptr<FrequenciesSet> pFS2(nestedReader.read(pWA->getNAlphabet(0), sAFS, data, false));
      map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_[st + "_" + it->first] = it->second;
      }
      pFS.reset(new WordFromUniqueFrequenciesSet(pWA, pFS2.release()));
    }
    else
    {
      if (args.find("frequency1") == args.end())
        throw Exception("BppOFrequenciesSetFormat::read. Missing argument 'frequency' or 'frequency1' for frequencies set 'Word'.");
      vector<string> v_sAFS;
      vector<FrequenciesSet*> v_AFS;
      unsigned int nbfreq = 1;

      while (args.find("frequency" + TextTools::toString(nbfreq)) != args.end())
      {
        v_sAFS.push_back(args["frequency" + TextTools::toString(nbfreq++)]);
      }

      if (v_sAFS.size() != pWA->getLength())
        throw Exception("BppOFrequenciesSetFormat::read. Number of frequencies (" + TextTools::toString(v_sAFS.size()) + ") does not match length of the words (" + TextTools::toString(pWA->getLength()) + ")");

      for (size_t i = 0; i < v_sAFS.size(); ++i)
      {
        BppOFrequenciesSetFormat nestedReader(alphabetCode_, false, warningLevel_);
        pFS.reset(nestedReader.read(pWA->getNAlphabet(i), v_sAFS[i], data, false));
        map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
        for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
        {
          unparsedArguments_[TextTools::toString(i + 1) + "_" + it->first] = it->second;
        }
        v_AFS.push_back(pFS.release());
      }

      pFS.reset(new WordFromIndependentFrequenciesSet(pWA, v_AFS));
    }
  }
  //From Model
  else if (freqName == "FromModel")
  {
    if (args.find("model") == args.end())
      throw Exception("Missing argument 'model' for frequencies " + freqName + ".");
    
    BppOSubstitutionModelFormat nestedReader(alphabetCode_, false, true, false, false, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_);

    SubstitutionModel* model=nestedReader.read(alphabet, args["model"], data, false);
    pFS.reset(new FromModelFrequenciesSet(model));
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedArguments_["FromModel." + it->first] = it->second;
    }
  }
  
  
  // INDEPENDENT CODON
  else if (freqName == "Codon")
  {
    if (!(alphabetCode_ & CODON))
      throw Exception("Codon alphabet not supported.");
    if (!AlphabetTools::isCodonAlphabet(alphabet))
      throw Exception("BppOFrequenciesSetFormat::read.\n\t Bad alphabet type "
                      + alphabet->getAlphabetType() + " for frequenciesset " + freqName + ".");
    if (!geneticCode_)
      throw Exception("BppOFrequenciesSetFormat::read(). No genetic code specified! Consider using 'setGeneticCode'.");

    const CodonAlphabet* pWA = dynamic_cast<const CodonAlphabet*>(alphabet);

    if (args.find("frequency") != args.end())
    {
      string sAFS = args["frequency"];

      BppOFrequenciesSetFormat nestedReader(alphabetCode_, false, warningLevel_);
      auto_ptr<FrequenciesSet> pFS2(nestedReader.read(pWA->getNAlphabet(0), sAFS, data, false));
      map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_["123_" + it->first] = it->second;
      }
      
      pFS.reset(new CodonFromUniqueFrequenciesSet(geneticCode_, pFS2.release(), "Codon"));
    }
    else
    {
      if (args.find("frequency1") == args.end())
        throw Exception("BppOFrequenciesSetFormat::read. Missing argument 'frequency' or 'frequency1' for frequencies set.");
      vector<string> v_sAFS;
      vector<FrequenciesSet*> v_AFS;
      unsigned int nbfreq = 1;

      while (args.find("frequency" + TextTools::toString(nbfreq)) != args.end())
      {
        v_sAFS.push_back(args["frequency" + TextTools::toString(nbfreq++)]);
      }

      if (v_sAFS.size() != 3)
        throw Exception("BppOFrequenciesSetFormat::read. Number of frequencies (" + TextTools::toString(v_sAFS.size()) + ") is not three");

      for (size_t i = 0; i < v_sAFS.size(); ++i)
      {
        BppOFrequenciesSetFormat nestedReader(alphabetCode_, false, warningLevel_);
        pFS.reset(nestedReader.read(pWA->getNAlphabet(i), v_sAFS[i], data, false));
        map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
        for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
        {
          unparsedArguments_[TextTools::toString(i + 1) + "_" + it->first] = it->second;
        }
        v_AFS.push_back(pFS.release());
      }

      if (!geneticCode_)
        throw Exception("BppOFrequenciesSetFormat::read(). No genetic code specified! Consider using 'setGeneticCode'.");
      pFS.reset(new CodonFromIndependentFrequenciesSet(geneticCode_, v_AFS, "Codon"));
    }
  }

  // CODON PER AA Frequencies
  else if (freqName == "FullPerAA")
  {
    if (!(alphabetCode_ & CODON))
      throw Exception("Codon alphabet not supported.");
    if (!AlphabetTools::isCodonAlphabet(alphabet))
      throw Exception("BppOFrequenciesSetFormat::read.\n\t Bad alphabet type "
                      + alphabet->getAlphabetType() + " for frequenciesset " + freqName + ".");

    if (!geneticCode_)
      throw Exception("BppOFrequenciesSetFormat::read(). No genetic code specified! Consider using 'setGeneticCode'.");
    
    const ProteicAlphabet* pPA = dynamic_cast<const ProteicAlphabet*>(geneticCode_->getTargetAlphabet());

    unsigned short method=1;
    if (args.find("method") != args.end()){
      if (args["method"]=="local")
        method=2;
      else
        if (args["method"]=="binary")
          method=3;
    }

    if (args.find("protein_frequencies") != args.end())
    {
      string sPFS = args["protein_frequencies"];
      BppOFrequenciesSetFormat nestedReader(alphabetCode_, false, warningLevel_);
      auto_ptr<ProteinFrequenciesSet> pPFS(dynamic_cast<ProteinFrequenciesSet*>(nestedReader.read(pPA, sPFS, data, false)));
      map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_["FullPerAA." + it->first] = it->second;
      }
      pFS.reset(new FullPerAACodonFrequenciesSet(geneticCode_, pPFS.release(), method));
    }
    else
      pFS.reset(new FullPerAACodonFrequenciesSet(geneticCode_, method));
  }

  // codeml frequencies syntax
  
  else if (AlphabetTools::isCodonAlphabet(alphabet))
  {
    if (!(alphabetCode_ & CODON))
      throw Exception("Codon alphabet not supported.");
    if (!geneticCode_)
      throw Exception("BppOFrequenciesSetFormat::read(). No genetic code specified! Consider using 'setGeneticCode'.");
    
    const CodonAlphabet* pWA = dynamic_cast<const CodonAlphabet*>(alphabet);

    if (args.find("genetic_code") != args.end()) {
      ApplicationTools::displayWarning("'genetic_code' argument is no longer supported inside model description, and has been supersided by a global 'genetic_code' option.");
      throw Exception("BppOFrequenciesSetFormat::read. Deprecated 'genetic_code' argument.");
    }
    short opt = -1;
    string mgmtStopCodon = "quadratic";
    
    if (freqName == "F0")
    {
      opt = CodonFrequenciesSet::F0;
    }
    else if (freqName == "F1X4")
    {
      opt = CodonFrequenciesSet::F1X4;
      
      if (args.find("mgmtStopCodons") != args.end()){
        mgmtStopCodon = args["mgmtStopCodons"];
        ApplicationTools::displayResult("StopCodon frequencies distribution ", mgmtStopCodon);
      }
      if (args.find("frequency") != args.end())
      {
        string sAFS = args["frequency"];
        
        BppOFrequenciesSetFormat nestedReader(alphabetCode_, false, warningLevel_);
        auto_ptr<FrequenciesSet> pFS2(nestedReader.read(pWA->getNAlphabet(0), sAFS, data, false));
        if (pFS2->getName()!="Full")
          throw Exception("BppOFrequenciesSetFormat::read. The frequency option in F1X4 can only be Full");
        
        map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
        
        for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
        {
          unparsedArguments_["123_" + it->first] = it->second;
        }
      }
      
      if (args.find("123_theta") != args.end())
          unparsedArguments_["123_Full.theta"] = args["123_theta"];
      if (args.find("123_theta1") != args.end())
        unparsedArguments_["123_Full.theta1"] = args["123_theta1"];
      if (args.find("123_theta2") != args.end())
        unparsedArguments_["123_Full.theta2"] = args["123_theta2"];
      if (args.find("theta") != args.end())
        unparsedArguments_["123_Full.theta"] = args["123_theta"];
      if (args.find("theta1") != args.end())
        unparsedArguments_["123_Full.theta1"] = args["123_theta1"];
      if (args.find("theta2") != args.end())
        unparsedArguments_["123_Full.theta2"] = args["123_theta2"];
    }
    else if (freqName == "F3X4")
    {
      opt = CodonFrequenciesSet::F3X4;

      if (args.find("mgmtStopCodons") != args.end()){
        mgmtStopCodon = args["mgmtStopCodons"];
        ApplicationTools::displayResult("StopCodon frequencies distribution ", mgmtStopCodon);
      }

      if (args.find("frequency1") != args.end() ||
          args.find("frequency2") != args.end() ||
          args.find("frequency3") != args.end())
        {
        vector<string> v_sAFS;

        for (unsigned int nbfreq = 1; nbfreq<=3; nbfreq++)
          if (args.find("frequency" + TextTools::toString(nbfreq)) != args.end())
            v_sAFS.push_back(args["frequency" + TextTools::toString(nbfreq)]);
          else
            v_sAFS.push_back("");
        
        for (unsigned i = 0; i < v_sAFS.size(); i++)
          {
            BppOFrequenciesSetFormat nestedReader(alphabetCode_, false, warningLevel_);
            if (v_sAFS[i]!=""){
              pFS.reset(nestedReader.read(pWA->getNAlphabet(i), v_sAFS[i], data, false));
              if (pFS->getName()!="Full")
                throw Exception("BppOFrequenciesSetFormat::read. The frequency options in F3X4 can only be Full");

              map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
              for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
                {
                  unparsedArguments_[TextTools::toString(i + 1) + "_" + it->first] = it->second;
                }
            }
          }
        }
      for (unsigned int i = 1 ; i <= 3; i++){
        
        if (args.find(TextTools::toString(i)+"_theta") != args.end())
          unparsedArguments_[TextTools::toString(i)+"_Full.theta"] = args[TextTools::toString(i)+"_theta"];
        if (args.find(TextTools::toString(i)+"_theta1") != args.end())
          unparsedArguments_[TextTools::toString(i)+"_Full.theta1"] = args[TextTools::toString(i)+"_theta1"];
        if (args.find(TextTools::toString(i)+"_theta2") != args.end())
          unparsedArguments_[TextTools::toString(i)+"_Full.theta2"] = args[TextTools::toString(i)+"_theta2"];
      }
    }
    else if (freqName == "F61")
    {
      opt = CodonFrequenciesSet::F61;
    }
    if (opt != -1){
      unsigned short method=1;
      if (args.find("method") != args.end()){
        if (args["method"]=="local")
          method=2;
        else
          if (args["method"]=="binary")
            method=3;
      }      
      pFS.reset(CodonFrequenciesSet::getFrequenciesSetForCodons(opt, geneticCode_, mgmtStopCodon, method));
    }
    else
      throw Exception("Unknown frequency option: " + freqName);
  }
  
  //MVAprotein freq set for COaLA model
  else if (freqName == "MVAprotein")
  {
	  pFS.reset(new MvaFrequenciesSet(dynamic_cast<const ProteicAlphabet*>(alphabet)));
	  dynamic_cast<MvaFrequenciesSet*>(pFS.get())->setParamValues(args);
  }
  else
    throw Exception("Unknown frequency option: " + freqName);

  // Update parameter args:
  vector<string> pnames = pFS->getParameters().getParameterNames();

  string pref = pFS->getNamespace();
  for (size_t i = 0; i < pnames.size(); i++)
  {
    string name = pFS->getParameterNameWithoutNamespace(pnames[i]);
    if (args.find(name) != args.end())
      unparsedArguments_[pref + name] = args[name];
  }

  // Now look if some parameters are aliased:
  ParameterList pl = pFS->getIndependentParameters();
  string pname, pval, pname2;
  for (size_t i = 0; i < pl.size(); i++)
  {
    pname = pFS->getParameterNameWithoutNamespace(pl[i].getName());

    if (args.find(pname) == args.end())
      continue;
    pval = args[pname];

    if (((pval.rfind("_") != string::npos) && (TextTools::isDecimalInteger(pval.substr(pval.rfind("_")+1,string::npos)))) ||
        (pval.find("(") != string::npos))
      continue;
    bool found = false;
    for (size_t j = 0; j < pl.size() && !found; j++)
    {
      pname2 = pFS->getParameterNameWithoutNamespace(pl[j].getName());

      // if (j == i || args.find(pname2) == args.end()) continue; Julien 03/03/2010: This extra condition prevents complicated (nested) models to work properly...
      if (j == i)
        continue;
      if (pval == pname2)
      {
        // This is an alias...
        // NB: this may throw an exception if uncorrect! We leave it as is for now :s
        pFS->aliasParameters(pname2, pname);
        if (verbose_)
          ApplicationTools::displayResult("Parameter alias found", pname + "->" + pname2);
        found = true;
      }
    }
    // if (!TextTools::isDecimalNumber(pval) && !found)
    //   throw Exception("Incorrect parameter syntax: parameter " + pval + " was not found and can't be used as a value for parameter " + pname + ".");
  }

  // Forward arguments:
  if (args.find("init") != args.end())
  {
    unparsedArguments_["init"] = args["init"];
    unparsedArguments_["initFreqs"] = args["init"];
  }
  if (args.find("init.observedPseudoCount") != args.end())
  {
    unparsedArguments_["init.observedPseudoCount"] = args["init.observedPseudoCount"];
    unparsedArguments_["initFreqs.observedPseudoCount"] = args["init.observedPseudoCount"];
  }
  if (args.find("values") != args.end())
  {
    unparsedArguments_["initFreqs"] = "values" + args["values"];
    unparsedArguments_["init"] = "values" + args["values"];
  }

  if (parseArguments)
    initialize_(*pFS, data);
  return pFS.release();
}

void BppOFrequenciesSetFormat::write(const FrequenciesSet* pfreqset,
                                     OutputStream& out,
                                     std::map<std::string, std::string>& globalAliases,
                                     std::vector<std::string>& writtenNames) const
{
  if (!pfreqset)
  {
    out << "None";
    return;
  }
  ParameterList pl = pfreqset->getParameters();

  int p = out.getPrecision();
  out.setPrecision(12);
  bool comma(false);
  string name = pfreqset->getName();
  out << name << "(";

  if ((name == "Fixed") || (name == "F0"))
  {
    vector<double> vf = pfreqset->getFrequencies();
    out << "values=(";
    for (size_t i = 0; i < vf.size(); i++)
    {
      if (i != 0)
        out << ", ";
      out << vf[i];
    }
    out << ")";
    out << ")";
    out.setPrecision(p);
    return;
  }

  
  // For Word or Codon FS : length mgmt
  const WordFromIndependentFrequenciesSet* pWFI = dynamic_cast<const WordFromIndependentFrequenciesSet*>(pfreqset);
  if (name != "F3X4" && pWFI != NULL)
  {
    for (size_t i = 0; i < pWFI->getLength(); i++)
    {
      if (i != 0)
        out << ", ";
      out << "frequency" << i + 1 << "=";
      write(&pWFI->getFrequenciesSetForLetter(i), out, globalAliases, writtenNames);
    }
      comma = true;
  }
  
  const WordFromUniqueFrequenciesSet* pWFU = dynamic_cast<const WordFromUniqueFrequenciesSet*>(pfreqset);
  if (name != "F1X4" && pWFU != NULL)
  {
    for (size_t i = 0; i < pWFU->getLength(); i++)
    {
      if (i != 0)
        out << ", ";
      out << "frequency=";
      write(&pWFU->getFrequenciesSetForLetter(i), out, globalAliases, writtenNames);
    }
    comma = true;
  }

  //FromModel

  const FromModelFrequenciesSet* pFMFS = dynamic_cast<const FromModelFrequenciesSet*>(pfreqset);
  if (pFMFS != NULL)
  {
    if (comma)
      out << ", ";
    out << "model=";
    
    BppOSubstitutionModelFormat bIO(BppOSubstitutionModelFormat::ALL, true, true, true, false, 0);

    bIO.write(*(pFMFS->getModel()), out, globalAliases, writtenNames);
    comma = true;
  }

  
  // FullPerAA
  const FullPerAACodonFrequenciesSet* pFPA=dynamic_cast<const FullPerAACodonFrequenciesSet*>(pfreqset);
  if (pFPA != NULL)
  {
    const ProteinFrequenciesSet* ppfs=pFPA->getProteinFrequenciesSet();
    out << "protein_frequencies=";
    
    write(ppfs, out, globalAliases, writtenNames);
    
    comma = true;
    
    unsigned short meth=pFPA->getMethod();
    if (meth>1){
      if (comma)
        out << ",";
      out << "method=" << ((meth==2)?"local":"binary");
      comma=true;
    }
  }


  // method 
  if (dynamic_cast<const FullProteinFrequenciesSet*>(pfreqset))
  {
    size_t meth=dynamic_cast<const FullProteinFrequenciesSet*>(pfreqset)->getMethod();
    if (meth>1){
      if (comma)
        out << ",";
      out << "method=" << ((meth==2)?"local":"binary");
      comma=true;
    }
  }
  
  const FullCodonFrequenciesSet* pF=dynamic_cast<const FullCodonFrequenciesSet*>(pfreqset);
  if (pF != NULL)
  {
    unsigned short meth=pF->getMethod();
    if (meth>1){
      if (comma)
        out << ",";
      out << "method=" << ((meth==2)?"local":"binary");
      comma=true;
    }
  }

  // mgmtStopCodon 
  const CodonFromUniqueFrequenciesSet* pCFU = dynamic_cast<const CodonFromUniqueFrequenciesSet*>(pfreqset);
  if (pCFU != NULL)
  {
    if (pCFU->getMgmtStopCodon()!="quadratic")
    {
      if (comma)
        out << ",";
      out << "mgmtStopCodons=" << pCFU->getMgmtStopCodon();
      comma = true;
    }
  }
    
  const CodonFromIndependentFrequenciesSet* pCFI = dynamic_cast<const CodonFromIndependentFrequenciesSet*>(pfreqset);
  if (pCFI != NULL)
  {
    if (pCFI->getMgmtStopCodon()!="quadratic")
    {
      if (comma)
        out << ",";
      out << "mgmtStopCodons=" << pCFI->getMgmtStopCodon();
      comma = true;
    }
  }

// All remaining parameters
  const BppOParametrizableFormat* bIO = new BppOParametrizableFormat();

  bIO->write(pfreqset, out, globalAliases, pl.getParameterNames(), writtenNames, true, comma);
  delete bIO;
  
  out << ")";
  out.setPrecision(p);
}

void BppOFrequenciesSetFormat::initialize_(FrequenciesSet& freqSet, const SiteContainer* data)
{
  if (unparsedArguments_.find("init") != unparsedArguments_.end())
  {
    // Initialization using the "init" option
    string init = unparsedArguments_["init"];
    if (init == "observed")
    {
      if (!data)
        throw Exception("Missing data for observed frequencies");
      unsigned int psc = 0;
      if (unparsedArguments_.find("observedPseudoCount") != unparsedArguments_.end())
        psc = TextTools::to<unsigned int>(unparsedArguments_["observedPseudoCount"]);

      map<int, double> freqs;
      SequenceContainerTools::getFrequencies(*data, freqs, psc);

      freqSet.setFrequenciesFromAlphabetStatesFrequencies(freqs);
    }
    else if (init.substr(0, 6) == "values")
    {
      // Initialization using the "values" argument
      vector<double> frequencies;
      string rf = init.substr(6);

      StringTokenizer strtok(rf.substr(1, rf.length() - 2), ",");
      while (strtok.hasMoreToken())
        frequencies.push_back(TextTools::toDouble(strtok.nextToken()));
      freqSet.setFrequencies(frequencies);
    }
    else if (init == "balanced")
    {
      // Nothing to do here, this is the default instanciation.
    }
    else
      throw Exception("Unknown init argument");

    unparsedArguments_.erase(unparsedArguments_.find("init"));
    unparsedArguments_.erase(unparsedArguments_.find("initFreqs"));    
  }

  // Explicit initialization of each parameter
  ParameterList pl = freqSet.getIndependentParameters();
      
  for (size_t i = 0; i < pl.size(); ++i)
  {
    AutoParameter ap(pl[i]);
    if (verbose_)
      ap.setMessageHandler(ApplicationTools::warning);
    pl.setParameter(i, ap);
  }

  for (size_t i = 0; i < pl.size(); ++i)
  {
    const string pName = pl[i].getName();
    
    try {
      double value = ApplicationTools::getDoubleParameter(pName, unparsedArguments_, pl[i].getValue(), "", true, warningLevel_);
      
      pl[i].setValue(value);
      
      if (unparsedArguments_.find(pName) != unparsedArguments_.end())
        unparsedArguments_.erase(unparsedArguments_.find(pName));
      
      if (verbose_)
        ApplicationTools::displayResult("Parameter found", pName + "=" + TextTools::toString(pl[i].getValue()));
    }
    catch (Exception& e) {}
  }
  
  freqSet.matchParametersValues(pl);
}

