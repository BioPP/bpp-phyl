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
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/KeyvalTools.h>
// #include <Bpp/Numeric/Prob.all>
#include <Bpp/Numeric/AutoParameter.h>

//From bpp-seq:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

// From the STL:
#include <iomanip>

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
      if (alphabetCode_ & CODON)
        pFS.reset(new FixedCodonFrequenciesSet(dynamic_cast<const CodonAlphabet*>(alphabet)));
      else
        throw Exception("Codon alphabet not supported.");
    }
    else
    {
      pFS.reset(new FixedFrequenciesSet(alphabet, alphabet->getSize()));
    }
  }
  else if (freqName == "Full")
  {
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
        pFS.reset(new FullProteinFrequenciesSet(dynamic_cast<const ProteicAlphabet*>(alphabet)));
      else
        throw Exception("Protein alphabet not supported.");
    }
    else if (AlphabetTools::isCodonAlphabet(alphabet))
    {
      if (alphabetCode_ & CODON)
        pFS.reset(new FullCodonFrequenciesSet(dynamic_cast<const CodonAlphabet*>(alphabet)));
      else
        throw Exception("Codon alphabet not supported.");
    }
    else
    {
      pFS.reset(new FullFrequenciesSet(alphabet));
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
    if (!(alphabetCode_ & WORD))
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

      BppOFrequenciesSetFormat nestedReader(alphabetCode_, false);
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

      for (unsigned i = 0; i < v_sAFS.size(); i++)
      {
        BppOFrequenciesSetFormat nestedReader(alphabetCode_, false);
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
  // INDEPENDENT CODON
  else if (freqName == "Codon")
  {
    if (!(alphabetCode_ & CODON))
      throw Exception("Codon alphabet not supported.");
    if (!AlphabetTools::isCodonAlphabet(alphabet))
      throw Exception("BppOFrequenciesSetFormat::read.\n\t Bad alphabet type "
                      + alphabet->getAlphabetType() + " for frequenciesset " + freqName + ".");

    const CodonAlphabet* pWA = dynamic_cast<const CodonAlphabet*>(alphabet);

    if (args.find("frequency") != args.end())
    {
      string sAFS = args["frequency"];

      BppOFrequenciesSetFormat nestedReader(alphabetCode_, false);
      auto_ptr<FrequenciesSet> pFS2(nestedReader.read(pWA->getNAlphabet(0), sAFS, data, false));
      map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_["123_" + it->first] = it->second;
      }
      pFS.reset(new CodonFromUniqueFrequenciesSet(pWA, pFS2.release(), "Codon"));
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

      for (unsigned i = 0; i < v_sAFS.size(); i++)
      {
        BppOFrequenciesSetFormat nestedReader(alphabetCode_, false);
        pFS.reset(nestedReader.read(pWA->getNAlphabet(i), v_sAFS[i], data, false));
        map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
        for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
        {
          unparsedArguments_[TextTools::toString(i + 1) + "_" + it->first] = it->second;
        }
        v_AFS.push_back(pFS.release());
      }

      pFS.reset(new CodonFromIndependentFrequenciesSet(pWA, v_AFS, "Codon"));
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

    const CodonAlphabet* pWA = dynamic_cast<const CodonAlphabet*>(alphabet);

    if (args.find("genetic_code") == args.end())
      args["genetic_code"] = pWA->getAlphabetType();

    GeneticCode* pgc = SequenceApplicationTools::getGeneticCode(dynamic_cast<const NucleicAlphabet*>(pWA->getNAlphabet(0)), args["genetic_code"]);
    if (pgc->getSourceAlphabet()->getAlphabetType() != pWA->getAlphabetType())
      throw Exception("Mismatch between genetic code and codon alphabet");

    const ProteicAlphabet* pPA = dynamic_cast<const ProteicAlphabet*>(pgc->getTargetAlphabet());

    auto_ptr<ProteinFrequenciesSet> pPFS;

    if (args.find("protein_frequencies") != args.end())
    {
      string sPFS = args["protein_frequencies"];
      BppOFrequenciesSetFormat nestedReader(alphabetCode_, false);
      pPFS.reset(dynamic_cast<ProteinFrequenciesSet*>(nestedReader.read(pPA, sPFS, data, false)));
      map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_["FullPerAA." + it->first] = it->second;
      }
      pFS.reset(new FullPerAACodonFrequenciesSet(pgc, pPFS.release()));
    }
    else
      pFS.reset(new FullPerAACodonFrequenciesSet(pgc));
  }

  // codeml frequencies syntax
  
  else if (AlphabetTools::isCodonAlphabet(alphabet))
  {
    if (!(alphabetCode_ & CODON))
      throw Exception("Codon alphabet not supported.");
    const CodonAlphabet* pWA = dynamic_cast<const CodonAlphabet*>(alphabet);

    short opt = -1;
    string mgmtStopCodon="quadratic";
    
    if (freqName == "F0")
    {
      opt = CodonFrequenciesSet::F0;
    }
    else if (freqName == "F1X4")
    {
      opt = CodonFrequenciesSet::F1X4;
      
      if (args.find("mgmtStopCodon") != args.end()){
        mgmtStopCodon = args["mgmtStopCodon"];
        ApplicationTools::displayResult("StopCodon frequencies distribution ", mgmtStopCodon);
      }
      if (args.find("frequency") != args.end())
        {
          string sAFS = args["frequency"];

          BppOFrequenciesSetFormat nestedReader(alphabetCode_, false);
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
    }
    else if (freqName == "F3X4")
    {
      opt = CodonFrequenciesSet::F3X4;

      if (args.find("mgmtStopCodon") != args.end()){
        mgmtStopCodon = args["mgmtStopCodon"];
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
            BppOFrequenciesSetFormat nestedReader(alphabetCode_, false);
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
    if (opt != -1)
      pFS.reset(CodonFrequenciesSet::getFrequenciesSetForCodons(opt, *dynamic_cast<const CodonAlphabet*>(alphabet), mgmtStopCodon));
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

  vector<string> pnames = pFS->getParameters().getParameterNames();

  string pref = pFS->getNamespace();
  for (unsigned int i = 0; i < pnames.size(); i++)
  {
    string name = pFS->getParameterNameWithoutNamespace(pnames[i]);
    if (args.find(name) != args.end())
      unparsedArguments_[pref + name] = args[name];
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
                                     std::vector<std::string>& writtenNames) const
{
  if (!pfreqset)
  {
    out << "None";
    return;
  }
  ParameterList pl = pfreqset->getParameters();
  unsigned int p = out.getPrecision();
  out.setPrecision(12);
  bool flag(false);
  string name = pfreqset->getName();
  out << name << "(";


  if ((name == "Fixed") || (name == "F0"))
  {
    vector<double> vf = pfreqset->getFrequencies();
    out << "values=(";
    for (unsigned int i = 0; i < vf.size(); i++)
    {
      if (i != 0)
        out << ", ";
      out << vf[i];
    }
    out << ")";
  }
  else
  {
    if (name != "F1X4" && name != "F3X4" && name != "F61")
    {
      const WordFromIndependentFrequenciesSet* pWFI = dynamic_cast<const WordFromIndependentFrequenciesSet*>(pfreqset);
      if (pWFI != NULL)
      {
        for (unsigned int i = 0; i < pWFI->getLength(); i++)
        {
          if (i != 0)
            out << ", ";
          out << "frequency" << i + 1 << "=";
          write(&pWFI->getFrequenciesSetForLetter(i), out, writtenNames);
        }
        flag = true;
      }
      const WordFromUniqueFrequenciesSet* pWFU = dynamic_cast<const WordFromUniqueFrequenciesSet*>(pfreqset);
      if (pWFU != NULL)
      {
        for (unsigned int i = 0; i < pWFU->getLength(); i++)
        {
          if (i != 0)
            out << ", ";
          out << "frequency=";
          write(&pWFU->getFrequenciesSetForLetter(i), out, writtenNames);
        }
        flag = true;
      }
      const FullPerAACodonFrequenciesSet* pFPA=dynamic_cast<const FullPerAACodonFrequenciesSet*>(pfreqset);
      if (pFPA != NULL)
        {
          const ProteinFrequenciesSet* ppfs=pFPA->getProteinFrequenciesSet();
          out << "protein_frequencies=";
          
          write(ppfs, out, writtenNames);
      
          flag = true;
        }
    }

    for (unsigned int i = 0; i < pl.size(); i++)
      {
        if (find(writtenNames.begin(), writtenNames.end(), pl[i].getName()) == writtenNames.end())
          {
            if (flag)
              out << ",";
            else
              flag = true;
            string pname = pfreqset->getParameterNameWithoutNamespace(pl[i].getName());
            (out << pname << "=").enableScientificNotation(false) << pl[i].getValue();
            writtenNames.push_back(pl[i].getName());
          }
      }
  }
  
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
        psc = TextTools::toInt(unparsedArguments_["observedPseudoCount"]);

      map<int, double> freqs;
      SequenceContainerTools::getFrequencies(*data, freqs, psc);

      freqSet.setFrequenciesFromMap(freqs);
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
  }
  else
  {
    // Explicit initialization of each parameter
    ParameterList pl = freqSet.getParameters();

    for (unsigned int i = 0; i < pl.size(); i++)
    {
      AutoParameter ap(pl[i]);
      if (verbose_)
        ap.setMessageHandler(ApplicationTools::warning);
      pl.setParameter(i, ap);
    }

    for (unsigned int i = 0; i < pl.size(); i++)
    {
      const string pName = pl[i].getName();
      double value = ApplicationTools::getDoubleParameter(pName, unparsedArguments_, pl[i].getValue());

      pl[i].setValue(value);
      if (verbose_)
        ApplicationTools::displayResult("Parameter found", pName + "=" + TextTools::toString(pl[i].getValue()));
    }

    freqSet.matchParametersValues(pl);
  }
}


