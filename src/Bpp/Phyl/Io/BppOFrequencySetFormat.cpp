//
// File: BppOFrequencySetFormat.cpp
// Authors:
//   Laurent GuÃÂ©guen
// Created: lundi 9 juillet 2012, ÃÂ  12h 56
//

/*
  Copyright or ÃÂ© or Copr. Bio++ Development Team, (November 16, 2004)
  
  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.
  
  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".
  
  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.
  
  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.
  
  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/


#include "../Model/FrequencySet/CodonFrequencySet.h"
#include "../Model/FrequencySet/MvaFrequencySet.h"
#include "../Model/FrequencySet/NucleotideFrequencySet.h"
#include "../Model/FrequencySet/ProteinFrequencySet.h"
#include "BppOFrequencySetFormat.h"

// From bpp-core:
#include <Bpp/Io/FileTools.h>
#include <Bpp/Io/BppOParametrizableFormat.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/AutoParameter.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

// From the STL:
#include <iomanip>

#include "BppOTransitionModelFormat.h"

using namespace std;

unsigned char BppOFrequencySetFormat::DNA = 1;
unsigned char BppOFrequencySetFormat::RNA = 2;
unsigned char BppOFrequencySetFormat::NUCLEOTIDE = 1 | 2;
unsigned char BppOFrequencySetFormat::PROTEIN = 4;
unsigned char BppOFrequencySetFormat::CODON = 8;
unsigned char BppOFrequencySetFormat::WORD = 16;
unsigned char BppOFrequencySetFormat::ALL = 1 | 2 | 4 | 8 | 16;


std::shared_ptr<FrequencySet> BppOFrequencySetFormat::readFrequencySet(const Alphabet* alphabet, const std::string& freqDescription, const AlignedValuesContainer* data, bool parseArguments)
{
  unparsedArguments_.clear();
  string freqName;
  map<string, string> args;
  KeyvalTools::parseProcedure(freqDescription, freqName, args);
  shared_ptr<FrequencySet> pFS;

  if (freqName == "Fixed")
  {
    if (AlphabetTools::isNucleicAlphabet(alphabet))
    {
      if (alphabetCode_ & NUCLEOTIDE)
        pFS.reset(new FixedNucleotideFrequencySet(dynamic_cast<const NucleicAlphabet*>(alphabet)));
      else
        throw Exception("Nucleotide alphabet not supported.");
    }
    else if (AlphabetTools::isProteicAlphabet(alphabet))
    {
      if (alphabetCode_ & PROTEIN)
        pFS.reset(new FixedProteinFrequencySet(dynamic_cast<const ProteicAlphabet*>(alphabet)));
      else
        throw Exception("Protein alphabet not supported.");
    }
    else if (AlphabetTools::isCodonAlphabet(alphabet))
    {
      if (alphabetCode_ & CODON)
      {
        if (!geneticCode_)
          throw Exception("BppOFrequencySetFormat::readFrequencySet(). No genetic code specified! Consider using 'setGeneticCode'.");
        pFS.reset(new FixedCodonFrequencySet(geneticCode_));
      }
      else
      {
        throw Exception("Codon alphabet not supported.");
      }
    }
    else
    {
      pFS.reset(new FixedFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false))));
    }
  }
  else if (freqName == "Full")
  {
    unsigned short method = 1;
    if (args.find("method") != args.end())
    {
      if (args["method"] == "local")
        method = 2;
      else if (args["method"] == "binary")
        method = 3;
    }

    if (AlphabetTools::isNucleicAlphabet(alphabet))
    {
      if (alphabetCode_ & NUCLEOTIDE)
        pFS.reset(new FullNucleotideFrequencySet(dynamic_cast<const NucleicAlphabet*>(alphabet)));
      else
        throw Exception("Nucleotide alphabet not supported.");
    }
    else if (AlphabetTools::isProteicAlphabet(alphabet))
    {
      if (alphabetCode_ & PROTEIN)
        pFS.reset(new FullProteinFrequencySet(dynamic_cast<const ProteicAlphabet*>(alphabet), false, method));
      else
        throw Exception("Protein alphabet not supported.");
    }
    else if (AlphabetTools::isCodonAlphabet(alphabet))
    {
      if (alphabetCode_ & CODON)
      {
        if (!geneticCode_)
          throw Exception("BppOFrequencySetFormat::readFrequencySet(). No genetic code specified! Consider using 'setGeneticCode'.");
        pFS.reset(new FullCodonFrequencySet(geneticCode_));
      }
      else
      {
        throw Exception("Codon alphabet not supported.");
      }
    }
    else
    {
      // NB: jdutheil 25/09/14 => gap models will not be supported before we add the appropriate option!
      pFS.reset(new FullFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), false, method));
    }
  }
  else if (freqName == "Empirical")
  {
    string fname = args["file"];
    if (TextTools::isEmpty(fname))
      throw Exception("'file' argument missing for user-defined frequencies.");

    size_t nCol = 1;
    if (args.find("col") != args.end())
      nCol = size_t(TextTools::toInt(args["col"]));

    if (AlphabetTools::isNucleicAlphabet(alphabet))
    {
      if (alphabetCode_ & NUCLEOTIDE)
        pFS.reset(new UserNucleotideFrequencySet(dynamic_cast<const NucleicAlphabet*>(alphabet), fname, nCol));
      else
        throw Exception("Nucleotide alphabet not supported.");
    }
    else if (AlphabetTools::isProteicAlphabet(alphabet))
    {
      if (alphabetCode_ & PROTEIN)
        pFS.reset(new UserProteinFrequencySet(dynamic_cast<const ProteicAlphabet*>(alphabet), fname, nCol));
      else
        throw Exception("Protein alphabet not supported.");
    }
    else if (AlphabetTools::isCodonAlphabet(alphabet))
    {
      if (alphabetCode_ & CODON)
      {
        if (!geneticCode_)
          throw Exception("BppOFrequencySetFormat::readFrequencySet(). No genetic code specified! Consider using 'setGeneticCode'.");
        pFS.reset(new UserCodonFrequencySet(geneticCode_, fname, nCol));
      }
      else
      {
        throw Exception("Codon alphabet not supported.");
      }
    }
    else
    {
      pFS.reset(new UserFrequencySet(std::shared_ptr<const StateMap>(new CanonicalStateMap(alphabet, false)), fname, nCol));
    }
  }
  else if (freqName == "GC")
  {
    if (!AlphabetTools::isNucleicAlphabet(alphabet))
      throw Exception("Error, unvalid frequencies " + freqName + " with non-nucleic alphabet.");
    if (alphabetCode_ & NUCLEOTIDE)
      pFS.reset(new GCFrequencySet(dynamic_cast<const NucleicAlphabet*>(alphabet)));
    else
      throw Exception("Nucleotide alphabet not supported.");
  }

  // INDEPENDENTWORD
  else if (freqName == "Word")
  {
    if (!(alphabetCode_ & WORD))
      throw Exception("Word alphabet not supported.");
    if (!AlphabetTools::isWordAlphabet(alphabet))
      throw Exception("BppOFrequencySetFormat::readFrequencySet(...).\n\t Bad alphabet type "
                      + alphabet->getAlphabetType() + " for frequencies set " + freqName + ".");

    const WordAlphabet* pWA = dynamic_cast<const WordAlphabet*>(alphabet);

    if (pWA == NULL)
      throw Exception("BppOFrequencySetFormat::read : Word freq name is from WordAlphabet.");

    if (args.find("frequency") != args.end())
    {
      string sAFS = args["frequency"];

      unsigned int nbfreq = pWA->getLength();
      string st = "";
      for (unsigned i = 0; i < nbfreq; i++)
      {
        st += TextTools::toString(i + 1);
      }

      BppOFrequencySetFormat nestedReader(alphabetCode_, false, warningLevel_);
      auto pFS2(nestedReader.readFrequencySet(pWA->getNAlphabet(0), sAFS, data, false));
      map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_[st + "_" + it->first] = it->second;
      }
      pFS.reset(new WordFromUniqueFrequencySet(pWA, pFS2));
    }
    else
    {
      if (args.find("frequency1") == args.end())
        throw Exception("BppOFrequencySetFormat::read. Missing argument 'frequency' or 'frequency1' for frequencies set 'Word'.");
      vector<string> v_sAFS;
      vector<std::shared_ptr<FrequencySet> > v_AFS;
      unsigned int nbfreq = 1;

      while (args.find("frequency" + TextTools::toString(nbfreq)) != args.end())
      {
        v_sAFS.push_back(args["frequency" + TextTools::toString(nbfreq++)]);
      }

      if (v_sAFS.size() != pWA->getLength())
        throw Exception("BppOFrequencySetFormat::read. Number of frequencies (" + TextTools::toString(v_sAFS.size()) + ") does not match length of the words (" + TextTools::toString(pWA->getLength()) + ")");

      for (size_t i = 0; i < v_sAFS.size(); ++i)
      {
        BppOFrequencySetFormat nestedReader(alphabetCode_, false, warningLevel_);
        pFS = nestedReader.readFrequencySet(pWA->getNAlphabet(i), v_sAFS[i], data, false);
        map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
        for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
        {
          unparsedArguments_[TextTools::toString(i + 1) + "_" + it->first] = it->second;
        }
        v_AFS.push_back(pFS);
      }

      pFS.reset(new WordFromIndependentFrequencySet(pWA, v_AFS));
    }
  }
  // From Model
  else if (freqName == "FromModel")
  {
    if (args.find("model") == args.end())
      throw Exception("Missing argument 'model' for frequencies " + freqName + ".");

    BppOTransitionModelFormat nestedReader(alphabetCode_, false, true, false, false, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_);

    TransitionModel* model = nestedReader.readTransitionModel(alphabet, args["model"], data, false);
    pFS.reset(new FromModelFrequencySet(model));
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
      throw Exception("BppOFrequencySetFormat::read.\n\t Bad alphabet type "
                      + alphabet->getAlphabetType() + " for frequenciesset " + freqName + ".");
    if (!geneticCode_)
      throw Exception("BppOFrequencySetFormat::readFrequencySet(). No genetic code specified! Consider using 'setGeneticCode'.");

    const CodonAlphabet* pWA = dynamic_cast<const CodonAlphabet*>(alphabet);

    if (args.find("frequency") != args.end())
    {
      string sAFS = args["frequency"];

      BppOFrequencySetFormat nestedReader(alphabetCode_, false, warningLevel_);
      auto pFS2(nestedReader.readFrequencySet(pWA->getNucleicAlphabet(), sAFS, data, false));
      map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_["123_" + it->first] = it->second;
      }

      pFS.reset(new CodonFromUniqueFrequencySet(geneticCode_, pFS2, "Codon"));
    }
    else
    {
      if (args.find("frequency1") == args.end())
        throw Exception("BppOFrequencySetFormat::read. Missing argument 'frequency' or 'frequency1' for frequencies set.");
      vector<string> v_sAFS;
      vector<std::shared_ptr<FrequencySet> > v_AFS;
      unsigned int nbfreq = 1;

      while (args.find("frequency" + TextTools::toString(nbfreq)) != args.end())
      {
        v_sAFS.push_back(args["frequency" + TextTools::toString(nbfreq++)]);
      }

      if (v_sAFS.size() != 3)
        throw Exception("BppOFrequencySetFormat::read. Number of frequencies (" + TextTools::toString(v_sAFS.size()) + ") is not three");

      for (size_t i = 0; i < v_sAFS.size(); ++i)
      {
        BppOFrequencySetFormat nestedReader(alphabetCode_, false, warningLevel_);
        pFS = nestedReader.readFrequencySet(pWA->getNucleicAlphabet(), v_sAFS[i], data, false);
        map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
        for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
        {
          unparsedArguments_[TextTools::toString(i + 1) + "_" + it->first] = it->second;
        }
        v_AFS.push_back(pFS);
      }

      if (!geneticCode_)
        throw Exception("BppOFrequencySetFormat::readFrequencySet(). No genetic code specified! Consider using 'setGeneticCode'.");
      pFS.reset(new CodonFromIndependentFrequencySet(geneticCode_, v_AFS, "Codon"));
    }
  }

  // CODON PER AA Frequencies
  else if (freqName == "FullPerAA")
  {
    if (!(alphabetCode_ & CODON))
      throw Exception("Codon alphabet not supported.");
    if (!AlphabetTools::isCodonAlphabet(alphabet))
      throw Exception("BppOFrequencySetFormat::read.\n\t Bad alphabet type "
                      + alphabet->getAlphabetType() + " for frequenciesset " + freqName + ".");

    if (!geneticCode_)
      throw Exception("BppOFrequencySetFormat::readFrequencySet(). No genetic code specified! Consider using 'setGeneticCode'.");

    const ProteicAlphabet* pPA = dynamic_cast<const ProteicAlphabet*>(geneticCode_->getTargetAlphabet());

    unsigned short method = 1;
    if (args.find("method") != args.end())
    {
      if (args["method"] == "local")
        method = 2;
      else if (args["method"] == "binary")
        method = 3;
    }

    if (args.find("protein_frequencies") != args.end())
    {
      string sPFS = args["protein_frequencies"];
      BppOFrequencySetFormat nestedReader(alphabetCode_, false, warningLevel_);
      auto pPFS(dynamic_pointer_cast<ProteinFrequencySet>(nestedReader.readFrequencySet(pPA, sPFS, data, false)));
      map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

      for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
      {
        unparsedArguments_["FullPerAA." + it->first] = it->second;
      }
      pFS.reset(new FullPerAACodonFrequencySet(geneticCode_, pPFS, method));
    }
    else
      pFS.reset(new FullPerAACodonFrequencySet(geneticCode_, method));
  }

  // codeml frequencies syntax

  else if (AlphabetTools::isCodonAlphabet(alphabet))
  {
    if (!(alphabetCode_ & CODON))
      throw Exception("Codon alphabet not supported.");
    if (!geneticCode_)
      throw Exception("BppOFrequencySetFormat::readFrequencySet(). No genetic code specified! Consider using 'setGeneticCode'.");

    const CodonAlphabet* pWA = dynamic_cast<const CodonAlphabet*>(alphabet);

    if (args.find("genetic_code") != args.end())
    {
      ApplicationTools::displayWarning("'genetic_code' argument is no longer supported inside model description, and has been supersided by a global 'genetic_code' option.");
      throw Exception("BppOFrequencySetFormat::read. Deprecated 'genetic_code' argument.");
    }
    short opt = -1;
    string mgmtStopCodon = "quadratic";

    if (freqName == "F0")
    {
      opt = CodonFrequencySet::F0;
    }
    else if (freqName == "F1X4")
    {
      opt = CodonFrequencySet::F1X4;

      if (args.find("mgmtStopCodons") != args.end())
        mgmtStopCodon = args["mgmtStopCodons"];
      else if (args.find("mgmtStopCodon") != args.end())
        mgmtStopCodon = args["mgmtStopCodon"];

      if (verbose_)
        ApplicationTools::displayResult("StopCodon frequencies distribution ", mgmtStopCodon);

      if (args.find("frequency") != args.end())
      {
        string sAFS = args["frequency"];

        BppOFrequencySetFormat nestedReader(alphabetCode_, false, warningLevel_);
        auto pFS2(nestedReader.readFrequencySet(pWA->getNucleicAlphabet(), sAFS, data, false));
        if (pFS2->getName() != "Full")
          throw Exception("BppOFrequencySetFormat::read. The frequency option in F1X4 can only be Full");

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
      opt = CodonFrequencySet::F3X4;

      if (args.find("mgmtStopCodons") != args.end())
        mgmtStopCodon = args["mgmtStopCodons"];
      else if (args.find("mgmtStopCodon") != args.end())
        mgmtStopCodon = args["mgmtStopCodon"];

      if (verbose_)
        ApplicationTools::displayResult("StopCodon frequencies distribution ", mgmtStopCodon);

      if (args.find("frequency1") != args.end() ||
          args.find("frequency2") != args.end() ||
          args.find("frequency3") != args.end())
      {
        vector<string> v_sAFS;

        for (unsigned int nbfreq = 1; nbfreq <= 3; nbfreq++)
        {
          if (args.find("frequency" + TextTools::toString(nbfreq)) != args.end())
            v_sAFS.push_back(args["frequency" + TextTools::toString(nbfreq)]);
          else
            v_sAFS.push_back("");
        }

        for (unsigned i = 0; i < v_sAFS.size(); i++)
        {
          BppOFrequencySetFormat nestedReader(alphabetCode_, false, warningLevel_);
          if (v_sAFS[i] != "")
          {
            pFS = nestedReader.readFrequencySet(pWA->getNucleicAlphabet(), v_sAFS[i], data, false);
            if (pFS->getName() != "Full")
              throw Exception("BppOFrequencySetFormat::read. The frequency options in F3X4 can only be Full");

            map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());
            for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
            {
              unparsedArguments_[TextTools::toString(i + 1) + "_" + it->first] = it->second;
            }
          }
        }
      }
      for (unsigned int i = 1; i <= 3; i++)
      {
        if (args.find(TextTools::toString(i) + "_theta") != args.end())
          unparsedArguments_[TextTools::toString(i) + "_Full.theta"] = args[TextTools::toString(i) + "_theta"];
        if (args.find(TextTools::toString(i) + "_theta1") != args.end())
          unparsedArguments_[TextTools::toString(i) + "_Full.theta1"] = args[TextTools::toString(i) + "_theta1"];
        if (args.find(TextTools::toString(i) + "_theta2") != args.end())
          unparsedArguments_[TextTools::toString(i) + "_Full.theta2"] = args[TextTools::toString(i) + "_theta2"];
      }
    }
    else if (freqName == "F61")
    {
      opt = CodonFrequencySet::F61;
    }
    if (opt != -1)
    {
      unsigned short method = 1;
      if (args.find("method") != args.end())
      {
        if (args["method"] == "local")
          method = 2;
        else if (args["method"] == "binary")
          method = 3;
      }
      pFS = CodonFrequencySet::getFrequencySetForCodons(opt, geneticCode_, mgmtStopCodon, method);
    }
    else
      throw Exception("Unknown frequency option: " + freqName);
  }

  // MVAprotein freq set for COaLA model
  else if (freqName == "MVAprotein")
  {
    pFS.reset(new MvaFrequencySet(dynamic_cast<const ProteicAlphabet*>(alphabet)));
    dynamic_pointer_cast<MvaFrequencySet>(pFS)->setParamValues(args);
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

    if (((pval.rfind("_") != string::npos) && (TextTools::isDecimalInteger(pval.substr(pval.rfind("_") + 1, string::npos)))) ||
        (pval.find("(") != string::npos))
      continue;
    bool found = false;
    for (size_t j = 0; j < pl.size() && !found; j++)
    {
      pname2 = pFS->getParameterNameWithoutNamespace(pl[j].getName());

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
    unparsedArguments_["init"] = args["init"];

  if (args.find("init.observedPseudoCount") != args.end())
    unparsedArguments_["init.observedPseudoCount"] = args["init.observedPseudoCount"];

  if (args.find("values") != args.end())
  {
    unparsedArguments_["init"] = "values" + args["values"];
    initialize_(*pFS, data);
  }
  else if (parseArguments)
    initialize_(*pFS, data);

  return pFS;
}

void BppOFrequencySetFormat::writeFrequencySet(
  const FrequencySet* pfreqset,
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

  bool oValues = false;

  if ((name == "Fixed") || (name == "F0") || (name == "Full"))
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

    vector<string> npl = pl.getParameterNames();
    writtenNames.insert(writtenNames.end(), npl.begin(), npl.end());

    comma = true;

    oValues = true;
  }

  auto ufs = dynamic_cast<const UserFrequencySet*>(pfreqset);
  if (ufs)
  {
    out << "file=" << ufs->getPath();
    size_t nCol = ufs->getColumnNumber();

    if (nCol != 1)
      out << ", col=" << nCol;

    oValues = true;
  }

  // For Word or Codon FS : length mgmt
  auto pWFI = dynamic_cast<const WordFromIndependentFrequencySet*>(pfreqset);
  if (name != "F3X4" && pWFI != NULL)
  {
    for (size_t i = 0; i < pWFI->getLength(); i++)
    {
      if (i != 0)
        out << ", ";
      out << "frequency" << i + 1 << "=";
      writeFrequencySet(pWFI->getFrequencySetForLetter(i).get(), out, globalAliases, writtenNames);
    }
    comma = true;
  }

  auto pWFU = dynamic_cast<const WordFromUniqueFrequencySet*>(pfreqset);
  if (name != "F1X4" && pWFU != NULL)
  {
    for (size_t i = 0; i < pWFU->getLength(); i++)
    {
      if (i != 0)
        out << ", ";
      out << "frequency=";
      writeFrequencySet(pWFU->getFrequencySetForLetter(i).get(), out, globalAliases, writtenNames);
    }
    comma = true;
  }

  // FromModel

  auto pFMFS = dynamic_cast<const FromModelFrequencySet*>(pfreqset);
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
  auto pFPA = dynamic_cast<const FullPerAACodonFrequencySet*>(pfreqset);
  if (pFPA != NULL)
  {
    auto ppfs = pFPA->getProteinFrequencySet();
    out << "protein_frequencies=";

    writeFrequencySet(ppfs.get(), out, globalAliases, writtenNames);

    comma = true;

    unsigned short meth = pFPA->getMethod();
    if (meth > 1)
    {
      if (comma)
        out << ",";
      out << "method=" << ((meth == 2) ? "local" : "binary");
      comma = true;
    }
  }


  // method
  if (dynamic_cast<const FullProteinFrequencySet*>(pfreqset))
  {
    size_t meth = dynamic_cast<const FullProteinFrequencySet*>(pfreqset)->getMethod();
    if (meth > 1)
    {
      if (comma)
        out << ",";
      out << "method=" << ((meth == 2) ? "local" : "binary");
      comma = true;
    }
  }

  auto pF = dynamic_cast<const FullCodonFrequencySet*>(pfreqset);
  if (pF != NULL)
  {
    unsigned short meth = pF->getMethod();
    if (meth > 1)
    {
      if (comma)
        out << ",";
      out << "method=" << ((meth == 2) ? "local" : "binary");
      comma = true;
    }
  }

  // mgmtStopCodon
  auto pCFU = dynamic_cast<const CodonFromUniqueFrequencySet*>(pfreqset);
  if (pCFU != NULL)
  {
    if (pCFU->getMgmtStopCodon() != "quadratic")
    {
      if (comma)
        out << ",";
      out << "mgmtStopCodons=" << pCFU->getMgmtStopCodon();
      comma = true;
    }
  }

  auto pCFI = dynamic_cast<const CodonFromIndependentFrequencySet*>(pfreqset);
  if (pCFI != NULL)
  {
    if (pCFI->getMgmtStopCodon() != "quadratic")
    {
      if (comma)
        out << ",";
      out << "mgmtStopCodons=" << pCFI->getMgmtStopCodon();
      comma = true;
    }
  }

// All remaining parameters
  if (!oValues)
  {
    const BppOParametrizableFormat* bIO = new BppOParametrizableFormat();

    bIO->write(pfreqset, out, globalAliases, pl.getParameterNames(), writtenNames, true, comma);
    delete bIO;
  }

  out << ")";
  out.setPrecision(p);
}

void BppOFrequencySetFormat::initialize_(FrequencySet& freqSet, const AlignedValuesContainer* data)
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
  }

  // Explicit initialization of each parameter
  ParameterList pl = freqSet.getIndependentParameters();

  for (size_t i = 0; i < pl.size(); ++i)
  {
    AutoParameter ap(pl[i]);
    if (verbose_)
      ap.setMessageHandler(ApplicationTools::warning.get());
    pl.setParameter(i, ap);
  }

  for (size_t i = 0; i < pl.size(); ++i)
  {
    const string pName = pl[i].getName();

    try
    {
      double value = ApplicationTools::getDoubleParameter(pName, unparsedArguments_, pl[i].getValue(), "", true, warningLevel_);

      pl[i].setValue(value);

      if (unparsedArguments_.find(pName) != unparsedArguments_.end())
        unparsedArguments_.erase(unparsedArguments_.find(pName));

      if (verbose_)
        ApplicationTools::displayResult("Parameter found", pName + "=" + TextTools::toString(pl[i].getValue()));
    }
    catch (Exception& e)
    {}
  }

  freqSet.matchParametersValues(pl);
}
