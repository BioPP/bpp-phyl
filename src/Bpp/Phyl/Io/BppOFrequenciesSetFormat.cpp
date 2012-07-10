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

#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/KeyvalTools.h>

#include "../Model/FrequenciesSet.all"
#include "../App/PhylogeneticsApplicationTools.h"

#include <Bpp/Seq/App/SequenceApplicationTools.h>

#include <Bpp/Numeric/Prob.all>

using namespace bpp;

// From the STL:
#include <iomanip>

using namespace std;

FrequenciesSet* BppOFrequenciesSetFormat::read(const Alphabet* alphabet,
                                               const std::string& freqDescription,
                                               std::map<std::string, std::string>& unparsedParameterValues)
{
  string freqName;
  map<string, string> args;
  KeyvalTools::parseProcedure(freqDescription, freqName, args);
  FrequenciesSet* pFS;

  if (freqName == "Fixed")
    {
      if (AlphabetTools::isNucleicAlphabet(alphabet))
        {
          pFS = new FixedNucleotideFrequenciesSet(dynamic_cast<const NucleicAlphabet*>(alphabet));
        }
      else if (AlphabetTools::isProteicAlphabet(alphabet))
        {
          pFS = new FixedProteinFrequenciesSet(dynamic_cast<const ProteicAlphabet*>(alphabet));
        }
      else if (AlphabetTools::isCodonAlphabet(alphabet))
        {
          pFS = new FixedCodonFrequenciesSet(dynamic_cast<const CodonAlphabet*>(alphabet));
        }
      else
        {
          pFS = new FixedFrequenciesSet(alphabet);
        }
    }
  else if (freqName == "Full")
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
    }
  else if (freqName == "GC")
    {
      if (!AlphabetTools::isNucleicAlphabet(alphabet))
        throw Exception("Error, unvalid frequencies " + freqName + " with non-nucleic alphabet.");
      
      pFS = new GCFrequenciesSet(dynamic_cast<const NucleicAlphabet*>(alphabet));
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
              st += TextTools::toString(i + 1);
            }
          
          map<string, string> unparsedParameterValuesNested;
          pFS2 = read(pWA->getNAlphabet(0), sAFS, unparsedParameterValuesNested);
          for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
            {
              unparsedParameterValues[st + "_" + it->first] = it->second;
            }
          pFS = new WordFromUniqueFrequenciesSet(pWA, pFS2);
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
          
          map<string, string> unparsedParameterValuesNested;
          for (unsigned i = 0; i < v_sAFS.size(); i++)
            {
              unparsedParameterValuesNested.clear();
              pFS = read(pWA->getNAlphabet(i), v_sAFS[i], unparsedParameterValuesNested);
              for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
                {
                  unparsedParameterValues[TextTools::toString(i + 1) + "_" + it->first] = it->second;
                }
              v_AFS.push_back(pFS);
            }
          
          pFS = new WordFromIndependentFrequenciesSet(pWA, v_AFS);
        }
    }
  // INDEPENDENT CODON
  else if (freqName == "Codon")
    {
      if (!AlphabetTools::isCodonAlphabet(alphabet))
        throw Exception("BppOFrequenciesSetFormat::read.\n\t Bad alphabet type "
                        + alphabet->getAlphabetType() + " for frequenciesset " + freqName + ".");
      
      const CodonAlphabet* pWA = dynamic_cast<const CodonAlphabet*>(alphabet);
      
      if (args.find("frequency") != args.end())
        {
          string sAFS = args["frequency"];
          
          unsigned int nbfreq = pWA->getLength();
          FrequenciesSet* pFS2;
          string st = "";
          for (unsigned i = 0; i < nbfreq; i++)
            {
              st += TextTools::toString(i + 1);
            }
          
          map<string, string> unparsedParameterValuesNested;
          pFS2 = read(pWA->getNAlphabet(0), sAFS, unparsedParameterValuesNested);
          
          for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
            {
              unparsedParameterValues[st + "_" + it->first] = it->second;
            }
          pFS = new CodonFromUniqueFrequenciesSet(pWA, pFS2, "Codon");
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
          
          map<string, string> unparsedParameterValuesNested;
          for (unsigned i = 0; i < v_sAFS.size(); i++)
            {
              unparsedParameterValuesNested.clear();
              pFS = read(pWA->getNAlphabet(i), v_sAFS[i], unparsedParameterValuesNested);
              for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
                {
                  unparsedParameterValues[TextTools::toString(i + 1) + "_" + it->first] = it->second;
                }
              v_AFS.push_back(pFS);
            }
          
          pFS = new CodonFromIndependentFrequenciesSet(pWA, v_AFS, "Codon");
        }
    }
  
  // CODON PER AA Frequencies 
  else if (freqName == "FullPerAA")
    {
      if (!AlphabetTools::isCodonAlphabet(alphabet))
        throw Exception("BppOFrequenciesSetFormat::read.\n\t Bad alphabet type "
                        + alphabet->getAlphabetType() + " for frequenciesset " + freqName + ".");
      
      const CodonAlphabet* pWA = dynamic_cast<const CodonAlphabet*>(alphabet);
      
      if (args.find("genetic_code") == args.end())
        args["genetic_code"] = pWA->getAlphabetType();

      GeneticCode* pgc = SequenceApplicationTools::getGeneticCode(dynamic_cast<const NucleicAlphabet*>(pWA->getNAlphabet(0)), args["genetic_code"]);
      if (pgc->getSourceAlphabet()->getAlphabetType() != pWA->getAlphabetType())
        throw Exception("Mismatch between genetic code and codon alphabet");

      const ProteicAlphabet* pPA=dynamic_cast<const ProteicAlphabet*>(pgc->getTargetAlphabet());

      ProteinFrequenciesSet* pPFS=0;
      
      map<string, string> unparsedParameterValuesNested;
      
      if (args.find("protein_frequencies") != args.end()){
        string sPFS= args["protein_frequencies"];
        pPFS= dynamic_cast<ProteinFrequenciesSet*>(read(pPA, sPFS, unparsedParameterValuesNested));

        for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
          {
            unparsedParameterValues["FullPerAA." + it->first] = it->second;
          }
        pFS = new FullPerAACodonFrequenciesSet(pgc, pPFS);
      }

      else
        pFS = new FullPerAACodonFrequenciesSet(pgc);
    }
  else if (AlphabetTools::isCodonAlphabet(alphabet))
    {
      short opt = -1;

      if (freqName == "F0"){
        opt = CodonFrequenciesSet::F0;
      }
      else if (freqName == "F1X4"){
        opt = CodonFrequenciesSet::F1X4;
      }
      else if (freqName == "F3X4"){
        opt = CodonFrequenciesSet::F3X4;
      }
      else if (freqName == "F61"){
        opt = CodonFrequenciesSet::F61;
      }
      if (opt != -1)
        pFS = CodonFrequenciesSet::getFrequenciesSetForCodons(opt, *dynamic_cast<const CodonAlphabet*>(alphabet));
      else
        throw Exception("Unknown frequency option: " + freqName);
    }
  else
    throw Exception("Unknown frequency option: " + freqName);

  // initial values
  
  if (args.find("values")!= args.end())
    {
      string rf= args["values"];
      // Initialization using the "values" argument
      vector<double> frequencies;
      StringTokenizer strtok(rf.substr(1, rf.length() - 2), ",");
      while (strtok.hasMoreToken())
        frequencies.push_back(TextTools::toDouble(strtok.nextToken()));
      pFS->setFrequencies(frequencies);
    }

  // Forward arguments:
  if (args.find("init") != args.end())
    {
      unparsedParameterValues["init"] = args["init"];
      unparsedParameterValues["initFreqs"] = args["init"];
    }
  if (args.find("init.observedPseudoCount") != args.end())
    {
      unparsedParameterValues["init.observedPseudoCount"] = args["init.observedPseudoCount"];
      unparsedParameterValues["initFreqs.observedPseudoCount"] = args["init.observedPseudoCount"];
    }

  vector<string> pnames = pFS->getParameters().getParameterNames();

  string pref=pFS->getNamespace();
  for (unsigned int i = 0; i < pnames.size(); i++)
    {
      string name = pFS->getParameterNameWithoutNamespace(pnames[i]);
      if (args.find(name) != args.end())
          unparsedParameterValues[pref + name] = args[name];
    }

  return pFS;
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
  out << pfreqset->getName() << "(";

  
  if ((pfreqset->getName()=="Fixed") || (pfreqset->getName()=="F0")){
    vector<double> vf=pfreqset->getFrequencies();
    out << "values=(" ;
    for (unsigned int i=0;i<vf.size();i++){
      if (i!=0)
        out << ", ";
      out << vf[i];
    }
    out << ")";
  }
  else {
    const WordFromIndependentFrequenciesSet* pWFI=dynamic_cast<const WordFromIndependentFrequenciesSet*>(pfreqset);
    if (pWFI!=NULL){
      for (unsigned int i=0; i< pWFI->getLength(); i++){
        if (i!=0)
          out << ", ";
        out << "frequency" << i+1 << "="; 
        write(&pWFI->getFrequenciesSetForLetter(i), out, writtenNames);
      }
      flag=true;
    }
    const WordFromUniqueFrequenciesSet* pWFU=dynamic_cast<const WordFromUniqueFrequenciesSet*>(pfreqset);
    if (pWFU!=NULL){
      for (unsigned int i=0; i< pWFU->getLength(); i++){
        if (i!=0)
          out << ", ";
        out << "frequency=";
        write(&pWFU->getFrequenciesSetForLetter(i), out, writtenNames);
      }
      flag=true;
    }

    for (unsigned int i = 0; i < pl.size(); i++)
      {
        if (find(writtenNames.begin(),writtenNames.end(),pl[i].getName())==writtenNames.end()){
          if (flag)
            out << ",";
          else
            flag=true;
          string pname = pfreqset->getParameterNameWithoutNamespace(pl[i].getName());
          (out << pname << "=").enableScientificNotation(false) << pl[i].getValue();
          writtenNames.push_back(pl[i].getName());
        }
      }
  }
  out << ")";
  out.setPrecision(p);
}
