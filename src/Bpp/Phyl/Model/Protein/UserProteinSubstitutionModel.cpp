//
// File: UserProteinSubstitutionModel.cpp
// Created by: Julien Dutheil
// Created on: Wed Aug 26 16:27 2005
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

#include "UserProteinSubstitutionModel.h"

#include <Bpp/Io/FileTools.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Numeric/VectorTools.h>

//From SeqLib:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

// From the STL:
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/******************************************************************************/

UserProteinSubstitutionModel::UserProteinSubstitutionModel(
    const ProteicAlphabet* alpha, const std::string& path, const std::string& prefix) : 
  AbstractParameterAliasable(prefix),
  AbstractReversibleProteinSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), prefix),
  path_(path),
  freqSet_(0)
{
  readFromFile();
  freqSet_ = new FixedProteinFrequenciesSet(alpha, freq_);
  updateMatrices();  
}

UserProteinSubstitutionModel::UserProteinSubstitutionModel(
    const ProteicAlphabet* alpha, const std::string& path,
    ProteinFrequenciesSet* freqSet, const std::string& prefix,
    bool initFreqs) : 
  AbstractParameterAliasable(prefix),
  AbstractReversibleProteinSubstitutionModel(alpha, new CanonicalStateMap(alpha, false), prefix),
  path_(path),
  freqSet_(freqSet)
{
  readFromFile();
  freqSet_->setNamespace(prefix+freqSet_->getNamespace());
  if (initFreqs) freqSet->setFrequencies(freq_);
  else freq_ = freqSet_->getFrequencies();
  addParameters_(freqSet_->getParameters());
  updateMatrices();  
}

/******************************************************************************/

std::string UserProteinSubstitutionModel::getName() const
{
  if (TextTools::hasSubstring(freqSet_->getNamespace(), "+F.") )
    return "Empirical+F"; 
  else  
    return "Empirical";
}

/******************************************************************************/

void UserProteinSubstitutionModel::readFromFile()
{
  ifstream in(path_.c_str(), ios::in);
  //Read exchangeability matrix:
  for (unsigned int i = 1; i < 20; i++)
  {
    string line = FileTools::getNextLine(in);
    StringTokenizer st(line);
    for(unsigned int j = 0; j < i; j++) {
      double s = TextTools::toDouble(st.nextToken());
      exchangeability_(i,j) = exchangeability_(j,i) = s;
    }
  }
  //Read frequencies:
  unsigned int fCount = 0;
  while (in && fCount < 20)
  {
    string line = FileTools::getNextLine(in);
    StringTokenizer st(line);
    while(st.hasMoreToken() && fCount < 20)
    {
      freq_[fCount] = TextTools::toDouble(st.nextToken());
      fCount++;
    }
  }
  double sf = VectorTools::sum(freq_);
  if (sf - 1 > 0.000001)
  {
    ApplicationTools::displayMessage("WARNING!!! Frequencies sum to " + TextTools::toString(sf) + ", frequencies have been scaled.");
    sf *= 1./sf;
  }

  //Now build diagonal of the exchangeability matrix:
  for (unsigned int i = 0; i < 20; i++)
  {
    double sum = 0;
    for(unsigned int j = 0; j < 20; j++)
    {
      if(j!=i) sum += exchangeability_(i,j);
    }
    exchangeability_(i,i) = -sum;
  }

  //Closing stream:
  in.close();
}

/******************************************************************************/

void UserProteinSubstitutionModel::setFreqFromData(const SequenceContainer& data, double pseudoCount)
{
  map<int, int> counts;
  SequenceContainerTools::getCounts(data, counts);
  double t = 0;
  for (int i = 0; i < static_cast<int>(size_); i++)
  {
    t += (counts[i] + pseudoCount);
  }
  for (size_t i = 0; i < size_; ++i) freq_[i] = (static_cast<double>(counts[static_cast<int>(i)]) + pseudoCount) / t;
  freqSet_->setFrequencies(freq_);
  //Update parameters and re-compute generator and eigen values:
  matchParametersValues(freqSet_->getParameters());
}

/******************************************************************************/

