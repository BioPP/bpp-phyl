// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/TextTools.h>

#include "UserProteinSubstitutionModel.h"

// From SeqLib:
#include <Bpp/Seq/Container/SequenceContainerTools.h>

using namespace bpp;

// From the STL:
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

/******************************************************************************/

UserProteinSubstitutionModel::UserProteinSubstitutionModel(
    std::shared_ptr<const ProteicAlphabet> alpha,
    const string& path,
    const string& prefix) :
  AbstractParameterAliasable(prefix),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), prefix),
  path_(path),
  freqSet_(nullptr)
{
  readFromFile();
  freqSet_.reset(new FixedProteinFrequencySet(alpha, freq_));
  updateMatrices_();
}

UserProteinSubstitutionModel::UserProteinSubstitutionModel(
    std::shared_ptr<const ProteicAlphabet> alpha,
    const string& path,
    unique_ptr<ProteinFrequencySetInterface> freqSet,
    const string& prefix,
    bool initFreqs) :
  AbstractParameterAliasable(prefix),
  AbstractReversibleProteinSubstitutionModel(alpha, make_shared<CanonicalStateMap>(alpha, false), prefix),
  path_(path),
  freqSet_(std::move(freqSet))
{
  readFromFile();
  freqSet_->setNamespace(prefix + freqSet_->getName() + ".");
  if (initFreqs)
    freqSet->setFrequencies(freq_);
  else
    freq_ = freqSet_->getFrequencies();
  addParameters_(freqSet_->getParameters());
  updateMatrices_();
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
  if (!FileTools::fileExists(path_.c_str()))
    throw Exception("UserProteinSubstitutionModel::readFromFile. Frequencies file not found : " +  path_);

  ifstream in(path_.c_str(), ios::in);
  // Read exchangeability matrix:
  for (unsigned int i = 1; i < 20; i++)
  {
    string line = FileTools::getNextLine(in);
    StringTokenizer st(line);
    for (unsigned int j = 0; j < i; j++)
    {
      double s = TextTools::toDouble(st.nextToken());
      exchangeability_(i, j) = exchangeability_(j, i) = s;
    }
  }
  // Read frequencies:
  unsigned int fCount = 0;
  while (in && fCount < 20)
  {
    string line = FileTools::getNextLine(in);
    StringTokenizer st(line);
    while (st.hasMoreToken() && fCount < 20)
    {
      freq_[fCount] = TextTools::toDouble(st.nextToken());
      fCount++;
    }
  }
  double sf = VectorTools::sum(freq_);
  if (sf - 1 > 0.000001)
  {
    ApplicationTools::displayMessage("WARNING!!! Frequencies sum to " + TextTools::toString(sf) + ", frequencies have been scaled.");
    sf *= 1. / sf;
  }

  // Now build diagonal of the exchangeability matrix:
  for (unsigned int i = 0; i < 20; i++)
  {
    double sum = 0;
    for (unsigned int j = 0; j < 20; j++)
    {
      if (j != i)
        sum += exchangeability_(i, j);
    }
    exchangeability_(i, i) = -sum;
  }

  // Closing stream:
  in.close();
}

/******************************************************************************/

void UserProteinSubstitutionModel::setFreqFromData(const SequenceDataInterface& data, double pseudoCount)
{
  map<int, double> counts;
  SequenceContainerTools::getFrequencies(data, counts, pseudoCount);
  for (auto i : counts)
  {
    freq_[(size_t)i.first] = i.second;
  }

  freqSet_->setFrequencies(freq_);
  // Update parameters and re-compute generator and eigen values:
  matchParametersValues(freqSet_->getParameters());
}

/******************************************************************************/
