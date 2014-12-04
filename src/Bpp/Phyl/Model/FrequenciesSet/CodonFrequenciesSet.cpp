//
// File: CodonFrequenciesSet.cpp
// Created by: Laurent Gueguen
// Created on: lundi 2 avril 2012, à 14h 15
//

/*
   Copyright or (c) or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "CodonFrequenciesSet.h"
#include "NucleotideFrequenciesSet.h"
#include "../StateMap.h"

// From bpp-core:
#include <Bpp/Numeric/Prob/Simplex.h>

using namespace bpp;
using namespace std;

// ////////////////////////////
// FullCodonFrequenciesSet

FullCodonFrequenciesSet::FullCodonFrequenciesSet(const GeneticCode* gCode, bool allowNullFreqs, unsigned short method, const string& name) :
  AbstractFrequenciesSet(
    new CanonicalStateMap(gCode->getSourceAlphabet(), false),
    "Full.",
    name),
  pgc_(gCode),
  sFreq_(gCode->getSourceAlphabet()->getSize()  - gCode->getNumberOfStopCodons(), method, allowNullFreqs, "Full.")
{
  vector<double> vd;
  double r = 1. / static_cast<double>(sFreq_.dimension());

  for (size_t i = 0; i < sFreq_.dimension(); i++)
  {
    vd.push_back(r);
  }

  sFreq_.setFrequencies(vd);
  addParameters_(sFreq_.getParameters());
  updateFreq_();
}

FullCodonFrequenciesSet::FullCodonFrequenciesSet(
  const GeneticCode* gCode,
  const vector<double>& initFreqs,
  bool allowNullFreqs,
  unsigned short method,
  const string& name) :
  AbstractFrequenciesSet(new CanonicalStateMap(gCode->getSourceAlphabet(), false), "Full.", name),
  pgc_(gCode),
  sFreq_(gCode->getSourceAlphabet()->getSize() - gCode->getNumberOfStopCodons(), method, allowNullFreqs, "Full.")
{
  if (initFreqs.size() != getAlphabet()->getSize())
    throw Exception("FullCodonFrequenciesSet(constructor). There must be " + TextTools::toString(gCode->getSourceAlphabet()->getSize()) + " frequencies.");

  double sum = 0;
  for (size_t i = 0; i < getAlphabet()->getSize(); i++)
  {
    if (!pgc_->isStop(static_cast<int>(i)))
      sum += initFreqs[i];
  }

  vector<double> vd;
  for (size_t i = 0; i < getAlphabet()->getSize(); i++)
  {
    if (!gCode->isStop(static_cast<int>(i)))
      vd.push_back(initFreqs[i] / sum);
  }

  sFreq_.setFrequencies(vd);
  addParameters_(sFreq_.getParameters());
  updateFreq_();
}

FullCodonFrequenciesSet::FullCodonFrequenciesSet(const FullCodonFrequenciesSet& fcfs) :
  AbstractFrequenciesSet(fcfs),
  pgc_(fcfs.pgc_),
  sFreq_(fcfs.sFreq_)
{}

FullCodonFrequenciesSet& FullCodonFrequenciesSet::operator=(const FullCodonFrequenciesSet& fcfs)
{
  AbstractFrequenciesSet::operator=(fcfs);
  pgc_ = fcfs.pgc_;
  sFreq_ = fcfs.sFreq_;

  return *this;
}

void FullCodonFrequenciesSet::setNamespace(const std::string& nameSpace)
{
  sFreq_.setNamespace(nameSpace);
  AbstractFrequenciesSet::setNamespace(nameSpace);
}

void FullCodonFrequenciesSet::setFrequencies(const vector<double>& frequencies)
{
  if (frequencies.size() != getAlphabet()->getSize())
    throw DimensionException("FullFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());

  double sum = 0;
  for (size_t i = 0; i < getAlphabet()->getSize(); i++)
  {
    if (!pgc_->isStop(static_cast<int>(i)))
      sum += frequencies[i];
  }

  vector<double> vd;
  for (size_t i = 0; i < getAlphabet()->getSize(); i++)
  {
    if (!pgc_->isStop(static_cast<int>(i)))
      vd.push_back(frequencies[i] / sum);
  }

  sFreq_.setFrequencies(vd);
  setParametersValues(sFreq_.getParameters());
  updateFreq_();
}

void FullCodonFrequenciesSet::fireParameterChanged(const ParameterList& parameters)
{
  AbstractFrequenciesSet::fireParameterChanged(parameters);
  
  sFreq_.matchParametersValues(parameters);
  updateFreq_();
}

void FullCodonFrequenciesSet::updateFreq_()
{
  size_t nbstop = 0;

  for (size_t j = 0; j < getAlphabet()->getSize(); j++)
  {
    if (pgc_->isStop(static_cast<int>(j)))
    {
      getFreq_(j) = 0;
      nbstop++;
    }
    else
      getFreq_(j) = sFreq_.prob(j - nbstop);
  }
}


// ////////////////////////////
// FullPerAACodonFrequenciesSet

FullPerAACodonFrequenciesSet::FullPerAACodonFrequenciesSet(
  const GeneticCode* gencode,
  ProteinFrequenciesSet* ppfs,
  unsigned short method) :
  AbstractFrequenciesSet(new CanonicalStateMap(gencode->getSourceAlphabet(), false), "FullPerAA.", "FullPerAA"),
  pgc_(gencode),
  ppfs_(ppfs),
  vS_()
{
  const ProteicAlphabet* ppa = dynamic_cast<const ProteicAlphabet*>(pgc_->getTargetAlphabet());
  const StateMap& aaStates = ppfs_->getStateMap();
  for (size_t i = 0; i < aaStates.getNumberOfModelStates(); ++i)
  {
    int aa = aaStates.getAlphabetStateAsInt(i);
    vector<int> vc = pgc_->getSynonymous(aa);
    vS_.push_back(Simplex(vc.size(), method, 0, ""));

    Simplex& si = vS_[i];

    si.setNamespace("FullPerAA." + ppa->getAbbr(aa) + "_");
    addParameters_(si.getParameters());
  }

  ppfs_->setNamespace("FullPerAA." + ppfs_->getName() + ".");
  addParameters_(ppfs_->getParameters());

  updateFrequencies();
}

FullPerAACodonFrequenciesSet::FullPerAACodonFrequenciesSet(const GeneticCode* gencode, unsigned short method) :
  AbstractFrequenciesSet(new CanonicalStateMap(gencode->getSourceAlphabet(), false), "FullPerAA.", "FullPerAA"),
  pgc_(gencode),
  ppfs_(new FixedProteinFrequenciesSet(dynamic_cast<const ProteicAlphabet*>(gencode->getTargetAlphabet()), "FullPerAA.")),
  vS_()
{
  const ProteicAlphabet* ppa = dynamic_cast<const ProteicAlphabet*>(pgc_->getTargetAlphabet());
  const StateMap& aaStates = ppfs_->getStateMap();
  for (size_t i = 0; i < aaStates.getNumberOfModelStates(); ++i)
  {
    int aa = aaStates.getAlphabetStateAsInt(i);
    vector<int> vc = pgc_->getSynonymous(aa);
    vS_.push_back(Simplex(vc.size(), method, 0, ""));

    Simplex& si = vS_[i];
    si.setNamespace("FullPerAA." + ppa->getAbbr(aa) + "_");
    addParameters_(si.getParameters());
  }

  updateFrequencies();
}

FullPerAACodonFrequenciesSet::FullPerAACodonFrequenciesSet(const FullPerAACodonFrequenciesSet& ffs) :
  CodonFrequenciesSet(ffs),
  AbstractFrequenciesSet(ffs),
  pgc_(ffs.pgc_),
  ppfs_(ffs.ppfs_->clone()),
  vS_(ffs.vS_)
{
  updateFrequencies();
}

FullPerAACodonFrequenciesSet& FullPerAACodonFrequenciesSet::operator=(const FullPerAACodonFrequenciesSet& ffs)
{
  CodonFrequenciesSet::operator=(ffs);
  AbstractFrequenciesSet::operator=(ffs);
  pgc_ = ffs.pgc_;
  ppfs_.reset(ffs.ppfs_->clone());
  vS_ = ffs.vS_;

  return *this;
}

void FullPerAACodonFrequenciesSet::fireParameterChanged(const ParameterList& parameters)
{
  AbstractFrequenciesSet::fireParameterChanged(parameters);

  if (dynamic_cast<AbstractFrequenciesSet*>(ppfs_.get()))
    (dynamic_cast<AbstractFrequenciesSet*>(ppfs_.get()))->matchParametersValues(parameters);
  for (size_t i = 0; i < vS_.size(); i++)
  {
    vS_[i].matchParametersValues(parameters);
  }
  updateFrequencies();
}

void FullPerAACodonFrequenciesSet::updateFrequencies()
{
  const StateMap& aaStates = ppfs_->getStateMap();
  for (size_t i = 0; i < aaStates.getNumberOfModelStates(); ++i)
  {
    int aa = aaStates.getAlphabetStateAsInt(i);
    std::vector<int> vc = pgc_->getSynonymous(aa);
    for (size_t j = 0; j < vc.size(); j++)
    {
      //NB: only one alphabet state per model state here, as it is a CodonFreqSet.
      getFreq_(getStateMap().getModelStates(vc[j])[0]) =
        static_cast<double>(vc.size()) * ppfs_->getFrequencies()[i] * vS_[i].prob(j);
    }
  }
  normalize();
}

void FullPerAACodonFrequenciesSet::setFrequencies(const vector<double>& frequencies)
{
  if (frequencies.size() != getAlphabet()->getSize())
    throw DimensionException("FullPerAAFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());

  double bigS = 0;
  Vdouble vaa;
  
  const StateMap& aaStates = ppfs_->getStateMap();
  for (size_t i = 0; i < aaStates.getNumberOfModelStates(); ++i)
  {
    int aa = aaStates.getAlphabetStateAsInt(i);
    std::vector<int> vc = pgc_->getSynonymous(aa);
    Vdouble vp;
    double s = 0;
    
    for (size_t j = 0; j < vc.size(); j++)
    {
      size_t index = pgc_->getSourceAlphabet()->getStateIndex(vc[j]);
      vp.push_back(frequencies[index-1]);
      s += frequencies[index-1];
    }
      
    vp /= s;
    vS_[i].setFrequencies(vp);

    matchParametersValues(vS_[i].getParameters());

    bigS += s / static_cast<double>(vc.size());
    vaa.push_back(s / static_cast<double>(vc.size()));
  }

  vaa /= bigS; // to avoid counting of stop codons
  ppfs_->setFrequencies(vaa);
  matchParametersValues(ppfs_->getParameters());
  updateFrequencies();
}

void FullPerAACodonFrequenciesSet::setNamespace(const std::string& prefix)
{
  const ProteicAlphabet* ppa = dynamic_cast<const ProteicAlphabet*>(pgc_->getTargetAlphabet());

  AbstractFrequenciesSet::setNamespace(prefix);
  ppfs_->setNamespace(prefix + ppfs_->getName() + ".");
  for (size_t i = 0; i < vS_.size(); i++)
  {
    vS_[i].setNamespace(prefix + ppa->getAbbr(static_cast<int>(i)) + "_");
  }
}


// ///////////////////////////////////////////
// / FixedCodonFrequenciesSet

FixedCodonFrequenciesSet::FixedCodonFrequenciesSet(
  const GeneticCode* gCode,
  const vector<double>& initFreqs,
  const string& name) :
  AbstractFrequenciesSet(new CanonicalStateMap(gCode->getSourceAlphabet(), false), "Fixed.", name),
  pgc_(gCode)
{
  setFrequencies(initFreqs);
}

FixedCodonFrequenciesSet::FixedCodonFrequenciesSet(const GeneticCode* gCode, const string& name) :
  AbstractFrequenciesSet(new CanonicalStateMap(gCode->getSourceAlphabet(), false), "Fixed.", name),
  pgc_(gCode)
{
  size_t size = gCode->getSourceAlphabet()->getSize() - gCode->getNumberOfStopCodons();

  for (size_t i = 0; i < gCode->getSourceAlphabet()->getSize(); i++)
  {
    getFreq_(i) = (gCode->isStop(static_cast<int>(i))) ? 0 : 1. / static_cast<double>(size);
  }
}

void FixedCodonFrequenciesSet::setFrequencies(const vector<double>& frequencies)
{
  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(getAlphabet());
  if (frequencies.size() != ca->getSize())
    throw DimensionException("FixedFrequenciesSet::setFrequencies", frequencies.size(), ca->getSize());
  double sum = 0.0;

  for (size_t i = 0; i < frequencies.size(); i++)
  {
    if (!(pgc_->isStop(static_cast<int>(i))))
      sum += frequencies[i];
  }

  for (size_t i = 0; i < ca->getSize(); i++)
  {
    getFreq_(i) = (pgc_->isStop(static_cast<int>(i))) ? 0 : frequencies[i] / sum;
  }
}


// ///////////////////////////////////////////////////////////////////
// // CodonFromIndependentFrequenciesSet


CodonFromIndependentFrequenciesSet::CodonFromIndependentFrequenciesSet(
  const GeneticCode* gCode,
  const std::vector<FrequenciesSet*>& freqvector,
  const string& name,
  const string& mgmtStopCodon) :
  WordFromIndependentFrequenciesSet(gCode->getSourceAlphabet(), freqvector, "", name),
  mStopNeigh_(),
  mgmtStopCodon_(2),
  pgc_(gCode)
{
  if (mgmtStopCodon == "uniform")
    mgmtStopCodon_ = 0;
  else if (mgmtStopCodon == "linear")
    mgmtStopCodon_ = 1;

  // fill the map of the stop codons

  vector<int> vspcod = gCode->getStopCodonsAsInt();
  for (size_t ispcod = 0; ispcod < vspcod.size(); ispcod++)
  {
    size_t pow = 1;
    int nspcod = vspcod[ispcod];
    for (size_t ph = 0; ph < 3; ph++)
    {
      size_t nspcod0 = static_cast<size_t>(nspcod) - pow * static_cast<size_t>(getAlphabet()->getNPosition(nspcod, 2 - ph));
      for (size_t dec = 0; dec < 4; dec++)
      {
        size_t vois = nspcod0 + pow * dec;
        if (!pgc_->isStop(static_cast<int>(vois)))
          mStopNeigh_[nspcod].push_back(static_cast<int>(vois));
      }
      pow *= 4;
    }
  }
  updateFrequencies();
}

const CodonAlphabet* CodonFromIndependentFrequenciesSet::getAlphabet() const
{
  return dynamic_cast<const CodonAlphabet*>(WordFromIndependentFrequenciesSet::getAlphabet());
}

CodonFromIndependentFrequenciesSet::CodonFromIndependentFrequenciesSet(const CodonFromIndependentFrequenciesSet& iwfs) :
  WordFromIndependentFrequenciesSet(iwfs),
  mStopNeigh_(iwfs.mStopNeigh_),
  mgmtStopCodon_(iwfs.mgmtStopCodon_),
  pgc_(iwfs.pgc_)
{
  updateFrequencies();
}

CodonFromIndependentFrequenciesSet& CodonFromIndependentFrequenciesSet::operator=(const CodonFromIndependentFrequenciesSet& iwfs)
{
  WordFromIndependentFrequenciesSet::operator=(iwfs);
  mStopNeigh_ = iwfs.mStopNeigh_;
  mgmtStopCodon_ = iwfs.mgmtStopCodon_;
  pgc_ = iwfs.pgc_;
  return *this;
}

void CodonFromIndependentFrequenciesSet::updateFrequencies()
{
  WordFromIndependentFrequenciesSet::updateFrequencies();

  size_t s = getAlphabet()->getSize();

  if (mgmtStopCodon_ != 0)
  {
    // The frequencies of the stop codons are distributed to all
    // neighbour non-stop codons
    double f[64];
    for (size_t i = 0; i < s; i++)
    {
      f[i] = 0;
    }

    std::map<int, Vint>::iterator mStopNeigh_it(mStopNeigh_.begin());
    while (mStopNeigh_it != mStopNeigh_.end())
    {
      int stNb = mStopNeigh_it->first;
      Vint vneigh = mStopNeigh_it->second;
      double sneifreq = 0;
      for (size_t vn = 0; vn < vneigh.size(); vn++)
      {
        sneifreq += pow(getFreq_(static_cast<size_t>(vneigh[vn])), mgmtStopCodon_);
      }
      double x = getFreq_(static_cast<size_t>(stNb)) / sneifreq;
      for (size_t vn = 0; vn < vneigh.size(); vn++)
      {
        f[vneigh[vn]] += pow(getFreq_(static_cast<size_t>(vneigh[vn])), mgmtStopCodon_) * x;
      }
      getFreq_(static_cast<size_t>(stNb)) = 0;
      mStopNeigh_it++;
    }

    for (size_t i = 0; i < s; i++)
    {
      getFreq_(i) += f[i];
    }
  }
  else
  {
    double sum = 0.;
    for (size_t i = 0; i < s; i++)
    {
      if (!pgc_->isStop(static_cast<int>(i)))
        sum += getFreq_(i);
    }

    for (size_t i = 0; i < s; i++)
    {
      if (pgc_->isStop(static_cast<int>(i)))
        getFreq_(i) = 0;
      else
        getFreq_(i) /= sum;
    }
  }
}

// ///////////////////////////////////////////////////////////////////
// // CodonFromUniqueFrequenciesSet


CodonFromUniqueFrequenciesSet::CodonFromUniqueFrequenciesSet(
  const GeneticCode* gCode,
  FrequenciesSet* pfreq,
  const string& name,
  const string& mgmtStopCodon) :
  WordFromUniqueFrequenciesSet(gCode->getSourceAlphabet(), pfreq, "", name),
  mStopNeigh_(),
  mgmtStopCodon_(2),
  pgc_(gCode)
{
  if (mgmtStopCodon == "uniform")
    mgmtStopCodon_ = 0;
  else if (mgmtStopCodon == "linear")
    mgmtStopCodon_ = 1;

  // fill the map of the stop codons

  vector<int> vspcod = gCode->getStopCodonsAsInt();
  for (size_t ispcod = 0; ispcod < vspcod.size(); ispcod++)
  {
    size_t pow = 1;
    int nspcod = vspcod[ispcod];
    for (int ph = 0; ph < 3; ph++)
    {
      size_t nspcod0 = static_cast<size_t>(nspcod) - pow * static_cast<size_t>(getAlphabet()->getNPosition(nspcod, static_cast<unsigned int>(2 - ph)));
      for (size_t dec = 0; dec < 4; dec++)
      {
        size_t vois = nspcod0 + pow * dec;
        if (!pgc_->isStop(static_cast<int>(vois)))
          mStopNeigh_[nspcod].push_back(static_cast<int>(vois));
      }
      pow *= 4;
    }
  }

  updateFrequencies();
}

const CodonAlphabet* CodonFromUniqueFrequenciesSet::getAlphabet() const
{
  return dynamic_cast<const CodonAlphabet*>(WordFromUniqueFrequenciesSet::getAlphabet());
}


CodonFromUniqueFrequenciesSet::CodonFromUniqueFrequenciesSet(const CodonFromUniqueFrequenciesSet& iwfs) :
  WordFromUniqueFrequenciesSet(iwfs),
  mStopNeigh_(iwfs.mStopNeigh_),
  mgmtStopCodon_(iwfs.mgmtStopCodon_),
  pgc_(iwfs.pgc_)
{
  updateFrequencies();
}

CodonFromUniqueFrequenciesSet& CodonFromUniqueFrequenciesSet::operator=(const CodonFromUniqueFrequenciesSet& iwfs)
{
  WordFromUniqueFrequenciesSet::operator=(iwfs);
  mStopNeigh_ = iwfs.mStopNeigh_;
  mgmtStopCodon_ = iwfs.mgmtStopCodon_;
  pgc_ = iwfs.pgc_;
  return *this;
}

void CodonFromUniqueFrequenciesSet::updateFrequencies()
{
  WordFromUniqueFrequenciesSet::updateFrequencies();

  size_t s = getAlphabet()->getSize();

  if (mgmtStopCodon_ != 0)
  {
    // The frequencies of the stop codons are distributed to all
    // neighbour non-stop codons
    double f[64];
    for (size_t i = 0; i < s; i++)
    {
      f[i] = 0;
    }

    std::map<int, Vint>::iterator mStopNeigh_it(mStopNeigh_.begin());
    while (mStopNeigh_it != mStopNeigh_.end())
    {
      int stNb = mStopNeigh_it->first;
      Vint vneigh = mStopNeigh_it->second;
      double sneifreq = 0;
      for (size_t vn = 0; vn < vneigh.size(); vn++)
      {
        sneifreq += pow(getFreq_(static_cast<size_t>(vneigh[vn])), mgmtStopCodon_);
      }
      double x = getFreq_(static_cast<size_t>(stNb)) / sneifreq;
      for (size_t vn = 0; vn < vneigh.size(); vn++)
      {
        f[vneigh[vn]] += pow(getFreq_(static_cast<size_t>(vneigh[vn])), mgmtStopCodon_) * x;
      }
      getFreq_(static_cast<size_t>(stNb)) = 0;
      mStopNeigh_it++;
    }

    for (size_t i = 0; i < s; i++)
    {
      getFreq_(i) += f[i];
    }
  }
  else
  {
    double sum = 0.;
    for (size_t i = 0; i < s; i++)
    {
      if (!pgc_->isStop(static_cast<int>(i)))
        sum += getFreq_(i);
    }

    for (unsigned int i = 0; i < s; i++)
    {
      if (pgc_->isStop(static_cast<int>(i)))
        getFreq_(i) = 0;
      else
        getFreq_(i) /= sum;
    }
  }
}

/*********************************************************************/

FrequenciesSet* CodonFrequenciesSet::getFrequenciesSetForCodons(short option, const GeneticCode* gCode, const string& mgmtStopCodon, unsigned short method)
{
  FrequenciesSet* codonFreqs;

  if (option == F0)
    codonFreqs = new FixedCodonFrequenciesSet(gCode, "F0");
  else if (option == F1X4)
    codonFreqs = new CodonFromUniqueFrequenciesSet(gCode, new FullNucleotideFrequenciesSet(gCode->getSourceAlphabet()->getNucleicAlphabet()), "F1X4", mgmtStopCodon);
  else if (option == F3X4)
  {
    vector<FrequenciesSet*> v_AFS(3);
    v_AFS[0] = new FullNucleotideFrequenciesSet(gCode->getSourceAlphabet()->getNucleicAlphabet());
    v_AFS[1] = new FullNucleotideFrequenciesSet(gCode->getSourceAlphabet()->getNucleicAlphabet());
    v_AFS[2] = new FullNucleotideFrequenciesSet(gCode->getSourceAlphabet()->getNucleicAlphabet());
    codonFreqs = new CodonFromIndependentFrequenciesSet(gCode, v_AFS, "F3X4", mgmtStopCodon);
  }
  else if (option == F61)
    codonFreqs = new FullCodonFrequenciesSet(gCode, false, method, "F61");
  else
    throw Exception("FrequenciesSet::getFrequencySetForCodons(). Unvalid codon frequency set argument.");

  return codonFreqs;
}

/******************************************************************************/

const short CodonFrequenciesSet::F0   = 0;
const short CodonFrequenciesSet::F1X4 = 1;
const short CodonFrequenciesSet::F3X4 = 2;
const short CodonFrequenciesSet::F61  = 3;

/******************************************************************************/
