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

using namespace bpp;
using namespace std;

// ////////////////////////////
// FullCodonFrequenciesSet

FullCodonFrequenciesSet::FullCodonFrequenciesSet(const CodonAlphabet* alphabet, bool allowNullFreqs, const string& name) :
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Full.", name)
{
  unsigned int size = alphabet->getSize() - alphabet->numberOfStopCodons();
  unsigned int j = 0;

  for (unsigned int i = 0; i < alphabet->getSize() - 1; i++)
  {
    if (alphabet->isStop(i))
    {
      getFreq_(i) = 0;
    }
    else
    {
      addParameter_(new Parameter(
        "Full.theta" + TextTools::toString(i + 1),
        1. / (size - j),
        allowNullFreqs ?
        &Parameter::PROP_CONSTRAINT_IN :
        &FrequenciesSet::FREQUENCE_CONSTRAINT_MILLI));
      getFreq_(i) = 1. / size;
      j++;
    }
  }
  unsigned int i = alphabet->getSize() - 1;
  getFreq_(i) = (alphabet->isStop(i)) ? 0 : 1. / size;
}


FullCodonFrequenciesSet::FullCodonFrequenciesSet(const CodonAlphabet* alphabet, const vector<double>& initFreqs, bool allowNullFreqs, const string& name) :
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Full.", name)
{
  if (initFreqs.size() != alphabet->getSize())
    throw Exception("FullCodonFrequenciesSet(constructor). There must be " + TextTools::toString(alphabet->getSize()) + " frequencies.");
  double sum = 0.0;

  for (unsigned int i = 0; i < initFreqs.size(); i++)
  {
    if (!alphabet->isStop(i))
    {
      sum += initFreqs[i];
    }
  }

  double y = 1;
  for (unsigned int i = 0; i < alphabet->getSize() - 1; i++)
  {
    if (alphabet->isStop(i))
    {
      getFreq_(i) = 0;
    }
    else
    {
      addParameter_(new Parameter(
        "Full.theta" + TextTools::toString(i + 1),
        initFreqs[i] / sum / y,
        allowNullFreqs ?
        &Parameter::PROP_CONSTRAINT_IN :
        &FrequenciesSet::FREQUENCE_CONSTRAINT_MILLI));
      getFreq_(i) = initFreqs[i] / sum;
      y -= initFreqs[i] / sum;
    }
  }
  unsigned int i = alphabet->getSize() - 1;
  getFreq_(i) = (alphabet->isStop(i)) ? 0 : initFreqs[i] / sum;
}

void FullCodonFrequenciesSet::setFrequencies(const vector<double>& frequencies)
{
  if (frequencies.size() != getAlphabet()->getSize()) throw DimensionException("FullFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());
  const CodonAlphabet* alphabet = getAlphabet();

  double sum = 0.0;
  unsigned int i;
  for (i = 0; i < frequencies.size(); i++)
  {
    if (!(alphabet->isStop(i)))
      sum += frequencies[i];
  }

  double y = 1;
  for (i = 0; i < alphabet->getSize() - 1; i++)
  {
    if (alphabet->isStop(i))
    {
      getFreq_(i) = 0;
    }
    else
    {
      getParameter_("theta" + TextTools::toString(i + 1)).setValue(frequencies[i] / sum / y);
      y -= frequencies[i] / sum;
      getFreq_(i) = frequencies[i] / sum;
    }
  }
  i = alphabet->getSize() - 1;
  getFreq_(i) = (alphabet->isStop(i)) ? 0 : frequencies[i] / sum;
}

void FullCodonFrequenciesSet::fireParameterChanged(const ParameterList& parameters)
{
  const CodonAlphabet* alphabet = getAlphabet();
  double y = 1;
  unsigned int i;
  for (i = 0; i < alphabet->getSize() - 1; i++)
  {
    if (!(alphabet->isStop(i)))
    {
      getFreq_(i) = getParameter_("theta" + TextTools::toString(i + 1)).getValue() * y;
      y *= 1 - getParameter_("theta" + TextTools::toString(i + 1)).getValue();
    }
  }

  i = alphabet->getSize() - 1;
  getFreq_(i) = (alphabet->isStop(i)) ? 0 : y;
}


// ////////////////////////////
// FullPerAACodonFrequenciesSet

FullPerAACodonFrequenciesSet::FullPerAACodonFrequenciesSet(const GeneticCode* gencode,
                                                           const ProteinFrequenciesSet* ppfs):
  AbstractFrequenciesSet(gencode->getSourceAlphabet()->getSize(), gencode->getSourceAlphabet(), "FullPerAA.", "FullPerAA"),
  pgc_(gencode),
  ppfs_(ppfs->clone()),
  vS_()
{
  const ProteicAlphabet* ppa=dynamic_cast<const ProteicAlphabet*>(pgc_->getTargetAlphabet());
  
  for (unsigned int i = 0; i < ppa->getSize(); i++)
    {
      vector<int> vc=pgc_->getSynonymous(i);
      vS_.push_back(Simplex(vc.size(),1, ""));

      Simplex& si=vS_[i];
      si.setNamespace("FullPerAA."+ppa->intToChar(i)+"_");
      addParameters_(si.getParameters());
    }

  ppfs_->setNamespace("FullPerAA.");

  addParameters_(ppfs_->getParameters());

  updateFrequencies();
}

FullPerAACodonFrequenciesSet::FullPerAACodonFrequenciesSet(const GeneticCode* gencode) :
  AbstractFrequenciesSet(gencode->getSourceAlphabet()->getSize(), gencode->getSourceAlphabet(), "FullPerAA.", "FullPerAA"),
  pgc_(gencode),
  ppfs_(new FixedProteinFrequenciesSet(dynamic_cast<const ProteicAlphabet*>(gencode->getTargetAlphabet()), "FullPerAA.")),
  vS_()
{
  const ProteicAlphabet* ppa=dynamic_cast<const ProteicAlphabet*>(pgc_->getTargetAlphabet());
  
  for (unsigned int i = 0; i < ppa->getSize(); i++)
    {
      vector<int> vc=pgc_->getSynonymous(i);
      vS_.push_back(Simplex(vc.size(),1, ""));

      Simplex& si=vS_[i];
      si.setNamespace("FullPerAA."+ppa->intToChar(i)+"_");
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
{}

FullPerAACodonFrequenciesSet::~FullPerAACodonFrequenciesSet()
{
  if (ppfs_)
    delete ppfs_;
  ppfs_=0;
}

FullPerAACodonFrequenciesSet& FullPerAACodonFrequenciesSet::operator=(const FullPerAACodonFrequenciesSet& ffs)
{
  if (ppfs_)
    delete ppfs_;
  
  CodonFrequenciesSet::operator=(ffs);
  AbstractFrequenciesSet::operator=(ffs);
  pgc_=ffs.pgc_;
  ppfs_=ffs.ppfs_->clone();
  vS_=ffs.vS_;

  return *this;
}

void FullPerAACodonFrequenciesSet::fireParameterChanged(const ParameterList& parameters)
{
  if (dynamic_cast<AbstractFrequenciesSet*>(ppfs_))
    (dynamic_cast<AbstractFrequenciesSet*>(ppfs_))->matchParametersValues(parameters);
  for (unsigned int i=0;i<vS_.size();i++)
    vS_[i].matchParametersValues(parameters);
  updateFrequencies();
}

void FullPerAACodonFrequenciesSet::updateFrequencies()
{
  const ProteicAlphabet* ppa = dynamic_cast<const ProteicAlphabet*>(pgc_->getTargetAlphabet());
  
  for (unsigned int i = 0; i < ppa->getSize(); i++)
    {
      std::vector<int> vc= pgc_->getSynonymous(i);
      for (unsigned int j=0;j<vc.size();j++)
        getFreq_(vc[j]) = (ppfs_->getFrequencies())[i]*vS_[i].prob(j);     
    }
}

void FullPerAACodonFrequenciesSet::setFrequencies(const vector<double>& frequencies)
{
  if (frequencies.size() != getAlphabet()->getSize()) throw DimensionException("FullParAAFrequenciesSet::setFrequencies", frequencies.size(), getAlphabet()->getSize());

  const ProteicAlphabet* ppa = dynamic_cast<const ProteicAlphabet*>(pgc_->getTargetAlphabet());

  vector<double> vaa;
  for (unsigned int i=0;i<ppa->getSize();i++){
    vector<double> vp;
    double s=0;
    std::vector<int> vc= pgc_->getSynonymous(i);
    for (unsigned int j=0;j<vc.size();j++){
      vp.push_back(frequencies[vc[j]]);
      s+=frequencies[vc[j]];
    }
    vaa.push_back(s);
    vp/=s;
    vS_[i].setFrequencies(vp);
    matchParametersValues(vS_[i].getParameters());
  }

  ppfs_->setFrequencies(vaa);
  matchParametersValues(ppfs_->getParameters());
  updateFrequencies();
}

void FullPerAACodonFrequenciesSet::setNamespace(const std::string& prefix)
{
  const ProteicAlphabet* ppa = dynamic_cast<const ProteicAlphabet*>(pgc_->getTargetAlphabet());

  AbstractFrequenciesSet::setNamespace(prefix);
  ppfs_->setNamespace(prefix);
  for (unsigned int i=0;i<vS_.size();i++)
    vS_[i].setNamespace(prefix+ppa->intToChar(i)+"_");
}


// ///////////////////////////////////////////
// / FixedCodonFrequenciesSet

FixedCodonFrequenciesSet::FixedCodonFrequenciesSet(const CodonAlphabet* alphabet, const vector<double>& initFreqs, const string& name) :
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Fixed.", name)
{
  setFrequencies(initFreqs);
}

FixedCodonFrequenciesSet::FixedCodonFrequenciesSet(const CodonAlphabet* alphabet, const string& name) :
  AbstractFrequenciesSet(alphabet->getSize(), alphabet, "Fixed.", name)
{
  unsigned int size = alphabet->getSize() - alphabet->numberOfStopCodons();

  for (unsigned int i = 0; i < alphabet->getSize(); i++)
  {
    getFreq_(i) = (alphabet->isStop(i)) ? 0 : 1. / size;
  }
}

void FixedCodonFrequenciesSet::setFrequencies(const vector<double>& frequencies)
{
  const CodonAlphabet* ca = dynamic_cast<const CodonAlphabet*>(getAlphabet());
  if (frequencies.size() != ca->getSize()) throw DimensionException("FixedFrequenciesSet::setFrequencies", frequencies.size(), ca->getSize());
  double sum = 0.0;
  unsigned int i;

  for (i = 0; i < frequencies.size(); i++)
  {
    if (!(ca->isStop(i)))
      sum += frequencies[i];
  }

  for (i = 0; i < ca->getSize(); i++)
  {
    getFreq_(i) = (ca->isStop(i)) ? 0 : frequencies[i] / sum;
  }
}


// ///////////////////////////////////////////////////////////////////
// // CodonFromIndependentFrequenciesSet


CodonFromIndependentFrequenciesSet::CodonFromIndependentFrequenciesSet(
                                                                       const CodonAlphabet* pCA,
                                                                       const std::vector<FrequenciesSet*>& freqvector,
                                                                       const string& name) :
  WordFromIndependentFrequenciesSet(pCA, freqvector, "", name)
{
}

const CodonAlphabet* CodonFromIndependentFrequenciesSet::getAlphabet() const
{
  return dynamic_cast<const CodonAlphabet*>(WordFromIndependentFrequenciesSet::getAlphabet());
}

CodonFromIndependentFrequenciesSet::CodonFromIndependentFrequenciesSet(const CodonFromIndependentFrequenciesSet& iwfs) :
  WordFromIndependentFrequenciesSet(iwfs)
{
}

CodonFromIndependentFrequenciesSet& CodonFromIndependentFrequenciesSet::operator=(const CodonFromIndependentFrequenciesSet& iwfs)
{
  WordFromIndependentFrequenciesSet::operator=(iwfs);
  return *this;
}

void CodonFromIndependentFrequenciesSet::updateFrequencies()
{
  WordFromIndependentFrequenciesSet::updateFrequencies();
  
  unsigned int s = getAlphabet()->getSize();
  double sum=0;
  for (unsigned int i = 0; i < s; i++)
    {
      if (getAlphabet()->isStop(i))
        getFreq_(i) = 0;
      else
        sum+=getFreq_(i);
    }
  
  for (unsigned int i = 0; i < s; i++)
    getFreq_(i)=getFreq_(i)/sum;
}

void CodonFromIndependentFrequenciesSet::setFrequencies(const vector<double>& frequencies) 
{
  unsigned int s = getAlphabet()->getSize();
  double sum=0;
  vector<double> freq;
  for (unsigned int i = 0; i < s; i++)
    if (!getAlphabet()->isStop(i))
      sum+=frequencies[i];

  for (unsigned int i = 0; i < s; i++)
    if (getAlphabet()->isStop(i))
      freq.push_back(0);
    else
      freq.push_back(frequencies[i]/sum);

  WordFromIndependentFrequenciesSet::setFrequencies(freq);
}


// ///////////////////////////////////////////////////////////////////
// // CodonFromUniqueFrequenciesSet


CodonFromUniqueFrequenciesSet::CodonFromUniqueFrequenciesSet(const CodonAlphabet* pCA, FrequenciesSet* pfreq, const string& name) :
  WordFromUniqueFrequenciesSet(pCA, pfreq, "", name)
{
}

const CodonAlphabet* CodonFromUniqueFrequenciesSet::getAlphabet() const
{
  return dynamic_cast<const CodonAlphabet*>(WordFromUniqueFrequenciesSet::getAlphabet());
}


CodonFromUniqueFrequenciesSet::CodonFromUniqueFrequenciesSet(const CodonFromUniqueFrequenciesSet& iwfs) :
  WordFromUniqueFrequenciesSet(iwfs)
{
}

CodonFromUniqueFrequenciesSet& CodonFromUniqueFrequenciesSet::operator=(const CodonFromUniqueFrequenciesSet& iwfs)
{
  WordFromUniqueFrequenciesSet::operator=(iwfs);
  return *this;
}

void CodonFromUniqueFrequenciesSet::updateFrequencies()
{
  WordFromUniqueFrequenciesSet::updateFrequencies();
  const CodonAlphabet* pCA=dynamic_cast<const CodonAlphabet*>(getAlphabet());
  
  unsigned int s = getAlphabet()->getSize();
  double sum=0;
  for (unsigned int i = 0; i < s; i++)
    {
      if (pCA->isStop(i))
        getFreq_(i) = 0;
      else
        sum+=getFreq_(i);
    }
  
  for (unsigned int i = 0; i < s; i++){
    getFreq_(i)=getFreq_(i)/sum;
  }
}

void CodonFromUniqueFrequenciesSet::setFrequencies(const vector<double>& frequencies) 
{
  const CodonAlphabet* pCA=dynamic_cast<const CodonAlphabet*>(getAlphabet());
  
  unsigned int s = getAlphabet()->getSize();
  double sum=0;
  vector<double> freq;
  for (unsigned int i = 0; i < s; i++)
    if (!pCA->isStop(i))
      sum+=frequencies[i];

  for (unsigned int i = 0; i < s; i++)
    if (pCA->isStop(i))
      freq.push_back(0);
    else
      freq.push_back(frequencies[i]/sum);

  WordFromUniqueFrequenciesSet::setFrequencies(freq);
}


/*********************************************************************/

FrequenciesSet* CodonFrequenciesSet::getFrequenciesSetForCodons(short option, const CodonAlphabet& CA)
{
  FrequenciesSet* codonFreqs;

  if (option == F0)
    codonFreqs = new FixedCodonFrequenciesSet(&CA, "F0");
  else if (option == F1X4)
    codonFreqs = new CodonFromUniqueFrequenciesSet(&CA, new FullNucleotideFrequenciesSet(CA.getNucleicAlphabet()), "F1X4");
  else if (option == F3X4)
  {
    vector<FrequenciesSet*> v_AFS(3);
    v_AFS[0] = new FullNucleotideFrequenciesSet(CA.getNucleicAlphabet());
    v_AFS[1] = new FullNucleotideFrequenciesSet(CA.getNucleicAlphabet());
    v_AFS[2] = new FullNucleotideFrequenciesSet(CA.getNucleicAlphabet());
    codonFreqs = new CodonFromIndependentFrequenciesSet(&CA,v_AFS, "F3X4");
  }
  else if (option == F61)
    codonFreqs = new FullCodonFrequenciesSet(&CA, "F61");
  else throw Exception("FrequenciesSet::getFrequencySetForCodons(). Unvalid codon frequency set argument.");

  return codonFreqs;
}

/******************************************************************************/

const short CodonFrequenciesSet::F0   = 0;
const short CodonFrequenciesSet::F1X4 = 1;
const short CodonFrequenciesSet::F3X4 = 2;
const short CodonFrequenciesSet::F61  = 3;

/******************************************************************************/

