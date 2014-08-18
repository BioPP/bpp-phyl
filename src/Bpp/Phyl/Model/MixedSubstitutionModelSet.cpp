//
// File: MixedSubstitutionModelSet.cpp
// Created by: Laurent Guéguen
// Created on: mercredi 25 mai 2011, à 22h 12
//

/*
   Copyright or <A9> or Copr. Bio++ Development Team, (November 16, 2004)

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

#include "MixedSubstitutionModelSet.h"
#include "MixedSubstitutionModel.h"
#include "MixtureOfASubstitutionModel.h"

using namespace bpp;
using namespace std;

MixedSubstitutionModelSet::MixedSubstitutionModelSet(const MixedSubstitutionModelSet& set) :
  SubstitutionModelSet(set),
  vpHyperNodes_()
{
  for (size_t i = 0; i < set.vpHyperNodes_.size(); i++)
  {
    vpHyperNodes_.push_back(new HyperNode(*set.vpHyperNodes_[i]));
  }
}

MixedSubstitutionModelSet::MixedSubstitutionModelSet(const SubstitutionModelSet& set) :
  SubstitutionModelSet(set),
  vpHyperNodes_()
{}

MixedSubstitutionModelSet::~MixedSubstitutionModelSet()
{
  for (size_t i = 0; i < vpHyperNodes_.size(); i++)
  {
    delete vpHyperNodes_[i];
  }
}

void MixedSubstitutionModelSet::clear()
{
  SubstitutionModelSet::clear();
  for (size_t i = 0; i < vpHyperNodes_.size(); i++)
  {
    delete vpHyperNodes_[i];
  }
}

MixedSubstitutionModelSet& MixedSubstitutionModelSet::operator=(const MixedSubstitutionModelSet& set)
{
  SubstitutionModelSet::operator=(set);
  for (size_t i = 0; i < vpHyperNodes_.size(); i++)
  {
    if (vpHyperNodes_[i] != NULL)
      delete vpHyperNodes_[i];
  }
  vpHyperNodes_.clear();

  for (size_t i = 0; i < set.vpHyperNodes_.size(); i++)
  {
    vpHyperNodes_.push_back(new HyperNode(*set.vpHyperNodes_[i]));
  }

  return *this;
}

void MixedSubstitutionModelSet::addEmptyHyperNode()
{
  vpHyperNodes_.push_back(new HyperNode(this));
}

void MixedSubstitutionModelSet::addHyperNode(const HyperNode& hn)
{
  vpHyperNodes_.push_back(new HyperNode(hn));
}

bool MixedSubstitutionModelSet::complete()
{
  MixedSubstitutionModelSet::HyperNode nhn(this);
  size_t i;
  for (i = 0; i < vpHyperNodes_.size(); i++)
  {
    nhn += *vpHyperNodes_[i];
  }

  size_t nbm = getNumberOfModels();
  for (i = 0; i < nbm; i++)
  {
    const MixedSubstitutionModel* pSM = dynamic_cast<const MixedSubstitutionModel*>(getModel(i));
    if (pSM)
    {
      if (nhn.getNode(i).size() != pSM->getNumberOfModels())
        break;
    }
  }

  if (i == nbm)
    return false;

  addEmptyHyperNode();
  for (i = 0; i < nbm; i++)
  {
    const MixedSubstitutionModel* pSM = dynamic_cast<const MixedSubstitutionModel*>(getModel(i));
    if (pSM)
    {
      const MixedSubstitutionModelSet::HyperNode::Node& nd = nhn.getNode(i);
      int snd = static_cast<int>(nd.size());
      int vs = static_cast<int>(pSM->getNumberOfModels());
      Vint an;

      int j(0), k(0);
      while (j < vs)
      {
        while ((k < snd) && (nd[static_cast<size_t>(k)] < j))
          k++;
        if ((k >= snd) || (nd[static_cast<size_t>(k)] > j))
          an.push_back(j);
        j++;
      }
      addToHyperNode(i, an);
    }
  }

  return true;
}

void MixedSubstitutionModelSet::addToHyperNode(size_t nM, const Vint& vnS, int nH)
{
  if (nH >= static_cast<int>(vpHyperNodes_.size()))
    throw BadIntegerException("MixedSubstitutionModelSet::addToHyperNode. Bad HyperNode number", nH);
  if (nH < 0)
    nH = static_cast<int>(vpHyperNodes_.size() - 1);

  if (nM >= getNumberOfModels())
    throw IndexOutOfBoundsException("MixedSubstitutionModelSet::addToHyperNode. Bad Mixed Model number", nM, 0, getNumberOfModels());

  vpHyperNodes_[static_cast<size_t>(nH)]->addToModel(nM, vnS);
}

bool MixedSubstitutionModelSet::hasExclusivePaths() const
{
  HyperNode tthn(this);

  size_t nhn = getNumberOfHyperNodes();
  for (size_t i = 0; i < nhn; i++)
  {
    if (tthn.intersects(getHyperNode(i)))
      return false;
    else
      tthn += getHyperNode(i);
  }

  return true;
}

void MixedSubstitutionModelSet::fireParameterChanged(const ParameterList& parameters)
{
  SubstitutionModelSet::fireParameterChanged(parameters);

  // should be restricted only when probability related parameters are changed
  computeHyperNodesProbabilities();
}


void MixedSubstitutionModelSet::computeHyperNodesProbabilities()
{
  size_t nbm = getNumberOfModels();

  // Looking for the first Mixed model

  size_t fmM = 0;

  MixedSubstitutionModel* pfSM = 0;
  for (fmM = 0; fmM < nbm; fmM++)
  {
    pfSM = dynamic_cast<MixedSubstitutionModel*>(getModel(fmM));
    if (pfSM != NULL)
      break;
  }
  if (fmM == nbm)
    return;

  // Compute the probabilities of the hypernodes from the first mixed
  // model

  size_t nbh = getNumberOfHyperNodes();


  for (size_t nh = 0; nh < nbh; nh++)
  {
    MixedSubstitutionModelSet::HyperNode& h = getHyperNode(nh);
    const MixedSubstitutionModelSet::HyperNode::Node& fnd = h.getNode(fmM);

    double fprob = 0;
    for (size_t i = 0; i < fnd.size(); i++)
    {
      fprob += pfSM->getNProbability(static_cast<size_t>(fnd[i]));
    }
    h.setProbability(fprob);
  }

  // Sets the new probabilities & rates of the mixmodels

  for (size_t iM = fmM + 1; iM < nbm; iM++)
  {
    pfSM = dynamic_cast<MixedSubstitutionModel*>(getModel(iM));
    if (pfSM != NULL)
    {
      for (size_t nh = 0; nh < nbh; nh++)
      {
        MixedSubstitutionModelSet::HyperNode& h = getHyperNode(nh);
        const MixedSubstitutionModelSet::HyperNode::Node& fnd = h.getNode(iM);
        double prob = 0;
        for (size_t j = 0; j < fnd.size(); j++)
        {
          prob += pfSM->getNProbability(static_cast<size_t>(fnd[j]));
        }

        // sets the real probabilities
        for (size_t j = 0; j < fnd.size(); j++)
        {
          pfSM->setNProbability(static_cast<size_t>(fnd[j]), h.getProbability() * pfSM->getNProbability(static_cast<size_t>(fnd[j])) / prob);
        }
      }

      // normalizes Vrates with the real probabilities

      pfSM->normalizeVRates();

      // sets the conditional probabilities

      for (size_t nh = 0; nh < nbh; nh++)
      {
        MixedSubstitutionModelSet::HyperNode& h = getHyperNode(nh);
        const MixedSubstitutionModelSet::HyperNode::Node& fnd = h.getNode(iM);
        for (size_t j = 0; j < fnd.size(); j++)
        {
          pfSM->setNProbability(static_cast<size_t>(fnd[j]), pfSM->getNProbability(static_cast<size_t>(fnd[j])) / h.getProbability());
        }
      }
    }
  }
}

double MixedSubstitutionModelSet::getHyperNodeProbability(const HyperNode& hn) const
{
  size_t nbm = getNumberOfModels();
  double fprob = 1;

  for (size_t fmM = 0; fmM < nbm; fmM++)
  {
    const MixedSubstitutionModelSet::HyperNode::Node& fnd = hn.getNode(fmM);
    const MixedSubstitutionModel* pfSM = dynamic_cast<const MixedSubstitutionModel*>(getModel(fmM));
    if (pfSM != NULL)
    {
      double x = 0;

      for (size_t i = 0; i < fnd.size(); i++)
      {
        x += pfSM->getNProbability(static_cast<size_t>(fnd[i]));
      }

      fprob *= x;
    }
  }

  return fprob;
}

/**********************************************************/
/*************** HYPERNODE ********************************/
/***********************************************************/


MixedSubstitutionModelSet::HyperNode::HyperNode(const MixedSubstitutionModelSet* pMSMS) : vNumbers_(pMSMS->getNumberOfModels()),
  vUnused_(),
  proba_(1)
{
  for (size_t i = 0; i < pMSMS->getNumberOfModels(); i++)
  {
    const MixedSubstitutionModel* pSM = dynamic_cast<const MixedSubstitutionModel*>(pMSMS->getModel(i));
    if (!pSM)
      vUnused_.push_back(static_cast<int>(i));
  }
}

MixedSubstitutionModelSet::HyperNode::HyperNode(const HyperNode& hn) : vNumbers_(hn.vNumbers_),
  vUnused_(hn.vUnused_),
  proba_(hn.proba_)
{}


MixedSubstitutionModelSet::HyperNode& MixedSubstitutionModelSet::HyperNode::operator=(const MixedSubstitutionModelSet::HyperNode& hn)
{
  vNumbers_.clear();
  vNumbers_.resize(hn.vNumbers_.size());
  for (size_t i = 0; i < hn.vNumbers_.size(); i++)
  {
    vNumbers_[i] = hn.vNumbers_[i];
  }
  vUnused_.clear();
  vUnused_.resize(hn.vUnused_.size());
  for (size_t i = 0; i < hn.vUnused_.size(); i++)
  {
    vUnused_[i] = hn.vUnused_[i];
  }

  proba_ = hn.proba_;

  return *this;
}

void MixedSubstitutionModelSet::HyperNode::addToModel(size_t nM, const Vint& vnS)
{
  if (nM >= vNumbers_.size())
    throw IndexOutOfBoundsException("MixedSubstitutionModelSet::HyperNode::addToModel. Bad Mixed model Number", nM, 0, vNumbers_.size());

  vNumbers_[nM].insertN(vnS);
}

void MixedSubstitutionModelSet::HyperNode::setModel(size_t nM, const Vint& vnS)
{
  if (nM >= vNumbers_.size())
    throw IndexOutOfBoundsException("MixedSubstitutionModelSet::HyperNode::setModel. Bad Mixed model Number", nM, 0, vNumbers_.size());

  vNumbers_[nM] = vnS;
}

bool MixedSubstitutionModelSet::HyperNode::isComplete() const
{
  size_t k;
  size_t vUs = vUnused_.size();
  for (size_t i = 0; i < vNumbers_.size(); ++i)
  {
    for (k = 0; k < vUs; k++)
    {
      if (vUnused_[k] == static_cast<int>(i))
        break;
    }
    if ((k == vUs) && vNumbers_[i].size() == 0)
      return false;
  }
  return true;
}

bool MixedSubstitutionModelSet::HyperNode::operator<=(const HyperNode& hn) const
{
  for (size_t i = 0; i < vNumbers_.size(); i++)
  {
    if (!( vNumbers_[i] <= hn.vNumbers_[i]))
      return false;
  }

  return true;
}

bool MixedSubstitutionModelSet::HyperNode::intersects(const HyperNode& hn) const
{
  for (size_t i = 0; i < vNumbers_.size(); i++)
  {
    if (vNumbers_[i].intersects(hn.vNumbers_[i]))
      return true;
  }

  return false;
}

bool MixedSubstitutionModelSet::HyperNode::operator>=(const HyperNode& hn) const
{
  return hn >= *this;
}

MixedSubstitutionModelSet::HyperNode& MixedSubstitutionModelSet::HyperNode::operator+=(const HyperNode& hn)
{
  for (size_t i = 0; i < vNumbers_.size(); i++)
  {
    vNumbers_[i] += hn.vNumbers_[i];
  }

  return *this;
}

/**********************************************************/
/******************** NODE ********************************/
/***********************************************************/

void MixedSubstitutionModelSet::HyperNode::Node::insertN(const Vint& vn)
{
  vector<int>::iterator it;
  vector<int>::const_iterator it2;

  for (it2 = vn.begin(); it2 != vn.end(); it2++)
  {
    for (it = vNumb_.begin(); it != vNumb_.end(); it++)
    {
      if (*it >= *it2)
        break;
    }
    if (it == vNumb_.end())
      vNumb_.push_back(*it2);
    else if (*it != *it2)
      vNumb_.insert(it, *it2);
  }
}

MixedSubstitutionModelSet::HyperNode::Node& MixedSubstitutionModelSet::HyperNode::Node::operator+=(const Node& n)
{
  insertN(n.vNumb_);

  return *this;
}

bool MixedSubstitutionModelSet::HyperNode::Node::operator<=(const Node& n) const
{
  vector<int>::const_iterator it(vNumb_.begin());
  vector<int>::const_iterator it2(n.vNumb_.begin());

  for ( ; it != vNumb_.end(); it++)
  {
    while (it2 != n.vNumb_.end()  && (*it2 < *it))
      it2++;
    if (it2 == n.vNumb_.end() || (*it2 > *it))
      return false;
    it++;
  }
  return true;
}

bool MixedSubstitutionModelSet::HyperNode::Node::operator>=(const Node& n) const
{
  return n <= *this;
}

bool MixedSubstitutionModelSet::HyperNode::Node::intersects(const Node& n) const
{
  vector<int>::const_iterator it(vNumb_.begin());
  vector<int>::const_iterator it2(n.vNumb_.begin());

  for ( ; it != vNumb_.end(); it++)
  {
    while (it2 != n.vNumb_.end()  && (*it2 < *it))
      it2++;

    if (it2 == n.vNumb_.end())
      return false;
    if (*it2 == *it)
      return true;
    it++;
  }
  return false;
}

