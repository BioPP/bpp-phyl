// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "../../Model/MixedTransitionModel.h"
#include "../../Model/MixtureOfASubstitutionModel.h"
#include "MixedSubstitutionModelSet.h"

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
  vpHyperNodes_.push_back(new HyperNode(shared_from_this()));
}

void MixedSubstitutionModelSet::addHyperNode(const HyperNode& hn)
{
  vpHyperNodes_.push_back(new HyperNode(hn));
}

bool MixedSubstitutionModelSet::complete()
{
  MixedSubstitutionModelSet::HyperNode nhn(shared_from_this());
  size_t i;
  for (auto nodei:vpHyperNodes_)
  {
    nhn += *nodei;
  }

  size_t nbm = getNumberOfModels();
  for (i = 0; i < nbm; i++)
  {
    try
    {
      auto& pSM = dynamic_cast<const MixedTransitionModelInterface&>(model(i));
      if (nhn.getNode(i).size() != pSM.getNumberOfModels())
        break;
    }
    catch (exception&)
    {}
  }

  if (i == nbm)
    return false;

  addEmptyHyperNode();
  for (i = 0; i < nbm; i++)
  {
    try
    {
      auto& pSM = dynamic_cast<const MixedTransitionModelInterface&>(model(i));
      const MixedSubstitutionModelSet::HyperNode::Node& nd = nhn.getNode(i);
      auto snd = nd.size();
      auto vs = pSM.getNumberOfModels();
      Vuint an;

      uint j(0), k(0);
      while (j < vs)
      {
        while ((k < snd) && ((uint)nd[k] < j))
          k++;
        if ((k >= snd) || ((uint)nd[k] > j))
          an.push_back(j);
        j++;
      }
      addToHyperNode(i, an);
    }
    catch (exception&)
    {}
  }

  return true;
}

void MixedSubstitutionModelSet::addToHyperNode(size_t nM, const Vuint& vnS, int nH)
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
  HyperNode tthn(shared_from_this());

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

  MixedTransitionModelInterface* pfSM = nullptr;
  for (fmM = 0; fmM < nbm; fmM++)
  {
    pfSM = dynamic_cast<MixedTransitionModelInterface*>(&model(fmM));
    if (pfSM)
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
    pfSM = dynamic_cast<MixedTransitionModelInterface*>(&model(iM));
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
        // sets the REAL probabilities
        for (size_t j = 0; j < fnd.size(); j++)
        {
          pfSM->setNProbability(static_cast<size_t>(fnd[j]), h.getProbability() * pfSM->getNProbability(static_cast<size_t>(fnd[j])) / prob);
        }
      }

      // normalizes Vrates with the real probabilities

      pfSM->normalizeVRates();

      // and now sets the CONDITIONAL probabilities

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
    try
    {
      auto& pfSM = dynamic_cast<const MixedTransitionModelInterface&>(model(fmM));
      double x = 0;

      for (size_t i = 0; i < fnd.size(); ++i)
      {
        x += pfSM.getNProbability(static_cast<size_t>(fnd[i]));
      }

      fprob *= x;
    }
    catch (exception&)
    {}
  }

  return fprob;
}

/**********************************************************/
/*************** HYPERNODE ********************************/
/***********************************************************/


MixedSubstitutionModelSet::HyperNode::HyperNode(std::shared_ptr<const MixedSubstitutionModelSet> pMSMS) :
  vNumbers_(pMSMS->getNumberOfModels()),
  vUnused_(),
  proba_(1)
{
  for (size_t i = 0; i < pMSMS->getNumberOfModels(); i++)
  {
    auto pSM = dynamic_cast<const MixedTransitionModelInterface*>(&pMSMS->model(i));
    if (!pSM)
      vUnused_.push_back(uint(i));
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

void MixedSubstitutionModelSet::HyperNode::addToModel(size_t nM, const Vuint& vnS)
{
  if (nM >= vNumbers_.size())
    throw IndexOutOfBoundsException("MixedSubstitutionModelSet::HyperNode::addToModel. Bad Mixed model Number", nM, 0, vNumbers_.size());

  vNumbers_[nM].insertN(vnS);
}

void MixedSubstitutionModelSet::HyperNode::setModel(size_t nM, const Vuint& vnS)
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
      if (vUnused_[k] == i)
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
  return hn <= *this;
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

void MixedSubstitutionModelSet::HyperNode::Node::insertN(const Vuint& vn)
{
  vector<uint>::iterator it;
  vector<uint>::const_iterator it2;

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
  vector<uint>::const_iterator it2(n.vNumb_.begin());

  for (const auto& it : vNumb_)
  {
    while (it2 != n.vNumb_.end()  && (*it2 < it))
      it2++;
    if (it2 == n.vNumb_.end() || (*it2 > it))
      return false;
  }
  return true;
}

bool MixedSubstitutionModelSet::HyperNode::Node::operator>=(const Node& n) const
{
  return n <= *this;
}

bool MixedSubstitutionModelSet::HyperNode::Node::intersects(const Node& n) const
{
  vector<uint>::const_iterator it2(n.vNumb_.begin());

  for (const auto& it : vNumb_)
  {
    while (it2 != n.vNumb_.end()  && (*it2 < it))
      it2++;

    if (it2 == n.vNumb_.end())
      return false;
    if (*it2 == it)
      return true;
  }
  return false;
}
