// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Numeric/VectorTools.h>
#include <algorithm>

#include "GivenDataSubstitutionProcessSiteSimulator.h"
#include "SimpleSubstitutionProcessSiteSimulator.h"
#include "SubstitutionProcessSequenceSimulator.h"

// From bpp-seq:
#include <Bpp/Seq/Container/VectorSiteContainer.h>

#include <Bpp/Phyl/Likelihood/PartitionSequenceEvolution.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

SubstitutionProcessSequenceSimulator::SubstitutionProcessSequenceSimulator(const SequenceEvolution& evol) :
  mProcess_(),
  vMap_(),
  seqNames_(),
  mvPosNames_()
{
  vector<size_t> nProc = evol.getSubstitutionProcessNumbers();

  vector<shared_ptr<PhyloNode>> vpn = evol.substitutionProcess(nProc[0]).getParametrizablePhyloTree()->getAllLeaves();

  // set ups seqnames for all processes
  for (auto& vi : vpn)
  {
    seqNames_.push_back(vi->getName());
  }

  for (size_t i = 0; i < nProc.size(); ++i)
  {
    const auto sp = evol.getSubstitutionProcess(nProc[i]);

    mProcess_[nProc[i]] = make_unique<SimpleSubstitutionProcessSiteSimulator>(sp);

    vector<string> seqNames2;

    vector<shared_ptr<PhyloNode>> vpn2 = sp->getParametrizablePhyloTree()->getAllLeaves();
    for (size_t i2 = 0; i2 < vpn2.size(); i2++)
    {
      seqNames2.push_back(vpn2[i2]->getName());
    }

    mvPosNames_[nProc[i]].resize(seqNames_.size());

    for (size_t j = 0; j < seqNames_.size(); j++)
    {
      if (!VectorTools::contains(seqNames2, seqNames_[j]))
        throw Exception("SubstitutionProcessSequenceSimulator, unknown sequence name: " +  seqNames_[j]);
      mvPosNames_[nProc[i]][j] = VectorTools::which(seqNames2, seqNames_[j]);
    }
  }

  // set up position specific processes

  auto pse = dynamic_cast<const PartitionSequenceEvolution*>(&evol);

  if (pse)
  {
    setMap(pse->getProcessNumbersPerSite());
  }
  else
    throw Exception("SubstitutionProcessSequenceSimulator::SubstitutionProcessSequenceSimulator(SequenceEvolution) not set for this type of process. Ask developers.");
}

/******************************************************************************/

SubstitutionProcessSequenceSimulator::SubstitutionProcessSequenceSimulator(const SubstitutionProcessSequenceSimulator& spss) :
  mProcess_(),
  vMap_(spss.vMap_),
  seqNames_(spss.seqNames_),
  mvPosNames_(spss.mvPosNames_)
{
  for (auto it : spss.mProcess_)
  {
    auto dit = std::dynamic_pointer_cast<SimpleSubstitutionProcessSiteSimulator>(it.second);
    auto git = std::dynamic_pointer_cast<GivenDataSubstitutionProcessSiteSimulator>(it.second);

    if (dit)
      mProcess_[it.first] = std::make_shared<SimpleSubstitutionProcessSiteSimulator>(*dit);
    else if (git)
      mProcess_[it.first] = std::make_shared<GivenDataSubstitutionProcessSiteSimulator>(*git);
    else
      throw Exception("SubstitutionProcessSequenceSimulator::SubstitutionProcessSequenceSimulator: unknown type of site simulator.");
  }
}

/******************************************************************************/

SubstitutionProcessSequenceSimulator& SubstitutionProcessSequenceSimulator::operator=(const SubstitutionProcessSequenceSimulator& spss)
{
  vMap_ = spss.vMap_;
  seqNames_ = spss.seqNames_;
  mvPosNames_ = spss.mvPosNames_;

  mProcess_.clear();

  for (auto it : spss.mProcess_)
  {
    auto dit = std::dynamic_pointer_cast<SimpleSubstitutionProcessSiteSimulator>(it.second);
    auto git = std::dynamic_pointer_cast<GivenDataSubstitutionProcessSiteSimulator>(it.second);

    if (dit)
      mProcess_[it.first] = std::make_shared<SimpleSubstitutionProcessSiteSimulator>(*dit);
    else if (git)
      mProcess_[it.first] = std::make_shared<GivenDataSubstitutionProcessSiteSimulator>(*git);
    else
      throw Exception("SubstitutionProcessSequenceSimulator::SubstitutionProcessSequenceSimulator: unknown type of site simulator.");
  }

  return *this;
}

/******************************************************************************/

void SubstitutionProcessSequenceSimulator::outputInternalSequences(bool yn)
{
  for (auto it : mProcess_)
  {
    it.second->outputInternalSites(yn);
  }
}

/******************************************************************************/

void SubstitutionProcessSequenceSimulator::setMap(std::vector<size_t> vMap)
{
  vMap_.clear();

  for (size_t i = 0; i < vMap.size(); i++)
  {
    if (mProcess_.find(vMap[i]) == mProcess_.end())
      throw Exception("SubstitutionProcessSequenceSimulator::setMap: unknown Process number" + TextTools::toString(vMap[i]));
    else
      vMap_.push_back(vMap[i]);
  }
}

/******************************************************************************/

unique_ptr<SiteContainerInterface> SubstitutionProcessSequenceSimulator::simulate(
    size_t numberOfSites) const
{
  if (numberOfSites > vMap_.size())
    throw BadIntegerException("SubstitutionProcessSequenceSimulator::simulate. Too many sites to simulate.", (int)numberOfSites);

  auto sites = make_unique<VectorSiteContainer>(seqNames_, getAlphabet());
  sites->setSequenceNames(seqNames_, true);

  Vint vval(seqNames_.size());

  for (size_t j = 0; j < numberOfSites; ++j)
  {
    auto site = mProcess_.find(vMap_[j])->second->simulateSite();

    const vector<size_t>& vPosNames = mvPosNames_.find(vMap_[j])->second;
    for (size_t vn = 0; vn < vPosNames.size(); vn++)
    {
      vval[vn] = site->getValue(vPosNames[vn]);
    }

    auto site2 = make_unique<Site>(vval, sites->getAlphabet(), static_cast<int>(j));
    sites->addSite(site2);
  }
  return sites;
}

/******************************************************************************/

unique_ptr<SiteContainerInterface> SubstitutionProcessSequenceSimulator::simulate(
    const vector<double>& rates) const
{
  size_t numberOfSites = rates.size();

  if (numberOfSites > vMap_.size())
    throw Exception("SubstitutionProcessSequenceSimulator::simulate : some sites do not have attributed process");

  auto sites = make_unique<VectorSiteContainer>(seqNames_, getAlphabet());
  sites->setSequenceNames(seqNames_, true);

  Vint vval(seqNames_.size());

  for (size_t j = 0; j < numberOfSites; ++j)
  {
    auto site = mProcess_.find(vMap_[j])->second->simulateSite(rates[j]);

    const vector<size_t>& vPosNames = mvPosNames_.find(vMap_[j])->second;
    for (size_t vn = 0; vn < vPosNames.size(); vn++)
    {
      vval[vn] = site->getValue(vPosNames[vn]);
    }

    auto site2 = make_unique<Site>(vval, sites->getAlphabet(), static_cast<int>(j));
    sites->addSite(site2);
  }
  return sites;
}

/******************************************************************************/

unique_ptr<SiteContainerInterface> SubstitutionProcessSequenceSimulator::simulate(
    const vector<size_t>& states) const
{
  size_t numberOfSites = states.size();

  auto sites = make_unique<VectorSiteContainer>(seqNames_, getAlphabet());
  sites->setSequenceNames(seqNames_, true);

  Vint vval(seqNames_.size());

  for (size_t j = 0; j < numberOfSites; ++j)
  {
    auto site = mProcess_.find(vMap_[j])->second->simulateSite(states[j]);

    const vector<size_t>& vPosNames = mvPosNames_.find(vMap_[j])->second;
    for (size_t vn = 0; vn < vPosNames.size(); vn++)
    {
      vval[vn] = site->getValue(vPosNames[vn]);
    }

    auto site2 = make_unique<Site>(vval, sites->getAlphabet(), static_cast<int>(j));
    sites->addSite(site2);
  }
  return sites;
}

/******************************************************************************/

unique_ptr<SiteContainerInterface> SubstitutionProcessSequenceSimulator::simulate(
    const vector<double>& rates,
    const vector<size_t>& states) const
{
  size_t numberOfSites = rates.size();
  if (states.size() != numberOfSites)
    throw Exception("SubstitutionProcessSequenceSimulator::simulate, 'rates' and 'states' must have the same length.");

  auto sites = make_unique<VectorSiteContainer>(seqNames_, getAlphabet());
  sites->setSequenceNames(seqNames_, true);

  Vint vval(seqNames_.size());

  for (size_t j = 0; j < numberOfSites; ++j)
  {
    auto site = mProcess_.find(vMap_[j])->second->simulateSite(states[j], rates[j]);

    const vector<size_t>& vPosNames = mvPosNames_.find(vMap_[j])->second;
    for (size_t vn = 0; vn < vPosNames.size(); vn++)
    {
      vval[vn] = site->getValue(vPosNames[vn]);
    }

    auto site2 = make_unique<Site>(vval, sites->getAlphabet(), static_cast<int>(j));
    sites->addSite(site2);
  }
  return sites;
}

/******************************************************************************/

shared_ptr<const Alphabet> SubstitutionProcessSequenceSimulator::getAlphabet() const
{
  if (mProcess_.size() == 0)
    return nullptr;

  return mProcess_.begin()->second->getAlphabet();
}

/******************************************************************************/

const Alphabet& SubstitutionProcessSequenceSimulator::alphabet() const
{
  if (mProcess_.size() == 0)
    throw NullPointerException("SubstitutionProcessSequenceSimulator::alphabet(). No process attached.");

  return mProcess_.begin()->second->alphabet();
}

/******************************************************************************/
