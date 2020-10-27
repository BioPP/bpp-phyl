//
// File: SubstitutionProcessSequenceSimulator.cpp
// Created by: Julien Dutheil
//             Bastien Boussau
//             Laurent Guéguen
// Created on: Wed Feb  4 16:30:51 2004
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

#include "SubstitutionProcessSequenceSimulator.h"

#include "SimpleSubstitutionProcessSiteSimulator.h"
#include "GivenDataSubstitutionProcessSiteSimulator.h"

#include <algorithm>

#include <Bpp/Numeric/VectorTools.h>

// From SeqLib:
#include <Bpp/Seq/Container/VectorSiteContainer.h>

#include <Bpp/Phyl/NewLikelihood/PartitionSequenceEvolution.h>

using namespace bpp;
using namespace std;

SubstitutionProcessSequenceSimulator::SubstitutionProcessSequenceSimulator(const SequenceEvolution& evol) :
  mProcess_(),
  vMap_(),
  seqNames_(),
  mvPosNames_()
{
  vector<size_t> nProc=evol.getSubstitutionProcessNumbers();
  
  vector<shared_ptr<PhyloNode> > vpn= evol.getSubstitutionProcess(nProc[0]).getParametrizablePhyloTree().getAllLeaves();

  // set ups seqnames for all processes
  for (size_t i=0;i<vpn.size();i++)
    seqNames_.push_back(vpn[i]->getName());

  for (size_t i=0; i< nProc.size(); i++)
  {
    const SubstitutionProcess& sp=evol.getSubstitutionProcess(nProc[i]);

    mProcess_[nProc[i]] = std::make_shared<SimpleSubstitutionProcessSiteSimulator>(sp);

    vector<string> seqNames2;
    
    vector<shared_ptr<PhyloNode> > vpn2= sp.getParametrizablePhyloTree().getAllLeaves();
    for (size_t i2=0;i2<vpn2.size();i2++)
      seqNames2.push_back(vpn2[i2]->getName());

    mvPosNames_[nProc[i]].resize(seqNames_.size());

    for (size_t j=0; j<seqNames_.size(); j++)
      mvPosNames_[nProc[i]][j]=VectorTools::which(seqNames2,seqNames_[j]);
  }
  
  // set up position specific processes

  auto pse = dynamic_cast<const PartitionSequenceEvolution*>(&evol);

  if (pse)
  {
    setMap(pse->getProcessNumbersPerSite());
  }
  else
    throw Exception("SubstitutionProcessSequenceSimulator::SubstitutionProcessSequenceSimulator(SequenceEvolution) not set for this type of process. Ask developpers.");
  
}
  
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
      mProcess_[it.first]=std::make_shared<SimpleSubstitutionProcessSiteSimulator>(*dit);
    else if (git)
      mProcess_[it.first]=std::make_shared<GivenDataSubstitutionProcessSiteSimulator>(*git);
    else
      throw Exception("SubstitutionProcessSequenceSimulator::SubstitutionProcessSequenceSimulator: unknown type of site simulator.");
  }
}

    
SubstitutionProcessSequenceSimulator& SubstitutionProcessSequenceSimulator::operator=(const SubstitutionProcessSequenceSimulator& spss)
{
  vMap_=spss.vMap_;
  seqNames_=spss.seqNames_;
  mvPosNames_=spss.mvPosNames_;
  
  mProcess_.clear();
  
  for (auto it : spss.mProcess_)
  {
    auto dit = std::dynamic_pointer_cast<SimpleSubstitutionProcessSiteSimulator>(it.second);
    auto git = std::dynamic_pointer_cast<GivenDataSubstitutionProcessSiteSimulator>(it.second);
    
    if (dit)
      mProcess_[it.first]=std::make_shared<SimpleSubstitutionProcessSiteSimulator>(*dit);
    else if (git)
      mProcess_[it.first]=std::make_shared<GivenDataSubstitutionProcessSiteSimulator>(*git);
    else
      throw Exception("SubstitutionProcessSequenceSimulator::SubstitutionProcessSequenceSimulator: unknown type of site simulator.");
  }
  
  return *this;
}


SubstitutionProcessSequenceSimulator::~SubstitutionProcessSequenceSimulator()
{
}

void SubstitutionProcessSequenceSimulator::outputInternalSequences(bool yn)
{
  for (auto it : mProcess_)
    it.second->outputInternalSites(yn);
}

void SubstitutionProcessSequenceSimulator::setMap(std::vector<size_t> vMap)
{
  vMap_.clear();

  for (size_t i=0; i<vMap.size(); i++)
    if (mProcess_.find(vMap[i])==mProcess_.end())
      throw Exception("SubstitutionProcessSequenceSimulator::setMap: unknown Process number" + TextTools::toString(vMap[i]));
    else
      vMap_.push_back(vMap[i]);
}


std::shared_ptr<SiteContainer> SubstitutionProcessSequenceSimulator::simulate(size_t numberOfSites) const
{
  auto sites = make_shared<VectorSiteContainer>(seqNames_, getAlphabet());
  sites->setSequencesNames(seqNames_);

  Vint vval(seqNames_.size());
  
  for (size_t j = 0; j < numberOfSites; j++)
  {
    Site* site=mProcess_.find(vMap_[j])->second->simulateSite();

    const vector<size_t>& vPosNames=mvPosNames_.find(vMap_[j])->second;
    for (size_t vn=0;vn<vPosNames.size(); vn++)
      vval[vn]=site->getValue(vPosNames[vn]);

    sites->addSite(*new Site(vval,sites->getAlphabet(),static_cast<int>(j)));
    delete site;
  }
  return sites;
}

std::shared_ptr<SiteContainer> SubstitutionProcessSequenceSimulator::simulate(const vector<double>& rates) const
{
  size_t numberOfSites=rates.size();
  
  if (numberOfSites>vMap_.size())
    throw Exception("SubstitutionProcessSequenceSimulator::simulate : some sites do not have attributed process");

  auto sites = make_shared<VectorSiteContainer>(seqNames_, getAlphabet());
  sites->setSequencesNames(seqNames_);

  Vint vval(seqNames_.size());
  
  for (size_t j = 0; j < numberOfSites; j++)
  {
    Site* site=mProcess_.find(vMap_[j])->second->simulateSite(rates[j]);

    const vector<size_t>& vPosNames=mvPosNames_.find(vMap_[j])->second;
    for (size_t vn=0;vn<vPosNames.size(); vn++)
      vval[vn]=site->getValue(vPosNames[vn]);
    
    sites->addSite(*new Site(vval,sites->getAlphabet(),static_cast<int>(j)));
    delete site;
  }
  return sites;
}

std::shared_ptr<SiteContainer> SubstitutionProcessSequenceSimulator::simulate(const vector<size_t>& states) const
{
  size_t numberOfSites=states.size();

  auto sites = make_shared<VectorSiteContainer>(seqNames_, getAlphabet());
  sites->setSequencesNames(seqNames_);

  Vint vval(seqNames_.size());
  
  for (size_t j = 0; j < numberOfSites; j++)
  {
    Site* site=mProcess_.find(vMap_[j])->second->simulateSite(states[j]);

    const vector<size_t>& vPosNames=mvPosNames_.find(vMap_[j])->second;
    for (size_t vn=0;vn<vPosNames.size(); vn++)
      vval[vn]=site->getValue(vPosNames[vn]);
    
    sites->addSite(*new Site(vval,sites->getAlphabet(),static_cast<int>(j)));
    delete site;
  }
  return sites;
}

std::shared_ptr<SiteContainer> SubstitutionProcessSequenceSimulator::simulate(const vector<double>& rates, const vector<size_t>& states) const
{
  size_t numberOfSites = rates.size();
  if (states.size() != numberOfSites)
    throw Exception("SubstitutionProcessSequenceSimulator::simulate, 'rates' and 'states' must have the same length.");

  auto sites = make_shared<VectorSiteContainer>(seqNames_, getAlphabet());
  sites->setSequencesNames(seqNames_);

  Vint vval(seqNames_.size());
  
  for (size_t j = 0; j < numberOfSites; j++)
  {
    Site* site=mProcess_.find(vMap_[j])->second->simulateSite(states[j], rates[j]);

    const vector<size_t>& vPosNames=mvPosNames_.find(vMap_[j])->second;
    for (size_t vn=0;vn<vPosNames.size(); vn++)
      vval[vn]=site->getValue(vPosNames[vn]);
    
    sites->addSite(*new Site(vval,sites->getAlphabet(),static_cast<int>(j)));
    delete site;
  }
  return sites;
}

  
const Alphabet* SubstitutionProcessSequenceSimulator::getAlphabet() const
{
  if (mProcess_.size()==0)
    return NULL;
  
  return mProcess_.begin()->second->getAlphabet();
}

/******************************************************************************/
    

