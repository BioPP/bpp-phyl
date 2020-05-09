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
#include <algorithm>

#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

#include <Bpp/Phyl/NewLikelihood/ProcessComputationTree.h>
// From SeqLib:
#include <Bpp/Seq/Container/VectorSiteContainer.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

SimpleSubstitutionProcessSequenceSimulator::SimpleSubstitutionProcessSequenceSimulator(const SubstitutionProcess& process) :
  process_(&process),
  alphabet_(process_->getStateMap().getAlphabet()),
  phyloTree_(&process_->getParametrizablePhyloTree()),
  tree_(ProcessComputationTree(*process_)),
  seqIndexes_(),
  seqNames_(),
  speciesNodes_(),
  nbNodes_(),
  nbClasses_(process_->getNumberOfClasses()),
  nbStates_(process_->getNumberOfStates()),
  continuousRates_(false),
  outputInternalSequences_(false)
{
  init();
}

/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::init()
{
  // Initialize sons & fathers of tree_ Nodes    
  // set sequence names

  if (outputInternalSequences_) {
    auto vCN= phyloTree_->getAllNodes();
    seqNames_.resize(vCN.size());    
    seqIndexes_.resize(vCN.size());    
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      auto index = phyloTree_->getNodeIndex(vCN[i]);
      seqNames_[i] = (phyloTree_->isLeaf(vCN[i]))?vCN[i]->getName():TextTools::toString(index);
      seqIndexes_[i] = index; 
    }
  }
  else {
    auto vCN= phyloTree_->getAllLeaves();
    seqNames_.resize(vCN.size());    
    seqIndexes_.resize(vCN.size());    
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      seqNames_[i] = vCN[i]->getName();
      seqIndexes_[i] = phyloTree_->getNodeIndex(vCN[i]);
    }
  }
  
  // Initialize cumulative pxy for edges that have models
  auto edges = tree_.getAllEdges();

  const auto dRate = process_->getRateDistribution();
  
  for (auto& edge : edges)
  {
    const auto model = edge->getModel();
    if (!model)
      continue;
    
    const auto transmodel = dynamic_cast<const TransitionModel*>(model);
    if (!transmodel)
      throw Exception("SubstitutionProcessSequenceSimulator::init : model "  + model->getName() + " on branch " + TextTools::toString(tree_.getEdgeIndex(edge)) + " is not a TransitionModel.");
    
    VVVdouble* cumpxy_node_ = &edge->cumpxy_;
    cumpxy_node_->resize(nbClasses_);
    
    for (size_t c = 0; c < nbClasses_; c++)
    {
      double brlen = dRate->getCategory(c) * phyloTree_->getEdge(edge->getSpeciesIndex())->getLength();
    
      VVdouble* cumpxy_node_c_ = &(*cumpxy_node_)[c];
    
      cumpxy_node_c_->resize(nbStates_);
    
      // process transition probabilities already consider rates &
      // branch length

      const Matrix<double>* P;
      
      const auto& vSub(edge->subModelNumbers());

      if (vSub.size()==0)
        P = &transmodel->getPij_t(brlen);
      else
      {
        if (vSub.size()>1)
          throw Exception("SubstitutionProcessSequenceSimulator::init : only 1 submodel can be used.");
        
        const auto* mmodel = dynamic_cast<const MixedTransitionModel*>(transmodel);
        
        const auto* model2 = mmodel->getNModel(vSub[0]);
        
        P = &model2->getPij_t(brlen);
      }
      
      for (size_t x = 0; x < nbStates_; x++)
      {
        Vdouble* cumpxy_node_c_x_ = &(*cumpxy_node_c_)[x];
        cumpxy_node_c_x_->resize(nbStates_);
        (*cumpxy_node_c_x_)[0] = (*P)(x, 0);
        for (size_t y = 1; y < nbStates_; y++)
        {
          (*cumpxy_node_c_x_)[y] = (*cumpxy_node_c_x_)[y - 1] + (*P)(x, y);
        }
      }
    }
  }

  // Initialize cumulative prob for mixture nodes
  auto nodes = tree_.getAllNodes();

  for (auto node:nodes)
  {
    if (node->isMixture()) // set probas to chose
    {
      auto outEdges = tree_.getOutgoingEdges(node);
      Vdouble vprob(0);

      for (auto edge : outEdges)
      {
        auto model = dynamic_cast<const MixedTransitionModel*>(edge->getModel());
        if (!model)
          throw Exception("SubstitutionProcessSequenceSimulator::init : model in edge " + TextTools::toString(tree_.getEdgeIndex(edge)) + " is not a mixture.");

        const auto& vNb(edge->subModelNumbers());

        double x=0.;
        for (auto nb:vNb)
          x += model->getNProbability(nb);

        vprob.push_back(x);
        node->sons_.push_back(tree_.getSon(edge));
      }

      vprob /= VectorTools::sum(vprob);

      node->cumProb_ = VectorTools::cumSum(vprob);
    }
  }

}

/******************************************************************************/

Site* SimpleSubstitutionProcessSequenceSimulator::simulateSite() const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  size_t initialStateIndex = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  const vector<double>& freqs = process_->getRootFrequencies();
  for (size_t i = 0; i < nbStates_; i++)
  {
    cumprob += freqs[i];
    if (r <= cumprob)
    {
      initialStateIndex = i;
      break;
    }
  }

  return simulateSite(initialStateIndex);
}


Site* SimpleSubstitutionProcessSequenceSimulator::simulateSite(double rate) const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  size_t initialStateIndex = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  const vector<double>& freqs = process_->getRootFrequencies();
  for (size_t i = 0; i < nbStates_; i++)
  {
    cumprob += freqs[i];
    if (r <= cumprob)
    {
      initialStateIndex = i;
      break;
    }
  }

  return simulateSite(initialStateIndex, rate);
}


Site* SimpleSubstitutionProcessSequenceSimulator::simulateSite(size_t ancestralStateIndex) const
{
  shared_ptr<SimProcessNode> root = tree_.getRoot();
  root->state_ = ancestralStateIndex;

  if (continuousRates_ && process_->getRateDistribution())
  {
    double rate = process_->getRateDistribution()->randC();
    evolveInternal(root, rate);
  }
  else
  {
    size_t rateClass = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(nbClasses_);
    evolveInternal(root, rateClass);
  }
  
  
  // Now create a Site object:
  Vint site(seqNames_.size());
  for (size_t i = 0; i < seqNames_.size(); ++i)
  {
    site[i] = process_->getStateMap().getAlphabetStateAsInt(speciesNodes_.at(seqIndexes_[i])->state_);
  }
  return new Site(site, alphabet_);
}


Site* SimpleSubstitutionProcessSequenceSimulator::simulateSite(size_t ancestralStateIndex, double rate) const
{
  shared_ptr<SimProcessNode> root = tree_.getRoot();
  root->state_ = ancestralStateIndex;

  evolveInternal(root, rate);
  
  // Now create a Site object:
  Vint site(seqNames_.size());
  for (size_t i = 0; i < seqNames_.size(); ++i)
  {
    site[i] = process_->getStateMap().getAlphabetStateAsInt(speciesNodes_.at(seqIndexes_[i])->state_);
  }
  return new Site(site, alphabet_);
}



/******************************************************************************/

SiteContainer* SimpleSubstitutionProcessSequenceSimulator::simulate(size_t numberOfSites) const
{
  VectorSiteContainer* sites = new VectorSiteContainer(seqNames_.size(), alphabet_);
  sites->setSequencesNames(seqNames_);
  for (size_t j = 0; j < numberOfSites; j++)
  {
    Site* site = simulateSite();
    site->setPosition(static_cast<int>(j));
    sites->addSite(*site);
    delete site;
  }
  return sites;
}

/******************************************************************************/

New_SiteSimulationResult* SimpleSubstitutionProcessSequenceSimulator::dSimulateSite() const
{
  // Draw an initial state randomly according to root frequencies:
  size_t ancestralStateIndex = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  const auto& freqs = process_->getRootFrequencies();
  for (size_t i = 0; i < nbStates_; i++)
  {
    cumprob += freqs[i];
    if (r <= cumprob)
    {
      ancestralStateIndex = i;
      break;
    }
  }

  return dSimulateSite(ancestralStateIndex);
}


New_SiteSimulationResult* SimpleSubstitutionProcessSequenceSimulator::dSimulateSite(size_t ancestralStateIndex) const
{
  shared_ptr<SimProcessNode> root = tree_.getRoot();
  root->state_ = ancestralStateIndex;

  New_SiteSimulationResult* ssr = new New_SiteSimulationResult(phyloTree_, &process_->getStateMap(), ancestralStateIndex);

  // Draw a random rate:
  if (continuousRates_ && process_->getRateDistribution())
  {
    // Draw a random rate:
    double rate = process_->getRateDistribution()->randC();
    evolveInternal(root, rate, ssr);
  }
  else
  {
    size_t rateClass = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(nbClasses_);
    evolveInternal(root, rateClass, ssr);
  }
  
  return ssr;
}

New_SiteSimulationResult* SimpleSubstitutionProcessSequenceSimulator::dSimulateSite(double rate) const
{
  // Draw an initial state randomly according to root frequencies:
  size_t ancestralStateIndex = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  const auto& freqs = process_->getRootFrequencies();
  for (size_t i = 0; i < nbStates_; i++)
  {
    cumprob += freqs[i];
    if (r <= cumprob)
    {
      ancestralStateIndex = i;
      break;
    }
  }

  shared_ptr<SimProcessNode> root = tree_.getRoot();
  root->state_ = ancestralStateIndex;

  New_SiteSimulationResult* ssr = new New_SiteSimulationResult(phyloTree_, &process_->getStateMap(), ancestralStateIndex);

  evolveInternal(root, rate, ssr);

  return ssr;
}

New_SiteSimulationResult* SimpleSubstitutionProcessSequenceSimulator::dSimulateSite(size_t ancestralStateIndex, double rate) const
{
  shared_ptr<SimProcessNode> root = tree_.getRoot();
  root->state_ = ancestralStateIndex;

  New_SiteSimulationResult* ssr = new New_SiteSimulationResult(phyloTree_, &process_->getStateMap(), ancestralStateIndex);

  evolveInternal(root, rate, ssr);

  return ssr;
}


/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::evolveInternal(std::shared_ptr<SimProcessNode> node, size_t rateClass, New_SiteSimulationResult * ssr) const
{
  if (node->isSpeciation())
  {
    auto vEdge = tree_.getOutgoingEdges(node);

    for (auto edge : vEdge)
    {
      auto son = tree_.getSon(edge);

      if (edge->getModel())
      {
        if (ssr) // Detailed simulation
        {
          auto tm = dynamic_cast<const SubstitutionModel*>(edge->getModel());

          if (!tm)
            throw Exception("SimpleSubstitutionProcessSequenceSimulator::EvolveInternal : detailed simulation not possible for non-markovian model on edge " + TextTools::toString(son->getSpeciesIndex()) + " for model " + edge->getModel()->getName());

          SimpleMutationProcess process(tm);

          double brlen = process_->getRateDistribution()->getCategory(rateClass) * phyloTree_->getEdge(edge->getSpeciesIndex())->getLength();

          MutationPath mp = process.detailedEvolve(node->state_, brlen);
  
          son->state_ = mp.getFinalState();

          // Now append infos in ssr:
          ssr->addNode(edge->getSpeciesIndex(), mp);
          
          speciesNodes_[son->getSpeciesIndex()] = son;
        }
        else
        {
          Vdouble* cumy = &edge->cumpxy_[rateClass][node->state_];
          
          double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
          for (size_t y = 0; y < nbStates_; y++)
          {
            if (rand < (*cumy)[y])
            {
              son->state_ = y;
              speciesNodes_[son->getSpeciesIndex()] = son;
              break;
            }
          }
        }
      }
      else
      {
        son->state_ = node->state_;
        speciesNodes_[son->getSpeciesIndex()] = son;
      }
      
      evolveInternal(son, rateClass, ssr);
    }
  }
  else
    if (node->isMixture())
    {
      const auto& cumProb = node->cumProb_;
      
      double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);

      for (size_t y = 0; y < cumProb.size(); y++)
      {
        if (rand < cumProb[y])
        {
          auto son = node->sons_[y];
          son->state_ = node->state_;
          speciesNodes_[son->getSpeciesIndex()] = son;
          evolveInternal(son, rateClass, ssr);
          break;
        }
      }
    }
    else
      throw Exception("SimpleSubstitutionProcessSequenceSimulator::evolveInternal : unknown property for node " + tree_.getNodeIndex(node));
}

/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::evolveInternal(std::shared_ptr<SimProcessNode> node, double rate, New_SiteSimulationResult * ssr) const
{
  if (node->isSpeciation())
  {
    auto vEdge = tree_.getOutgoingEdges(node);

    for (auto edge : vEdge)
    {
      auto son = tree_.getSon(edge);

      if (edge->getModel())
      {
        auto tm = dynamic_cast<const TransitionModel*>(edge->getModel());
        
        double brlen = rate * phyloTree_->getEdge(edge->getSpeciesIndex())->getLength(); 
        if (ssr) // Detailed simulation
        {
          auto sm = dynamic_cast<const SubstitutionModel*>(edge->getModel());

          if (!sm)
            throw Exception("SimpleSubstitutionProcessSequenceSimulator::EvolveInternal : detailed simulation not possible for non-markovian model on edge " + TextTools::toString(son->getSpeciesIndex()) + " for model " + tm->getName());

          SimpleMutationProcess process(sm);

          MutationPath mp = process.detailedEvolve(node->state_, brlen);
  
          son->state_ = mp.getFinalState();

          // Now append infos in ssr:
          ssr->addNode(edge->getSpeciesIndex(), mp);
          
          speciesNodes_[son->getSpeciesIndex()] = son;
        }
        else
        {
          // process transition probabilities already consider rates &
          // branch length

          const Matrix<double>* P;
          
          const auto& vSub(edge->subModelNumbers());
          
          if (vSub.size()==0)
            P = &tm->getPij_t(brlen);
          else
          {
            if (vSub.size()>1)
              throw Exception("SubstitutionProcessSequenceSimulator::init : only 1 submodel can be used.");
            
            const auto* mmodel = dynamic_cast<const MixedTransitionModel*>(tm);
            const auto* model = mmodel->getNModel(vSub[0]);
            
            P = &model->getPij_t(brlen);
          }

          double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
          for (size_t y = 0; y < nbStates_; y++)
          {
            rand -= (*P)(node->state_,y);
            if (rand <= 0)
            {
              son->state_ = y;
              speciesNodes_[son->getSpeciesIndex()] = son;
              break;
            }
          }
        }
      }
      else
      {
        son->state_ = node->state_;
        speciesNodes_[son->getSpeciesIndex()] = son;
      }
      
      evolveInternal(son, rate, ssr);
    }
  }
  else
    if (node->isMixture())
    {
      const auto& cumProb = node->cumProb_;
      
      double rand2 = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);

      for (size_t y = 0; y < cumProb.size(); y++)
      {
        if (rand2 < cumProb[y])
        {
          auto son = node->sons_[y];
          son->state_ = node->state_;
          speciesNodes_[son->getSpeciesIndex()] = son;
          evolveInternal(son, rate, ssr);
          break;
        }
      }
    }
    else
      throw Exception("SimpleSubstitutionProcessSequenceSimulator::evolveInternal : unknown property for node " + tree_.getNodeIndex(node));
}


void SimpleSubstitutionProcessSequenceSimulator::outputInternalSequences(bool yn)
{
  outputInternalSequences_ = yn;
  
  if (outputInternalSequences_) {
    auto vCN= phyloTree_->getAllNodes();
    seqNames_.resize(vCN.size());    
    seqIndexes_.resize(vCN.size());    
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      auto index = phyloTree_->getNodeIndex(vCN[i]);
      seqNames_[i] = (phyloTree_->isLeaf(vCN[i]))?vCN[i]->getName():TextTools::toString(index);
      seqIndexes_[i] = index; 
    }
  }
  else {
    auto vCN= phyloTree_->getAllLeaves();
    seqNames_.resize(vCN.size());    
    seqIndexes_.resize(vCN.size());    
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      seqNames_[i] = vCN[i]->getName();
      seqIndexes_[i] = phyloTree_->getNodeIndex(vCN[i]);
    }
  }
}


/******************************************************************************/
/******************************************************************************/
/******************    SubstitutionProcessSequenceSimulator       *************/
/******************************************************************************/


SubstitutionProcessSequenceSimulator::SubstitutionProcessSequenceSimulator(const SequenceEvolution& evol) :
  mProcess_(),
  vMap_(),
  seqNames_(),
  mvPosNames_()
{
  vector<size_t> nProc=evol.getSubstitutionProcessNumbers();
  
  vector<shared_ptr<PhyloNode> > vpn= evol.getSubstitutionProcess(nProc[0]).getParametrizablePhyloTree().getAllLeaves();
  
  for (size_t i=0;i<vpn.size();i++)
    seqNames_.push_back(vpn[i]->getName());

  for (size_t i=0; i< nProc.size(); i++)
  {
    const SubstitutionProcess& sp=evol.getSubstitutionProcess(nProc[i]);
    
    mProcess_[nProc[i]]=new SimpleSubstitutionProcessSequenceSimulator(sp);

    vector<string> seqNames2;
    
    vector<shared_ptr<PhyloNode> > vpn2= sp.getParametrizablePhyloTree().getAllLeaves();
    for (size_t i2=0;i2<vpn2.size();i2++)
      seqNames2.push_back(vpn2[i2]->getName());

    mvPosNames_[nProc[i]].resize(seqNames_.size());

    for (size_t j=0; j<seqNames_.size(); j++)
      mvPosNames_[nProc[i]][j]=VectorTools::which(seqNames2,seqNames_[j]);
  }
}
  
SubstitutionProcessSequenceSimulator::SubstitutionProcessSequenceSimulator(const std::map<size_t, const SubstitutionProcess&>& mSP) :
  mProcess_(),
  vMap_(),
  seqNames_(),
  mvPosNames_()
{
  if (mSP.size()!=0)    
  {
    vector<shared_ptr<PhyloNode> > vpn= mSP.begin()->second.getParametrizablePhyloTree().getAllLeaves();
    for (size_t i=0;i<vpn.size();i++)
      seqNames_.push_back(vpn[i]->getName());
  }
  

  for (std::map<size_t, const SubstitutionProcess&>::const_iterator  it=mSP.begin(); it != mSP.end(); it++)
  {
    mProcess_[it->first]=new SimpleSubstitutionProcessSequenceSimulator(it->second);

    vector<string> seqNames2;
    
    vector<shared_ptr<PhyloNode> > vpn2= it->second.getParametrizablePhyloTree().getAllLeaves();
    for (size_t i=0;i<vpn2.size();i++)
      seqNames2.push_back(vpn2[i]->getName());

    mvPosNames_[it->first].resize(seqNames_.size());

    for (size_t i=0; i<seqNames_.size(); i++)
      mvPosNames_[it->first][i]=VectorTools::which(seqNames2,seqNames_[i]);
  }
}
  
SubstitutionProcessSequenceSimulator::SubstitutionProcessSequenceSimulator(const SubstitutionProcessCollection& spc) :
  mProcess_(),
  vMap_(),
  seqNames_(),
  mvPosNames_()
{
  vector<size_t> procN=spc.getSubstitutionProcessNumbers();

  if (procN.size()==0)
    return;
  
  vector<shared_ptr<PhyloNode> > vpn= spc.getSubstitutionProcess(procN[0]).getParametrizablePhyloTree().getAllLeaves();
  for (size_t i=0;i<vpn.size();i++)
    seqNames_.push_back(vpn[i]->getName());
  
  for (size_t i=0; i<procN.size(); i++)
  {
    const SubstitutionProcess& sp=spc.getSubstitutionProcess(procN[i]);
    
    mProcess_[procN[i]]=new SimpleSubstitutionProcessSequenceSimulator(sp);

    vector<string> seqNames2;
    
    vector<shared_ptr<PhyloNode> > vpn2= sp.getParametrizablePhyloTree().getAllLeaves();
    for (size_t i2=0;i2<vpn2.size();i2++)
      seqNames2.push_back(vpn2[i2]->getName());

    mvPosNames_[procN[i]].resize(seqNames_.size());

    for (size_t i2=0; i2<seqNames_.size(); i2++){
      mvPosNames_[procN[i]][i2]=VectorTools::which(seqNames2,seqNames_[i2]);
    }
    
  }
}

SubstitutionProcessSequenceSimulator::SubstitutionProcessSequenceSimulator(const SubstitutionProcessSequenceSimulator& spss) :
  mProcess_(),
  vMap_(spss.vMap_),
  seqNames_(spss.seqNames_),
  mvPosNames_(spss.mvPosNames_)
{
  for (std::map<size_t, SimpleSubstitutionProcessSequenceSimulator*>::const_iterator  it=spss.mProcess_.begin(); it != spss.mProcess_.end(); it++)
    mProcess_[it->first]=new SimpleSubstitutionProcessSequenceSimulator(*it->second);
}

    
SubstitutionProcessSequenceSimulator& SubstitutionProcessSequenceSimulator::operator=(const SubstitutionProcessSequenceSimulator& spss)
{
  vMap_=spss.vMap_;
  seqNames_=spss.seqNames_;
  mvPosNames_=spss.mvPosNames_;
  
  for (std::map<size_t, SimpleSubstitutionProcessSequenceSimulator*>::const_iterator  it=mProcess_.begin(); it != mProcess_.end(); it++)
    delete it->second;

  mProcess_.clear();
  
  for (std::map<size_t, SimpleSubstitutionProcessSequenceSimulator*>::const_iterator  it=spss.mProcess_.begin(); it != spss.mProcess_.end(); it++)
    mProcess_[it->first]=new SimpleSubstitutionProcessSequenceSimulator(*it->second);

  return *this;
}


SubstitutionProcessSequenceSimulator::~SubstitutionProcessSequenceSimulator()
{
  for (std::map<size_t, SimpleSubstitutionProcessSequenceSimulator*>::const_iterator  it=mProcess_.begin(); it != mProcess_.end(); it++)
    delete it->second;
}

void SubstitutionProcessSequenceSimulator::outputInternalSequences(bool yn)
{
  for (std::map<size_t, SimpleSubstitutionProcessSequenceSimulator*>::iterator  it=mProcess_.begin(); it != mProcess_.end(); it++)
    it->second->outputInternalSequences(yn);
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


SiteContainer* SubstitutionProcessSequenceSimulator::simulate(size_t numberOfSites) const
{
  VectorSiteContainer* sites = new VectorSiteContainer(seqNames_.size(), getAlphabet());
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

SiteContainer* SubstitutionProcessSequenceSimulator::simulate(const vector<double>& rates) const
{
  size_t numberOfSites=rates.size();
  
  if (numberOfSites>vMap_.size())
    throw Exception("SubstitutionProcessSequenceSimulator::simulate : some sites do not have attributed process");

  VectorSiteContainer* sites = new VectorSiteContainer(seqNames_.size(), getAlphabet());
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

SiteContainer* SubstitutionProcessSequenceSimulator::simulate(const vector<size_t>& states) const
{
  size_t numberOfSites=states.size();

  VectorSiteContainer* sites = new VectorSiteContainer(seqNames_.size(), getAlphabet());
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

SiteContainer* SubstitutionProcessSequenceSimulator::simulate(const vector<double>& rates, const vector<size_t>& states) const
{
  size_t numberOfSites = rates.size();
  if (states.size() != numberOfSites)
    throw Exception("SubstitutionProcessSequenceSimulator::simulate, 'rates' and 'states' must have the same length.");

  VectorSiteContainer* sites = new VectorSiteContainer(seqNames_.size(), getAlphabet());
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
    

