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

// From SeqLib:
#include <Bpp/Seq/Container/VectorSiteContainer.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

SimpleSubstitutionProcessSequenceSimulator::SimpleSubstitutionProcessSequenceSimulator(
  const SubstitutionProcess& process) throw (Exception) :
  process_(&process),
  alphabet_(process_->getTransitionModel(process_->getParametrizablePhyloTree().getNodeIndex(process_->getParametrizablePhyloTree().getOutgoingNeighbors(process_->getParametrizablePhyloTree().getRoot())[0]),0).getAlphabet()),
  supportedStates_(process_->getTransitionModel(process_->getParametrizablePhyloTree().getNodeIndex(process_->getParametrizablePhyloTree().getOutgoingNeighbors(process_->getParametrizablePhyloTree().getRoot())[0]),0).getAlphabetStates()),
  phyloTree_(&process_->getParametrizablePhyloTree()),
  tree_(process_->getParametrizablePhyloTree()),
  leaves_(tree_.getAllLeaves()),
  seqNames_(),
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
  std::vector<std::shared_ptr<SimProcessNode> > vCN=tree_.getAllNodes();
    
  for (size_t j=0; j<vCN.size(); j++)
    vCN[j]->updateTree(&tree_, tree_.getNodeIndex(vCN[j]));

  // set sequence names

  if (outputInternalSequences_) {
    seqNames_.resize(vCN.size());
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      if (vCN[i]->isLeaf()) {
        seqNames_[i] = vCN[i]->getName();
      }
      else {
        seqNames_[i] = TextTools::toString(vCN[i]->getId() );
      }
    }
  }
  else {
    seqNames_.resize(leaves_.size());
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      seqNames_[i] = leaves_[i]->getName();
    }
  }
  
  // Initialize cumulative pxy:
  vector<shared_ptr<SimProcessNode> > nodes = tree_.getAllNodes();
  nodes.erase(std::find(nodes.begin(), nodes.end(), tree_.getRoot()));
  
  nbNodes_ = nodes.size();

  for (size_t i = 0; i < nodes.size(); i++)
  {
    shared_ptr<SimProcessNode> node = nodes[i];

    node->process_ = process_;
    VVVdouble* cumpxy_node_ = &node->cumpxy;
    cumpxy_node_->resize(nbClasses_);
    for (size_t c = 0; c < nbClasses_; c++)
    {
      VVdouble* cumpxy_node_c_ = &(*cumpxy_node_)[c];

      cumpxy_node_c_->resize(nbStates_);

      // process transition probabilities already consider rates &
      // branch length
      
      const RowMatrix<double>& P = process_->getTransitionProbabilities(node->getId(),c);

      for (size_t x = 0; x < nbStates_; x++)
      {
        Vdouble* cumpxy_node_c_x_ = &(*cumpxy_node_c_)[x];
        cumpxy_node_c_x_->resize(nbStates_);
        (*cumpxy_node_c_x_)[0] = P(x, 0);
        for (size_t y = 1; y < nbStates_; y++)
        {
          (*cumpxy_node_c_x_)[y] = (*cumpxy_node_c_x_)[y - 1] + P(x, y);
        }
      }
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

/******************************************************************************/

Site* SimpleSubstitutionProcessSequenceSimulator::simulateSite(size_t ancestralStateIndex) const
{
  if (continuousRates_ && process_->getRateDistribution())
  {
    // Draw a random rate:
    double rate = process_->getRateDistribution()->randC();
    // Make this state evolve:
    return simulateSite(ancestralStateIndex, rate);
  }
  else
  {
    size_t rateClass = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(nbClasses_);

    // Launch recursion:
    shared_ptr<SimProcessNode> root = tree_.getRoot();
    root->state = ancestralStateIndex;
    for (size_t i = 0; i < root->getNumberOfSons(); ++i)
    {
      evolveInternal(root->getSon(i), rateClass);
    }
    // Now create a Site object:
    Vint site(leaves_.size());
    for (size_t i = 0; i < leaves_.size(); ++i)
    {
      site[i] = process_->getTransitionModel(leaves_[i]->getId(), rateClass).getAlphabetStateAsInt(leaves_[i]->state);
    }
    return new Site(site, alphabet_);
  }
}


/******************************************************************************/

Site* SimpleSubstitutionProcessSequenceSimulator::simulateSite(size_t ancestralStateIndex, double rate) const
{
  size_t rateClass = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(nbClasses_);

  // Launch recursion:
  shared_ptr<SimProcessNode> root = tree_.getRoot();
  root->state = ancestralStateIndex;
  for (size_t i = 0; i < root->getNumberOfSons(); i++)
  {
    evolveInternal(root->getSon(i), rateClass, rate);
  }
  // Now create a Site object:
  size_t n = nbNodes_ + 1 ;
  if (! outputInternalSequences_) {
    n = leaves_.size() ;
  }
  Vint site(n);
  
  if (outputInternalSequences_) {
    vector<shared_ptr<SimProcessNode> > nodes = tree_.getAllNodes();
    for (size_t i = 0; i < n; i++)
    {
      size_t i2 = (i==n-1)?i-1:i; //  because the root has no model
      site[i] = process_->getTransitionModel(nodes[i2]->getId(),rateClass).getAlphabetStateAsInt(nodes[i2]->state);
    }
  }
  else {
    for (size_t i = 0; i < n; i++)
    {
      site[i] = process_->getTransitionModel(leaves_[i]->getId(),rateClass).getAlphabetStateAsInt(leaves_[i]->state);
    }
  }
  
  return new Site(site, alphabet_);
}

/******************************************************************************/

Site* SimpleSubstitutionProcessSequenceSimulator::simulateSite(double rate) const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  size_t ancestralStateIndex = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  vector<double> freqs = process_->getRootFrequencies();
  for (size_t i = 0; i < nbStates_; i++)
  {
    cumprob += freqs[i];
    if (r <= cumprob)
    {
      ancestralStateIndex = i;
      break;
    }
  }
  // Make this state evolve:
  return simulateSite(ancestralStateIndex, rate);
}

/******************************************************************************/

SiteContainer* SimpleSubstitutionProcessSequenceSimulator::simulate(size_t numberOfSites) const
{
  vector<size_t> ancestralStateIndices(numberOfSites, 0);
  for (size_t j = 0; j < numberOfSites; j++)
  {
    double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
    double cumprob = 0;
    const vector<double>& freqs = process_->getRootFrequencies();
    for (size_t i = 0; i < nbStates_; i++)
    {
      cumprob += freqs[i];
      if (r <= cumprob)
      {
        ancestralStateIndices[j] = i;
        break;
      }
    }
  }
  if (continuousRates_)
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
  else
  {
    // More efficient to do site this way:
    // Draw random rates:
    vector<size_t> rateClasses(numberOfSites);
    size_t nCat = nbClasses_;
    for (size_t j = 0; j < numberOfSites; j++)
    {
      rateClasses[j] = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(nCat);
    }
    // Make these states evolve:
    SiteContainer* sites = multipleEvolve(ancestralStateIndices, rateClasses);
    return sites;
  }
}

/******************************************************************************/

New_SiteSimulationResult* SimpleSubstitutionProcessSequenceSimulator::dSimulateSite() const
{
  // Draw an initial state randomly according to root frequencies:
  size_t ancestralStateIndex = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  vector<double> freqs = process_->getRootFrequencies();
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

/******************************************************************************/

New_SiteSimulationResult* SimpleSubstitutionProcessSequenceSimulator::dSimulateSite(size_t ancestralStateIndex) const
{
  // Draw a random rate:
  if (continuousRates_ && process_->getRateDistribution())
  {
    // Draw a random rate:
    double rate = process_->getRateDistribution()->randC();
    return dSimulateSite(ancestralStateIndex, rate);
  }
  else
  {
    size_t rateClass = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(nbClasses_);
    New_SiteSimulationResult* ssr = new New_SiteSimulationResult(phyloTree_, alphabet_, ancestralStateIndex);
    dEvolve(ancestralStateIndex, rateClass, *ssr);
    return ssr;
    // NB: this is more efficient than dSimulate(initialState, rDist_->rand())
  }
}

/******************************************************************************/

New_SiteSimulationResult* SimpleSubstitutionProcessSequenceSimulator::dSimulateSite(size_t ancestralStateIndex, double rate) const
{
  size_t rateClass = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(nbClasses_);

  // Make this state evolve:
  New_SiteSimulationResult* ssr = new New_SiteSimulationResult(phyloTree_, alphabet_, ancestralStateIndex);
  dEvolve(ancestralStateIndex, rateClass, rate, *ssr);
  return ssr;
}

/******************************************************************************/

New_SiteSimulationResult* SimpleSubstitutionProcessSequenceSimulator::dSimulateSite(double rate) const
{
  // Draw an initial state randomly according to equilibrum frequencies:
  size_t ancestralStateIndex = 0;
  double r = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double cumprob = 0;
  const vector<double>& freqs = process_->getRootFrequencies();
  for (size_t i = 0; i < nbStates_; i++)
  {
    cumprob += freqs[i];
    if (r <= cumprob)
    {
      ancestralStateIndex = i;
      break;
    }
  }
  return dSimulateSite(ancestralStateIndex, rate);
}




/******************************************************************************/

size_t SimpleSubstitutionProcessSequenceSimulator::evolve(const SimProcessNode* node, size_t initialStateIndex, size_t rateClass) const
{
  const Vdouble* cumpxy_node_c_x_ = &node->cumpxy[rateClass][initialStateIndex];
  double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  for (size_t y = 0; y < nbStates_; y++)
  {
    if (rand < (*cumpxy_node_c_x_)[y]) return y;
  }
  throw Exception("SimpleSubstitutionProcessSequenceSimulator::evolve. The impossible happened! rand = " + TextTools::toString(rand) + ".");
}

/******************************************************************************/

size_t SimpleSubstitutionProcessSequenceSimulator::evolve(const SimProcessNode* node, size_t initialStateIndex, size_t rateClass, double rate) const
{
  double cumpxy = 0;
  double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double l = rate * node->getDistanceToFather();
  
  const TransitionModel* model = &node->process_->getTransitionModel(node->getId(), rateClass);
  
  for (size_t y = 0; y < nbStates_; y++)
  {
    cumpxy += model->Pij_t(initialStateIndex, y, l);
    if (rand < cumpxy) return y;
  }
  MatrixTools::print(model->getPij_t(l));
  throw Exception("SimpleSubstitutionProcessSequenceSimulator::evolve. The impossible happened! rand = " + TextTools::toString(rand) + ".");
}

/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::multipleEvolve(
    const SimProcessNode* node,
    const std::vector<size_t>& initialStateIndices,
    const vector<size_t>& rateClasses,
    std::vector<size_t>& finalStateIndices) const
{
  const VVVdouble* cumpxy_node_ = &node->cumpxy;
  
  for (size_t i = 0; i < initialStateIndices.size(); i++)
  {
    const Vdouble* cumpxy_node_c_x_ = &(*cumpxy_node_)[rateClasses[i]][initialStateIndices[i]];
    double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
    for (size_t y = 0; y < nbStates_; y++)
    {
      if (rand < (*cumpxy_node_c_x_)[y])
      {
        finalStateIndices[i] = y;
        break;
      }
    }
  }
}

/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::evolveInternal(SimProcessNode* node, size_t rateClass) const
{
  if (!node->hasFather())
  {
    cerr << "DEBUG: SimpleSubstitutionProcessSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  node->state = evolve(node, node->getFather()->state, rateClass);
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    evolveInternal(node->getSon(i), rateClass);
  }
}

/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::evolveInternal(SimProcessNode* node, size_t rateClass, double rate) const
{
  if (!node->hasFather())
  {
    cerr << "DEBUG: SimpleSubstitutionProcessSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  node->state = evolve(node, node->getFather()->state, rateClass, rate);
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    evolveInternal(node->getSon(i), rateClass, rate);
  }
}

/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::multipleEvolveInternal(SimProcessNode* node, const vector<size_t>& rateClasses) const
{
  if (!node->hasFather())
  {
    cerr << "DEBUG: SimpleSubstitutionProcessSequenceSimulator::multipleEvolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  const vector<size_t>* initialStates = &node->getFather()->states;
  size_t n = initialStates->size();
  node->states.resize(n); // allocation.
  multipleEvolve(node, node->getFather()->states, rateClasses, node->states);
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    multipleEvolveInternal(node->getSon(i), rateClasses);
  }
}

/******************************************************************************/

SiteContainer* SimpleSubstitutionProcessSequenceSimulator::multipleEvolve(
    const std::vector<size_t>& initialStateIndices,
    const vector<size_t>& rateClasses) const
{
  // Launch recursion:
  shared_ptr<SimProcessNode> root = tree_.getRoot();

  root->states = initialStateIndices;
  for (size_t i = 0; i < root->getNumberOfSons(); i++)
  {
    multipleEvolveInternal(root->getSon(i), rateClasses);
  }

  // Now create a SiteContainer object:
  AlignedSequenceContainer* sites = new AlignedSequenceContainer(alphabet_);
  size_t n = leaves_.size();
  size_t nbSites = initialStateIndices.size();
  const TransitionModel* model = 0;

  if (outputInternalSequences_) {
    vector<shared_ptr<SimProcessNode> > nodes = tree_.getAllNodes();
    size_t nn = nbNodes_ + 1 ;
    for (size_t i = 0; i < nn; i++)
    {
      vector<int> content(nbSites);
      vector<size_t>& states = nodes[i]->states;

      size_t i2 = (i==nn-1)?i-1:i; // at the root, there is no model, so we take the model of node n-1.
      model = &nodes[i2]->process_->getTransitionModel(nodes[i2]->getId(), rateClasses[0]);

      for (size_t j = 0; j < nbSites; j++)
      {
        content[j] = model->getAlphabetStateAsInt(states[j]);
      }
      if (nodes[i]->isLeaf())
        sites->addSequence(BasicSequence(nodes[i]->getName(), content, alphabet_), false);
      else
        sites->addSequence(BasicSequence(TextTools::toString(nodes[i]->getId()), content, alphabet_), false);
    }
  }
  else
  {
    for (size_t i = 0; i < n; i++)
    {
      vector<int> content(nbSites);
      vector<size_t>& states = leaves_[i]->states;
      model = &leaves_[i]->process_->getTransitionModel(leaves_[i]->getId(), rateClasses[0]);
      
      for (size_t j = 0; j < nbSites; j++)
      {
        content[j] = model->getAlphabetStateAsInt(states[j]);
      }
      sites->addSequence(BasicSequence(leaves_[i]->getName(), content, alphabet_), false);
    }
  }
  
  return sites;
}

/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::dEvolve(size_t initialState, size_t rateClass, double rate, New_SiteSimulationResult& ssr) const
{
  // Launch recursion:
  shared_ptr<SimProcessNode> root = tree_.getRoot();
  root->state = initialState;
  for (size_t i = 0; i < root->getNumberOfSons(); i++)
  {
    dEvolveInternal(root->getSon(i), rateClass, rate, ssr);
  }
}


/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::dEvolve(size_t initialState, size_t rateClass, New_SiteSimulationResult& ssr) const
{
  // Launch recursion:
  shared_ptr<SimProcessNode> root = tree_.getRoot();
  root->state = initialState;
  for (size_t i = 0; i < root->getNumberOfSons(); i++)
  {
    dEvolveInternal(root->getSon(i), rateClass, ssr);
  }
}

/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::dEvolveInternal(SimProcessNode* node, size_t rateClass, double rate, New_SiteSimulationResult& ssr) const
{
  if (!node->hasFather())
  {
    cerr << "DEBUG: SimpleSubstitutionProcessSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  
  const TransitionModel* tm=&node->process_->getTransitionModel(node->getId(), rateClass);
  
  if (dynamic_cast<const SubstitutionModel*>(tm)==0)
    throw Exception("SimpleSubstitutionProcessSequenceSimulator::dEvolveInternal : detailed simulation not possible for non-markovian model");

  SimpleMutationProcess process(dynamic_cast<const SubstitutionModel*>(tm));
                                
  MutationPath mp = process.detailedEvolve(node->getFather()->state, node->getDistanceToFather() * rate);
  node->state = mp.getFinalState();

  // Now append infos in ssr:
  ssr.addNode(node->getId(), mp);

  // Now jump to son nodes:
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    dEvolveInternal(node->getSon(i), rateClass, rate, ssr);
  }
}

/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::dEvolveInternal(SimProcessNode* node, size_t rateClass, New_SiteSimulationResult& ssr) const
{
  if (!node->hasFather())
  {
    cerr << "DEBUG: SimpleSubstitutionProcessSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl;
    return;
  }

  const TransitionModel* tm=&node->process_->getTransitionModel(node->getId(), rateClass);
  
  if (dynamic_cast<const SubstitutionModel*>(tm)==0)
    throw Exception("SimpleSubstitutionProcessSequenceSimulator::dEvolveInternal : detailed simulation not possible for non-markovian model");

  SimpleMutationProcess process(dynamic_cast<const SubstitutionModel*>(tm));

  MutationPath mp = process.detailedEvolve(node->getFather()->state, node->getDistanceToFather());
  node->state = mp.getFinalState();

  // Now append infos in ssr:
  ssr.addNode(node->getId(), mp);

  // Now jump to son nodes:
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    dEvolveInternal(node->getSon(i), rateClass, ssr);
  }
}


void SimpleSubstitutionProcessSequenceSimulator::outputInternalSequences(bool yn) {
  outputInternalSequences_ = yn;
  if (outputInternalSequences_) {
    vector<shared_ptr<SimProcessNode> > nodes = tree_.getAllNodes();
    seqNames_.resize(nodes.size());
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      if (nodes[i]->isLeaf()) {
        seqNames_[i] = nodes[i]->getName();
      }
      else {
        seqNames_[i] = TextTools::toString( nodes[i]->getId() );
      }
    }
  }
  else {
    seqNames_.resize(leaves_.size());
    for (size_t i = 0; i < seqNames_.size(); i++)
    {
      seqNames_[i] = leaves_[i]->getName();
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
  resetSiteSimulators(numberOfSites);
  
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
  resetSiteSimulators(numberOfSites);
  
  if (numberOfSites>vMap_.size())
    throw Exception("SubstitutionProcessSequenceSimulator::simulate some sites do not have attributed process");

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
  resetSiteSimulators(numberOfSites);

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

  resetSiteSimulators(numberOfSites);

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
    

