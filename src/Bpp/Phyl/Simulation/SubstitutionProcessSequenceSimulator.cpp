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
  alphabet_(process_->getSubstitutionModel(0,0).getAlphabet()),
  supportedStates_(process_->getSubstitutionModel(0,0).getAlphabetStates()),
  templateTree_(&process_->getTree()),
  tree_(*&process_->getTree()),
  leaves_(tree_.getLeaves()),
  seqNames_(),
  nbNodes_(),
  nbClasses_(process_->getNumberOfClasses()),
  nbStates_(process_->getNumberOfStates()),
  continuousRates_(false)
{
  init();
}

/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::init()
{
  seqNames_.resize(leaves_.size());
  for (size_t i = 0; i < seqNames_.size(); i++)
  {
    seqNames_[i] = leaves_[i]->getName();
  }

  // Initialize cumulative pxy:
  vector<SPNode*> nodes = tree_.getNodes();
  nodes.pop_back(); // remove root
  nbNodes_ = nodes.size();

  for (size_t i = 0; i < nodes.size(); i++)
  {
    SPNode* node = nodes[i];
    node->getInfos().process_ = process_;
    VVVdouble* cumpxy_node_ = &node->getInfos().cumpxy;
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
    SPNode* root = tree_.getRootNode();
    root->getInfos().state = ancestralStateIndex;
    for (size_t i = 0; i < root->getNumberOfSons(); ++i)
    {
      evolveInternal(root->getSon(i), rateClass);
    }
    // Now create a Site object:
    Vint site(leaves_.size());
    for (size_t i = 0; i < leaves_.size(); ++i)
    {
      site[i] = process_->getSubstitutionModel(leaves_[i]->getId(), rateClass).getAlphabetStateAsInt(leaves_[i]->getInfos().state);
    }
    return new Site(site, alphabet_);
  }
}


/******************************************************************************/

Site* SimpleSubstitutionProcessSequenceSimulator::simulateSite(size_t ancestralStateIndex, double rate) const
{
  size_t rateClass = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(nbClasses_);

  // Launch recursion:
  SPNode* root = tree_.getRootNode();
  root->getInfos().state = ancestralStateIndex;
  for (size_t i = 0; i < root->getNumberOfSons(); i++)
  {
    evolveInternal(root->getSon(i), rateClass, rate);
  }
  // Now create a Site object:
  Vint site(leaves_.size());
  for (size_t i = 0; i < leaves_.size(); i++)
  {
    site[i] = process_->getSubstitutionModel(leaves_[i]->getId(),rateClass).getAlphabetStateAsInt(leaves_[i]->getInfos().state);
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

SiteSimulationResult* SimpleSubstitutionProcessSequenceSimulator::dSimulateSite() const
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

SiteSimulationResult* SimpleSubstitutionProcessSequenceSimulator::dSimulateSite(size_t ancestralStateIndex) const
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
    SiteSimulationResult* ssr = new SiteSimulationResult(templateTree_, alphabet_, ancestralStateIndex);
    dEvolve(ancestralStateIndex, rateClass, *ssr);
    return ssr;
    // NB: this is more efficient than dSimulate(initialState, rDist_->rand())
  }
}

/******************************************************************************/

SiteSimulationResult* SimpleSubstitutionProcessSequenceSimulator::dSimulateSite(size_t ancestralStateIndex, double rate) const
{
  size_t rateClass = RandomTools::giveIntRandomNumberBetweenZeroAndEntry<size_t>(nbClasses_);

  // Make this state evolve:
  SiteSimulationResult* ssr = new SiteSimulationResult(templateTree_, alphabet_, ancestralStateIndex);
  dEvolve(ancestralStateIndex, rateClass, rate, *ssr);
  return ssr;
}

/******************************************************************************/

SiteSimulationResult* SimpleSubstitutionProcessSequenceSimulator::dSimulateSite(double rate) const
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

size_t SimpleSubstitutionProcessSequenceSimulator::evolve(const SPNode* node, size_t initialStateIndex, size_t rateClass) const
{
  const Vdouble* cumpxy_node_c_x_ = &node->getInfos().cumpxy[rateClass][initialStateIndex];
  double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  for (size_t y = 0; y < nbStates_; y++)
  {
    if (rand < (*cumpxy_node_c_x_)[y]) return y;
  }
  throw Exception("SimpleSubstitutionProcessSequenceSimulator::evolve. The impossible happened! rand = " + TextTools::toString(rand) + ".");
}

/******************************************************************************/

size_t SimpleSubstitutionProcessSequenceSimulator::evolve(const SPNode* node, size_t initialStateIndex, size_t rateClass, double rate) const
{
  double cumpxy = 0;
  double rand = RandomTools::giveRandomNumberBetweenZeroAndEntry(1.);
  double l = rate * node->getDistanceToFather();
  
  const SubstitutionModel* model = &node->getInfos().process_->getSubstitutionModel(node->getId(), rateClass);
  
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
    const SPNode* node,
    const std::vector<size_t>& initialStateIndices,
    const vector<size_t>& rateClasses,
    std::vector<size_t>& finalStateIndices) const
{
  const VVVdouble* cumpxy_node_ = &node->getInfos().cumpxy;
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

void SimpleSubstitutionProcessSequenceSimulator::evolveInternal(SPNode* node, size_t rateClass) const
{
  if (!node->hasFather())
  {
    cerr << "DEBUG: SimpleSubstitutionProcessSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  node->getInfos().state = evolve(node, node->getFather()->getInfos().state, rateClass);
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    evolveInternal(node->getSon(i), rateClass);
  }
}

/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::evolveInternal(SPNode* node, size_t rateClass, double rate) const
{
  if (!node->hasFather())
  {
    cerr << "DEBUG: SimpleSubstitutionProcessSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  node->getInfos().state = evolve(node, node->getFather()->getInfos().state, rateClass, rate);
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    evolveInternal(node->getSon(i), rateClass, rate);
  }
}

/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::multipleEvolveInternal(SPNode* node, const vector<size_t>& rateClasses) const
{
  if (!node->hasFather())
  {
    cerr << "DEBUG: SimpleSubstitutionProcessSequenceSimulator::multipleEvolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  const vector<size_t>* initialStates = &node->getFather()->getInfos().states;
  size_t n = initialStates->size();
  node->getInfos().states.resize(n); // allocation.
  multipleEvolve(node, node->getFather()->getInfos().states, rateClasses, node->getInfos().states);
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
  SPNode* root = tree_.getRootNode();
  root->getInfos().states = initialStateIndices;
  for (size_t i = 0; i < root->getNumberOfSons(); i++)
  {
    multipleEvolveInternal(root->getSon(i), rateClasses);
  }
  // Now create a SiteContainer object:
  AlignedSequenceContainer* sites = new AlignedSequenceContainer(alphabet_);
  size_t n = leaves_.size();
  size_t nbSites = initialStateIndices.size();
  const SubstitutionModel* model = 0;
  for (size_t i = 0; i < n; i++)
  {
    vector<int> content(nbSites);
    vector<size_t>& states = leaves_[i]->getInfos().states;
    model = &leaves_[i]->getInfos().process_->getSubstitutionModel(leaves_[i]->getId(), rateClasses[0]);

    for (size_t j = 0; j < nbSites; j++)
    {
      content[j] = model->getAlphabetStateAsInt(states[j]);
    }
    sites->addSequence(BasicSequence(leaves_[i]->getName(), content, alphabet_), false);
  }
  return sites;
}

/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::dEvolve(size_t initialState, size_t rateClass, double rate, SiteSimulationResult& ssr) const
{
  // Launch recursion:
  SPNode* root = tree_.getRootNode();
  root->getInfos().state = initialState;
  for (size_t i = 0; i < root->getNumberOfSons(); i++)
  {
    dEvolveInternal(root->getSon(i), rateClass, rate, ssr);
  }
}


/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::dEvolve(size_t initialState, size_t rateClass, SiteSimulationResult& ssr) const
{
  // Launch recursion:
  SPNode* root = tree_.getRootNode();
  root->getInfos().state = initialState;
  for (size_t i = 0; i < root->getNumberOfSons(); i++)
  {
    dEvolveInternal(root->getSon(i), rateClass, ssr);
  }
}

/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::dEvolveInternal(SPNode* node, size_t rateClass, double rate, SiteSimulationResult& ssr) const
{
  if (!node->hasFather())
  {
    cerr << "DEBUG: SimpleSubstitutionProcessSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  
  SimpleMutationProcess process(&node->getInfos().process_->getSubstitutionModel(node->getId(), rateClass));
                                
  MutationPath mp = process.detailedEvolve(node->getFather()->getInfos().state, node->getDistanceToFather() * rate);
  node->getInfos().state = mp.getFinalState();

  // Now append infos in ssr:
  ssr.addNode(node->getId(), mp);

  // Now jump to son nodes:
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    dEvolveInternal(node->getSon(i), rateClass, rate, ssr);
  }
}

/******************************************************************************/

void SimpleSubstitutionProcessSequenceSimulator::dEvolveInternal(SPNode* node, size_t rateClass, SiteSimulationResult& ssr) const
{
  if (!node->hasFather())
  {
    cerr << "DEBUG: SimpleSubstitutionProcessSequenceSimulator::evolveInternal. Forbidden call of method on root node." << endl;
    return;
  }
  
  SimpleMutationProcess process(&node->getInfos().process_->getSubstitutionModel(node->getId(), rateClass));
                                
  MutationPath mp = process.detailedEvolve(node->getFather()->getInfos().state, node->getDistanceToFather());
  node->getInfos().state = mp.getFinalState();

  // Now append infos in ssr:
  ssr.addNode(node->getId(), mp);

  // Now jump to son nodes:
  for (size_t i = 0; i < node->getNumberOfSons(); i++)
  {
    dEvolveInternal(node->getSon(i), rateClass, ssr);
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
  
  seqNames_=evol.getSubstitutionProcess(nProc[0]).getTree().getLeavesNames();

  for (size_t i=0; i< nProc.size(); i++)
  {
    const SubstitutionProcess& sp=evol.getSubstitutionProcess(nProc[i]);
    
    mProcess_[nProc[i]]=new SimpleSubstitutionProcessSequenceSimulator(sp);

    const vector<string>& seqNames2=sp.getTree().getLeavesNames();
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
    seqNames_=mSP.begin()->second.getTree().getLeavesNames();
  

  for (std::map<size_t, const SubstitutionProcess&>::const_iterator  it=mSP.begin(); it != mSP.end(); it++)
  {
    mProcess_[it->first]=new SimpleSubstitutionProcessSequenceSimulator(it->second);

    const vector<string>& seqNames2=it->second.getTree().getLeavesNames();
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
  
  seqNames_=spc.getSubstitutionProcess(procN[0]).getTree().getLeavesNames();

  for (size_t i=0; i<procN.size(); i++)
  {
    const SubstitutionProcess& sp=spc.getSubstitutionProcess(procN[i]);
    
    mProcess_[procN[i]]=new SimpleSubstitutionProcessSequenceSimulator(sp);

    const vector<string>& seqNames2=sp.getTree().getLeavesNames();
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

    

