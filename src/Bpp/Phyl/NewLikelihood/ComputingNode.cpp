
 // File: ComputingNode.cpp
 // Created by: Laurent Guéguen
 // Created on: mercredi 3 juillet 2013, à 00h 10


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

#include "ComputingNode.h"
#include <Bpp/Numeric/Constraints.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>

using namespace bpp;
using namespace std;

ComputingNode::ComputingNode(const SubstitutionModel* model) :
  Node(),
  AbstractParametrizable(""),
  model_(model),
  nbStates_(model->getNumberOfStates()),
  scale_(1),
  probabilities_(),
  probabilitiesD1_(),
  probabilitiesD2_(),
  computeProbabilities_(true),
  computeProbabilitiesD1_(true),
  computeProbabilitiesD2_(true),
  vLogStates_(model->getNumberOfStates())
{
  addParameter_(new Parameter("scale", 1, &Parameter::R_PLUS_STAR));
}

ComputingNode::ComputingNode(int num, string st):
  Node(num, st),
  AbstractParametrizable(""),
  model_(0),
  nbStates_(0),
  scale_(1),
  probabilities_(),
  probabilitiesD1_(),
  probabilitiesD2_(),
  computeProbabilities_(true),
  computeProbabilitiesD1_(true),
  computeProbabilitiesD2_(true),
  vLogStates_()
{
  addParameter_(new Parameter("scale", 1, &Parameter::R_PLUS_STAR));
}

ComputingNode::ComputingNode():
  Node(),
  AbstractParametrizable(""),
  model_(0),
  nbStates_(0),
  scale_(1),
  probabilities_(),
  probabilitiesD1_(),
  probabilitiesD2_(),
  computeProbabilities_(true),
  computeProbabilitiesD1_(true),
  computeProbabilitiesD2_(true),
  vLogStates_()
{
  addParameter_(new Parameter("scale", 1, &Parameter::R_PLUS_STAR));
}

ComputingNode::ComputingNode(const Node& cn) :
  Node(cn),
  AbstractParametrizable(""),
  model_(0),
  nbStates_(0),
  scale_(1),
  probabilities_(),
  probabilitiesD1_(),
  probabilitiesD2_(),
  computeProbabilities_(true),
  computeProbabilitiesD1_(true),
  computeProbabilitiesD2_(true),
  vLogStates_()
{
  addParameter_(new Parameter("scale", 1, &Parameter::R_PLUS_STAR));
}

ComputingNode::ComputingNode(const ComputingNode& cn) :
  Node(cn),
  AbstractParametrizable(cn),
  model_(cn.model_),
  nbStates_(cn.nbStates_),
  scale_(cn.scale_),
  probabilities_(cn.probabilities_),
  probabilitiesD1_(cn.probabilitiesD1_),
  probabilitiesD2_(cn.probabilitiesD2_),
  computeProbabilities_(cn.computeProbabilities_),
  computeProbabilitiesD1_(cn.computeProbabilitiesD1_),
  computeProbabilitiesD2_(cn.computeProbabilitiesD2_),
  vLogStates_(cn.vLogStates_)
{
}

ComputingNode& ComputingNode::operator=(const ComputingNode& cn)
{
  Node::operator=(cn);
  AbstractParametrizable::operator=(cn);

  model_=cn.model_;
  nbStates_=cn.nbStates_;
  
  scale_=cn.scale_;
  probabilities_=cn.probabilities_;
  probabilitiesD1_=cn.probabilitiesD1_;
  probabilitiesD2_=cn.probabilitiesD2_;
  computeProbabilities_=cn.computeProbabilities_;
  computeProbabilitiesD1_=cn.computeProbabilitiesD1_;
  computeProbabilitiesD2_=cn.computeProbabilitiesD2_;
  vLogStates_=cn.vLogStates_;

  return *this;
}

void ComputingNode::setSubstitutionModel(const SubstitutionModel* pSM)
{
  model_=pSM;

  if (model_)
    nbStates_=model_->getNumberOfStates();
  
  computeProbabilities_=true;
  computeProbabilitiesD1_=true;
  computeProbabilitiesD2_=true;
  vLogStates_.resize(nbStates_);
}

void ComputingNode::fireParameterChanged(const ParameterList& pl)
{
  scale_=getParameterValue("scale");
  
  computeProbabilities_=true;
  computeProbabilitiesD1_=true;
  computeProbabilitiesD2_=true;
}  

void ComputingNode::update(bool flag)
{
  computeProbabilities_=flag;
  computeProbabilitiesD1_=flag;
  computeProbabilitiesD2_=flag;
}  
  
void ComputingNode::computeTransitionProbabilities() const
{
  if (computeProbabilities_) {
    computeProbabilities_ = false; //We record that we did this computation.
    probabilities_ = model_->getPij_t(scale_*getDistanceToFather());
  }
}

void ComputingNode::computeTransitionProbabilitiesD1() const
{
  if (computeProbabilitiesD1_) {
    computeProbabilitiesD1_ = false; //We record that we did this computation.
    //The transition matrix was never computed before. We therefore have to compute it first:
    probabilitiesD1_ = model_->getdPij_dt(scale_*getDistanceToFather());
    MatrixTools::scale(probabilitiesD1_,scale_);
  }
}

void ComputingNode::computeTransitionProbabilitiesD2() const
{
  if (computeProbabilitiesD2_) {
    computeProbabilitiesD2_ = false; //We record that we did this computation.
    //The transition matrix was never computed before. We therefore have to compute it first:
    probabilitiesD2_ = model_->getd2Pij_dt2(scale_*getDistanceToFather());
    MatrixTools::scale(probabilitiesD2_,scale_*scale_);
  }
}

void ComputingNode::updatedSubTreeNodes(Vint& lId) const
{
  if (!isUp2dateTransitionProbabilities())
    lId.push_back(getId());

  if (!isLeaf()){
    size_t nS=getNumberOfSons();
    for (size_t i=0; i<nS; i++)
      getSon(i)->updatedSubTreeNodes(lId);
  }
}
