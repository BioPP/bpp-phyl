//
// File: MixedSubstitutionModelSet.cpp
// Created by: Laurent Guéguen
// Created on: mercredi 25 mai 2011, à 22h 12
//

/*
   Copyright or <A9> or Copr. CNRS, (November 16, 2004)

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

using namespace bpp;
using namespace std;

MixedSubstitutionModelSet::MixedSubstitutionModelSet(const MixedSubstitutionModelSet& set) :
  SubstitutionModelSet(set), vpHyperNodes_()
{
  for (unsigned int i=0;i<set.vpHyperNodes_.size();i++)
    vpHyperNodes_.push_back(new HyperNode(*set.vpHyperNodes_[i]));  
}

MixedSubstitutionModelSet::~MixedSubstitutionModelSet() 
{
  for (unsigned int i=0;i<vpHyperNodes_.size();i++)
    delete vpHyperNodes_[i];
}

void MixedSubstitutionModelSet::clear() 
{
  SubstitutionModelSet::clear();
  for (unsigned int i=0;i<vpHyperNodes_.size();i++)
    delete vpHyperNodes_[i];
}

MixedSubstitutionModelSet& MixedSubstitutionModelSet::operator=(const MixedSubstitutionModelSet& set)
{
  SubstitutionModelSet::operator=(set);
  for (unsigned int i=0;i<vpHyperNodes_.size();i++)
    if (vpHyperNodes_[i]!=NULL)
      delete vpHyperNodes_[i];
  vpHyperNodes_.clear();
  
  for (unsigned int i=0;i<set.vpHyperNodes_.size();i++)
    vpHyperNodes_.push_back(new HyperNode(*set.vpHyperNodes_[i]));
  
  return *this;
}

void MixedSubstitutionModelSet::addHyperNode()
{
  vpHyperNodes_.push_back(new HyperNode(this));
}

void MixedSubstitutionModelSet::addToHyperNode(int nM, const Vint& vnS, int nH)
{
  if (nH>=(int)vpHyperNodes_.size())
    throw BadIntegerException("Bad HyperNode number",nH);
  if (nH<0)
    nH=vpHyperNodes_.size()-1;
  
  if (nM<0 || nM>=(int)getNumberOfModels())
    throw BadIntegerException("Bad Mixed Model number",nM);

  vpHyperNodes_[nH]->addToModel(nM,vnS);
}

/**********************************************************/
/*************** HYPERNODE ********************************/
/***********************************************************/


MixedSubstitutionModelSet::HyperNode::HyperNode(const MixedSubstitutionModelSet *pMSMS): vNumbers_(pMSMS->getNumberOfModels())
{}


MixedSubstitutionModelSet::HyperNode::HyperNode(const HyperNode& hn): vNumbers_(hn.vNumbers_)
{}

MixedSubstitutionModelSet::HyperNode& MixedSubstitutionModelSet::HyperNode::operator=(const MixedSubstitutionModelSet::HyperNode& hn)
{
  vNumbers_.clear();
  vNumbers_.resize(hn.vNumbers_.size());
  for (unsigned int i=0;i<hn.vNumbers_.size();i++)
    vNumbers_[i]=hn.vNumbers_[i];
  
  return *this;
}

void MixedSubstitutionModelSet::HyperNode::addToModel(int nM, const Vint& vnS)
{
  if ((nM<0) || (nM>=(int)vNumbers_.size()))
    throw BadIntegerException("Bad Mixed model",nM);
      
  vNumbers_[nM].insertN(vnS);
}

bool MixedSubstitutionModelSet::HyperNode::operator<=(const HyperNode& hn) const
{
  for (unsigned int i=0;i<vNumbers_.size();i++){
    if (!( vNumbers_[i]<=hn.vNumbers_[i]))
      return false;
  }

  return true;
}

bool MixedSubstitutionModelSet::HyperNode::intersects(const HyperNode& hn) const
{
  for (unsigned int i=0;i<vNumbers_.size();i++){
    if (vNumbers_[i].intersects(hn.vNumbers_[i]))
      return true;
  }

  return false;
}

bool MixedSubstitutionModelSet::HyperNode::operator>=(const HyperNode& hn) const
{
  return hn>=*this;
}

/**********************************************************/
/******************** NODE ********************************/
/***********************************************************/

void MixedSubstitutionModelSet::HyperNode::Node::insertN(const Vint& vn)
{
  vector<int>::iterator it;
  vector<int>::const_iterator it2;
  
  for (it2=vn.begin();it2!=vn.end();it2++){
    for (it=vNumb_.begin();it!=vNumb_.end();it++)
      if (*it>=*it2)
        break;
    if (it==vNumb_.end())
        vNumb_.push_back(*it2);
    else if (*it!=*it2)
      vNumb_.insert(it,*it2);
  }
}

bool MixedSubstitutionModelSet::HyperNode::Node::operator<=(const Node& n) const
{
  vector<int>::const_iterator it(vNumb_.begin());
  vector<int>::const_iterator it2(n.vNumb_.begin());
  
  for (;it!=vNumb_.end();it++){
    while (it2!=n.vNumb_.end()  && (*it2<*it))
      it2++;
    if (it2==n.vNumb_.end() || (*it2>*it))
      return false;
    it++;
  }
  return true;
}

bool MixedSubstitutionModelSet::HyperNode::Node::operator>=(const Node& n) const
{
  return n<=*this;
}

bool MixedSubstitutionModelSet::HyperNode::Node::intersects(const Node& n) const
{
  vector<int>::const_iterator it(vNumb_.begin());
  vector<int>::const_iterator it2(n.vNumb_.begin());
  
  for (;it!=vNumb_.end();it++){
    while (it2!=n.vNumb_.end()  && (*it2<*it))
      it2++;
    if (*it2==*it)
      return true;
    if (it2==n.vNumb_.end())
      return false;
    it++;
  }
  return false;
}

