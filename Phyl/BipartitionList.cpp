//
// File: BipartitionList.cpp
// Created by: Nicolas Galtier and Julien Dutheil
// Created on: Tue Apr 13 15:09 2007
//

/*
Copyright or © or Copr. CNRS, (November 16, 2004)

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

#include "BipartitionList.h"
#include "BipartitionTools.h"

#include "TreeTemplate.h"
#include "iotree"

// From Utils:
#include <Utils/exceptions>
#include <Utils/TextTools.h>
#include <Utils/FileTools.h>

// From SeqLib:
#include <Seq/alphabets>
#include <Seq/containers>
#include <Seq/ioseq>

// From NumCalc: 
#include <NumCalc/random>

using namespace bpp;

// From the STL:
#include <iostream>
#include <climits> //defines CHAR_BIT

using namespace std;

/****************************************************************/
/* utilitary classes required for sorting elements/bipartitions */
/****************************************************************/

class StringAndInt
{
  public:
    unsigned int ind;
    string str;
};

bool operator < (StringAndInt sai1, StringAndInt sai2)
{
  if(sai1.str < sai2.str) return true;
  return false;
}

/******************************************************************************/

class IntAndInt
{
  public:
    unsigned int ind;
    int val;
};

bool operator < (IntAndInt iai1, IntAndInt iai2)
{
  if(iai1.val < iai2.val) return true;
  return false;
}



/******************************************************************************/

BipartitionList::BipartitionList(const Tree& tr, bool sorted, vector<int>* index):
  sorted_(sorted)
{
  unsigned int nbbip;

  elements_ = tr.getLeavesNames();

  if(tr.isRooted())
    nbbip = tr.getNumberOfNodes() - 2;
  else
    nbbip = tr.getNumberOfNodes() - 1;

  if(sorted) std::sort(elements_.begin(), elements_.end());

  unsigned int lword  = BipartitionTools::LWORD;
  unsigned int nbword = (elements_.size() + lword - 1) / lword;
  unsigned int nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  for(unsigned int i = 0; i < nbbip; i++)
  {
    bitBipartitionList_.push_back(new int[nbint]);
    for(unsigned int j = 0; j < nbint; j++)
    {
      bitBipartitionList_[i][j] = 0;
    }
  }

  unsigned int cpt = 0;
  vector<string> underlyingNames;
  const Tree* tree = &tr;
  const TreeTemplate<Node>* ttree = dynamic_cast<const TreeTemplate<Node> *>(tree);
  if (ttree)
  {
    //Gain some time...
    buildBitBipartitions(ttree->getRootNode(), bitBipartitionList_, elements_, &cpt, index);
  }
  else
  {
    TreeTemplate<Node> tmp(tr);
    buildBitBipartitions(tmp.getRootNode(), bitBipartitionList_, elements_, &cpt, index);
  }
}

/******************************************************************************/

BipartitionList::BipartitionList(const vector<string>& elements, const vector<int*>& bitBipL)
{
  unsigned int lword  = BipartitionTools::LWORD;
  unsigned int nbword = (elements.size() + lword - 1) / lword;
  unsigned int nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  for(unsigned int i = 0; i < bitBipL.size(); i++)
  {
    bitBipartitionList_.push_back(new int[nbint]);
    for(unsigned int j = 0; j < nbint; j++)
    {
      bitBipartitionList_[i][j] = bitBipL[i][j];
    }
  }

  elements_ = elements;

  vector<string> cpelements_ = elements;
  std::sort(cpelements_.begin(), cpelements_.end());
  if(cpelements_ == elements) sorted_=true; else sorted_=false;
}

/******************************************************************************/

BipartitionList::BipartitionList(const BipartitionList& bipL)
{

  unsigned int lword  = BipartitionTools::LWORD;
  unsigned int nbword = (bipL.getNumberOfElements() + lword - 1) / lword;
  unsigned int nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  bitBipartitionList_.resize(bipL.getNumberOfBipartitions());
  vector<int*> bitBipL = bipL.getBitBipartitionList();
  for(unsigned int i = 0; i < bipL.getNumberOfBipartitions(); i++)
  {
    bitBipartitionList_[i] = new int[nbint];
    for(unsigned int j = 0; j < nbint; j++)
    {
      bitBipartitionList_[i][j] = bitBipL[i][j];
    }
  }

  elements_ = bipL.elements_;
  sorted_ = bipL.sorted_;
}

/******************************************************************************/

BipartitionList & BipartitionList::operator=(const BipartitionList& bipL)
{
  unsigned int lword  = BipartitionTools::LWORD;
  unsigned int nbword = (bipL.getNumberOfElements() + lword - 1) / lword;
  unsigned int nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  for(unsigned int i = 0; i < bitBipartitionList_.size(); i++)
    delete[] bitBipartitionList_[i];
  bitBipartitionList_.resize(bipL.getNumberOfBipartitions());
  vector<int*> bitBipL = bipL.getBitBipartitionList();
  for(unsigned int i = 0; i < bipL.getNumberOfBipartitions(); i++)
  {
    bitBipartitionList_[i] = new int[nbint];
    for(unsigned int j = 0; j < nbint; j++)
    {
      bitBipartitionList_[i][j] = bitBipL[i][j];
    }
  }

  elements_ = bipL.elements_;
  sorted_ = bipL.sorted_;
  return *this;
}

/******************************************************************************/

BipartitionList::~BipartitionList()
{
  for(unsigned int i = 0; i < bitBipartitionList_.size(); i++)
    delete[] bitBipartitionList_[i];
}

/******************************************************************************/

map<string, bool> BipartitionList::getBipartition(unsigned int i) const throw (Exception)
{
  map<string, bool> bip;

  if(i>=bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  for(unsigned int j = 0; j < elements_.size(); j++)
  {
    if(BipartitionTools::testBit(bitBipartitionList_[i], j))
      bip[elements_[j]] = true;
    else
      bip[elements_[j]] = false;
  }
  return bip;
}

/******************************************************************************/

int* BipartitionList::getBitBipartition(unsigned int i) throw (Exception)
{
  if(i >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  return bitBipartitionList_[i];
}

/******************************************************************************/

bool BipartitionList::haveSameElementsThan(map<string, bool>& bipart) const
{
  vector<string> elements = elements_;
  vector<string> keys;

  map<string,bool>::iterator it;

  for(it = bipart.begin(); it != bipart.end(); it++)
  {
    keys.push_back(it->first);
  }

  std::sort(elements.begin(), elements.end());
  std::sort(keys.begin(), keys.end());

  if(elements == keys) return true;
  return false;
}

/******************************************************************************/

void BipartitionList::addBipartition(map<string, bool>& bipart, bool checkElements) throw (Exception)
{
  if(checkElements && !BipartitionList::haveSameElementsThan(bipart))
    throw Exception("Distinct bipartition element sets");

  unsigned int lword  = BipartitionTools::LWORD;
  unsigned int nbword = (elements_.size() + lword - 1) / lword;
  unsigned int nbint  = nbword * lword / (CHAR_BIT * sizeof(int));
  bitBipartitionList_.push_back(new int[nbint]);
  unsigned int ind    = bitBipartitionList_.size() - 1;
  for(unsigned int j = 0; j < nbint; j++)
    bitBipartitionList_[ind][j] = 0;

  for(unsigned int i = 0; i < elements_.size(); i++)
  {
    if(bipart[elements_[i]] == true)
      BipartitionTools::bit1(bitBipartitionList_[ind], i);
    else
      BipartitionTools::bit0(bitBipartitionList_[ind], i);
  }
}

/******************************************************************************/

void BipartitionList::deleteBipartition(unsigned int i) throw(Exception)
{
  if(i >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  delete[] bitBipartitionList_[i];
  bitBipartitionList_.erase(bitBipartitionList_.begin()+i);
}

/******************************************************************************/

bool BipartitionList::containsBipartition(map<string, bool>& bipart, bool checkElements) const throw (Exception)
{
  unsigned int i, j;
  bool dac, padac;

  if(checkElements && !BipartitionList::haveSameElementsThan(bipart))
    throw Exception("Distinct bipartition element sets");

  for(i = 0; i < bitBipartitionList_.size(); i++)
  {
    dac = padac = false;
    for(j = 0; j < elements_.size(); j++)
    {
      if(BipartitionTools::testBit(bitBipartitionList_[i], j))
      {
        if(bipart[elements_[j]]) dac = true;
        else padac = true;
      }
      else
      {
        if(bipart[elements_[j]]) padac = true;
        else dac = true;
	    }
      if(dac && padac) break;
    }
    if(j == elements_.size())
      return true;
  }
  return false;
}

/******************************************************************************/

bool BipartitionList::areIdentical(unsigned int k1, unsigned int k2) const throw (Exception)
{
  bool dac, padac;

  if(k1 >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");
  if(k2 >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  dac = padac = false;
  for(unsigned int j = 0; j < elements_.size(); j++)
  {
    if(BipartitionTools::testBit(bitBipartitionList_[k1], j))
    {
      if(BipartitionTools::testBit(bitBipartitionList_[k2], j)) dac = true;
      else padac = true;
    }
    else
    {
      if(BipartitionTools::testBit(bitBipartitionList_[k2], j)) padac = true;
      else dac = true;
    }
    if(dac && padac) return false;
  }
  return true;
}

/******************************************************************************/

bool BipartitionList::areCompatible(unsigned int k1, unsigned int k2) const throw (Exception)
{
  bool uu, uz, zu, zz;

  if(k1 >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");
  if(k2 >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  uu = uz = zu = zz = false;

  for(unsigned int j = 0; j < elements_.size(); j++)
  {
    if(BipartitionTools::testBit(bitBipartitionList_[k1], j))
    {
      if(BipartitionTools::testBit(bitBipartitionList_[k2], j)) uu = true;
      else uz = true;
    }
    else
    {
      if(BipartitionTools::testBit(bitBipartitionList_[k2], j)) zu = true;
      else zz = true;
    }
    if(uu && uz && zu && zz) return false;
  }

  return true;
}

/******************************************************************************/

bool BipartitionList::areAllCompatible() const
{
  for(unsigned int i = 0; i < bitBipartitionList_.size(); i++)
  {
    for(unsigned int j = i + 1; j < bitBipartitionList_.size(); j++)
    {
      if(!BipartitionList::areCompatible(i, j))
        return false;
    }
  }
  return true;
}

/******************************************************************************/

bool BipartitionList::areAllCompatibleWith(map<string, bool> & bipart, bool checkElements) const throw (Exception)
{
  if(checkElements && !haveSameElementsThan(bipart))
    throw Exception("Distinct bipartition element sets");
  unsigned int nbBip = bitBipartitionList_.size();
  const_cast<BipartitionList *>(this)->addBipartition(bipart, false);

  for(unsigned int i = 0; i < nbBip; i++)
  {
    if(!areCompatible(i, nbBip))
    {
      const_cast<BipartitionList *>(this)->deleteBipartition(nbBip);
      return false;
    }
  }
  const_cast<BipartitionList *>(this)->deleteBipartition(nbBip);
  return true;
}

/******************************************************************************/

void BipartitionList::sortElements()
{
  vector<StringAndInt> relements_;
  StringAndInt sai;
  unsigned int nbbip;

  for(unsigned int i = 0; i < elements_.size(); i++)
  {
    sai.str = elements_[i];
    sai.ind = i;
    relements_.push_back(sai);
  }

  std::sort(relements_.begin(), relements_.end());

  for(unsigned int i = 0; i < elements_.size(); i++)
    elements_[i] = relements_[i].str;

  nbbip = bitBipartitionList_.size();
  bitBipartitionList_.resize(2 * nbbip);
  unsigned int lword  = BipartitionTools::LWORD;
  unsigned int nbword = (elements_.size() + lword - 1) / lword;
  unsigned int nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  for(unsigned int j = nbbip; j < 2 * nbbip; j++)
  {
	  bitBipartitionList_[j] = new int[nbint];
	  for(unsigned int k = 0; k < nbint; k++)
	    bitBipartitionList_[j][k] = 0;
    for(unsigned int i = 0; i < elements_.size(); i++)
    {
	    if(BipartitionTools::testBit(bitBipartitionList_[j - nbbip], relements_[i].ind))
	      BipartitionTools::bit1(bitBipartitionList_[j], i);
	    else
	      BipartitionTools::bit0(bitBipartitionList_[j], i);
	  }
  }

  for(unsigned int j = 0; j < nbbip; j++)
    delete[] bitBipartitionList_[j];

  bitBipartitionList_.erase(bitBipartitionList_.begin(), bitBipartitionList_.begin() + nbbip);
  sorted_ = true;
}

/******************************************************************************/

unsigned int BipartitionList::getPartitionSize(unsigned int k) const throw (Exception)
{
  unsigned int size = 0;
  if(k >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  for(unsigned int i = 0; i < elements_.size(); i++)
    if(BipartitionTools::testBit(bitBipartitionList_[k], i))
      size++;

  if(size<=elements_.size()/2) return size;
  else return elements_.size() - size;
}

/******************************************************************************/

void BipartitionList::removeTrivialBipartitions()
{
  unsigned int size = bitBipartitionList_.size();
  for(unsigned int i = size; i > 0; i--)
    if(BipartitionList::getPartitionSize(i - 1) < 2)
      BipartitionList::deleteBipartition(i - 1);
}

/******************************************************************************/

void BipartitionList::addTrivialBipartitions(bool checkExisting)
{
  map<string, bool> bip;

  for(unsigned int i = 0; i < elements_.size(); i++)
    bip[elements_[i]] = false;
  for(unsigned int i = 0; i < elements_.size(); i++)
  {
    bip[elements_[i]] = true;
    if(checkExisting && BipartitionList::containsBipartition(bip, false))
      continue;
    BipartitionList::addBipartition(bip, false);
    bip[elements_[i]] = false;
  }
}

/******************************************************************************/

void BipartitionList::sortByPartitionSize()
{
  vector<int*> sortedBitBipL;
  vector<IntAndInt> iaiVec;
  IntAndInt iai;

  for(unsigned int i = 0; i < bitBipartitionList_.size(); i++)
  {
    iai.ind = i;
    iai.val = BipartitionList::getPartitionSize(i);
    iaiVec.push_back(iai);
  }

  std::sort(iaiVec.begin(), iaiVec.end());

  for(unsigned int i = 0; i < bitBipartitionList_.size(); i++)
    sortedBitBipL.push_back(bitBipartitionList_[iaiVec[i].ind]);

  bitBipartitionList_=sortedBitBipL;
}

/******************************************************************************/

void BipartitionList::flip(unsigned int k) throw (Exception)
{
  if(k >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");
  unsigned int lword = BipartitionTools::LWORD;
  unsigned int nbword= (elements_.size() + lword - 1) / lword;
  unsigned int nbint = nbword * lword / (CHAR_BIT * sizeof(int));
  int* flipbip = new int[nbint];
  for(unsigned int i = 0; i < nbint; i++) flipbip[i] = 0;
  BipartitionTools::bitNot(flipbip, bitBipartitionList_[k], nbint);
  delete[] bitBipartitionList_[k];
  bitBipartitionList_[k] = flipbip;
}

/******************************************************************************/

void BipartitionList::removeRedundantBipartitions()
{
  bool deletion = true;

  while(deletion)
  {
	  deletion = false;
    for(unsigned int i = 0; i < bitBipartitionList_.size(); i++)
    {
      for(unsigned int j = i + 1; j < bitBipartitionList_.size(); j++)
      {
        if(BipartitionList::areIdentical(i, j))
        {
          BipartitionList::deleteBipartition(j);
          deletion = true;
          break;
        }
      }
      if(deletion) break;
    }
  }
}

/******************************************************************************/

TreeTemplate<Node>* BipartitionList::toTree() const throw (Exception)
{
  BipartitionList* sortedBipL;
  vector<int*> sortedBitBipL;
  int* bip;
  vector<Node*> vecNd, sonNd;
  vector<bool> alive;
  unsigned int lword, nbword, nbint, ii;

  /* check, copy and prepare bipartition list */

  if(!BipartitionList::areAllCompatible())
    throw Exception("Trying to build a tree from incompatible bipartitions");

  sortedBipL = dynamic_cast<BipartitionList *>(clone());
  for(unsigned int i = 0; i < sortedBipL->getNumberOfBipartitions(); i++)
  {
    if(sortedBipL->getPartitionSize(i) > sortedBipL->getNumberOfElements() / 2)
      sortedBipL->flip(i);
  }
  sortedBipL->sortByPartitionSize();
  sortedBipL->removeRedundantBipartitions();
  sortedBitBipL = sortedBipL->getBitBipartitionList();

  for(unsigned int i = 0; i < sortedBipL->getNumberOfBipartitions(); i++)
    alive.push_back(true);
  vecNd.resize(sortedBipL->getNumberOfBipartitions() + 1);
  lword  = BipartitionTools::LWORD;
  nbword = (elements_.size() + lword - 1) / lword;
  nbint  = nbword * lword / (CHAR_BIT * sizeof(int));
  bip    = new int[1]; bip[0]=0;

  /* main loop: create one node per bipartition */
  for(unsigned int i = 0; i < sortedBipL->getNumberOfBipartitions(); i++)
  {
    if(sortedBipL->getPartitionSize(i) == 1)
    { //terminal
      for(unsigned int j = 0; j < sortedBipL->getNumberOfElements(); j++)
      {
        if(BipartitionTools::testBit(sortedBitBipL[i], j))
        {
          vecNd[i]=new Node(elements_[j]);
          break;
        }
      }
    }
    else
    { //internal
	    sonNd.clear();
      for(unsigned int j = 0; j < i; j++)
      {
        if(alive[j])
        {
          for(ii = 0; ii < nbint; ii++)
          {
			      BipartitionTools::bitOr(bip, sortedBitBipL[j] + ii, sortedBitBipL[i] + ii, 1);
            if(bip[0] != sortedBitBipL[i][ii]) break;
          }
          if(ii == nbint)
          {
            sonNd.push_back(vecNd[j]);
            alive[j] = false;
          }
        }
      }
      vecNd[i] = new Node();
      for(unsigned int k = 0; k < sonNd.size(); k++)
        vecNd[i]->addSon(sonNd[k]);
    }
  }
 
  /* create last node, which joins alive bipartitions = fatherless nodes */
  Node* rootNd = new Node();
  for(unsigned int i = 0; i < sortedBipL->getNumberOfBipartitions(); i++)
    if(alive[i])
      rootNd->addSon(vecNd[i]);

  /* construct tree and return */
  TreeTemplate<Node>* tr = new TreeTemplate<Node>(rootNd);
  tr->resetNodesId();
  delete sortedBipL;
  return tr;
}

/******************************************************************************/

vector<string> BipartitionList::buildBitBipartitions(const Node* nd, vector<int*>& bitbip, const vector<string>& elements, unsigned int* cpt, vector<int>* index) const
{
  vector<string> underelements_, retelements_;

  if(nd->getNumberOfSons() == 0)
    underelements_.push_back(nd->getName());

  for(unsigned int i = 0; i < nd->getNumberOfSons(); i++)
  {
    retelements_ = BipartitionList::buildBitBipartitions(nd->getSon(i), bitbip, elements, cpt, index);
    for(unsigned int j = 0; j < retelements_.size(); j++)
      underelements_.push_back(retelements_[j]);
  }

  if(!nd->hasFather()) return underelements_; // root node

  if(!nd->getFather()->hasFather())
  {
    unsigned int nbrootson = nd->getFather()->getNumberOfSons();
    if(nbrootson == 2 && nd == nd->getFather()->getSon(1))
    return underelements_; //son 2 of root node when root node has 2 sons
  }

  bool ones;
  if(underelements_.size() <= elements.size() / 2) ones = true;
  else ones = false;

  for(unsigned int i = 0; i < elements.size(); i++)
  {
    if(ones) BipartitionTools::bit0(bitbip[*cpt], i);
    else BipartitionTools::bit1(bitbip[*cpt], i);
  }

  for(unsigned int i = 0; i < underelements_.size(); i++)
  {
    unsigned int taxa_ind = 0;
    while(underelements_[i] != elements[taxa_ind])
      taxa_ind++;
	  if(ones) BipartitionTools::bit1(bitbip[*cpt], taxa_ind);
	  else BipartitionTools::bit0(bitbip[*cpt], taxa_ind);
  }

  (*cpt)++;

  if(index) index->push_back(nd->getId());

  return underelements_;
}

/******************************************************************************/

RowMatrix<int> BipartitionList::toMatrix() const
{
  vector< map<string, bool> > bipl;
  for(unsigned int i = 0; i < getNumberOfBipartitions(); i++)
  {
    bipl.push_back(getBipartition(i));
  }

  vector<string> el = getElementNames();
  
  RowMatrix<int> mat(el.size(), getNumberOfBipartitions());

  for(unsigned int j = 0; j < el.size(); j++)
  {
    for(unsigned int i = 0; i < getNumberOfBipartitions(); i++)
    {
      mat(j, i) = bipl[i][el[j]];
    }
  }
  return mat;
}

/******************************************************************************/

