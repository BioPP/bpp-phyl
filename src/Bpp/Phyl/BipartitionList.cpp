//
// File: BipartitionList.cpp
// Created by: Nicolas Galtier and Julien Dutheil
// Created on: Tue Apr 13 15:09 2007
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

#include "BipartitionList.h"
#include "BipartitionTools.h"

#include "TreeTemplate.h"

#include <Bpp/Exceptions.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Numeric/Random/RandomTools.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>


using namespace bpp;

// From the STL:
#include <iostream>
#include <climits> // defines CHAR_BIT

using namespace std;

/****************************************************************/
/* utilitary classes required for sorting elements/bipartitions */
/****************************************************************/

class StringAndInt
{
public:
  int ind;
  string str;

public:
  StringAndInt() : ind(0),
    str() {}
};

bool operator<(StringAndInt sai1, StringAndInt sai2)
{
  if (sai1.str < sai2.str)
    return true;
  return false;
}

/******************************************************************************/

class IntAndInt
{
public:
  size_t ind;
  int val;
};

bool operator<(IntAndInt iai1, IntAndInt iai2)
{
  if (iai1.val < iai2.val)
    return true;
  return false;
}


/******************************************************************************/

BipartitionList::BipartitionList(const Tree& tr, bool sorted, std::vector<int>* index) :
  bitBipartitionList_(),
  elements_(),
  sorted_(sorted)
{
  size_t nbbip;

  elements_ = tr.getLeavesNames();

  if (tr.isRooted())
    nbbip = tr.getNumberOfNodes() - 2;
  else
    nbbip = tr.getNumberOfNodes() - 1;

  if (sorted)
    std::sort(elements_.begin(), elements_.end());

  size_t lword  = static_cast<size_t>(BipartitionTools::LWORD);
  size_t nbword = (elements_.size() + lword - 1) / lword;
  size_t nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  for (size_t i = 0; i < nbbip; i++)
  {
    bitBipartitionList_.push_back(new int[nbint]);
    for (size_t j = 0; j < nbint; j++)
    {
      bitBipartitionList_[i][j] = 0;
    }
  }

  size_t cpt = 0;
  vector<string> underlyingNames;
  const Tree* tree = &tr;
  const TreeTemplate<Node>* ttree = dynamic_cast<const TreeTemplate<Node>*>(tree);
  if (ttree)
  {
    // Gain some time...
    buildBitBipartitions(ttree->getRootNode(), bitBipartitionList_, elements_, &cpt, index);
  }
  else
  {
    TreeTemplate<Node> tmp(tr);
    buildBitBipartitions(tmp.getRootNode(), bitBipartitionList_, elements_, &cpt, index);
  }
}

/******************************************************************************/

BipartitionList::BipartitionList(
  const std::vector<std::string>& elements,
  const std::vector<int*>& bitBipL) :
  bitBipartitionList_(),
  elements_(elements),
  sorted_()
{
  size_t lword  = static_cast<size_t>(BipartitionTools::LWORD);
  size_t nbword = (elements.size() + lword - 1) / lword;
  size_t nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  for (size_t i = 0; i < bitBipL.size(); i++)
  {
    bitBipartitionList_.push_back(new int[nbint]);
    for (size_t j = 0; j < nbint; j++)
    {
      bitBipartitionList_[i][j] = bitBipL[i][j];
    }
  }

  vector<string> cpelements_ = elements;
  std::sort(cpelements_.begin(), cpelements_.end());
  if (cpelements_ == elements)
    sorted_ = true;
  else
    sorted_ = false;
}

/******************************************************************************/

BipartitionList::BipartitionList(const BipartitionList& bipL) :
  bitBipartitionList_(),
  elements_(bipL.elements_),
  sorted_(bipL.sorted_)
{
  size_t lword  = static_cast<size_t>(BipartitionTools::LWORD);
  size_t nbword = (bipL.getNumberOfElements() + lword - 1) / lword;
  size_t nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  bitBipartitionList_.resize(bipL.getNumberOfBipartitions());
  vector<int*> bitBipL = bipL.getBitBipartitionList();
  for (size_t i = 0; i < bipL.getNumberOfBipartitions(); i++)
  {
    bitBipartitionList_[i] = new int[nbint];
    for (size_t j = 0; j < nbint; j++)
    {
      bitBipartitionList_[i][j] = bitBipL[i][j];
    }
  }
}

/******************************************************************************/

BipartitionList& BipartitionList::operator=(const BipartitionList& bipL)
{
  size_t lword  = static_cast<size_t>(BipartitionTools::LWORD);
  size_t nbword = (bipL.getNumberOfElements() + lword - 1) / lword;
  size_t nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  for (size_t i = 0; i < bitBipartitionList_.size(); i++)
  {
    delete[] bitBipartitionList_[i];
  }
  bitBipartitionList_.resize(bipL.getNumberOfBipartitions());
  vector<int*> bitBipL = bipL.getBitBipartitionList();
  for (size_t i = 0; i < bipL.getNumberOfBipartitions(); i++)
  {
    bitBipartitionList_[i] = new int[nbint];
    for (size_t j = 0; j < nbint; j++)
    {
      bitBipartitionList_[i][j] = bitBipL[i][j];
    }
  }

  elements_ = bipL.elements_;
  sorted_   = bipL.sorted_;
  return *this;
}

/******************************************************************************/

BipartitionList::~BipartitionList()
{
  for (size_t i = 0; i < bitBipartitionList_.size(); i++)
  {
    delete[] bitBipartitionList_[i];
  }
}

/******************************************************************************/

map<string, bool> BipartitionList::getBipartition(size_t i) const throw (Exception)
{
  map<string, bool> bip;

  if (i >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  for (size_t j = 0; j < elements_.size(); j++)
  {
    if (BipartitionTools::testBit(bitBipartitionList_[i], static_cast<int>(j)))
      bip[elements_[j]] = true;
    else
      bip[elements_[j]] = false;
  }
  return bip;
}

/******************************************************************************/

int* BipartitionList::getBitBipartition(size_t i) throw (Exception)
{
  if (i >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  return bitBipartitionList_[i];
}

/******************************************************************************/

bool BipartitionList::haveSameElementsThan(map<string, bool>& bipart) const
{
  vector<string> elements = elements_;
  vector<string> keys;

  map<string, bool>::iterator it;

  for (it = bipart.begin(); it != bipart.end(); it++)
  {
    keys.push_back(it->first);
  }

  std::sort(elements.begin(), elements.end());
  std::sort(keys.begin(), keys.end());

  if (elements == keys)
    return true;
  return false;
}

/******************************************************************************/

void BipartitionList::addBipartition(map<string, bool>& bipart, bool checkElements) throw (Exception)
{
  if (checkElements && !BipartitionList::haveSameElementsThan(bipart))
    throw Exception("Distinct bipartition element sets");

  size_t lword  = static_cast<size_t>(BipartitionTools::LWORD);
  size_t nbword = (elements_.size() + lword - 1) / lword;
  size_t nbint  = nbword * lword / (CHAR_BIT * sizeof(int));
  bitBipartitionList_.push_back(new int[nbint]);
  size_t ind    = bitBipartitionList_.size() - 1;
  for (size_t j = 0; j < nbint; j++)
  {
    bitBipartitionList_[ind][j] = 0;
  }

  for (size_t i = 0; i < elements_.size(); i++)
  {
    if (bipart[elements_[i]] == true)
      BipartitionTools::bit1(bitBipartitionList_[ind], static_cast<int>(i));
    else
      BipartitionTools::bit0(bitBipartitionList_[ind], static_cast<int>(i));
  }
}

/******************************************************************************/

void BipartitionList::deleteBipartition(size_t i) throw (Exception)
{
  if (i >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  delete[] bitBipartitionList_[i];
  bitBipartitionList_.erase(bitBipartitionList_.begin() + static_cast<ptrdiff_t>(i));
}

/******************************************************************************/

bool BipartitionList::containsBipartition(map<string, bool>& bipart, bool checkElements) const throw (Exception)
{
  size_t i, j;
  bool dac, padac;

  if (checkElements && !BipartitionList::haveSameElementsThan(bipart))
    throw Exception("Distinct bipartition element sets");

  for (i = 0; i < bitBipartitionList_.size(); i++)
  {
    dac = padac = false;
    for (j = 0; j < elements_.size(); j++)
    {
      if (BipartitionTools::testBit(bitBipartitionList_[i], static_cast<int>(j)))
      {
        if (bipart[elements_[j]])
          dac = true;
        else
          padac = true;
      }
      else
      {
        if (bipart[elements_[j]])
          padac = true;
        else
          dac = true;
      }
      if (dac && padac)
        break;
    }
    if (j == elements_.size())
      return true;
  }
  return false;
}

/******************************************************************************/

bool BipartitionList::areIdentical(size_t k1, size_t k2) const throw (Exception)
{
  bool dac, padac;

  if (k1 >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");
  if (k2 >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  dac = padac = false;
  for (size_t j = 0; j < elements_.size(); j++)
  {
    if (BipartitionTools::testBit(bitBipartitionList_[k1], static_cast<int>(j)))
    {
      if (BipartitionTools::testBit(bitBipartitionList_[k2], static_cast<int>(j)))
        dac = true;
      else
        padac = true;
    }
    else
    {
      if (BipartitionTools::testBit(bitBipartitionList_[k2], static_cast<int>(j)))
        padac = true;
      else
        dac = true;
    }
    if (dac && padac)
      return false;
  }
  return true;
}

/******************************************************************************/

bool BipartitionList::areCompatible(size_t k1, size_t k2) const throw (Exception)
{
  bool uu, uz, zu, zz;

  if (k1 >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");
  if (k2 >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  uu = uz = zu = zz = false;

  for (size_t j = 0; j < elements_.size(); j++)
  {
    if (BipartitionTools::testBit(bitBipartitionList_[k1], static_cast<int>(j)))
    {
      if (BipartitionTools::testBit(bitBipartitionList_[k2], static_cast<int>(j)))
        uu = true;
      else
        uz = true;
    }
    else
    {
      if (BipartitionTools::testBit(bitBipartitionList_[k2], static_cast<int>(j)))
        zu = true;
      else
        zz = true;
    }
    if (uu && uz && zu && zz)
      return false;
  }

  return true;
}

/******************************************************************************/

bool BipartitionList::areAllCompatible() const
{
  for (size_t i = 0; i < bitBipartitionList_.size(); i++)
  {
    for (size_t j = i + 1; j < bitBipartitionList_.size(); j++)
    {
      if (!BipartitionList::areCompatible(i, j))
        return false;
    }
  }
  return true;
}

/******************************************************************************/

bool BipartitionList::areAllCompatibleWith(map<string, bool>& bipart, bool checkElements) const throw (Exception)
{
  if (checkElements && !haveSameElementsThan(bipart))
    throw Exception("Distinct bipartition element sets");
  size_t nbBip = bitBipartitionList_.size();
  const_cast<BipartitionList*>(this)->addBipartition(bipart, false);

  for (size_t i = 0; i < nbBip; i++)
  {
    if (!areCompatible(i, nbBip))
    {
      const_cast<BipartitionList*>(this)->deleteBipartition(nbBip);
      return false;
    }
  }
  const_cast<BipartitionList*>(this)->deleteBipartition(nbBip);
  return true;
}

/******************************************************************************/

void BipartitionList::sortElements()
{
  vector<StringAndInt> relements_;
  StringAndInt sai;
  size_t nbbip;

  for (size_t i = 0; i < elements_.size(); i++)
  {
    sai.str = elements_[i];
    sai.ind = static_cast<int>(i);
    relements_.push_back(sai);
  }

  std::sort(relements_.begin(), relements_.end());

  for (size_t i = 0; i < elements_.size(); i++)
  {
    elements_[i] = relements_[i].str;
  }

  nbbip = bitBipartitionList_.size();
  bitBipartitionList_.resize(2 * nbbip);
  size_t lword  = static_cast<size_t>(BipartitionTools::LWORD);
  size_t nbword = (elements_.size() + lword - 1) / lword;
  size_t nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  for (size_t j = nbbip; j < 2 * nbbip; j++)
  {
    bitBipartitionList_[j] = new int[nbint];
    for (size_t k = 0; k < nbint; k++)
    {
      bitBipartitionList_[j][k] = 0;
    }
    for (size_t i = 0; i < elements_.size(); i++)
    {
      if (BipartitionTools::testBit(bitBipartitionList_[j - nbbip], relements_[i].ind))
        BipartitionTools::bit1(bitBipartitionList_[j], static_cast<int>(i));
      else
        BipartitionTools::bit0(bitBipartitionList_[j], static_cast<int>(i));
    }
  }

  for (size_t j = 0; j < nbbip; j++)
  {
    delete[] bitBipartitionList_[j];
  }

  bitBipartitionList_.erase(bitBipartitionList_.begin(), bitBipartitionList_.begin() + static_cast<ptrdiff_t>(nbbip));
  sorted_ = true;
}

/******************************************************************************/

size_t BipartitionList::getPartitionSize(size_t k) const throw (Exception)
{
  size_t size = 0;
  if (k >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  for (size_t i = 0; i < elements_.size(); i++)
  {
    if (BipartitionTools::testBit(bitBipartitionList_[k], static_cast<int>(i)))
      size++;
  }

  if (size <= elements_.size() / 2)
    return size;
  else
    return elements_.size() - size;
}

/******************************************************************************/

void BipartitionList::removeTrivialBipartitions()
{
  size_t size = bitBipartitionList_.size();
  for (size_t i = size; i > 0; i--)
  {
    if (BipartitionList::getPartitionSize(i - 1) < 2)
      BipartitionList::deleteBipartition(i - 1);
  }
}

/******************************************************************************/

void BipartitionList::addTrivialBipartitions(bool checkExisting)
{
  map<string, bool> bip;

  for (size_t i = 0; i < elements_.size(); i++)
  {
    bip[elements_[i]] = false;
  }
  for (size_t i = 0; i < elements_.size(); i++)
  {
    bip[elements_[i]] = true;
    if (checkExisting && BipartitionList::containsBipartition(bip, false))
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

  for (size_t i = 0; i < bitBipartitionList_.size(); i++)
  {
    iai.ind = i;
    iai.val = static_cast<int>(BipartitionList::getPartitionSize(i));
    iaiVec.push_back(iai);
  }

  std::sort(iaiVec.begin(), iaiVec.end());

  for (size_t i = 0; i < bitBipartitionList_.size(); i++)
  {
    sortedBitBipL.push_back(bitBipartitionList_[iaiVec[i].ind]);
  }

  bitBipartitionList_ = sortedBitBipL;
}

/******************************************************************************/

void BipartitionList::flip(size_t k) throw (Exception)
{
  if (k >= bitBipartitionList_.size())
    throw Exception("Bipartition index exceeds BipartitionList size");
  size_t lword = static_cast<size_t>(BipartitionTools::LWORD);
  size_t nbword = (elements_.size() + lword - 1) / lword;
  size_t nbint = nbword * lword / (CHAR_BIT * sizeof(int));
  int* flipbip = new int[nbint];
  for (size_t i = 0; i < nbint; i++)
  {
    flipbip[i] = 0;
  }
  BipartitionTools::bitNot(flipbip, bitBipartitionList_[k], nbint);
  delete[] bitBipartitionList_[k];
  bitBipartitionList_[k] = flipbip;
}

/******************************************************************************/

void BipartitionList::removeRedundantBipartitions()
{
  bool deletion = true;

  while (deletion)
  {
    deletion = false;
    for (size_t i = 0; i < bitBipartitionList_.size(); i++)
    {
      for (size_t j = i + 1; j < bitBipartitionList_.size(); j++)
      {
        if (BipartitionList::areIdentical(i, j))
        {
          BipartitionList::deleteBipartition(j);
          deletion = true;
          break;
        }
      }
      if (deletion)
        break;
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
  size_t lword, nbword, nbint, ii;

  /* check, copy and prepare bipartition list */

  if (!BipartitionList::areAllCompatible())
    throw Exception("Trying to build a tree from incompatible bipartitions");

  sortedBipL = dynamic_cast<BipartitionList*>(clone());
  for (size_t i = 0; i < sortedBipL->getNumberOfBipartitions(); i++)
  {
    if (sortedBipL->getPartitionSize(i) > sortedBipL->getNumberOfElements() / 2)
      sortedBipL->flip(i);
  }
  sortedBipL->sortByPartitionSize();
  sortedBipL->removeRedundantBipartitions();
  sortedBitBipL = sortedBipL->getBitBipartitionList();

  for (size_t i = 0; i < sortedBipL->getNumberOfBipartitions(); i++)
  {
    alive.push_back(true);
  }
  vecNd.resize(sortedBipL->getNumberOfBipartitions() + 1);
  lword  = static_cast<size_t>(BipartitionTools::LWORD);
  nbword = (elements_.size() + lword - 1) / lword;
  nbint  = nbword * lword / (CHAR_BIT * sizeof(int));
  bip    = new int[1]; bip[0] = 0;

  /* main loop: create one node per bipartition */
  for (size_t i = 0; i < sortedBipL->getNumberOfBipartitions(); i++)
  {
    if (sortedBipL->getPartitionSize(i) == 1)
    { // terminal
      for (size_t j = 0; j < sortedBipL->getNumberOfElements(); j++)
      {
        if (BipartitionTools::testBit(sortedBitBipL[i], static_cast<int>(j)))
        {
          vecNd[i] = new Node(elements_[j]);
          break;
        }
      }
    }
    else
    { // internal
      sonNd.clear();
      for (size_t j = 0; j < i; j++)
      {
        if (alive[j])
        {
          for (ii = 0; ii < nbint; ii++)
          {
            BipartitionTools::bitOr(bip, sortedBitBipL[j] + ii, sortedBitBipL[i] + ii, 1);
            if (bip[0] != sortedBitBipL[i][ii])
              break;
          }
          if (ii == nbint)
          {
            sonNd.push_back(vecNd[j]);
            alive[j] = false;
          }
        }
      }
      vecNd[i] = new Node();
      for (size_t k = 0; k < sonNd.size(); k++)
      {
        vecNd[i]->addSon(sonNd[k]);
      }
    }
  }

  /* create last node, which joins alive bipartitions = fatherless nodes */
  Node* rootNd = new Node();
  for (size_t i = 0; i < sortedBipL->getNumberOfBipartitions(); i++)
  {
    if (alive[i])
      rootNd->addSon(vecNd[i]);
  }

  /* construct tree and return */
  TreeTemplate<Node>* tr = new TreeTemplate<Node>(rootNd);
  tr->resetNodesId();
  delete sortedBipL;
  return tr;
}

/******************************************************************************/

vector<string> BipartitionList::buildBitBipartitions(const Node* nd, vector<int*>& bitbip, const vector<string>& elements, size_t* cpt, vector<int>* index) const
{
  vector<string> underelements_, retelements_;

  if (nd->getNumberOfSons() == 0)
    underelements_.push_back(nd->getName());

  for (size_t i = 0; i < nd->getNumberOfSons(); i++)
  {
    retelements_ = BipartitionList::buildBitBipartitions(nd->getSon(i), bitbip, elements, cpt, index);
    for (size_t j = 0; j < retelements_.size(); j++)
    {
      underelements_.push_back(retelements_[j]);
    }
  }

  if (!nd->hasFather())
    return underelements_;  // root node

  if (!nd->getFather()->hasFather())
  {
    size_t nbrootson = nd->getFather()->getNumberOfSons();
    if (nbrootson == 2 && nd == nd->getFather()->getSon(1))
      return underelements_;  // son 2 of root node when root node has 2 sons
  }

  bool ones;
  if (underelements_.size() <= elements.size() / 2)
    ones = true;
  else
    ones = false;

  for (size_t i = 0; i < elements.size(); i++)
  {
    if (ones)
      BipartitionTools::bit0(bitbip[*cpt], static_cast<int>(i));
    else
      BipartitionTools::bit1(bitbip[*cpt], static_cast<int>(i));
  }

  for (size_t i = 0; i < underelements_.size(); i++)
  {
    size_t taxa_ind = 0;
    while (underelements_[i] != elements[taxa_ind])
      taxa_ind++;
    if (ones)
      BipartitionTools::bit1(bitbip[*cpt], static_cast<int>(taxa_ind));
    else
      BipartitionTools::bit0(bitbip[*cpt], static_cast<int>(taxa_ind));
  }

  (*cpt)++;

  if (index)
    index->push_back(nd->getId());

  return underelements_;
}

/******************************************************************************/

RowMatrix<int> BipartitionList::toMatrix() const
{
  vector< map<string, bool> > bipl;
  for (size_t i = 0; i < getNumberOfBipartitions(); i++)
  {
    bipl.push_back(getBipartition(i));
  }

  vector<string> el = getElementNames();

  RowMatrix<int> mat(el.size(), getNumberOfBipartitions());

  for (size_t j = 0; j < el.size(); j++)
  {
    for (size_t i = 0; i < getNumberOfBipartitions(); i++)
    {
      mat(j, i) = bipl[i][el[j]];
    }
  }
  return mat;
}

/******************************************************************************/

