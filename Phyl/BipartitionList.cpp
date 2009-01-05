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
#include <limits.h> //defines CHAR_BIT

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

BipartitionList::BipartitionList(const Tree & tr, bool sorted, vector<int> * index):
  _sorted(sorted)
{
  unsigned int nbbip;

  _elements = tr.getLeavesNames();

  if(tr.isRooted())
    nbbip = tr.getNumberOfNodes() - 2;
  else
    nbbip = tr.getNumberOfNodes() - 1;

  if(sorted) std::sort(_elements.begin(), _elements.end());

  unsigned int lword  = BipartitionTools::LWORD;
  unsigned int nbword = (_elements.size() + lword - 1) / lword;
  unsigned int nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  for(unsigned int i = 0; i < nbbip; i++)
  {
    _bitBipartitionList.push_back(new int[nbint]);
    for(unsigned int j = 0; j < nbint; j++)
    {
      _bitBipartitionList[i][j] = 0;
    }
  }

  unsigned int cpt = 0;
  vector<string> underlyingNames;
  try
  {
    //Gain some time...
    buildBitBipartitions(dynamic_cast<const TreeTemplate<Node> &>(tr).getRootNode(), _bitBipartitionList, _elements, &cpt, index);
  }
  catch(exception & e)
  {
    TreeTemplate<Node> tmp(tr);
    buildBitBipartitions(tmp.getRootNode(), _bitBipartitionList, _elements, &cpt, index);
  }
}

/******************************************************************************/

BipartitionList::BipartitionList(const vector<string> & elements, const vector<int*> & bitBipL)
{
  unsigned int lword  = BipartitionTools::LWORD;
  unsigned int nbword = (elements.size() + lword - 1) / lword;
  unsigned int nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  for(unsigned int i = 0; i < bitBipL.size(); i++)
  {
    _bitBipartitionList.push_back(new int[nbint]);
    for(unsigned int j = 0; j < nbint; j++)
    {
      _bitBipartitionList[i][j] = bitBipL[i][j];
    }
  }

  _elements = elements;

  vector<string> cp_elements = elements;
  std::sort(cp_elements.begin(), cp_elements.end());
  if(cp_elements == elements) _sorted=true; else _sorted=false;
}

/******************************************************************************/

BipartitionList::BipartitionList(const BipartitionList & bipL)
{

  unsigned int lword  = BipartitionTools::LWORD;
  unsigned int nbword = (bipL.getNumberOfElements() + lword - 1) / lword;
  unsigned int nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  _bitBipartitionList.resize(bipL.getNumberOfBipartitions());
  vector<int*> bitBipL = bipL.getBitBipartitionList();
  for(unsigned int i = 0; i < bipL.getNumberOfBipartitions(); i++)
  {
    _bitBipartitionList[i] = new int[nbint];
    for(unsigned int j = 0; j < nbint; j++)
    {
      _bitBipartitionList[i][j] = bitBipL[i][j];
    }
  }

  _elements = bipL._elements;
  _sorted = bipL._sorted;
}

/******************************************************************************/

BipartitionList & BipartitionList::operator=(const BipartitionList & bipL)
{
  unsigned int lword  = BipartitionTools::LWORD;
  unsigned int nbword = (bipL.getNumberOfElements() + lword - 1) / lword;
  unsigned int nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  for(unsigned int i = 0; i < _bitBipartitionList.size(); i++)
    delete[] _bitBipartitionList[i];
  _bitBipartitionList.resize(bipL.getNumberOfBipartitions());
  vector<int*> bitBipL = bipL.getBitBipartitionList();
  for(unsigned int i = 0; i < bipL.getNumberOfBipartitions(); i++)
  {
    _bitBipartitionList[i] = new int[nbint];
    for(unsigned int j = 0; j < nbint; j++)
    {
      _bitBipartitionList[i][j] = bitBipL[i][j];
    }
  }

  _elements = bipL._elements;
  _sorted = bipL._sorted;
  return *this;
}

/******************************************************************************/

BipartitionList::~BipartitionList()
{
  for(unsigned int i = 0; i < _bitBipartitionList.size(); i++)
    delete[] _bitBipartitionList[i];
}

/******************************************************************************/

map<string, bool> BipartitionList::getBipartition(unsigned int i) const throw (Exception)
{
  map<string, bool> bip;

  if(i>=_bitBipartitionList.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  for(unsigned int j = 0; j < _elements.size(); j++)
  {
    if(BipartitionTools::testBit(_bitBipartitionList[i], j))
      bip[_elements[j]] = true;
    else
      bip[_elements[j]] = false;
  }
  return bip;
}

/******************************************************************************/

int* BipartitionList::getBitBipartition(unsigned int i) throw (Exception)
{
  if(i >= _bitBipartitionList.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  return _bitBipartitionList[i];
}

/******************************************************************************/

bool BipartitionList::haveSameElementsThan(map <string, bool> bipart) const
{
  vector<string> elements = _elements;
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

void BipartitionList::addBipartition(map<string, bool> & bipart, bool checkElements) throw (Exception)
{
  if(checkElements && !BipartitionList::haveSameElementsThan(bipart))
    throw Exception("Distinct bipartition element sets");

  unsigned int lword  = BipartitionTools::LWORD;
  unsigned int nbword = (_elements.size() + lword - 1) / lword;
  unsigned int nbint  = nbword * lword / (CHAR_BIT * sizeof(int));
  _bitBipartitionList.push_back(new int[nbint]);
  unsigned int ind    = _bitBipartitionList.size() - 1;
  for(unsigned int j = 0; j < nbint; j++)
    _bitBipartitionList[ind][j] = 0;

  for(unsigned int i = 0; i < _elements.size(); i++)
  {
    if(bipart[_elements[i]] == true)
      BipartitionTools::bit1(_bitBipartitionList[ind], i);
    else
      BipartitionTools::bit0(_bitBipartitionList[ind], i);
  }
}

/******************************************************************************/

void BipartitionList::deleteBipartition(unsigned int i) throw(Exception)
{
  if(i >= _bitBipartitionList.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  delete[] _bitBipartitionList[i];
  _bitBipartitionList.erase(_bitBipartitionList.begin()+i);
}

/******************************************************************************/

bool BipartitionList::containsBipartition(map<string, bool> & bipart, bool checkElements) const throw (Exception)
{
  unsigned int i, j;
  bool dac, padac;

  if(checkElements && !BipartitionList::haveSameElementsThan(bipart))
    throw Exception("Distinct bipartition element sets");

  for(i = 0; i < _bitBipartitionList.size(); i++)
  {
    dac = padac = false;
    for(j = 0; j < _elements.size(); j++)
    {
      if(BipartitionTools::testBit(_bitBipartitionList[i], j))
      {
        if(bipart[_elements[j]]) dac = true;
        else padac = true;
      }
      else
      {
        if(bipart[_elements[j]]) padac = true;
        else dac = true;
	    }
      if(dac && padac) break;
    }
    if(j == _elements.size())
      return true;
  }
  return false;
}

/******************************************************************************/

bool BipartitionList::areIdentical(unsigned int k1, unsigned int k2) const throw (Exception)
{
  bool dac, padac;

  if(k1 >= _bitBipartitionList.size())
    throw Exception("Bipartition index exceeds BipartitionList size");
  if(k2 >= _bitBipartitionList.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  dac = padac = false;
  for(unsigned int j = 0; j < _elements.size(); j++)
  {
    if(BipartitionTools::testBit(_bitBipartitionList[k1], j))
    {
      if(BipartitionTools::testBit(_bitBipartitionList[k2], j)) dac = true;
      else padac = true;
    }
    else
    {
      if(BipartitionTools::testBit(_bitBipartitionList[k2], j)) padac = true;
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

  if(k1 >= _bitBipartitionList.size())
    throw Exception("Bipartition index exceeds BipartitionList size");
  if(k2 >= _bitBipartitionList.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  uu = uz = zu = zz = false;

  for(unsigned int j = 0; j < _elements.size(); j++)
  {
    if(BipartitionTools::testBit(_bitBipartitionList[k1], j))
    {
      if(BipartitionTools::testBit(_bitBipartitionList[k2], j)) uu = true;
      else uz = true;
    }
    else
    {
      if(BipartitionTools::testBit(_bitBipartitionList[k2], j)) zu = true;
      else zz = true;
    }
    if(uu && uz && zu && zz) return false;
  }

  return true;
}

/******************************************************************************/

bool BipartitionList::areAllCompatible() const
{
  for(unsigned int i = 0; i < _bitBipartitionList.size(); i++)
  {
    for(unsigned int j = i + 1; j < _bitBipartitionList.size(); j++)
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
  unsigned int nbBip = _bitBipartitionList.size();
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
  vector<StringAndInt> r_elements;
  StringAndInt sai;
  unsigned int nbbip;

  for(unsigned int i = 0; i < _elements.size(); i++)
  {
    sai.str = _elements[i];
    sai.ind = i;
    r_elements.push_back(sai);
  }

  std::sort(r_elements.begin(), r_elements.end());

  for(unsigned int i = 0; i < _elements.size(); i++)
    _elements[i] = r_elements[i].str;

  nbbip = _bitBipartitionList.size();
  _bitBipartitionList.resize(2 * nbbip);
  unsigned int lword  = BipartitionTools::LWORD;
  unsigned int nbword = (_elements.size() + lword - 1) / lword;
  unsigned int nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  for(unsigned int j = nbbip; j < 2 * nbbip; j++)
  {
	  _bitBipartitionList[j] = new int[nbint];
	  for(unsigned int k = 0; k < nbint; k++)
	    _bitBipartitionList[j][k] = 0;
    for(unsigned int i = 0; i < _elements.size(); i++)
    {
	    if(BipartitionTools::testBit(_bitBipartitionList[j - nbbip], r_elements[i].ind))
	      BipartitionTools::bit1(_bitBipartitionList[j], i);
	    else
	      BipartitionTools::bit0(_bitBipartitionList[j], i);
	  }
  }

  for(unsigned int j = 0; j < nbbip; j++)
    delete[] _bitBipartitionList[j];

  _bitBipartitionList.erase(_bitBipartitionList.begin(), _bitBipartitionList.begin() + nbbip);
  _sorted = true;
}

/******************************************************************************/

unsigned int BipartitionList::getPartitionSize(unsigned int k) const throw (Exception)
{
  unsigned int size = 0;
  if(k >= _bitBipartitionList.size())
    throw Exception("Bipartition index exceeds BipartitionList size");

  for(unsigned int i = 0; i < _elements.size(); i++)
    if(BipartitionTools::testBit(_bitBipartitionList[k], i))
      size++;

  if(size<=_elements.size()/2) return size;
  else return _elements.size() - size;
}

/******************************************************************************/

void BipartitionList::removeTrivialBipartitions()
{
  unsigned int size = _bitBipartitionList.size();
  for(unsigned int i = size; i > 0; i--)
    if(BipartitionList::getPartitionSize(i - 1) < 2)
      BipartitionList::deleteBipartition(i - 1);
}

/******************************************************************************/

void BipartitionList::addTrivialBipartitions(bool checkExisting)
{
  map<string, bool> bip;

  for(unsigned int i = 0; i < _elements.size(); i++)
    bip[_elements[i]] = false;
  for(unsigned int i = 0; i < _elements.size(); i++)
  {
    bip[_elements[i]] = true;
    if(checkExisting && BipartitionList::containsBipartition(bip, false))
      continue;
    BipartitionList::addBipartition(bip, false);
    bip[_elements[i]] = false;
  }
}

/******************************************************************************/

void BipartitionList::sortByPartitionSize()
{
  vector<int*> sortedBitBipL;
  vector<IntAndInt> iaiVec;
  IntAndInt iai;

  for(unsigned int i = 0; i < _bitBipartitionList.size(); i++)
  {
    iai.ind = i;
    iai.val = BipartitionList::getPartitionSize(i);
    iaiVec.push_back(iai);
  }

  std::sort(iaiVec.begin(), iaiVec.end());

  for(unsigned int i = 0; i < _bitBipartitionList.size(); i++)
    sortedBitBipL.push_back(_bitBipartitionList[iaiVec[i].ind]);

  _bitBipartitionList=sortedBitBipL;
}

/******************************************************************************/

void BipartitionList::flip(unsigned int k) throw (Exception)
{
  if(k >= _bitBipartitionList.size())
    throw Exception("Bipartition index exceeds BipartitionList size");
  unsigned int lword = BipartitionTools::LWORD;
  unsigned int nbword= (_elements.size() + lword - 1) / lword;
  unsigned int nbint = nbword * lword / (CHAR_BIT * sizeof(int));
  int* flipbip = new int[nbint];
  for(unsigned int i = 0; i < nbint; i++) flipbip[i] = 0;
  BipartitionTools::bitNot(flipbip, _bitBipartitionList[k], nbint);
  delete[] _bitBipartitionList[k];
  _bitBipartitionList[k] = flipbip;
}

/******************************************************************************/

void BipartitionList::removeRedundantBipartitions()
{
  bool deletion = true;

  while(deletion)
  {
	  deletion = false;
    for(unsigned int i = 0; i < _bitBipartitionList.size(); i++)
    {
      for(unsigned int j = i + 1; j < _bitBipartitionList.size(); j++)
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
  nbword = (_elements.size() + lword - 1) / lword;
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
          vecNd[i]=new Node(_elements[j]);
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
        vecNd[i]->addSon(*sonNd[k]);
    }
  }
 
  /* create last node, which joins alive bipartitions = fatherless nodes */
  Node* rootNd = new Node();
  for(unsigned int i = 0; i < sortedBipL->getNumberOfBipartitions(); i++)
    if(alive[i])
      rootNd->addSon(*vecNd[i]);

  /* construct tree and return */
  TreeTemplate<Node>* tr = new TreeTemplate<Node>(*rootNd);
  tr->resetNodesId();
  delete sortedBipL;
  return tr;
}

/******************************************************************************/

vector<string> BipartitionList::buildBitBipartitions(const Node * nd, vector<int*> & bitbip, const vector<string> & elements, unsigned int* cpt, vector<int> * index) const
{
  vector<string> under_elements, ret_elements;

  if(nd->getNumberOfSons() == 0)
    under_elements.push_back(nd->getName());

  for(unsigned int i = 0; i < nd->getNumberOfSons(); i++)
  {
    ret_elements = BipartitionList::buildBitBipartitions(nd->getSon(i), bitbip, elements, cpt, index);
    for(unsigned int j = 0; j < ret_elements.size(); j++)
      under_elements.push_back(ret_elements[j]);
  }

  if(!nd->hasFather()) return under_elements; // root node

  if(!nd->getFather()->hasFather())
  {
    unsigned int nbrootson = nd->getFather()->getNumberOfSons();
    if(nbrootson == 2 && nd == nd->getFather()->getSon(1))
    return under_elements; //son 2 of root node when root node has 2 sons
  }

  bool ones;
  if(under_elements.size() <= elements.size() / 2) ones = true;
  else ones = false;

  for(unsigned int i = 0; i < elements.size(); i++)
  {
    if(ones) BipartitionTools::bit0(bitbip[*cpt], i);
    else BipartitionTools::bit1(bitbip[*cpt], i);
  }

  for(unsigned int i = 0; i < under_elements.size(); i++)
  {
    unsigned int taxa_ind = 0;
    while(under_elements[i] != elements[taxa_ind])
      taxa_ind++;
	  if(ones) BipartitionTools::bit1(bitbip[*cpt], taxa_ind);
	  else BipartitionTools::bit0(bitbip[*cpt], taxa_ind);
  }

  (*cpt)++;

  if(index) index->push_back(nd->getId());

  return under_elements;
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

