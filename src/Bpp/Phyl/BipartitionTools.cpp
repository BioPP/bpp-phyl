//
// File: BipartitionTools.cpp
// Created by: Nicolas Galtier & Julien Dutheil
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

// From SeqLib
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>

using namespace bpp;

// From STL
#include <iostream>
#include <algorithm>
#include <limits.h> // defines CHAR_BIT

using namespace std;

/******************************************************************************/

int BipartitionTools::LWORD = static_cast<int>(CHAR_BIT * sizeof(int));

/******************************************************************************/

/* functions dealing with int* seen as arrays of bits */
/* (provided by Manolo Gouy) */

void BipartitionTools::bit1(int* plist, int num)
{
  // num--;
  plist += (num / LWORD);
  *plist |= (1 << (num % LWORD));
}

/******************************************************************************/

void BipartitionTools::bit0(int* plist, int num)
{
  // num--;
  plist += (num / LWORD);
  *plist &=  ~(1 << (num % LWORD));
}

/******************************************************************************/

void BipartitionTools::bitAnd(int* listet, int* list1, int* list2, size_t len)
{
  for (size_t i = 0; i < len; i++)
  {
    listet[i] = list1[i] & list2[i];
  }
}

/******************************************************************************/

void BipartitionTools::bitOr(int* listou, int* list1, int* list2, size_t len)
{
  for (size_t i = 0; i < len; i++)
  {
    listou[i] = list1[i] | list2[i];
  }
}

/******************************************************************************/

void BipartitionTools::bitNot(int* listnon, int* list, size_t len)
{
  for (size_t i = 0; i < len; i++)
  {
    listnon[i] = ~list[i];
  }
}

/******************************************************************************/
bool BipartitionTools::testBit(int* plist, int num)
{
  // num--;
  plist += (num / LWORD);
  return (*plist) & (1 << (num % LWORD));
}

/******************************************************************************/

BipartitionList* BipartitionTools::buildBipartitionPair(
  const BipartitionList& bipartL1, size_t i1,
  const BipartitionList& bipartL2, size_t i2,
  bool checkElements) throw (Exception)
{
  vector<int*> bitBipL1, bitBipL2, twoBitBipL;
  vector<string> elements;

  if (i1 >= bipartL1.getNumberOfBipartitions())
    throw Exception("Bipartition index exceeds BipartitionList size");
  if (i2 >= bipartL2.getNumberOfBipartitions())
    throw Exception("Bipartition index exceeds BipartitionList size");

  if (checkElements && !VectorTools::haveSameElements(bipartL1.getElementNames(), bipartL2.getElementNames()))
    throw Exception("Distinct bipartition element sets");

  /* get sorted bit bipartition lists */
  /* (if input is sorted: easy; otherwise: first copy, then sort) */

  if (bipartL1.isSorted())
  {
    elements = bipartL1.getElementNames();
    bitBipL1 = bipartL1.getBitBipartitionList();
  }
  else
  {
    BipartitionList provBipartL(bipartL1.getElementNames(), bipartL1.getBitBipartitionList());
    provBipartL.sortElements();
    elements = provBipartL.getElementNames();
    bitBipL1 = provBipartL.getBitBipartitionList();
  }

  if (bipartL2.isSorted())
  {
    bitBipL2 = bipartL2.getBitBipartitionList();
  }
  else
  {
    BipartitionList provBipartL(bipartL2.getElementNames(), bipartL2.getBitBipartitionList());
    provBipartL.sortElements();
    bitBipL2 = provBipartL.getBitBipartitionList();
  }

  /* create a new BipartitionList with just the two focal bipartitions */

  twoBitBipL.push_back(bitBipL1[i1]);
  twoBitBipL.push_back(bitBipL2[i2]);
  BipartitionList* twoBipL = new BipartitionList(elements, twoBitBipL);
  return twoBipL;
}

/******************************************************************************/

bool BipartitionTools::areIdentical(
  const BipartitionList& bipartL1, size_t i1,
  const BipartitionList& bipartL2, size_t i2,
  bool checkElements)
{
  BipartitionList* twoBipL = buildBipartitionPair(bipartL1, i1, bipartL2, i2, checkElements);
  bool test = twoBipL->areIdentical(0, 1);
  delete twoBipL;
  return test;
}

/******************************************************************************/

bool BipartitionTools::areCompatible(
  const BipartitionList& bipartL1, size_t i1,
  const BipartitionList& bipartL2, size_t i2,
  bool checkElements)
{
  BipartitionList* twoBipL = buildBipartitionPair(bipartL1, i1, bipartL2, i2, checkElements);
  bool test = twoBipL->areCompatible(0, 1);
  delete twoBipL;
  return test;
}

/******************************************************************************/

BipartitionList* BipartitionTools::mergeBipartitionLists(
  const vector<BipartitionList*>& vecBipartL,
  bool checkElements) throw (Exception)
{
  vector<string> elements;
  vector<int*> mergedBitBipL;
  int* provBitBip;
  BipartitionList* mergedBipL;

  if (vecBipartL.size() == 0)
    throw Exception("Empty vector passed");

  if (checkElements)
  {
    for (size_t i = 1; i < vecBipartL.size(); ++i)
    {
      if (!VectorTools::haveSameElements(vecBipartL[0]->getElementNames(), vecBipartL[0]->getElementNames()))
        throw Exception("BipartitionTools::mergeBipartitionLists. Distinct bipartition element sets");
    }
  }

  size_t lword  = static_cast<size_t>(BipartitionTools::LWORD);
  size_t nbword = (vecBipartL[0]->getElementNames().size() + lword - 1) / lword;
  size_t nbint  = nbword * lword / (CHAR_BIT * sizeof(int));

  elements = vecBipartL[0]->getElementNames();
  if (!vecBipartL[0]->isSorted())
    std::sort(elements.begin(), elements.end());

  for (size_t i = 0; i < vecBipartL.size(); i++)
  {
    vector<int*> bitBipL;
    if (vecBipartL[i]->isSorted())
    {
      bitBipL = vecBipartL[i]->getBitBipartitionList();
    }
    else
    {
      // We don't need the extra recopy here, do we?
      // BipartitionList provBipartL(BipartitionList(vecBipartL[i]->getElementNames(), vecBipartL[i]->getBitBipartitionList()));
      BipartitionList provBipartL(vecBipartL[i]->getElementNames(), vecBipartL[i]->getBitBipartitionList());
      provBipartL.sortElements();
      bitBipL = provBipartL.getBitBipartitionList();
    }
    for (size_t j = 0; j < bitBipL.size(); j++)
    {
      provBitBip = new int[nbint];
      for (size_t k = 0; k < nbint; k++)
      {
        provBitBip[k] = bitBipL[j][k];
      }
      mergedBitBipL.push_back(provBitBip);
    }
  }

  mergedBipL = new BipartitionList(elements, mergedBitBipL);
  return mergedBipL;
}

/******************************************************************************/

VectorSiteContainer* BipartitionTools::MRPEncode(
  const vector<BipartitionList*>& vecBipartL) throw (Exception)
{
  vector<string> all_elements;
  map<string, bool> bip;
  vector<string> bip_elements;
  const DNA* alpha = &AlphabetTools::DNA_ALPHABET;
  vector<string> sequences;

  if (vecBipartL.size() == 0)
    throw Exception("Empty vector passed");

  vector< vector<string> > vecElementLists;
  for (size_t i = 0; i < vecBipartL.size(); i++)
  {
    vecElementLists.push_back(vecBipartL[i]->getElementNames());
  }

  all_elements = VectorTools::vectorUnion(vecElementLists);

  sequences.resize(all_elements.size());

  for (size_t i = 0; i < vecBipartL.size(); i++)
  {
    for (size_t j = 0; j < vecBipartL[i]->getNumberOfBipartitions(); j++)
    {
      bip = vecBipartL[i]->getBipartition(j);
      bip_elements = MapTools::getKeys(bip);

      for (size_t k = 0; k < all_elements.size(); k++)
      {
        if (VectorTools::contains(bip_elements, all_elements[k]))
        {
          if (bip[all_elements[k]])
            sequences[k].push_back('C');
          else
            sequences[k].push_back('A');
        }
        else
          sequences[k].push_back('N');
      }
    }
  }

  vector<const Sequence*> vec_sequences;
  for (size_t i = 0; i < all_elements.size(); i++)
  {
    const Sequence* seq = new BasicSequence(all_elements[i], sequences[i], alpha);
    vec_sequences.push_back(seq);
  }

  VectorSequenceContainer vec_seq_cont(vec_sequences, alpha);
  for (size_t i = 0; i < all_elements.size(); i++)
  {
    delete vec_sequences[i];
  }

  VectorSiteContainer* vec_site_cont = new VectorSiteContainer(vec_seq_cont);

  return vec_site_cont;
}

/******************************************************************************/

VectorSiteContainer* BipartitionTools::MRPEncodeMultilabel(
                                                 const vector<BipartitionList*>& vecBipartL) throw (Exception)
{
    vector<string> all_elements;
    map<string, bool> bip;
    vector<string> bip_elements;
    const DNA* alpha = &AlphabetTools::DNA_ALPHABET;
    vector<string> sequences;
    
    if (vecBipartL.size() == 0)
        throw Exception("Empty vector passed");
    
    vector< vector<string> > vecElementLists;
    for (size_t i = 0; i < vecBipartL.size(); i++)
    {
        vecElementLists.push_back(vecBipartL[i]->getElementNames());
    }
    
    all_elements = VectorTools::vectorUnion(vecElementLists);
    
    sequences.resize(all_elements.size());
    
    for (size_t i = 0; i < vecBipartL.size(); i++)
    {
        for (size_t j = 0; j < vecBipartL[i]->getNumberOfBipartitions(); j++)
        {
            bip = vecBipartL[i]->getBipartition(j);
            bip_elements = MapTools::getKeys(bip);
            //Check for multilabel trees: if a taxa found on both sides, do not consider the entire bipartition
            vector< string > zeroes;
            vector< string > ones;
            for (size_t k = 0; k < all_elements.size(); k++)
            {
                if (VectorTools::contains(bip_elements, all_elements[k]))
                {
                    if (bip[all_elements[k]])
                        ones.push_back(all_elements[k]);
                    else
                        zeroes.push_back(all_elements[k]);
                }
            }
            vector<string> inter = VectorTools::vectorIntersection(ones, zeroes);
            if (inter.size() != 0) { //some taxa found on both sides of the bipartition
                for (size_t k = 0; k < all_elements.size(); k++)
                {
                    sequences[k].push_back('N');
                }
            }
            else
            {
                for (size_t k = 0; k < all_elements.size(); k++)
                {
                    if (VectorTools::contains(bip_elements, all_elements[k]))
                    {
                        if (bip[all_elements[k]])
                            sequences[k].push_back('C');
                        else
                            sequences[k].push_back('A');
                    }
                    else
                        sequences[k].push_back('N');
                }
            }
        }
    }
    
    vector<const Sequence*> vec_sequences;
    for (size_t i = 0; i < all_elements.size(); i++)
    {
        const Sequence* seq = new BasicSequence(all_elements[i], sequences[i], alpha);
        vec_sequences.push_back(seq);
    }
    
    VectorSequenceContainer vec_seq_cont(vec_sequences, alpha);
    for (size_t i = 0; i < all_elements.size(); i++)
    {
        delete vec_sequences[i];
    }
    
    VectorSiteContainer* vec_site_cont = new VectorSiteContainer(vec_seq_cont);
    
    return vec_site_cont;
}

/******************************************************************************/


