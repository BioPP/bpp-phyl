// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Exceptions.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>

#include "BipartitionList.h"
#include "BipartitionTools.h"
#include "TreeTemplate.h"

// From SeqLib
#include <Bpp/Seq/Alphabet/DNA.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSequenceContainer.h>

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

unique_ptr<BipartitionList> BipartitionTools::buildBipartitionPair(
    const BipartitionList& bipartL1, size_t i1,
    const BipartitionList& bipartL2, size_t i2,
    bool checkElements)
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
  return make_unique<BipartitionList>(elements, twoBitBipL);
}

/******************************************************************************/

bool BipartitionTools::areIdentical(
    const BipartitionList& bipartL1, size_t i1,
    const BipartitionList& bipartL2, size_t i2,
    bool checkElements)
{
  auto twoBipL = buildBipartitionPair(bipartL1, i1, bipartL2, i2, checkElements);
  return twoBipL->areIdentical(0, 1);
}

/******************************************************************************/

bool BipartitionTools::areCompatible(
    const BipartitionList& bipartL1, size_t i1,
    const BipartitionList& bipartL2, size_t i2,
    bool checkElements)
{
  auto twoBipL = buildBipartitionPair(bipartL1, i1, bipartL2, i2, checkElements);
  return twoBipL->areCompatible(0, 1);
}

/******************************************************************************/

std::unique_ptr<BipartitionList> BipartitionTools::mergeBipartitionLists(
    const vector<unique_ptr<BipartitionList>>& vecBipartL,
    bool checkElements)
{
  vector<string> elements;
  vector<int*> mergedBitBipL;
  int* provBitBip;
  unique_ptr<BipartitionList> mergedBipL;

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

  mergedBipL = make_unique<BipartitionList>(elements, mergedBitBipL);
  return mergedBipL;
}

/******************************************************************************/

unique_ptr<VectorSiteContainer> BipartitionTools::MRPEncode(
    const vector<unique_ptr<BipartitionList>>& vecBipartL)
{
  vector<string> all_elements;
  map<string, bool> bip;
  vector<string> bip_elements;
  shared_ptr<const Alphabet> alpha = AlphabetTools::DNA_ALPHABET;
  vector<string> sequences;

  if (vecBipartL.size() == 0)
    throw Exception("Empty vector passed");

  vector< vector<string>> vecElementLists;
  for (auto& vi : vecBipartL)
  {
    vecElementLists.push_back(vi->getElementNames());
  }

  all_elements = VectorTools::vectorUnion(vecElementLists);

  sequences.resize(all_elements.size());

  for (auto& vi : vecBipartL)
  {
    for (size_t j = 0; j < vi->getNumberOfBipartitions(); j++)
    {
      bip = vi->getBipartition(j);
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

  vector<unique_ptr<Sequence>> vec_sequences;
  for (size_t i = 0; i < all_elements.size(); i++)
  {
    vec_sequences.push_back(make_unique<Sequence>(all_elements[i], sequences[i], alpha));
  }

  VectorSequenceContainer vec_seq_cont(alpha, vec_sequences);

  auto vec_site_cont = make_unique<VectorSiteContainer>(vec_seq_cont);

  return vec_site_cont;
}

/******************************************************************************/

unique_ptr<VectorSiteContainer> BipartitionTools::MRPEncodeMultilabel(
    const vector<unique_ptr<BipartitionList>>& vecBipartL)
{
  vector<string> all_elements;
  map<string, bool> bip;
  vector<string> bip_elements;
  shared_ptr<const Alphabet> alpha = AlphabetTools::DNA_ALPHABET;
  vector<string> sequences;

  if (vecBipartL.size() == 0)
    throw Exception("Empty vector passed");

  vector< vector<string>> vecElementLists;
  for (auto& vi : vecBipartL)
  {
    vecElementLists.push_back(vi->getElementNames());
  }

  all_elements = VectorTools::vectorUnion(vecElementLists);

  sequences.resize(all_elements.size());

  for (auto& vi : vecBipartL)
  {
    for (size_t j = 0; j < vi->getNumberOfBipartitions(); j++)
    {
      bip = vi->getBipartition(j);
      bip_elements = MapTools::getKeys(bip);
      // Check for multilabel trees: if a taxa found on both sides, do not consider the entire bipartition
      vector< string > zeroes;
      vector< string > ones;
      for (auto& element : all_elements)
      {
        if (VectorTools::contains(bip_elements, element))
        {
          if (bip[element])
            ones.push_back(element);
          else
            zeroes.push_back(element);
        }
      }
      vector<string> inter = VectorTools::vectorIntersection(ones, zeroes);
      if (inter.size() != 0)  // some taxa found on both sides of the bipartition
      {
        for (auto& sk : sequences)
        {
          sk.push_back('N');
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

  vector<unique_ptr<Sequence>> vec_sequences;
  for (size_t i = 0; i < all_elements.size(); i++)
  {
    vec_sequences.push_back(make_unique<Sequence>(all_elements[i], sequences[i], alpha));
  }

  VectorSequenceContainer vec_seq_cont(alpha, vec_sequences);

  return make_unique<VectorSiteContainer>(vec_seq_cont);
}

/******************************************************************************/
