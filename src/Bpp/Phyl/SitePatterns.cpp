// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "SitePatterns.h"

// From the bpp-seq library:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Container/AlignmentData.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>

using namespace bpp;
using namespace std;

/******************************************************************************/

SitePatterns::SitePatterns(
    const AlignmentDataInterface& sequences,
    std::vector<std::string> names) :
  names_(),
  sites_(),
  weights_(),
  indices_(),
  alpha_(sequences.getAlphabet())
{
  names_ = sequences.getSequenceNames();
  if (names.size() != 0)
    names_ = VectorTools::vectorIntersection(names_, names);
  init_(sequences, names_);
}

void SitePatterns::init_(
    const AlignmentDataInterface& sequences,
    std::vector<std::string> names)
{
  // positions of the names in sequences list
  std::vector<size_t> posseq;
  for (const auto& n : names)
  {
    posseq.push_back(sequences.getSequencePosition(n));
  }

  int nbSeq = static_cast<int>(sequences.getNumberOfSequences());
  std::vector<size_t> posnseq;

  std::stable_sort(posseq.begin(), posseq.end());
  for (int i = nbSeq - 1; i >= 0; i--)
  {
    if (!std::binary_search(posseq.begin(), posseq.end(), i))
      posnseq.push_back(static_cast<size_t>(i));
  }

  // Then build Sortable sites with correct sequences
  size_t nbSites = sequences.getNumberOfSites();

  vector<SortableSite> ss(nbSites);
  for (size_t i = 0; i < nbSites; i++)
  {
    CoreSiteInterface* currentSite = sequences.site(i).clone();
    for (auto pos : posnseq)
    {
      currentSite->deleteElement(pos);
    }

    SortableSite* ssi = &ss[i];
    ssi->siteS = currentSite->toString();
    ssi->siteP = currentSite;
    ssi->originalPosition = i;
  }

  if (nbSites > 0)
  {
    // Quick sort according to site contents:
    sort(ss.begin(), ss.end());

    // Now build patterns:
    SortableSite* ss0 = &ss[0];
    auto previousSite = ss0->siteP;
    indices_.resize(Eigen::Index(nbSites));
    indices_[Eigen::Index(ss0->originalPosition)] = 0;
    sites_.push_back(shared_ptr<const CoreSiteInterface>(previousSite));
    weights_.push_back(1);

    size_t currentPos = 0;
    for (size_t i = 1; i < nbSites; ++i)
    {
      SortableSite* ssi = &ss[i];
      auto currentSite = ssi->siteP;

      bool siteExists = SymbolListTools::areSymbolListsIdentical(*currentSite, *previousSite);
      if (siteExists)
      {
        weights_[currentPos]++;
        delete currentSite;
      }
      else
      {
        sites_.push_back(shared_ptr<const CoreSiteInterface>(currentSite));
        weights_.push_back(1);
        currentPos++;
        previousSite = currentSite;
      }
      indices_[Eigen::Index(ssi->originalPosition)] = currentPos;
    }
  }
}

/******************************************************************************/

unique_ptr<AlignmentDataInterface> SitePatterns::getSites() const
{
  if (sites_.size() == 0)
    throw Exception("SitePatterns::getSites : empty set.");

  unique_ptr<AlignmentDataInterface> sites;

  if (dynamic_pointer_cast<const Site>(sites_[0]))
  {
    // Copy the sites
    vector<unique_ptr<Site>> vSites;
    for (auto& s : sites_)
    {
      auto ptr = unique_ptr<Site>(dynamic_cast<Site*>(s->clone()));
      vSites.push_back(std::move(ptr));
    }
    sites.reset(new VectorSiteContainer(vSites, alpha_));
    sites->setSequenceNames(names_, true);
    return sites;
  }

  if (dynamic_pointer_cast<const ProbabilisticSite>(sites_[0]))
  {
    // Copy the sites
    vector<unique_ptr<ProbabilisticSite>> vSites;
    for (auto& s : sites_)
    {
      vSites.push_back(unique_ptr<ProbabilisticSite>(dynamic_cast<ProbabilisticSite*>(s->clone())));
    }
    sites.reset(new ProbabilisticVectorSiteContainer(vSites, alpha_));
    sites->setSequenceNames(names_, true);
    return sites;
  }

  throw Exception("SitePatterns::getSites(). Unsupported site type.");
}

/******************************************************************************/
