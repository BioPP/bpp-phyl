//
// File: SimpleSubstitutionProcessSequenceSimulator.h
// Created by: Laurent Guéguen
// Created on: dimanche 24 mai 2020, à 07h 30
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

#ifndef _SIMPLE_SUBSTITUTION_PROCESS_SEQUENCE_SIMULATOR_H_
#define _SIMPLE_SUBSTITUTION_PROCESS_SEQUENCE_SIMULATOR_H_

#include "SiteSimulator.h"
#include "SimpleSubstitutionProcessSiteSimulator.h"
#include "GivenDataSubstitutionProcessSiteSimulator.h"

#include "SequenceSimulator.h"

#include "../NewLikelihood/SequenceEvolution.h"

namespace bpp
{
/**
 * @brief Sequences simulation under a unique substitution process.
 *
 */

class SimpleSubstitutionProcessSequenceSimulator :
  public virtual SequenceSimulator
{
private:
  std::shared_ptr<SiteSimulator> siteSim_;

public:
  SimpleSubstitutionProcessSequenceSimulator(const SubstitutionProcess& process) :
    siteSim_(std::make_shared<SimpleSubstitutionProcessSiteSimulator>(process))
  {}

  /*
   * @brief A posterior simulation, from a position in an alignment.
   *
   */

  SimpleSubstitutionProcessSequenceSimulator(std::shared_ptr<LikelihoodCalculationSingleProcess> calcul, size_t pos) :
    siteSim_(std::make_shared<GivenDataSubstitutionProcessSiteSimulator>(calcul, pos))
  {}

  SimpleSubstitutionProcessSequenceSimulator(std::shared_ptr<SiteSimulator> simul) :
    siteSim_(simul) {}


  virtual ~SimpleSubstitutionProcessSequenceSimulator()
  {}

  SimpleSubstitutionProcessSequenceSimulator(const SimpleSubstitutionProcessSequenceSimulator& nhss) :
    siteSim_(nhss.siteSim_)
  {}

  SimpleSubstitutionProcessSequenceSimulator* clone() const { return new SimpleSubstitutionProcessSequenceSimulator(*this); }

public:
  /**
   * @name The SequenceSimulator interface
   *
   */

  std::shared_ptr<SiteContainer> simulate(size_t numberOfSites) const;


  const SiteSimulator& getSiteSimulator(size_t pos) const
  {
    return *siteSim_;
  }

  std::vector<std::string> getSequencesNames() const
  {
    return siteSim_->getSequencesNames();
  }

  /**
   * @name SiteSimulator and SequenceSimulator interface
   *
   * @{
   */
  const Alphabet* getAlphabet() const { return siteSim_->getAlphabet(); }
  /** @} */

  /**
   * @brief Sets whether we will output the internal sequences or not.
   *
   *
   * @param yn Tell if we should output internal sequences.
   */
  void outputInternalSequences(bool yn)
  {
    siteSim_->outputInternalSites(yn);
  }
};
} // end of namespace bpp.

#endif//_SIMPLESUBSTITUTIONPROCESSSEQUENCESIMULATOR_H_
