//
// File: SubstitutionProcessSequenceSimulator.h
// Created by: Laurent Guéguen
// Created on: jeudi 9 avril 2015, à 17h 17
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

#ifndef _SUBSTITUTION_PROCESS_SEQUENCE_SIMULATOR_H_
#define _SUBSTITUTION_PROCESS_SEQUENCE_SIMULATOR_H_

#include "SiteSimulator.h"
#include "SequenceSimulator.h"

#include <Bpp/Phyl/NewLikelihood/SequenceEvolution.h>

namespace bpp
{
/**
 * @brief Sequences simulation under position specific substitution process.
 *
 */

class SubstitutionProcessSequenceSimulator :
  public SequenceSimulator
{
protected:
  /**
   * @brief the map of the process simulators.
   *
   */

  std::map<size_t, std::shared_ptr<SiteSimulator> > mProcess_;

  /**
   * @brief The vector of the site specific process in mProcess_;
   * is mutable because can be changed for each simulation (for ex
   * in case of HMM).
   */

  mutable std::vector<size_t> vMap_;

  /**
   * @brief all processes trees must have at least the same sequence
   * names as the first process of the map.
   *
   */

  std::vector<std::string> seqNames_;

  /**
   * @brief correspondance map of seqNames positions of the several trees.
   * Reference is the tree of the first process of the map.
   *
   * mvPosNames[process id][i] is the position in the id_th tree
   * leaves names of the i_th name of seqName_.
   */

  std::map<size_t, std::vector<size_t> > mvPosNames_;

public:
  SubstitutionProcessSequenceSimulator(const SequenceEvolution& evol);

  SubstitutionProcessSequenceSimulator(const SubstitutionProcessSequenceSimulator&);

  SubstitutionProcessSequenceSimulator& operator=(const SubstitutionProcessSequenceSimulator&);

  SubstitutionProcessSequenceSimulator* clone() const { return new SubstitutionProcessSequenceSimulator(*this); }

  ~SubstitutionProcessSequenceSimulator();

  const SiteSimulator& getSiteSimulator(size_t pos) const
  {
    if (pos > vMap_.size())
      throw BadIntegerException("Out of range position for SubstitutionProcessSequenceSimulator", (int)pos);
    return *mProcess_.at(vMap_[pos]);
  }

  std::vector<std::string> getSequencesNames() const
  {
    return seqNames_;
  }

  /**
   * @brief Sets whether we will output the internal sequences or not.
   *
   *
   * @param yn Tell if we should output internal sequences.
   */

  void outputInternalSequences(bool yn);

  /**
   * @brief reset the set of processes.
   *
   */

  void setMap(std::vector<size_t> vMap);

  std::shared_ptr<SiteContainer> simulate(size_t numberOfSites) const;

  std::shared_ptr<SiteContainer> simulate(const std::vector<double>& rates) const;

  std::shared_ptr<SiteContainer> simulate(const std::vector<size_t>& states) const;

  std::shared_ptr<SiteContainer> simulate(const std::vector<double>& rates, const std::vector<size_t>& states) const;

  const Alphabet* getAlphabet() const;
};
} // end of namespace bpp.

#endif//_SUBSTITUTIONPROCESSSEQUENCESIMULATOR_H_
