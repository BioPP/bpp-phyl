// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_PROTEIN_COALA_H
#define BPP_PHYL_MODEL_PROTEIN_COALA_H


#include "../AbstractSubstitutionModel.h"
#include "CoalaCore.h"
#include "ProteinSubstitutionModel.h"

// From bpp-seq:
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>

namespace bpp
{
/**
 * @brief The Coala branch-heterogeneous amino-acid substitution model.
 *
 * This branch-heterogeneous model allows each branch to have its own set of amino acid equilibrium frequencies. It makes use of a Correspondence Analysis to reduce the number of parameters to be
 * optimized through Maximum Likelihood, focusing on most of the compositional variation observed in the data. The same COA is used for all branches.
 * An empirical exchangeability matrix is used, common to all branches. The choice of this matrix is up to the user. A user-defined empirical matrix can also be employed.
 * The model may also be used as a branch-homogeneous but non-stationary model, where the same estimated equilibrium frequencies are used for all branches, but different equilibrium frequencies
 * (that are optimized) are used on the root. Finally, the model may be employed as a branch-homogeneous and stationary model, where the frequencies at the root are similar to the ones on branches.
 *
 * @author Mathieu Groussin
 * @param alpha The alphabet (Protein)
 * @param nbAxes The number of principal axes of the COA that have to be taken into account to optimize the 20 branch-specific equilibrium frequencies. This number is common to all branches, as
 * well as on the root, where frequencies are optimized with a MVAprotein object (See the ProteinFrequencySet class).
 * @param exch The exchangeability matrix. The matrices currently available are DSO78, JTT92, WAG01 or LG08. A user-defined matrix can be specified with the 'file' argument.
 * @param file [optional] Used only to specify the file containing the user-defined exchangeabilities, written in PAML format.
 */

class Coala :
  public AbstractReversibleProteinSubstitutionModel,
  public CoalaCore
{
protected:
  bool init_;
  unsigned int nbrOfAxes_;
  std::string exch_;
  std::string file_;
  bool param_;

public:
  Coala(std::shared_ptr<const ProteicAlphabet> alpha,
        const ProteinSubstitutionModelInterface& model,
        unsigned int nbAxes = 0,
        bool param = true);

  virtual ~Coala() {}

  Coala* clone() const override { return new Coala(*this); }

public:
  std::string getName() const override { return "Coala"; }
  std::string getExch() const { return exch_; }
  void setFreqFromData(const SequenceDataInterface& data, double pseudoCount = 0) override;
  std::string getEmpiricalMatrixFile() const { return file_; }

protected:
  void readFromFile(std::string& file);
  void computeEquilibriumFrequencies();
  void updateMatrices_() override;
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_COALA_H
