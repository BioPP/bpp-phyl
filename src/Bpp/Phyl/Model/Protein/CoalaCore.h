// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_PROTEIN_COALACORE_H
#define BPP_PHYL_MODEL_PROTEIN_COALACORE_H



// From bpp-core:
#include <Bpp/Numeric/Matrix/Matrix.h>
#include <Bpp/Numeric/ParameterList.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/ProteicAlphabet.h>
#include <Bpp/Seq/Container/SequenceData.h>


namespace bpp
{
/**
 * @brief This class is the core class inherited by the Coala class. COaLA is a branch-heterogeneous amino-acid substitution model.
 *
 * This class allows to compute the COA from the alignment, to define the parameter (axis positions), and implements a function used to compute the equilibrium frequencies from a set of
 * coordinates along the principal axes of the COA.
 *
 * @author Mathieu Groussin
 * @param nbAxes The number of principal axes of the COA that have to be taken into account to optimize the 20 branch-specific equilibrium frequencies. This number is common to all branches, as
 * well as on the root, where frequencies are optimized with a MVAprotein object (See the ProteinFrequencySet class).
 * @param exch The exchangeability matrix. The matrices currently available are DSO78, JTT92, WAG01 or LG08. A user-defined matrix can be specified with the 'file' argument.
 */
class CoalaCore
{
protected:
  bool init_;
  size_t nbrOfAxes_;
  std::string exch_;
  RowMatrix<double> P_;
  RowMatrix<double> R_;
  std::vector<double> colWeights_;
  std::map<std::string, std::string> paramValues_;

public:
  CoalaCore(size_t nbAxes = 0, const std::string& exch = "LG08");

  virtual ~CoalaCore() {}

  CoalaCore* clone() const { return new CoalaCore(*this); }

public:
  size_t getNbrOfAxes() const { return nbrOfAxes_; }
  const RowMatrix<double>& getTppalAxesMatrix() const { return P_; }
  const RowMatrix<double>& getRowCoordinates() const { return R_; }
  const std::vector<double>& getColumnWeights() const { return colWeights_; }
  void setParamValues(const std::map<std::string, std::string>& valuesSettings) { paramValues_ = valuesSettings; }

protected:
  ParameterList computeCOA(const SequenceDataInterface& data, bool param = true);

  std::vector<double> prodMatrixVector(RowMatrix<double>& P, std::vector<double>& V);
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_PROTEIN_COALACORE_H
