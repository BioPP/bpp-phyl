// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#ifndef BPP_PHYL_MODEL_FREQUENCYSET_MVAFREQUENCYSET_H
#define BPP_PHYL_MODEL_FREQUENCYSET_MVAFREQUENCYSET_H


#include "../Protein/Coala.h"
#include "ProteinFrequencySet.h"

namespace bpp
{
/**
 * @brief A frequencies set used to estimate frequencies at the root with the COaLA model.
 * Frequencies at the root are optimized in the same way than the equlibrium frequencies on branches.
 * Hyperparameters are used, which represent positions along the principal axes obtained from a preliminary Correspondence Analysis.
 * From the optimized positions, the 20 frequencies are calculated.
 * @author Mathieu Groussin
 */
class MvaFrequencySet :
  public virtual ProteinFrequencySetInterface,
  public AbstractFrequencySet
{
public:
  /**
   * @brief Constructor
   */
  MvaFrequencySet(std::shared_ptr<const ProteicAlphabet> alpha);

  MvaFrequencySet* clone() const override { return new MvaFrequencySet(*this); }

  /*
     MvaFrequencySet& operator=(const MvaFrequencySet& mfs)
     {
      AbstractFrequencySet::operator=(mfs);
      tPpalAxes_ = mfs.tPpalAxes_;
      rowCoords_ = mfs.rowCoords_;
      nbrOfAxes_ = mfs.nbrOfAxes_;
      model_ = mfs.model_;
      columnWeights_ = mfs.columnWeights_;
      paramValues_ = mfs.paramValues_;
      return *this;
     }
   */

protected:
  RowMatrix<double> tPpalAxes_;
  RowMatrix<double> rowCoords_;
  size_t nbrOfAxes_;
  std::string model_;
  std::vector<double> columnWeights_;
  std::map<std::string, std::string> paramValues_;

public:
  std::shared_ptr<const ProteicAlphabet> getProteicAlphabet() const override
  { 
    return std::dynamic_pointer_cast<const ProteicAlphabet>(getAlphabet());
  }

  void setTransposeMatrixOfPpalAxes(const RowMatrix<double>& matrix) { tPpalAxes_ = matrix; }
  void setMatrixOfRowCoords(const RowMatrix<double>& matrix) { rowCoords_ = matrix; }
  void setNbrOfAxes(const size_t& nAxes) { nbrOfAxes_ = nAxes; }
  void setModelName(const std::string& modelName) { model_ = modelName; }
  void setVectorOfColumnWeights(const std::vector<double>& cw) { columnWeights_ = cw; }
  void setParamValues(std::map<std::string, std::string>& valuesSettings) {paramValues_ = valuesSettings;}

  void setFrequencies(const std::vector<double>& frequencies) override;

  void defineParameters();
  void fireParameterChanged(const ParameterList& parameters) override;
  void updateFrequencies();

  void initSet(const CoalaCore& coala);

  void computeReversePCA(const std::vector<double>& positions, std::vector<double>& tmpFreqs);
  void computeCoordsFirstSpacePCA(std::vector<double>& tmpFreqs, std::vector<double>& freqs);
  void computeReverseCOA(const std::vector<double>& positions, std::vector<double>& tmpFreqs);
  void computeCoordsFirstSpaceCOA(std::vector<double>& tmpFreqs, std::vector<double>& freqs);
};
} // end of namespace bpp.
#endif // BPP_PHYL_MODEL_FREQUENCYSET_MVAFREQUENCYSET_H
