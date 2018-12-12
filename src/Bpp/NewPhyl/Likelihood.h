//
// File: Likelihood.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-03
// Last modified: 2017-05-03
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 16, 2004)

  This software is a computer program whose purpose is to provide classes
  for phylogenetic data analysis.

  This software is governed by the CeCILL license under French law and
  abiding by the rules of distribution of free software. You can use,
  modify and/ or redistribute the software under the terms of the CeCILL
  license as circulated by CEA, CNRS and INRIA at the following URL
  "http://www.cecill.info".

  As a counterpart to the access to the source code and rights to copy,
  modify and redistribute granted by the license, users are provided only
  with a limited warranty and the software's author, the holder of the
  economic rights, and the successive licensors have only limited
  liability.

  In this respect, the user's attention is drawn to the risks associated
  with loading, using, modifying and/or developing or reproducing the
  software by the user in light of its specific status of free software,
  that may mean that it is complicated to manipulate, and that also
  therefore means that it is reserved for developers and experienced
  professionals having in-depth computer knowledge. Users are therefore
  encouraged to load and test the software's suitability as regards their
  requirements in conditions enabling the security of their systems and/or
  data to be ensured and, more generally, to use and operate it in the
  same conditions as regards security.

  The fact that you are presently reading this means that you have had
  knowledge of the CeCILL license and that you accept its terms.
*/

#pragma once
#ifndef BPP_NEWPHYL_LIKELIHOOD_H
#define BPP_NEWPHYL_LIKELIHOOD_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowCWise.h>
#include <functional>
#include <unordered_map>

namespace bpp {
  /** Conditional likelihoods are stored in a matrix of sizes (nbState, nbSite).
   * Rows represents states (nucleotides, proteins or codon).
   * Columns represents sites (one site for each column).
   * Conditional likelihood is thus accessed by m(state,site) for an eigen matrix.
   * Eigen defaults to ColMajor matrices, meaning data is stored column by column.
   * So with the current layout, values for a site are grouped together for locality.
   *
   * A Transition matrix is a (nbState,nbState) matrix.
   * tm(fromState, toState) = probability of going to toState from fromState.
   * This matches the convention from TransitionModel::getPij_t().
   *
   * Equilibrium frequencies are stored as a RowVector(nbState) : matrix with 1 row and n columns.
   * This choice allows to reuse the MatrixProduct numeric node directly.
   *
   * Initial conditional likelihood for leaves (sequences on the tree) should be computed outside of the
   * dataflow graph, and provided as NumericConstant<MatrixXd>.
   */
  inline MatrixDimension conditionalLikelihoodDimension (std::size_t nbState, std::size_t nbSite) {
    return {Eigen::Index (nbState), Eigen::Index (nbSite)};
  }
  inline MatrixDimension transitionMatrixDimension (std::size_t nbState) {
    return {Eigen::Index (nbState), Eigen::Index (nbState)};
  }
  inline MatrixDimension equilibriumFrequenciesDimension (std::size_t nbState) {
    return rowVectorDimension (Eigen::Index (nbState));
  }

  // Dataflow nodes for likelihood computation.
  namespace dataflow {
    
    /** @brief likelihood = f(equilibriumFrequencies, rootConditionalLikelihood).
     * - likelihood: RowVector(site).
     * - equilibriumFrequencies: RowVector(state).
     * - rootConditionalLikelihood: Matrix(state, site).
     *
     * likelihood(site) = sum_state equFreqs(state) * rootConditionalLikelihood(state, site).
     * Using matrix multiply: likelihood = equilibriumFrequencies * rootConditionalLikelihood.
     */

    using LikelihoodFromRootConditional =
      MatrixProduct<Eigen::RowVectorXd, Eigen::RowVectorXd, Eigen::MatrixXd>;

    /** @brief totalLikelihood = product_site likelihood(site).
     * - likelihood: RowVector (site).
     * - totalLikelihood: Extended float.
     */

    using TotalLogLikelihood = SumOfLogarithms<Eigen::RowVectorXd>;

    // /** @brief totalLikelihood = log(likelihood)
    //  * - likelihood: double.
    //  * - totalLogLikelihood: double.
    //  */

    // using LogLikelihood = CWiseLog<double>;

  } // namespace dataflow
} // namespace bpp

#endif // BPP_NEWPHYL_LIKELIHOOD_H
