//
// File: Model.h
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-05-29
// Last modified: 2017-05-29
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
#ifndef BPP_NEWPHYL_MODEL_H
#define BPP_NEWPHYL_MODEL_H

#include <Bpp/NewPhyl/DataFlow.h>
#include <Bpp/NewPhyl/DataFlowTemplates.h>
#include <Bpp/NewPhyl/LinearAlgebraFwd.h>
#include <Bpp/NewPhyl/PhylogenyTypes.h>
#include <map>
#include <memory>
#include <string>

namespace bpp {
class TransitionModel;

/** Interface for accessing model parameters by their names.
 * Associates Value<double> nodes to names.
 */
class ModelParameterAccessByName {
public:
	virtual ~ModelParameterAccessByName () = default;
	virtual DF::ValueRef<double> getModelParameter (const std::string & name) const = 0;
};

/** Impl for ModelParameterAccessByName: map of DF::Mutable<double>.
 * Represents the set of parameters for one TransitionModel.
 * Mutable nodes are created for each TransitionModel parameter.
 * They are registered under their non-namespaced names.
 * They can be changed to point to other Mutable<double> objects.
 */
class ModelParameterMap : public ModelParameterAccessByName {
public:
	/** Creates a new ModelParameterMap.
	 * One Mutable<double> object is created for each model parameter.
	 * It is registered in the map with its non-namespaced name.
	 * Mutable<double> nodes are initialized with values from the TransitionModel parameters.
	 * The reference to the model is only used in this constructor, not stored.
	 */
	ModelParameterMap (const TransitionModel & model);

	/// Impl of ModelParameterAccessByName, returns the parameter or throws.
	DF::ValueRef<double> getModelParameter (const std::string & name) const override;

	/** Access Mutable directly, or throws if not found (non-const version).
	 * Can change which Mutable is associated to this name.
	 * This is useful to have multiple models share a subset of parameter values.
	 */
	DF::MutableRef<double> & operator[] (const std::string & name);

	/** Access Mutable directly or throws if not found (const version).
	 *  Mutable cannot be changed, but its value can.
	 */
	const DF::MutableRef<double> & operator[] (const std::string & name) const;

	/// Direct map access (for info, debug, iteration).
	const std::map<std::string, DF::MutableRef<double>> & getMap () const {
		return mutableNodeByName_;
	}

private:
	std::map<std::string, DF::MutableRef<double>> mutableNodeByName_;
};

namespace DF {
	// Compute nodes

	// (model) -> Vector of freqs
	class EquilibriumFrequenciesFromModel;
	template <> struct Builder<EquilibriumFrequenciesFromModel> {
		// Dim == nbStates FIXME add a dim type for Equ Freq ?
		static ValueRef<VectorDouble> make (NodeRefVec && deps, const Dimension<VectorDouble> & dim);
	};

	// (model, branch length) -> transition matrix
	class TransitionMatrixFromModel;
	template <> struct Builder<TransitionMatrixFromModel> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, const TransitionMatrixDimension & dim);
	};

	// (model, branch length) -> d(transition matrix)
	class TransitionMatrixFromModelBrlenDerivative;
	template <> struct Builder<TransitionMatrixFromModelBrlenDerivative> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, const TransitionMatrixDimension & dim);
	};

	// (model, branch length) -> d2(transition matrix)
	class TransitionMatrixFromModelBrlenSecondDerivative;
	template <> struct Builder<TransitionMatrixFromModelBrlenSecondDerivative> {
		static ValueRef<MatrixDouble> make (NodeRefVec && deps, const TransitionMatrixDimension & dim);
	};
} // namespace DF
} // namespace bpp

#endif // BPP_NEWPHYL_MODEL_H
