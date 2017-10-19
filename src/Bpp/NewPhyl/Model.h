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
#include <Bpp/NewPhyl/DataFlowNumeric.h> // Matrix types
#include <Bpp/NewPhyl/DataFlowTemplates.h>
#include <Bpp/NewPhyl/Signed.h>
#include <memory>
#include <string>

namespace bpp {
class SubstitutionModel;

namespace Phyl {
	namespace DF {
		using namespace bpp::DF;

		// DF Node representing a model with its parameter (wrapper to master code).
		class Model : public Value<const SubstitutionModel *> {
		public:
			Model (std::unique_ptr<SubstitutionModel> model);
			~Model ();

			SizeType nbParameters () const noexcept;
			ParameterRef<double> getParameter (SizeType index);
			ParameterRef<double> getParameter (const std::string & name);
			const std::string & getParameterName (SizeType index);

			void compute () override final;
			std::string description () const override final;
			std::string debugInfo () const override final;

			static std::shared_ptr<Model> create (std::unique_ptr<SubstitutionModel> model);

		private:
			std::unique_ptr<SubstitutionModel> model_;
		};

		// Compute nodes

		struct EquilibriumFrequenciesFromModel : public Value<VectorDouble> {
			// -> vector of equilibrium frequencies by state
			using Dependencies = FunctionOfValues<const SubstitutionModel *>;
			EquilibriumFrequenciesFromModel (NodeRefVec && deps, SizeType nbStates);
			void compute () override final;
			std::string debugInfo () const override final;
			NodeRef derive (const Node & node) override final;
			static std::shared_ptr<EquilibriumFrequenciesFromModel> create (NodeRefVec && deps,
			                                                                SizeType nbStates);
			static std::shared_ptr<EquilibriumFrequenciesFromModel>
			create (ValueRef<const SubstitutionModel *> model, SizeType nbStates);
		};

		struct TransitionMatrixFromModel : public Value<MatrixDouble> {
			// (model, branch length) -> transition matrix
			using Dependencies = FunctionOfValues<const SubstitutionModel *, double>;
			TransitionMatrixFromModel (NodeRefVec && deps, SizeType nbStates);
			void compute () override final;
			std::string debugInfo () const override final;
			NodeRef derive (const Node & node) override final;
			static std::shared_ptr<TransitionMatrixFromModel> create (NodeRefVec && deps,
			                                                          SizeType nbStates);
			static std::shared_ptr<TransitionMatrixFromModel>
			create (ValueRef<const SubstitutionModel *> model, ValueRef<double> brlen, SizeType nbStates);
		};

		struct TransitionMatrixFromModelBrlenDerivative : public Value<MatrixDouble> {
			using Dependencies = FunctionOfValues<const SubstitutionModel *, double>;
			TransitionMatrixFromModelBrlenDerivative (NodeRefVec && deps, SizeType nbStates);
			void compute () override final;
			std::string debugInfo () const override final;
			static std::shared_ptr<TransitionMatrixFromModelBrlenDerivative> create (NodeRefVec && deps,
			                                                                         SizeType nbStates);
		};
	} // namespace DF
} // namespace Phyl
} // namespace bpp

#endif // BPP_NEWPHYL_MODEL_H
