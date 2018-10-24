//
// File: new_likelihood.cpp
// Authors:
//   Francois Gindraud (2017)
// Created: 2017-04-19 00:00:00
// Last modified: 2017-06-08
//

/*
  Copyright or © or Copr. Bio++ Development Team, (November 17, 2004)

  This software is a computer program whose purpose is to provide classes
  for numerical calculus. This file is part of the Bio++ project.

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

#define DOCTEST_CONFIG_IMPLEMENT
#include "doctest.h"

#include <Bpp/NewPhyl/FrequenciesSet.h>
#include <Bpp/NewPhyl/Model.h>
#include <Bpp/NewPhyl/Parametrizable.h>
#include <Bpp/NewPhyl/Likelihood.h>
#include <Bpp/Phyl/Model/Nucleotide/T92.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>

static bool enableDotOutput = false;

static void dotOutput(const std::string& testName, const std::vector<const bpp::dataflow::Node*>& nodes)
{
  if (enableDotOutput)
  {
    using bpp::dataflow::DotOptions;
    writeGraphToDot(
      "debug_" + testName + ".dot", nodes, DotOptions::DetailedNodeInfo | DotOptions::ShowDependencyIndex);
  }
}

using namespace bpp::dataflow;
using bpp::Dimension;
using bpp::MatrixDimension;

TEST_CASE("model")
{
  // Bpp model
  const bpp::NucleicAlphabet& alphabet = bpp::AlphabetTools::DNA_ALPHABET;
  auto model = std::unique_ptr<bpp::TransitionModel>(new bpp::T92(&alphabet, 3.));

  Context c;

  // Create set of raw parameters for the model
  auto modelParameterNodes = createParameterMap(c, *model);
  auto kappa = modelParameterNodes["kappa"];

  // Create model node
  auto modelNodeDeps = createDependencyVector(
    *model, [&modelParameterNodes](const std::string& name) { return modelParameterNodes[name]; });
  auto modelNode = ConfiguredModel::create(c, std::move(modelNodeDeps), std::move(model));

  // Create nodes generating numeric values from model.
  auto ef = EquilibriumFrequenciesFromModel::create(c, {modelNode}, bpp::rowVectorDimension(alphabet.getSize()));

  // Setup numerical derivation
  auto delta = NumericMutable<double>::create(c, 0.001);
  modelNode->config.delta = delta;
  modelNode->config.type = NumericalDerivativeType::ThreePoints;

  auto def_dkappa = ef->deriveAsValue(c, *kappa);

  // T92 : equilibrium frequencies do not depend on kappa
  CHECK (def_dkappa->getValue().isZero());

  dotOutput("model", {ef.get(), def_dkappa.get()});
}

int main(int argc, char** argv)
{
  const std::string keyword = "dot_output";
  for (int i = 1; i < argc; ++i)
  {
    if (argv[i] == keyword)
    {
      enableDotOutput = true;
    }
  }
  doctest::Context context;
  context.applyCommandLine(argc, argv);
  return context.run();
}
