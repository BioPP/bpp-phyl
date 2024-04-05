// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Io/BppODiscreteDistributionFormat.h>
#include <Bpp/Io/BppOParametrizableFormat.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Io/OutputStream.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/TextTools.h>

#include "../App/PhylogeneticsApplicationTools.h"
#include "../Model/MultinomialFromTransitionModel.h"
#include "../Model/TransitionFromTransitionModel.h"
#include "BppOBranchModelFormat.h"
#include "BppOFrequencySetFormat.h"

// From Numeric

#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/AutoParameter.h>


#include <Bpp/Text/StringTokenizer.h>

using namespace bpp;

// From the STL:
#include <iomanip>
#include <set>

using namespace std;

std::unique_ptr<BranchModelInterface> BppOBranchModelFormat::readBranchModel(
    shared_ptr<const Alphabet> alphabet,
    const string& modelDescription,
    const AlignmentDataInterface& data,
    bool parseArguments)
{
  unparsedArguments_.clear();
  unique_ptr<BranchModelInterface> model;
  string modelName = "";
  map<string, string> args;
  KeyvalTools::parseProcedure(modelDescription, modelName, args);

  // //////////////////////////////////
  // / MIXED MODELS
  // ////////////////////////////////

  if (modelName == "Multinomial")
  {
    // We have to parse the nested model first:
    if (args.find("model") == args.end())
      throw Exception("BppOBranchModelFormat::read. Missing argument 'model' for model " + modelName);
    string nestedModelDescription = args["model"];
    BppOTransitionModelFormat nestedReader(ALL, false, allowMixed_, allowGaps_, verbose_, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_);

    std::shared_ptr<TransitionModelInterface> nestedModel = nestedReader.readTransitionModel(alphabet, nestedModelDescription, data, false);
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    model.reset(new MultinomialFromTransitionModel(nestedModel));
  }

  if (!model)
    model = readTransitionModel(alphabet, modelDescription, data, parseArguments);
  else
  {
    if (verbose_)
      ApplicationTools::displayResult("Branch model", modelName);

    updateParameters_(*model, args);

    if (parseArguments)
      initialize_(*model, data);
  }

  return model;
}
