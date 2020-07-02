//
// File: BppOBranchModelFormat.cpp
// Created by: Laurent Guéguen
// Created on: mercredi 11 mars 2020, à 14h 21
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

#include "BppOBranchModelFormat.h"

#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/KeyvalTools.h>

#include "../Model/MultinomialFromTransitionModel.h"
#include "../Model/TransitionFromTransitionModel.h"

#include "../App/PhylogeneticsApplicationTools.h"

#include "BppOFrequencySetFormat.h"

#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>

#include <Bpp/Io/OutputStream.h>
#include <Bpp/Io/BppOParametrizableFormat.h>
#include <Bpp/Io/BppODiscreteDistributionFormat.h>

//From Numeric

#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/AutoParameter.h>


#include <Bpp/Text/StringTokenizer.h>

using namespace bpp;

// From the STL:
#include <iomanip>
#include <set>

using namespace std;

BranchModel* BppOBranchModelFormat::readBranchModel(
  const Alphabet* alphabet,
    const std::string& modelDescription,
    const AlignedValuesContainer* data,
    bool parseArguments)
{
  unparsedArguments_.clear();
  unique_ptr<BranchModel> model;
  string modelName = "";
  map<string, string> args;
  KeyvalTools::parseProcedure(modelDescription, modelName, args);

  // //////////////////////////////////
  // / MIXED MODELS
  // ////////////////////////////////

  if (modelName == "Multinomial")
  {
    // We have to parse the nested model first:
    if (args.find("model")==args.end())
      throw Exception("BppOBranchModelFormat::read. Missing argument 'model' for model " + modelName);
    string nestedModelDescription = args["model"];
    BppOTransitionModelFormat nestedReader(ALL, false, allowMixed_, allowGaps_, verbose_, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_);
    
    TransitionModel* nestedModel=nestedReader.readTransitionModel(alphabet, nestedModelDescription, data, false);
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    model.reset(new MultinomialFromTransitionModel(*nestedModel));
  }
  
  if (!model)
    model.reset(readTransitionModel(alphabet, modelDescription, data, parseArguments));
  else
  {
    if (verbose_)
      ApplicationTools::displayResult("Branch model", modelName);
  
    updateParameters_(model.get(), args);
  
    if (parseArguments)
      initialize_(*model, data);
  }
  
  return model.release();
}
