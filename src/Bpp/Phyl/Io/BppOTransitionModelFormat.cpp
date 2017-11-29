//
// File: BppOTransitionModelFormat.cpp
// Created by: Laurent Guéguen
// Created on: mercredi 4 juillet 2012, à 13h 58
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

#include "BppOTransitionModelFormat.h"
#include "BppORateDistributionFormat.h"

#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/Text/KeyvalTools.h>

#include "../Model/OneChangeTransitionModel.h"
#include "../Model/OneChangeRegisterTransitionModel.h"

#include "../App/PhylogeneticsApplicationTools.h"

#include "BppOFrequenciesSetFormat.h"

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

TransitionModel* BppOTransitionModelFormat::readTransitionModel(
  const Alphabet* alphabet,
    const std::string& modelDescription,
    const AlignedValuesContainer* data,
    bool parseArguments)
{
  unparsedArguments_.clear();
  unique_ptr<TransitionModel> model;
  string modelName = "";
  map<string, string> args;
  KeyvalTools::parseProcedure(modelDescription, modelName, args);

  if (modelName == "OneChange")
  {
    // We have to parse the nested model first:
    if (args.find("model")==args.end())
      throw Exception("BppOTransitionModelFormat::read. Missing argument 'model' for model 'OneChange'.");
    string nestedModelDescription = args["model"];
    BppOSubstitutionModelFormat nestedReader(ALL, false, allowMixed_, allowGaps_, verbose_, warningLevel_);
    if (geneticCode_)
      nestedReader.setGeneticCode(geneticCode_);
    
    SubstitutionModel* nestedModel=nestedReader.read(alphabet, nestedModelDescription, data, false);
    map<string, string> unparsedParameterValuesNested(nestedReader.getUnparsedArguments());

    // We look for the register:
    if (args.find("register")==args.end())
      model.reset(new OneChangeTransitionModel(*nestedModel));
    else
    {
      string registerDescription = args["register"];
      unique_ptr<SubstitutionRegister> reg(PhylogeneticsApplicationTools::getSubstitutionRegister(registerDescription, nestedModel->getStateMap(), geneticCode_));

      if (args.find("numReg") == args.end())
        throw Exception("Missing argument 'numReg' (number of event for register in model " + modelName);

      vector<size_t> vNumRegs;
      
      StringTokenizer nst(args["numReg"], "+");

      bool out=true;
      
      while (nst.hasMoreToken())
      {
        size_t n=TextTools::to<size_t>(nst.nextToken());
        vNumRegs.push_back(n);
        if (verbose_)
        {
          ApplicationTools::displayResult(out?"Register types":"", reg->getTypeName(n));
          out=false;
        }
      }

      model.reset(new OneChangeRegisterTransitionModel(*nestedModel, *reg, vNumRegs));
    }

    // Then we update the parameter set:
    for (map<string, string>::iterator it = unparsedParameterValuesNested.begin(); it != unparsedParameterValuesNested.end(); it++)
    {
      unparsedArguments_["OneChange." + it->first] = it->second;
    }

    delete nestedModel;
  }
  else
    model.reset(BppOSubstitutionModelFormat::read(alphabet, modelDescription, data, parseArguments));

  updateParameters_(model.get(), args);
  if (parseArguments)
    initialize_(*model, data);

  return model.release();
}
