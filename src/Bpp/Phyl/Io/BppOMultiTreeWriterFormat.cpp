// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Text/KeyvalTools.h>
#include <memory>
#include <string>

#include "BppOMultiTreeWriterFormat.h"
#include "Newick.h"
#include "NexusIoTree.h"
#include "Nhx.h"

using namespace bpp;
using namespace std;

OMultiTree* BppOMultiTreeWriterFormat::readOMultiTree(const std::string& description)
{
  unparsedArguments_.clear();
  string format = "";
  KeyvalTools::parseProcedure(description, format, unparsedArguments_);
  unique_ptr<OMultiTree> oTrees;
  if (format == "Newick")
  {
    bool allowComments = ApplicationTools::getBooleanParameter("allow_comments", unparsedArguments_, false, "", true, warningLevel_);
    oTrees.reset(new Newick(allowComments));
  }
  else if (format == "Nhx")
  {
    oTrees.reset(new Nhx());
  }
  else if (format == "Nexus")
  {
    oTrees.reset(new NexusIOTree());
  }
  else
  {
    throw Exception("Trees format '" + format + "' unknown.");
  }

  return oTrees.release();
}
