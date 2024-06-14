// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Text/KeyvalTools.h>
#include <memory>
#include <string>

#include "BppOTreeWriterFormat.h"
#include "Newick.h"
#include "NexusIoTree.h"
#include "Nhx.h"

using namespace bpp;
using namespace std;

OTree* BppOTreeWriterFormat::readOTree(const std::string& description)
{
  unparsedArguments_.clear();
  string format = "";
  KeyvalTools::parseProcedure(description, format, unparsedArguments_);
  unique_ptr<OTree> oTree;
  if (format == "Newick")
  {
    bool allowComments = ApplicationTools::getBooleanParameter("allow_comments", unparsedArguments_, false, "", true, warningLevel_);
    oTree.reset(new Newick(allowComments));
  }
  else if (format == "Nhx")
  {
    oTree.reset(new Nhx());
  }
  else if (format == "Nexus")
  {
    oTree.reset(new NexusIOTree());
  }
  else
  {
    throw Exception("Tree format '" + format + "' unknown.");
  }

  return oTree.release();
}
