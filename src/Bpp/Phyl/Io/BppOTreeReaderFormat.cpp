// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Text/KeyvalTools.h>
#include <memory>
#include <string>

#include "BppOTreeReaderFormat.h"
#include "Newick.h"
#include "NexusIoTree.h"
#include "Nhx.h"

using namespace bpp;
using namespace std;

unique_ptr<ITree> BppOTreeReaderFormat::readITree(const std::string& description)
{
  unparsedArguments_.clear();
  string format = "";
  KeyvalTools::parseProcedure(description, format, unparsedArguments_);
  unique_ptr<ITree> iTree;
  if (format == "Newick")
  {
    bool allowComments = ApplicationTools::getBooleanParameter("allow_comments", unparsedArguments_, false, "", true, warningLevel_);
    iTree.reset(new Newick(allowComments));
  }
  else if (format == "Nhx")
  {
    iTree.reset(new Nhx());
  }
  else if (format == "Nexus")
  {
    iTree.reset(new NexusIOTree());
  }
  else
  {
    throw Exception("Tree format '" + format + "' unknown.");
  }

  return iTree;
}
