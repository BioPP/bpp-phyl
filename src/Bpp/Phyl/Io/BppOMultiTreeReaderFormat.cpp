// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include <Bpp/Text/KeyvalTools.h>
#include <memory>
#include <string>

#include "BppOMultiTreeReaderFormat.h"
#include "Newick.h"
#include "NexusIoTree.h"
#include "Nhx.h"

using namespace bpp;
using namespace std;

IMultiTree* BppOMultiTreeReaderFormat::readIMultiTree(const std::string& description)
{
  unparsedArguments_.clear();
  string format = "";
  KeyvalTools::parseProcedure(description, format, unparsedArguments_);
  unique_ptr<IMultiTree> iTrees;
  if (format == "Newick")
  {
    bool allowComments = ApplicationTools::getBooleanParameter("allow_comments", unparsedArguments_, false, "", true, warningLevel_);
    iTrees.reset(new Newick(allowComments));
  }
  else if (format == "Nhx")
  {
    iTrees.reset(new Nhx());
  }
  else if (format == "Nexus")
  {
    iTrees.reset(new NexusIOTree());
  }
  else
  {
    throw Exception("Trees format '" + format + "' unknown.");
  }

  return iTrees.release();
}
