// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "IoTreeFactory.h"
#include "Newick.h"
#include "NexusIoTree.h"
#include "Nhx.h"

using namespace bpp;

const std::string IOTreeFactory::NEWICK_FORMAT = "Newick";
const std::string IOTreeFactory::NEXUS_FORMAT = "Nexus";
const std::string IOTreeFactory::NHX_FORMAT = "Nhx";

ITree* IOTreeFactory::createReader(const std::string& format)
{
  if (format == NEWICK_FORMAT)
    return new Newick();
  else if (format == NEXUS_FORMAT)
    return new NexusIOTree();
  else if (format == NHX_FORMAT)
    return new Nhx();
  else
    throw Exception("Format " + format + " is not supported for input.");
}

OTree* IOTreeFactory::createWriter(const std::string& format)
{
  if (format == NEWICK_FORMAT)
    return new Newick();
  else if (format == NEXUS_FORMAT)
    return new NexusIOTree();
  else if (format == NHX_FORMAT)
    return new Nhx();
  else
    throw Exception("Format " + format + " is not supported for output.");
}
