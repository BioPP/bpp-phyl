// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

#include "PhylipDistanceMatrixFormat.h"
#include "IoDistanceMatrixFactory.h"

using namespace bpp;

const std::string IODistanceMatrixFactory::PHYLIP_FORMAT = "Phylip";

IDistanceMatrix* IODistanceMatrixFactory::createReader(const std::string& format, bool extended)
{
  if (format == PHYLIP_FORMAT)
    return new PhylipDistanceMatrixFormat(extended);
  else
    throw Exception("Format " + format + " is not supported for input.");
}

ODistanceMatrix* IODistanceMatrixFactory::createWriter(const std::string& format, bool extended)
{
  if (format == PHYLIP_FORMAT)
    return new PhylipDistanceMatrixFormat(extended);
  else
    throw Exception("Format " + format + " is not supported for output.");
}
