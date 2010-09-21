#! /usr/bin/python

import os;
from os.path import *;

for root, dirs, files in os.walk('src/Bpp'):
    print "Creating generic include file " + str(root) + ".all"
    f = open(str(root) + ".all", 'w')
    head, tail = split(root)
    for file in files:
        base, ext = splitext(file)
        if (ext == ".h"):
          f.write("#include \"" + tail + "/" + file + "\"\n")
    f.close()
    if '.svn' in dirs:
        dirs.remove('.svn')  # don't visit SVN directories
