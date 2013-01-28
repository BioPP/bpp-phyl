#! /usr/bin/python

import os;
from os.path import *;
import sys;

for root, dirs, files in os.walk(sys.argv[1]):
    print "-- Creating generic include file: " + str(root) + ".all"
    f = open(str(root) + ".all", 'w')
    head, tail = split(root)
    files.sort()
    for file in files:
        base, ext = splitext(file)
        if (ext == ".h") or (ext == ".all"):
          f.write("#include \"" + tail + "/" + file + "\"\n")
    f.close()
