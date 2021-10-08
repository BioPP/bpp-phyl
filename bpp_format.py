#!/usr/bin/env python

# Script that can be used to format the file headers of bpp.

import os
import sys
import pathlib
import collections
import itertools
import enum
import datetime
import dateutil.parser
import shutil

# Token classes have a match method.
# This method returns unmatched text on match success.
# On match failure, it returns None.
class EmptyToken:
    ''' Detects an empty buf. Unmatched text is special, just return "". '''
    def match (self, buf):
        return "" if buf == "" else None
class LToken:
    ''' Matches text on the left. '''
    def __init__ (self, text):
        self.text = text.casefold ()
    def match (self, buf):
        return buf[len (self.text):] if len (buf) >= len (self.text) and buf[:len(self.text)].casefold () == self.text else None
class RToken:
    ''' Matches text on the right. '''
    def __init__ (self, text):
        self.text = text.casefold ()
    def match (self, buf):
        return buf[:-len (self.text)] if len (buf) >= len (self.text) and buf[-len(self.text):].casefold () == self.text else None

# Some combinators for tokens
class LStripToken:
    ''' Combinator: Strip left whitespace and matches subtoken. '''
    def __init__ (self, next_tok):
        self.next_tok = next_tok
    def match (self, buf):
        return self.next_tok.match (buf.lstrip ())
class RStripToken:
    ''' Combinator: Strip right whitespace and matches subtoken. '''
    def __init__ (self, next_tok):
        self.next_tok = next_tok
    def match (self, buf):
        return self.next_tok.match (buf.rstrip ())

class NotToken:
    ''' Combinator: Returns buf only if subtoken did not match, allow negative lookahead. '''
    def __init__ (self, token):
        self.token = token
    def match (self, buf):
        return None if self.token.match (buf) is not None else buf

class LSkipAnyUntil:
    def __init__ (self, token):
        self.token = token
    def match (self, buf):
        while len(buf) > 0:
            r = self.token.match (buf)
            if r is not None:
                return r
            buf = buf[1:]
        return None
class RSkipAnyUntil:
    def __init__ (self, token):
        self.token = token
    def match (self, buf):
        while len(buf) > 0:
            r = self.token.match (buf)
            if r is not None:
                return r
            buf = buf[:-1]
        return None

# Line parser allows to define a chain of tokens that functions as a sort of regexp.
class LineParser:
    '''
    Parse a line, by tokenizing (with whitespace) or regexp.
    Methods return new LineParser instances with remaining text (or a match failure value).
    LineParser can be evaluated as bool, tests if there was a match failure or not.
    '''
    def __init__ (self, text = None):
        self.text = text
    def __bool__ (self):
        return self.text is not None

    def __enter__ (self):
        return self
    def __exit__ (self, exc_type, exc_value, traceback):
        return None

    def match (self, tok):
        if bool (self):
            return LineParser (tok.match (self.text))
        return LineParser ()

    def match_eol (self):
        return self.match (LStripToken (EmptyToken ()))
    def match_tok (self, token):
        return self.match (LStripToken (LToken (token)))
    def match_toks (self, *tokens):
        ''' Match a sequence of tokens (chainable) '''
        if bool (self) and len (tokens) > 0:
            return self.match_tok (tokens[0]).match_toks (*tokens[1:])
        return self

    def match_or (self, *tests):
        ''' Branch: test two possibilities (chainable). Use lambdas '''
        if bool (self):
            for test in tests:
                r = test (self)
                if r:
                    return r
        return LineParser ()

    def as_tokens (self):
        ''' Returns the text as a token iterable (not chainable, fails if badmatch) '''
        return self.text.split ()
    def as_date (self):
        ''' Returns the text as a datetime obkect (not chainable, fails if badmatch) '''
        try:
          return dateutil.parser.parse (self.text.strip ())
        except:
          return self.text.strip ()

# Parse file line by line, store bounds and what has been parsed last
class FileParser:
    ''' Parse a file line by line (tracking our position), from both ends '''
    def __init__ (self, file_path):
        self.lines = file_path.read_text (encoding="latin-1").splitlines ()
        self.next_forward_line = 0 # next line to be read in fwd order
        self.next_backward_line_after = len (self.lines) # next line to be read (+1) in bwd order
        self.last_parsed = None

    def get_unparsed_lines (self):
        return self.lines[self.next_forward_line:self.next_backward_line_after]

    def parse_line (self):
        if self.next_forward_line < self.next_backward_line_after:
            index = self.next_forward_line
            self.next_forward_line += 1
            return LineParser (self.lines[index])
        else:
            raise Exception ("Reached end of forward parsing")
    def unparse_line (self):
        assert self.next_forward_line > 0
        self.next_forward_line -= 1

    def backward_parse_line (self):
        if self.next_forward_line < self.next_backward_line_after:
            index = self.next_backward_line_after - 1
            self.next_backward_line_after -= 1
            return LineParser (self.lines[index])
        else:
            raise Exception ("Reached end of backward parsing")
    def backward_unparse_line (self):
        assert self.next_backward_line_after < len (self.lines)
        self.next_backward_line_after += 1

# Enum of elements we can parse.
Element = enum.Enum ("Element", [
    "Space", "File", "Author", "CreationDate", "LastModificationDate", "License",
    "HeaderGuardTest", "HeaderGuardDef", "Include", "HeaderGuardEndif", "HeaderGuardOnce"
    ])

class BppFile:
    def __init__ (self, file_path):
        # File
        self.file_path = file_path.resolve ()
        print ("File: {}".format (self.file_path))
        # Get path to Bpp/
        try:
            self.bpp_dir = next ((d for d in self.file_path.parents if d.name == "Bpp"))
        except StopIteration:
            print ("File not in Bpp dir (test ?)")
            self.bpp_dir = None
        # Set info variables to default
        self.authors = []
        self.license_lines = None
        self.creation_date = None
        self.last_modification_date = None
        self.rel_includes = []
        self.abs_includes = []
        self.code_lines = []
        self.has_doctest = False

    def add_author (self, author_first_last_name):
        # Remove author from list (exclude suffixes)
        auth_list = [a for a in self.authors if a.split ()[:2] != author_first_last_name.split()[:2]]
        # Add it at end
        auth_list.append (author_first_last_name)
        self.authors = auth_list
        print ("Adding author: {}".format (author_first_last_name))

    # Parser setters
    def add_authors_from_elements (self, elements):
        # Split elements by "and"
        authors = (" ".join (g) for k, g in itertools.groupby(elements, lambda s: s == "and") if not k)
        for a in authors:
            self.authors.append (a)
            print ("Found author: {}".format (a))
    def add_license_line_from_elements (self, elements):
        self.license_lines = self.license_lines if self.license_lines is not None else []
        self.license_lines.append (" ".join (elements))
    def add_absolute_include_from_elements (self, elements, comment = None):
        include = " ".join (elements)
        self.abs_includes.append ((include, comment))
        print ("Found absolute include: {}".format (include))
    def add_relative_include_from_elements (self, elements, comment = None):
        relative_original_path = pathlib.Path (" ".join (elements))
        if relative_original_path.name == "doctest.h":
            self.has_doctest = True
            print ("Found doctest include")
        else:
            self.rel_includes.append ((relative_original_path, comment))
            print ("Found relative include: {}".format (relative_original_path))

    # Common parsing functions
    def parse_file_header (self, file_parser):
        while True: # Parsing file header
            line = file_parser.parse_line ()
            if line.match_eol ():
                file_parser.last_parsed = Element.Space
                continue # Empty line
            with line.match_tok ("//") as comment:
                if comment:
                    if comment.match_eol ():
                        file_parser.last_parsed = Element.Space
                        continue # Empty comment
                    if comment.match_toks ("File"):
                        file_parser.last_parsed = Element.File
                        continue # ignore filename
                    with comment.match_or (
                            lambda p: p.match_tok ("Author").match (NotToken (LToken ("s"))),
                            lambda p: p.match_tok ("Authors"),
                            lambda p: p.match_toks ("Created", "by")
                            ).match_tok (":") as author_field:
                        if author_field:
                            self.add_authors_from_elements (author_field.as_tokens ())
                            file_parser.last_parsed = Element.Author
                            continue
                    with comment.match_or (
                            lambda p: p.match_toks ("Created", "on"),
                            lambda p: p.match_tok ("Created")
                            ).match_tok (":") as creation_date:
                        if creation_date:
                            self.creation_date = creation_date.as_date ()
                            print ("Found creation date: {}".format (self.creation_date))
                            file_parser.last_parsed = Element.CreationDate
                            continue
                    with comment.match_tok ("Last").match_or (
                            lambda p: p.match_tok ("modification"),
                            lambda p: p.match_tok ("modified")
                            ).match_tok (":") as modif_date:
                        if modif_date:
                            self.last_modification_date = modif_date.as_date ()
                            print ("Found modification date: {}".format (self.last_modification_date))
                            file_parser.last_parsed = Element.LastModificationDate
                            continue
                    if file_parser.last_parsed == Element.Author:
                        self.add_authors_from_elements (comment.as_tokens ())
                        continue
                    continue
                else:
                  if line.match_tok ("/*").match_eol () or line.match_tok ("#include") or line.match_tok ("#ifndef") or line.match_tok ("namespace") or line.match_tok ("#define"):
                    file_parser.unparse_line ()
                    break # Found stuff of next steps, get out
                print ("Unexpected line (parsing file header): {}".format (line.text))
                file_parser.unparse_line ()
                break # Found unexpected line, let other pass try
          
    def parse_file_license (self, file_parser):
        while True: # Parsing license
            line = file_parser.parse_line ()
            if file_parser.last_parsed != Element.License:
                if line.match_eol ():
                    file_parser.last_parsed = Element.Space
                    continue # Empty line
                if line.match_tok ("/*").match_eol ():
                    file_parser.last_parsed = Element.License
                    continue # Entering license block
                print ("Unexpected line (parsing license): {}".format (line.text))
                file_parser.unparse_line ()
                break # Failed to find start of license
            else:
                if line.match_tok ("*/").match_eol ():
                    file_parser.last_parsed = Element.Space
                    break # End of license
                self.add_license_line_from_elements (line.as_tokens ())
                
    def parse_include_line (self, line):
        # Match a single include line ; returns LineParser with given state
        with line.match_tok ("#include") as include:
            if include:
                # Strip potential trailing comment
                comment = None
                stripped_of_comment = include.match (RSkipAnyUntil (RToken ("//")))
                if stripped_of_comment:
                    comment_elems = include.match (LSkipAnyUntil (LToken ("//"))).as_tokens ()
                    comment = " ".join (comment_elems)
                    include = stripped_of_comment
                # Continue processing include (either a quote or a system include)
                with include.match (LStripToken (LToken ("<"))).match (RStripToken (RToken (">"))) as absolute_inc:
                    if absolute_inc:
                        self.add_absolute_include_from_elements (absolute_inc.as_tokens (), comment)
                        return absolute_inc
                with include.match (LStripToken (LToken ("\""))).match (RStripToken (RToken ("\""))) as relative_inc:
                    if relative_inc:
                        self.add_relative_include_from_elements (relative_inc.as_tokens (), comment)
                        return relative_inc
        return LineParser ()
    def parse_file_code (self, file_parser):
        self.code_lines = file_parser.get_unparsed_lines ()

    # Common writing functions
    def write_file_header (self, f):
        f.write ("//\n")
        f.write ("// File: {}\n".format (self.file_path.name))
        f.write ("// Authors:\n")
        f.writelines ("//   {}\n".format (a) for a in self.authors)
        if self.creation_date is not None:
            f.write ("// Created: {}\n".format (self.creation_date))
        if self.last_modification_date is not None:
            f.write ("// Last modified: {}\n".format (self.last_modification_date))
        f.write ("//\n")
        f.write ("\n")
        
    def write_file_license (self, f):
        if self.license_lines:
            f.write ("/*\n")
            f.writelines ("  {}\n".format (l) for l in self.license_lines)
            f.write ("*/\n")
            f.write ("\n")

    def write_file_includes (self, f):
        # Sorted include files
        self.rel_includes.sort ()
        self.abs_includes.sort ()
        def include_line (inc, enclose_path):
            path, comment = inc
            enclosed_path = enclose_path (path)
            if comment:
                return "#include {} // {}\n".format (enclosed_path, comment)
            else:
                return "#include {}\n".format (enclosed_path)
        f.writelines (include_line (inc, lambda p : "<{}>".format (p)) for inc in self.abs_includes)
        f.write ("\n")
        f.writelines (include_line (inc, lambda p : "\"{}\"".format (p)) for inc in self.rel_includes)
        f.write ("\n")
    def write_file_code (self, f):
        f.writelines ("{}\n".format (l) for l in self.code_lines)

class BppHeader (BppFile):
    def __init__ (self, file_path, parse = True):
        BppFile.__init__ (self, file_path)
        # Then build header guard : BPP_DIR_DIR_FILE_H
        if self.bpp_dir:
            header_elements = ("Bpp",) + self.file_path.parent.relative_to (self.bpp_dir).parts + (self.file_path.stem, "h")
            self.header_guard = "_".join (n.upper() for n in header_elements)
            print ("Header guard: {}".format (self.header_guard))
        else:
            raise Exception ("Unable to generate header guard if file is not in bpp dir")
        # Get header info from file
        if parse:
            self.parse_file ()

    def parse_file (self):
        file_parser = FileParser (self.file_path)
        self.parse_file_header (file_parser)
        self.parse_file_license (file_parser)
        while True: # Parsing header inclusion
            line = file_parser.parse_line ()
            if line.match_eol ():
                file_parser.last_parsed = Element.Space
                continue # Empty line
            if line.match_tok ("#ifndef"):
                file_parser.last_parsed = Element.HeaderGuardTest
                continue # Remove current header guard
            if file_parser.last_parsed == Element.HeaderGuardTest and line.match_tok ("#define"):
                file_parser.last_parsed = Element.HeaderGuardDef
                continue # Remove current header guard part 2 (must follow part 1)
            if line.match_toks ("#pragma", "once"):
                file_parser.last_parsed = Element.HeaderGuardOnce
                continue # Remove pragma once
            if self.parse_include_line (line):
                file_parser.last_parsed = Element.Include
                continue # Include line upgraded
            if line.match_toks ("namespace"):
                file_parser.unparse_line ()
                break # Start of code
            print ("Unexpected line (parsing includes): {}".format (line.text))
            file_parser.unparse_line ()
            break # Try catching start of code anyway
        while True: # Parsing bottom of file
            line = file_parser.backward_parse_line ()
            if line.match_eol ():
                file_parser.last_parsed = Element.Space
                continue # Empty line
            if line.match_tok ("#endif"):
                file_parser.last_parsed = Element.HeaderGuardEndif
                continue # Rm end guard
            if line.match_tok ("}"):
                file_parser.backward_unparse_line ()
                break # Reached code
            print ("Unexpected line (parsing file footer): {}".format (line.text))
            file_parser.backward_unparse_line ()
            break # Try stopping nicely
        self.parse_file_code (file_parser)

    def write_file (self, file_path):
        with file_path.open ("w") as f:
            self.write_file_header (f)
            self.write_file_license (f)
            # f.write ("#pragma once\n") If we want to use it...
            f.write ("#ifndef {}\n".format (self.header_guard))
            f.write ("#define {}\n".format (self.header_guard))
            f.write ("\n")
            self.write_file_includes (f)
            self.write_file_code (f)
            f.write ("#endif // {}\n".format (self.header_guard))

class BppCpp (BppFile):
    def __init__ (self, file_path, parse = True):
        BppFile.__init__ (self, file_path)
        if parse:
            self.parse_file ()
    
    def parse_file (self):
        file_parser = FileParser (self.file_path)
        self.parse_file_header (file_parser)
        self.parse_file_license (file_parser)
        while True: # Parsing header inclusion
            line = file_parser.parse_line ()
            if line.match_eol ():
                file_parser.last_parsed = Element.Space
                continue # Empty line
            if self.parse_include_line (line):
                file_parser.last_parsed = Element.Include
                continue # Include line upgraded
            if line.match_toks ("#define", "DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN"):
                self.has_doctest = True
                continue # Doctest define
            if line.match_toks ("namespace"):
                file_parser.unparse_line ()
                break # Start of code
            print ("Unexpected line (parsing includes): {}".format (line.text))
            file_parser.unparse_line ()
            break # Try catching start of code anyway
        self.parse_file_code (file_parser)

    def write_file (self, file_path):
        with file_path.open ("w") as f:
            self.write_file_header (f)
            self.write_file_license (f)
            if self.has_doctest:
                f.write ("#define DOCTEST_CONFIG_IMPLEMENT_WITH_MAIN\n")
                f.write ("#include \"doctest.h\"\n")
                f.write ("\n")
            self.write_file_includes (f)
            self.write_file_code (f)

def create_bpp_file_object (path, parse_it):
    if path.suffix == ".h":
        return BppHeader (path, parse_it)
    elif path.suffix == ".cpp":
        return BppCpp (path, parse_it)
    else:
        raise Exception ("Unsupported file extension: {}".format (path))

if __name__ == "__main__":
    # Parse command line
    authors_to_add = []
    overwrite = False
    backup = False
    overwrite_checked = False
    update_modif_date = True
    create_from = None
    args = sys.argv[1:]
    while True:
        if len (args) >= 2:
            if args[0] in ("author", "a"):
                authors_to_add.append (args[1])
                args = args[2:]
                continue
            if args[0] in ("create_from", "c"):
                create_from = args[1]
                args = args[2:]
                continue
        if len (args) >= 1:
            if args[0] in ("overwrite", "O"):
                overwrite = True
                overwrite_checked = False
                args = args[1:]
                continue
            if args[0] in ("backup", "b"):
                backup = True
                args = args[1:]
                continue
            if args[0] in ("overwrite_checked", "o"):
                overwrite = False
                overwrite_checked = True
                args = args[1:]
                continue
            if args[0] == "keep_date":
                update_modif_date = False
                args = args[1:]
                continue
        if len (args) == 1 and args[0] not in ("help", "-h", "--help"):
            break
        print ("bpp_format.py [args] <file>")
        print ("\t<file> : file to reformat")
        print ("\tauthor [name] : adds [name] as recent author of the file")
        print ("\tbackup : copy <file> to <file>.bup.<file suffix>")
        print ("\toverwrite : write directly to <file> (default is writing in <file>.new.<file suffix>)")
        print ("\toverwrite_checked : same as overwrite, but shows a diff and ask for confirmation")
        print ("\tkeep_date : do not update the last modification date")
        print ("\tcreate_from [file] : create a new file from [file] (license, date, guard)")
        sys.exit ()

    # Create file object.
    # In create_from mode, do not parse underlying file, just create an empty file.
    file_path = pathlib.Path (args[0])
    
    lbpp_file=[]  
    if os.path.isdir(file_path):
      all_file = [f for f in pathlib.Path(args[0]).rglob("*.h")] + [f for f in pathlib.Path(args[0]).rglob("*.cpp")] 
    else:
      all_file=[file_path]


    lbpp_file=[]
    for f in all_file:
        lbpp_file.append(create_bpp_file_object (f, parse_it = create_from is None))

    # Update last modification date
    if update_modif_date:
      for bpp_file in lbpp_file:
        bpp_file.last_modification_date = datetime.date.today ()
    # Add authors
    for a in authors_to_add:
      for bpp_file in lbpp_file:
        bpp_file.add_author (a)
    # Fill info in create mode, by taking an example file
    if create_from is not None:
      for bpp_file in lbpp_file:
        example_file = create_bpp_file_object (pathlib.Path (create_from), True)
        bpp_file.license_lines = example_file.license_lines
        bpp_file.creation_date = datetime.date.today ()
        bpp_file.has_doctest = example_file.has_doctest
    # Write back

    for i in range(len(all_file)):
      file_path=all_file[i]
      bpp_file=lbpp_file[i]
      
      out_file_path = file_path if overwrite else file_path.with_suffix (".new" + file_path.suffix)
      if backup:
        print ("Writing backup to {}".format (file_path.with_suffix(".bup" + file_path.suffix)))
        shutil.copyfile(file_path, file_path.with_suffix(".bup" + file_path.suffix))
        
      print ("Writing to {}".format (out_file_path))
      bpp_file.write_file (out_file_path)
      # If overwrite_checked, show diff, ask before overwriting
      if overwrite_checked:
        import subprocess
        diff_cmd = subprocess.Popen (
                ["diff", "--color=always", "-U3", file_path.as_posix (), out_file_path.as_posix ()],
                stdout=subprocess.PIPE)
        less_cmd = subprocess.Popen (["less"], stdin=diff_cmd.stdout)
        less_cmd.wait ()
        diff_cmd.wait ()
        print ("Overwrite ? [y/N]")
        ans = input ()
        if ans in ("y", "Y"):
            print ("Renaming {} to {}".format (out_file_path, file_path))
            out_file_path.rename (file_path)
        else:
            print ("Delete {} ? [y/N]".format (out_file_path))
            ans2 = input ()
            if ans2 in ("y", "Y"):
                print ("Deleting {}".format (out_file_path))
                out_file_path.unlink ()

