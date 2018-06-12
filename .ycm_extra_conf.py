import os
import ycm_core

def DirectoryOfThisScript():
    return os.path.dirname(os.path.abspath(__file__))

project_root = DirectoryOfThisScript()

if project_root:
    database = ycm_core.CompilationDatabase(project_root)
else:
    database = None

SOURCE_EXTENSIONS = ['.cpp']

def IsHeaderFile( filename ):
  extension = os.path.splitext( filename )[ 1 ]
  return extension == '.h'

def FindCorrespondingSourceFile( filename ):
  if IsHeaderFile( filename ):
    basename = os.path.splitext( filename )[ 0 ]
    for extension in SOURCE_EXTENSIONS:
      replacement_file = basename + extension
      if os.path.exists( replacement_file ):
        return replacement_file
  return filename

def FlagsForFile( filename, **kwargs ):
  # If the file is a header, try to find the corresponding source file and
  # retrieve its flags from the compilation database if using one. This is
  # necessary since compilation databases don't have entries for header files.
  # In addition, use this source file as the translation unit. This makes it
  # possible to jump from a declaration in the header file to its definition in
  # the corresponding source file.
  filename = FindCorrespondingSourceFile( filename )

  compilation_info = database.GetCompilationInfoForFile( filename )
  if not compilation_info.compiler_flags_:
    return None

  # Bear in mind that compilation_info.compiler_flags_ does NOT return a
  # python list, but a "list-like" StringVec object.
  final_flags = list( compilation_info.compiler_flags_ )

  return {
    'flags': final_flags,
    'include_paths_relative_to_dir': compilation_info.compiler_working_dir_,
    'override_filename': filename
  }
