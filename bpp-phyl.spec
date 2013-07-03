%define _basename bpp-phyl
%define _version 2.1.0
%define _release 1
%define _prefix /usr

URL: http://biopp.univ-montp2.fr/

Name: %{_basename}
Version: %{_version}
Release: %{_release}
License: CECILL-2.0
Vendor: The Bio++ Project
Source: http://biopp.univ-montp2.fr/repos/sources/%{_basename}-%{_version}.tar.gz
Summary: Bio++ Phylogenetics library
Group: Development/Libraries/C and C++
Requires: bpp-core = %{_version}
Requires: bpp-seq = %{_version}

BuildRoot: %{_builddir}/%{_basename}-root
BuildRequires: cmake >= 2.6.0
BuildRequires: gcc-c++ >= 4.0.0
BuildRequires: libbpp-core2 = %{_version}
BuildRequires: libbpp-core-devel = %{_version}
BuildRequires: libbpp-seq9 = %{_version}
BuildRequires: libbpp-seq-devel = %{_version}

AutoReq: yes
AutoProv: yes

%description
This library contains utilitary and classes for phylogenetics and molecular evolution analysis.
It is part of the Bio++ project.

%package -n libbpp-phyl9
Summary: Bio++ Phylogenetics library
Group: Development/Libraries/C and C++

%description -n libbpp-phyl9
This library contains utilitary and classes for phylogenetics and molecular evolution analysis.
It is part of the Bio++ project.


%package -n libbpp-phyl-devel
Summary: Libraries, includes to develop applications with %{_basename}
Group: Development/Libraries/C and C++
Requires: libbpp-phyl9 = %{_version}
Requires: libbpp-seq9 = %{_version}
Requires: libbpp-seq-devel = %{_version}
Requires: libbpp-core2 = %{_version}
Requires: libbpp-core-devel = %{_version}

%description -n libbpp-phyl-devel
The libbpp-phyl-devel package contains the header files and static libraries for
building applications which use %{_basename}.

%prep
%setup -q

%build
CFLAGS="$RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix} -DBUILD_TESTING=OFF"
if [ %{_lib} == 'lib64' ] ; then
  CMAKE_FLAGS="$CMAKE_FLAGS -DLIB_SUFFIX=64"
fi
cmake $CMAKE_FLAGS .
make

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -n libbpp-phyl9 -p /sbin/ldconfig

%post -n libbpp-phyl-devel
createGeneric() {
  echo "-- Creating generic include file: $1.all"
  #Make sure we run into subdirectories first:
  dirs=()
  for file in "$1"/*
  do
    if [ -d "$file" ]
    then
      # Recursion:
      dirs+=( "$file" )
    fi
  done
  for dir in ${dirs[@]}
  do
    createGeneric $dir
  done
  #Now list all files, including newly created .all files:
  if [ -f $1.all ]
  then
    rm $1.all
  fi
  dir=`basename $1`
  for file in "$1"/*
  do
    if [ -f "$file" ] && ( [ "${file##*.}" == "h" ] || [ "${file##*.}" == "all" ] )
    then
      file=`basename $file`
      echo "#include \"$dir/$file\"" >> $1.all
    fi
  done;
}
# Actualize .all files
createGeneric %{_prefix}/include/Bpp
exit 0

%preun -n libbpp-phyl-devel
removeGeneric() {
  if [ -f $1.all ]
  then
    echo "-- Remove generic include file: $1.all"
    rm $1.all
  fi
  for file in "$1"/*
  do
    if [ -d "$file" ]
    then
      # Recursion:
      removeGeneric $file
    fi
  done
}
# Actualize .all files
removeGeneric %{_prefix}/include/Bpp
exit 0

%postun -n libbpp-phyl9 -p /sbin/ldconfig

%postun -n libbpp-phyl-devel
createGeneric() {
  echo "-- Creating generic include file: $1.all"
  #Make sure we run into subdirectories first:
  dirs=()
  for file in "$1"/*
  do
    if [ -d "$file" ]
    then
      # Recursion:
      dirs+=( "$file" )
    fi
  done
  for dir in ${dirs[@]}
  do
    createGeneric $dir
  done
  #Now list all files, including newly created .all files:
  if [ -f $1.all ]
  then
    rm $1.all
  fi
  dir=`basename $1`
  for file in "$1"/*
  do
    if [ -f "$file" ] && ( [ "${file##*.}" == "h" ] || [ "${file##*.}" == "all" ] )
    then
      file=`basename $file`
      echo "#include \"$dir/$file\"" >> $1.all
    fi
  done;
}
# Actualize .all files
createGeneric %{_prefix}/include/Bpp
exit 0

%files -n libbpp-phyl9
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/%{_lib}/lib*.so.*

%files -n libbpp-phyl-devel
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/%{_lib}/lib*.so
%{_prefix}/%{_lib}/lib*.a
%{_prefix}/include/*

%changelog
* Thu Mar 07 2013 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.1.0-1
- New RateDistribution classes
- New models for protein sequences (COaLA)
- Support for gaps in parsimony score
- Improved and extended support for BppO
- Several bugs fixed and warnings removed
* Thu Feb 09 2012 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.0.3-1
- Reorganized model hierarchy
- New pairwise models
- Several bugs fixed
* Thu Jun 09 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.0.2-1
- New substitution models, new substitution mapping tools.
* Mon Feb 28 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.0.1-1
* Mon Feb 07 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr> 2.0.0-1
* Thu Mar 25 2010 Julien Dutheil <julien.dutheil@univ-montp2.fr> 1.9.0-1
* Wed Jun 10 2009 Julien Dutheil <jdutheil@birc.au.dk> 1.8.0-1
* Thu Dec 11 2008 Julien Dutheil <jdutheil@birc.au.dk> 1.7.0-1
* Wed Sep 24 2008 Julien Dutheil <jdutheil@birc.au.dk> 1.6.0-1
* Mon Jul 21 2008 Julien Dutheil <jdutheil@birc.au.dk> 1.5.1-1
* Fri Jan 18 2008 Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 1.5.0-1
* Sat Jul 06 2007 Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 1.4.0-1
* Sat Jan 27 2007 Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 1.3.1-1
* Fri Jan 19 2007 Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 1.3.0-1
* Mon Aug 28 2006 Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 1.2.0-1
* Tue Apr 18 2006 Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 1.1.0-1
* Fri Nov 16 2005 Julien Dutheil <Julien.Dutheil@univ-montp2.fr> 1.0.0-1
- First draft of the spec file

