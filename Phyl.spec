%define name bpp-phyl
%define version 2.0.2
%define release 1
%define _prefix /usr

Summary: The Bio++ PhylLib library.
Name: %{name}
Version: %{version}
Release: %{release}
Vendor: The Bio++ Project
Source: http://biopp.univ-montp2.fr/Repositories/sources/%{name}-%{version}.tar.gz
License: CeCILL 2
Group: System Environment/Libraries
BuildRoot: %{_builddir}/%{name}-root
Packager: Julien Dutheil
AutoReq: yes
AutoProv: yes

%description
This library contains utilitary and classes for phylogenetics and molecular evolution analysis.
It is part of the Bio++ project.

%package devel
Summary: Libraries, includes to develop applications with %{name}.
Group: Development/Libraries
Requires: %{name} = %{version}

%description devel
The %{name}-devel package contains the header files and static libraries for
building applications which use %{name}.

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
rm -rf $RPM_BUILD_ROOT
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/%{_lib}/lib*.so
%{_prefix}/%{_lib}/lib*.so.*

%files devel
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/%{_lib}/lib*.a
%{_prefix}/include/*

%changelog
* Thu Jun 09 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr>
- Version 2.0.2. New substitution models, new substitution mapping tools.
* Mon Feb 28 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr>
- Version 2.0.1
* Mon Feb 07 2011 Julien Dutheil <julien.dutheil@univ-montp2.fr>
- Version 2.0.0
* Thu Mar 25 2010 Julien Dutheil <julien.dutheil@univ-montp2.fr>
- Version 1.9.0
* Wed Jun 10 2009 Julien Dutheil <jdutheil@birc.au.dk>
- Version 1.8.0
* Thu Dec 11 2008 Julien Dutheil <jdutheil@birc.au.dk>
- Version 1.7.0
* Wed Sep 24 2008 Julien Dutheil <jdutheil@birc.au.dk>
- Version 1.6.0
* Mon Jul 21 2008 Julien Dutheil <jdutheil@birc.au.dk>
- Version 1.5.1
* Fri Jan 18 2008 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- Version 1.5.0
* Sat Jul 06 2007 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- Version 1.4.0
i Sat Jan 27 2007 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- Version 1.3.1
* Fri Jan 19 2007 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- Version 1.3.0
* Mon Aug 28 2006 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- Version 1.2.0
* Tue Apr 18 2006 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- Version 1.1.0
* Fri Nov 16 2005 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- First draft of the spec file

