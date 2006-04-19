%define name Bpp-Phyl
%define version 1.1.0
%define release 1
%define _prefix /usr/local

Summary: The Bio++ PhylLib library.
Name: %{name}
Version: %{version}
Release: %{release}
Vendor: The Bio++ Project
Source: http://kimura.univ-montp2.fr/BioPP/Download/Sources/%{name}-%{version}.tar.gz
License: CeCILL 2
Group: System Environment/Libraries
BuildRoot: %{_builddir}/%{name}-root
Packager: Julien Dutheil
AutoReqProv: no
Requires: libstdc++6
Requires: Bpp-Utils >= 1.0.0
Requires: Bpp-NumCalc >= 1.0.0
Requires: Bpp-Seq >= 1.1.0

%description
This library contains utilitary and classes for phylogenetics and molecular evolution analysis.
It is part of the Bio++ project.

%package devel
Summary: Libraries, includes to develop applications with %{name}.
Group: Development/Libraries
Requires: %{name} = %{version}
Requires: Bpp-Utils-devel >= 1.0.0 
Requires: Bpp-NumCalc-devel >= 1.0.0
Requires: Bpp-Seq-devel >= 1.1.0

%description devel
The %{name}-devel package contains the header files and static libraries for
building applications which use %{name}.

%prep
%setup -q

%build
CFLAGS="$RPM_OPT_FLAGS" ./configure --prefix=%{_prefix}
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
%doc AUTHORS COPYING INSTALL NEWS README ChangeLog
%{_prefix}/lib/lib*.so
%{_prefix}/lib/lib*.so.*

%files devel
%defattr(-,root,root)
%doc AUTHORS COPYING INSTALL NEWS README ChangeLog
%{_prefix}/lib/lib*.a
%{_prefix}/lib/lib*.la
%{_prefix}/include/*

%changelog
* Tue Apr 18 2006 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- Version 1.1.0
* Fri Nov 16 2005 Julien Dutheil <Julien.Dutheil@univ-montp2.fr>
- First draft of the spec file

