%define _basename traitrelax
%define _version 1.0.0
%define _release 1
%define _prefix /usr

Name: %{_basename}
Version: %{_version}
Release: %{_release}
License: CECILL-2.0
Vendor: Julien Dutheil
Source: https://github.com/halabikeren/TraitRELAX/%{_basename}-%{_version}.tar.gz
Summary: The TraitRELAX program
Group: Productivity/Scientific/Other

Requires: libbpp-phyl9 = 2.4.1
Requires: libbpp-seq9 = 2.4.1
Requires: libbpp-core2 = 2.4.1

BuildRoot: %{_builddir}/%{_basename}-root
BuildRequires: cmake >= 2.6.0
BuildRequires: gcc-c++ >= 4.0.0
BuildRequires: groff
BuildRequires: texinfo >= 4.0.0
BuildRequires: libbpp-core2 = 2.4.1
BuildRequires: libbpp-core-devel = 2.4.1
BuildRequires: libbpp-seq9 = 2.4.1
BuildRequires: libbpp-seq-devel = 2.4.1
BuildRequires: libbpp-phyl9 = 2.4.1
BuildRequires: libbpp-phyl-devel = 2.4.1


AutoReq: yes
AutoProv: yes
%if 0%{?mdkversion} >= 201100 || %{?distribution} == "Mageia"
BuildRequires: xz
%define compress_program xz
%else
%if 0%{?mdkversion}
BuildRequires: lzma
%define compress_program lzma
%else
BuildRequires: gzip
%define compress_program gzip
%endif
%endif

%description
TraitRELAX - program for detection of  association of changes in phenotype traits with changes in selection intensity at the codon level across a phylogeny

%prep
%setup -q

%build
CFLAGS="$RPM_OPT_FLAGS"
CMAKE_FLAGS="-DCMAKE_INSTALL_PREFIX=%{_prefix} -DCOMPRESS_PROGRAM=%{compress_program}"
cmake $CMAKE_FLAGS .
make

%install
make DESTDIR=$RPM_BUILD_ROOT install

%clean
rm -rf $RPM_BUILD_ROOT

%post -p /sbin/ldconfig

%postun -p /sbin/ldconfig

%files
%defattr(-,root,root)
%doc AUTHORS.txt COPYING.txt INSTALL.txt ChangeLog
%{_prefix}/bin/*
%{_prefix}/share/info/*.info*
%{_prefix}/share/man/man1/*.1*

%changelog
* Thur Apr 30 2020 Keren Halabi <halabikeren@mail.tau.ac.il> 1.0.0
- First upload.

