INCLUDES = -I$(top_srcdir)/src \
	-I$(top_srcdir)/src/liblib \
	$(all_includes)
METASOURCES = AUTO
AM_CXXFLAGS = -DLOCALSTATEDIR=$(localstatedir)
bin_SCRIPTS = makepro.sh ssp.pl inssp.pl ssp2.pl inssp2.pl
bin_PROGRAMS = makepro select
makepro_SOURCES = makepro.cpp
makepro_LDADD = \
	$(top_builddir)/src/libpro/libpro.a \
	$(top_builddir)/src/libHDP/libHDP.a \
	$(top_builddir)/src/libseg/libseg.a \
	$(top_builddir)/src/libhsh/libhsh.a \
	$(top_builddir)/src/liblib/liblib.a \
	$(top_builddir)/src/sort/libsort.a \
	$(top_builddir)/src/ext/libpsl.a
select_SOURCES = select.cpp
select_LDADD = \
	$(top_builddir)/src/libpro/libpro.a \
	$(top_builddir)/src/libHDP/libHDP.a \
	$(top_builddir)/src/libseg/libseg.a \
	$(top_builddir)/src/libhsh/libhsh.a \
	$(top_builddir)/src/liblib/liblib.a \
	$(top_builddir)/src/sort/libsort.a \
	$(top_builddir)/src/ext/libpsl.a
