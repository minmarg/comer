INCLUDES =-I$(top_srcdir)/src  \
	-I$(top_srcdir)/src/liblib \
	-I$(top_srcdir)/src/lmpi \
	$(all_includes)
METASOURCES = AUTO
AM_CXXFLAGS = -DLOCALSTATEDIR=$(localstatedir)
bin_SCRIPTS = mpiscaler.sh
bin_PROGRAMS = mpiscaler
mpiscaler_SOURCES = \
	PScaler.cpp \
	ParallelUniversalScoreMatrix.cpp \
	PrlLSOUniversalScoreMatrix.cpp \
	PrlHDPUniversalScoreMatrix.cpp \
	mpiscaler.cpp
mpiscaler_LDADD = \
	$(top_builddir)/src/libpro/libpro.a \
	$(top_builddir)/src/libHDP/libHDP.a \
	$(top_builddir)/src/libseg/libseg.a \
	$(top_builddir)/src/libhsh/libhsh.a \
	$(top_builddir)/src/liblib/liblib.a \
	$(top_builddir)/src/ext/libpsl.a \
	$(top_builddir)/src/lmpi/libpmpi.a \
	-lmpi_cxx -lpthread
#	-lmpich -lpthread
