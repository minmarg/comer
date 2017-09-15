INCLUDES = -I$(top_srcdir)/src \
	-I$(top_srcdir)/src/ext \
	-I$(top_srcdir)/src/liblib \
	-I$(top_srcdir)/src/lmpi \
	-I$(top_srcdir)/src/libHDP \
	$(all_includes)
METASOURCES = AUTO
AM_CXXFLAGS = -DLOCALSTATEDIR=$(localstatedir)
bin_SCRIPTS = hdpclust.sh hdpclust.hyd.sh
bin_PROGRAMS = fsample hdpsim hdpclust hdprior
hdprior_SOURCES = \
	hdprior.cpp \
	HDPCalcPrior.cpp
hdpclust_SOURCES =  \
	hdpclust.cpp \
	HDPsampler.cpp \
	HDPsampler_rw.cpp \
	HDPsampler_JNMCMC.cpp \
	HDPsampler_SMCMC.cpp \
	HDPsampler_master.cpp \
	HDPsampler_slave.cpp
hdpsim_SOURCES = \
	hdpsim.cpp \
	HDPSmplGenerator.cpp
fsample_SOURCES = \
	ProFreqWriter.cpp \
	fsample.cpp
hdprior_LDADD = \
	$(top_builddir)/src/libHDP/libHDP.a \
	$(top_builddir)/src/libpro/libpro.a \
	$(top_builddir)/src/libhsh/libhsh.a \
	$(top_builddir)/src/liblib/liblib.a \
	$(top_builddir)/src/ext/libpsl.a 
hdpclust_LDADD = \
	$(top_builddir)/src/libHDP/libHDP.a \
	$(top_builddir)/src/liblib/liblib.a \
	$(top_builddir)/src/root/libroot.a \
	$(top_builddir)/src/smpl/libsmpl.a \
	$(top_builddir)/src/ext/libpsl.a \
        $(top_builddir)/src/lmpi/libpmpi.a \
	-lmpi_cxx -lpthread
#	-lmpich -lpthread
hdpsim_LDADD = \
	$(top_builddir)/src/libHDP/libHDP.a \
	$(top_builddir)/src/libpro/libpro.a \
	$(top_builddir)/src/libhsh/libhsh.a \
	$(top_builddir)/src/liblib/liblib.a \
	$(top_builddir)/src/ext/libpsl.a
fsample_LDADD = \
	$(top_builddir)/src/libpro/libpro.a \
	$(top_builddir)/src/libHDP/libHDP.a \
	$(top_builddir)/src/libhsh/libhsh.a \
	$(top_builddir)/src/liblib/liblib.a \
	$(top_builddir)/src/ext/libpsl.a \
        $(top_builddir)/src/lmpi/libpmpi.a \
	-lmpi_cxx -lpthread
#	-lmpich -lpthread
