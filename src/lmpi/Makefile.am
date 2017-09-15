INCLUDES = -I$(top_srcdir)/src \
	-I$(top_srcdir)/src/liblib \
	$(all_includes)
METASOURCES = AUTO
AM_CXXFLAGS = -DLOCALSTATEDIR=$(localstatedir)
noinst_LIBRARIES = libpmpi.a
libpmpi_a_SOURCES = MessageDispatcher.cpp
libpmpi_a_LIBADD = \
	$(top_builddir)/src/liblib/liblib.a
