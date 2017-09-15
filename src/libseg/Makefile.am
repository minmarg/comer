INCLUDES = -I$(top_srcdir)/src \
	-I$(top_srcdir)/src/liblib \
	$(all_includes)
METASOURCES = AUTO
noinst_LIBRARIES = libseg.a
AM_CXXFLAGS = -DLOCALSTATEDIR=$(localstatedir)
libseg_a_SOURCES = \
	segdata.cpp \
	SegmentStructure.cpp \
	SEGAbstract.cpp \
	SEGSequence.cpp
