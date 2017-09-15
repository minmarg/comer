INCLUDES = -I$(top_srcdir)/src \
	-I$(top_srcdir)/src/ext \
	-I$(top_srcdir)/src/liblib \
	$(all_includes)
METASOURCES = AUTO
noinst_LIBRARIES = libHDP.a
AM_CXXFLAGS = -DLOCALSTATEDIR=$(localstatedir)
libHDP_a_SOURCES = \
	basin.cpp \
	dish.cpp \
	menu.cpp \
	table.cpp \
	restaurant.cpp \
	rntchain.cpp \
	HDPscores.cpp \
	HDP_SSSScores.cpp \
	iHDP_SSSScores.cpp \
	HDPbase.cpp \
	HDPbase_rw.cpp \
	HDPbase_probvec.cpp \
	HDPbase_probmtx.cpp 
