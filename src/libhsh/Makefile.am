INCLUDES = -I$(top_srcdir)/src \
	-I$(top_srcdir)/src/liblib \
	$(all_includes)
METASOURCES = AUTO
noinst_LIBRARIES = libhsh.a
AM_CXXFLAGS = -DLOCALSTATEDIR=$(localstatedir)
libhsh_a_SOURCES = \
	CRCHashing.cpp \
	HashFunctions.cpp \
	HashTable.cpp \
	MD5Hashing.cpp \
	RJHashing.cpp \
	SBoxHashing.cpp
