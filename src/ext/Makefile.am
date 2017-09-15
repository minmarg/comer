INCLUDES = -I$(top_srcdir)/src \
	$(all_includes)
METASOURCES = AUTO
noinst_LIBRARIES = libpsl.a
AUTOMAKE_OPTIONS = subdir-objects
AM_CXXFLAGS = -DLOCALSTATEDIR=$(localstatedir)
libpsl_a_SOURCES = \
	cheb.cpp \
	log.cpp \
	exp.cpp \
	erf.cpp \
	expint.cpp \
	zeta.cpp \
	gamma.cpp \
	gammainc.cpp \
	beta.cpp \
	betainc.cpp \
	digamma.cpp \
	rng.cpp \
	\
	rv/rvnorm.cpp \
	rv/rvexp.cpp \
	rv/rvgamma.cpp \
	rv/rvchi2.cpp \
	rv/rvbeta.cpp \
	rv/rvt.cpp \
	rv/rvdisc.cpp \
	rv/rmvnorm.cpp \
	rv/rmvt.cpp \
	rv/rmvdir.cpp \
	\
	ivector.cpp \
	pslvector.cpp \
	pslmatrix.cpp \
	nsmatrix.cpp \
	spdmatrix.cpp
noinst_HEADERS = pslcodes.h
