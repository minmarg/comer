INCLUDES = -I$(top_srcdir)/src \
	-I$(top_srcdir)/src/ext \
	-I$(top_srcdir)/src/liblib \
	$(all_includes)
METASOURCES = AUTO
noinst_LIBRARIES = libpro.a
AUTOMAKE_OPTIONS = subdir-objects
AM_CXXFLAGS = -DLOCALSTATEDIR=$(localstatedir)
libpro_a_SOURCES = \
	srcpro/MOptions.cpp \
	srcpro/Configuration.cpp \
	srcpro/Database.cpp \
	srcpro/DBProfileProbs.cpp \
	srcpro/DescriptionVector.cpp \
	srcpro/DistributionMatrix.cpp \
	srcpro/FrequencyStore.cpp \
	srcpro/GapScheme.cpp \
	srcpro/SEGProfile.cpp \
	srcpro/VirtScores.cpp \
	srcpro/CVS2Scores.cpp \
	srcpro/SSSScores.cpp \
	srcpro/Serializer.cpp \
	srcpro/datapro.cpp \
	\
	srcsco/AbstractScoreMatrix.cpp \
	srcsco/AbstractUniversalScoreMatrix.cpp \
	srcsco/AttributableScores.cpp \
	srcsco/AttributableScoresFPI.cpp \
	srcsco/AttributableScoresII.cpp \
	srcsco/ProfileMatrix.cpp \
	srcsco/ScoringMatrix.cpp \
	srcsco/LSOSMatrix.cpp \
	srcsco/HDP0ScoreMatrix.cpp \
	srcsco/HDP0CtxScoreMatrix.cpp \
	srcsco/AdjustedScoreMatrix.cpp \
	srcsco/UniversalScoreMatrix.cpp \
	\
	srcaln/AlignmentSimulation.cpp \
	srcaln/AlignmentProbe.cpp \
	srcaln/FastAlignment.cpp \
	srcaln/InputMultipleAlignment.cpp \
	srcaln/IMAClusters.cpp \
	srcaln/IMACountFiles.cpp \
	srcaln/IMACounts.cpp \
	srcaln/ProfileShuffler.cpp \
	srcaln/ProfileAlignment.cpp \
	srcaln/HybridAlignment.cpp \
	srcaln/MAAlignment.cpp \
	srcaln/ProfileSearching.cpp \
	\
	tfopt/TargetFreqOptimizerH.cpp

