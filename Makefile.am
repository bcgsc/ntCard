bin_PROGRAMS = ntcard

ntcard_CPPFLAGS = -I$(top_srcdir)/Common

ntcard_SOURCES = \
        Common/Uncompress.cpp \
        Common/Uncompress.h \
        Common/SignalHandler.cpp \
        Common/SignalHandler.h \
        Common/StringUtil.h \
        Common/Fcontrol.cpp \
        Common/Fcontrol.h \
        Common/nthash.hpp \
        Common/ntHashIterator.hpp \
        ntcard.cpp


dist_doc_DATA = \
	ChangeLog \
	CITATION.bib \
	LICENSE \
	README.md

EXTRA_DIST=autogen.sh
