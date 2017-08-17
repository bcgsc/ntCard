bin_PROGRAMS = ntcard nthll

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

nthll_CPPFLAGS = -I$(top_srcdir)/Common

nthll_SOURCES = \
        Common/Uncompress.cpp \
        Common/Uncompress.h \
        Common/SignalHandler.cpp \
        Common/SignalHandler.h \
        Common/Fcontrol.cpp \
        Common/Fcontrol.h \
        Common/nthash.hpp \
        Common/ntHashIterator.hpp \
        nthll.cpp

dist_doc_DATA = \
	ChangeLog \
	CITATION.bib \
	LICENSE \
	README.md

EXTRA_DIST=autogen.sh
