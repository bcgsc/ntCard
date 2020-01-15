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

test-data.tar.gz:
	wget https://www.bcgsc.ca/sites/default/files/bioinformatics/software/abyss/releases/1.3.4/test-data.tar.gz

test_k12.hist: test-data.tar.gz
	./ntcard -k 12 -p test test-data.tar.gz

check: test_k12.hist
	diff -q test_k12.hist data/test_k12.hist.good
