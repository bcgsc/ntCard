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

EXTRA_DIST= autogen.sh

SUBDIRS = \
        vendor \
        data

test-reads.fq.gz:
	wget https://www.bcgsc.ca/downloads/btl/ntcard/test-reads.fq.gz

test-rna-reads.fq.gz:
	wget https://www.bcgsc.ca/downloads/btl/ntcard/test-rna-reads.fq.gz

test_k12.hist: test-reads.fq.gz
	./ntcard -k 12 -p test test-reads.fq.gz

test-rna_k12.hist: test-rna-reads.fq.gz
	./ntcard -k 12 -p test-rna test-rna-reads.fq.gz

check-dna: test_k12.hist
	diff -q test_k12.hist $(srcdir)/data/test_k12.hist.good

check-rna: test-rna_k12.hist
	diff -q test-rna_k12.hist $(srcdir)/data/test_k12.hist.good

check: check-dna check-rna

# Check the C++ source code for white-space errors with clang-format.
clang-format:
	for i in ntcard.cpp; do clang-format -style=file $$i >$$i.fixed; done
	for i in ntcard.cpp; do diff -su $$i $$i.fixed && rm -f $$i.fixed; done
	if ls *.fixed; then exit 1; fi

clean-local:
        test-reads.fq.gz
        test_k12.hist
        test-rna-reads.fq.gz
        test-rna_k12.hist
