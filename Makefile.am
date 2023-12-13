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
        ntcard.cpp

nthll_CPPFLAGS = -I$(top_srcdir)/Common

nthll_SOURCES = \
        Common/Uncompress.cpp \
        Common/Uncompress.h \
        Common/SignalHandler.cpp \
        Common/SignalHandler.h \
        Common/Fcontrol.cpp \
        Common/Fcontrol.h \
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

test-reads.fa: test-reads.fq.gz
	zcat test-reads.fq.gz |sed -n '1~4s/^@/>/p;2~4p' > test-reads.fa

test_k12.hist: test-reads.fq.gz
	./ntcard -k 12 -p test test-reads.fq.gz

test-rna_k12.hist: test-rna-reads.fq.gz
	./ntcard -k 12 -p test-rna test-rna-reads.fq.gz

test-fa_k12.hist: test-reads.fa
	./ntcard -k 12 -p test-fa test-reads.fa

test-gap_k12.hist: test-reads.fq.gz
	./ntcard -k 12 -p test-gap -g 2 test-reads.fq.gz

test-rna-gap_k12.hist: test-rna-reads.fq.gz
	./ntcard -k 12 -p test-rna-gap -g 2 test-rna-reads.fq.gz

test-fa-gap_k12.hist: test-reads.fa
	./ntcard -k 12 -p test-fa-gap -g 2 test-reads.fa

check-dna: test_k12.hist
	diff -q test_k12.hist $(srcdir)/data/test_k12.hist.good

check-dna-gap: test-gap_k12.hist
	diff -q test-gap_k12.hist $(srcdir)/data/test-gap_k12.hist.good

check-rna: test-rna_k12.hist
	diff -q test-rna_k12.hist $(srcdir)/data/test_k12.hist.good

check-rna-gap: test-rna-gap_k12.hist
	diff -q test-rna-gap_k12.hist $(srcdir)/data/test-gap_k12.hist.good

check-fa: test-fa_k12.hist
	diff -q test-fa_k12.hist $(srcdir)/data/test_k12.hist.good

check-fa-gap: test-fa-gap_k12.hist
	diff -q test-fa-gap_k12.hist $(srcdir)/data/test-gap_k12.hist.good

check: check-dna check-rna check-fa check-dna-gap check-rna-gap check-fa-gap

# Check the C++ source code for white-space errors with clang-format.
clang-format:
	for i in ntcard.cpp; do clang-format -style=file $$i >$$i.fixed; done
	for i in ntcard.cpp; do diff -su $$i $$i.fixed && rm -f $$i.fixed; done
	if ls *.fixed; then exit 1; fi

clean-local:
	$(RM) test-reads.fq.gz
	$(RM) test_k12.hist
	$(RM) test-gap_k12.hist
	$(RM) test-rna-reads.fq.gz
	$(RM) test-rna_k12.hist
	$(RM) test-rna-gap_k12.hist
	$(RM) test-reads.fa
	$(RM) test-fa_k12.hist
	$(RM) test-fa-gap_k12.hist
