
#include <btllib/seq_reader.hpp>
#include <fstream>
#include <omp.h>

#include "args.hpp"
#include "ntcard/ntcard.hpp"

int
main(int argc, char* argv[])
{
	auto args = ProgramArguments(argc, argv);

	omp_set_num_threads(args.num_threads);

	unsigned seq_reader_flags;
	if (args.seq_reader_long_mode) {
		seq_reader_flags = btllib::SeqReader::Flag::LONG_MODE;
	} else {
		seq_reader_flags = btllib::SeqReader::Flag::SHORT_MODE;
	}

	ntcard::NtCard* ntc;
	if (args.spaced_seed == "") {
		ntc = new ntcard::NtCard(args.kmer_length, args.left_bits, args.right_bits);
	} else {
		ntc = new ntcard::SeedNtCard(args.spaced_seed, args.left_bits, args.right_bits);
	}

#pragma omp parallel for
	for (const auto& file : args.input_files) {
		btllib::SeqReader reader(file, seq_reader_flags);
		for (const auto& record : reader) {
			if (record.seq.size() >= args.kmer_length) {
				ntc->process(record.seq);
			}
		}
	}
	auto hist = ntc->get_histogram(args.max_coverage);

	std::ofstream out(args.output_path);
	if (!args.output_jellyfish) {
		out << "F1\t" << hist[0] << std::endl;
		out << "F0\t" << hist[1] << std::endl;
	}
	for (size_t i = 2; i < hist.size(); i++) {
		out << i - 1 << "\t" << hist[i] << std::endl;
	}

	return 0;
}