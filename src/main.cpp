
#include <btllib/seq_reader.hpp>
#include <fstream>
#include <omp.h>

#include "args.hpp"
#include "ntcard/ntcard.hpp"

int
main(int argc, char* argv[])
{
	auto args = ProgramArguments(argc, argv);

	if (args.verbose) {
		args.print();
	}

	std::cout << "Initializing... " << std::flush;
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
	std::cout << "DONE" << std::endl;

	std::cout << "Processing data... " << std::endl;
#pragma omp parallel for
	for (auto& path : args.input_files) {
		std::string file_name = path.substr(path.find_last_of("/") + 1);
		btllib::SeqReader reader(path, seq_reader_flags);
		for (const auto& record : reader) {
			if (record.seq.size() >= args.kmer_length) {
				ntc->process(record.seq);
				if (args.verbose) {
#pragma omp critical
					std::cout << "[" << file_name << "] DONE: " << record.id << " ("
					          << record.seq.size() << "bps)" << std::endl;
				}
			}
		}
		if (!args.verbose) {
#pragma omp critical
			std::cout << "  - DONE: " << file_name << std::endl;
		}
	}

	std::cout << "Computing histogram... " << std::flush;
	auto hist = ntc->get_histogram(args.max_coverage);
	std::cout << "DONE" << std::endl;

	std::cout << "Writing output... " << std::flush;
	std::ofstream out(args.output_path);
	if (!args.output_jellyfish) {
		out << "F1\t" << hist[0] << std::endl;
		out << "F0\t" << hist[1] << std::endl;
	}
	for (size_t i = 2; i < hist.size(); i++) {
		out << i - 1 << "\t" << hist[i] << std::endl;
	}
	std::cout << "DONE" << std::endl;

	return 0;
}