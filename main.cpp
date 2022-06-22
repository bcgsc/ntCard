#include "ntcard.hpp"
#include <argparse/argparse.hpp>
#include <btllib/seq_reader.hpp>

#ifdef _OPENMP
#include <omp.h>
#endif

#define PROGRAM_NAME "ntCard"
#define PROGRAM_VERSION "1.2.2"
#define PROGRAM_DESCRIPTION "Estimates k-mer coverage histogram in input files"

int
main(int argc, char** argv)
{
	argparse::ArgumentParser args(PROGRAM_NAME, PROGRAM_VERSION);
	args.add_description(PROGRAM_DESCRIPTION);

	args.add_argument("-k", "--kmer-length")
	    .help("Length of the k-mers, ignored if using a spaced seed (see -s)")
	    .scan<'u', unsigned>();

	args.add_argument("-s", "--spaced-seed")
	    .help("Spaced seed pattern with 1's as cares and 0's as don't cares");

	args.add_argument("-t", "--num-threads")
	    .help("Number of threads")
	    .default_value((int)1)
	    .scan<'i', int>();

	args.add_argument("-c", "--max-coverage")
	    .help("Maximum coverage to estimate")
	    .required()
	    .scan<'u', unsigned>();

	args.add_argument("-l", "--left-bits")
	    .help("Number of bits to take from the left for sampling")
	    .default_value(11U)
	    .scan<'u', unsigned>();

	args.add_argument("-r", "--right-bits")
	    .help("Number of bits to take from the right as k-mer representations")
	    .default_value(27U)
	    .scan<'u', unsigned>();

	args.add_argument("--long-mode")
	    .help("Optimize sequence reader for long sequences (>5kb)")
	    .default_value(false)
	    .implicit_value(true);

	args.add_argument("files").remaining().required();

	try {
		args.parse_args(argc, argv);
	} catch (const std::runtime_error& err) {
		std::cerr << err.what() << std::endl;
		std::cerr << args;
		std::exit(1);
	}

	std::vector<std::string> files;
	try {
		files = args.get<std::vector<std::string> >("files");
	} catch (std::logic_error& e) {
		std::cout << "No files provided" << std::endl;
		std::exit(1);
	}

	auto num_threads = args.get<int>("-t");
#ifdef _OPENMP
	omp_set_num_threads(std::max(num_threads, (int)files.size()));
#endif

	auto max_coverage = args.get<unsigned>("-c");
	auto left_bits = args.get<unsigned>("-l");
	auto right_bits = args.get<unsigned>("-r");

	ntcard::NtCard* ntc;
	if (args.is_used("-s")) {
		ntc = new ntcard::SeedNtCard(args.get("-s"), left_bits, right_bits);
	} else if (args.is_used("-k")) {
		ntc = new ntcard::NtCard(args.get<unsigned>("-k"), left_bits, right_bits);
	} else {
		std::cerr << "Please specify k-mer length (-k) or spaced seed pattern (-s)" << std::endl;
		std::cerr << args;
		std::exit(1);
	}

	unsigned seq_reader_flags;
	if (args.get<bool>("--long-mode")) {
		seq_reader_flags = btllib::SeqReader::Flag::LONG_MODE;
	} else {
		seq_reader_flags = btllib::SeqReader::Flag::SHORT_MODE;
	}

	for (const auto& file : files) {
		btllib::SeqReader reader(file, seq_reader_flags);
		for (const auto& record : reader) {
			if ((args.is_used("-s") && record.seq.size() >= args.get("-s").size()) ||
			    args.is_used("-k") && record.seq.size() >= args.get<unsigned>("-k")) {
				ntc->process(record.seq);
			}
		}
	}

	auto hist = ntc->get_histogram(max_coverage);
	std::cout << "F1\t" << hist[0] << std::endl;
	std::cout << "F0\t" << hist[1] << std::endl;
	for (size_t i = 2; i < hist.size(); i++) {
		std::cout << i - 1 << "\t" << hist[i] << std::endl;
	}

	return 0;
}