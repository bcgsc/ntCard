#ifndef ARGS_HPP
#define ARGS_HPP

#include <argparse/argparse.hpp>
#include <string>
#include <vector>

#include "ntcard/version.hpp"

struct ProgramArguments
{
	unsigned kmer_length, num_threads, max_coverage, left_bits, right_bits;
	std::string spaced_seed, output_path;
	bool seq_reader_long_mode, output_jellyfish;
	std::vector<std::string> input_files;

	ProgramArguments(int argc, char* argv[])
	{
		argparse::ArgumentParser parser("ntcard", ntcard::VERSION);
		parser.add_description("Tool for estimating the k-mer coverage histogram of genomic data");

		parser.add_argument("-k", "--kmer-length")
		    .help("Length of the k-mers, ignored if using a spaced seed (see -s)")
		    .scan<'u', unsigned>();

		parser.add_argument("-s", "--spaced-seed")
		    .help("Spaced seed pattern with 1's as cares and 0's as don't cares");

		parser.add_argument("-c", "--max-coverage")
		    .help("Maximum coverage to estimate")
		    .required()
		    .scan<'u', unsigned>();

		parser.add_argument("-o", "--out-path").help("Path to output file").required();

		parser.add_argument("-t", "--num-threads")
		    .help("Number of threads")
		    .default_value((int)1)
		    .scan<'i', int>();

		parser.add_argument("-l", "--left-bits")
		    .help("Number of bits to take from the left for sampling")
		    .default_value(7U)
		    .scan<'u', unsigned>();

		parser.add_argument("-r", "--right-bits")
		    .help("Number of bits to take from the right as k-mer representations")
		    .default_value(27U)
		    .scan<'u', unsigned>();

		parser.add_argument("--long-mode")
		    .help("Optimize file reader for long sequences (>5kb)")
		    .default_value(false)
		    .implicit_value(true);

		parser.add_argument("--jellyfish")
		    .help("Output in jellyfish format")
		    .default_value(false)
		    .implicit_value(true);

		parser.add_argument("files").remaining().required();

		try {
			parser.parse_args(argc, argv);
		} catch (const std::runtime_error& err) {
			std::cerr << err.what() << std::endl;
			std::cerr << parser;
			std::exit(1);
		}

		if (parser.is_used("-s")) {
			spaced_seed = parser.get("-s");
			kmer_length = spaced_seed.size();
		} else if (parser.is_used("-k")) {
			spaced_seed = "";
			kmer_length = parser.get<unsigned>("-k");
		} else {
			std::cerr << "Please specify k-mer length (-k) or spaced seed pattern (-s)"
			          << std::endl;
			std::cerr << parser;
			std::exit(1);
		}

		max_coverage = parser.get<unsigned>("-c");
		output_path = parser.get("-o");
		num_threads = parser.get<int>("-t");
		left_bits = parser.get<unsigned>("-l");
		right_bits = parser.get<unsigned>("-r");
		seq_reader_long_mode = parser.get<bool>("--long-mode");
		output_jellyfish = parser.get<bool>("--jellyfish");

		try {
			input_files = parser.get<std::vector<std::string> >("files");
		} catch (std::logic_error& e) {
			std::cout << "No files provided" << std::endl;
			std::exit(1);
		}
	}
};

#endif