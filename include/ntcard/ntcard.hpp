#ifndef NTCARD_NTCARD_HPP
#define NTCARD_NTCARD_HPP

#include <btllib/nthash.hpp>
#include <string>
#include <vector>

#include "version.hpp"

namespace ntcard {

class NtCard
{
  private:
	unsigned kmer_size;
	unsigned left_bits, right_bits, left_mask;
	std::vector<uint16_t> counts;
	unsigned num_total;
	std::vector<unsigned> histogram;

	void update_estimations()
	{
		double* p = new double[histogram.size()];
		double* f = new double[histogram.size()];
		std::memset(p, 0, histogram.size() * sizeof(double));
		std::memset(f, 0, histogram.size() * sizeof(double));
		for (const auto& t_i : counts) {
			++p[t_i];
		}
		double log_p0 = log(p[0]) - (double)right_bits * log(2);
		f[0] = -log_p0 * ((uint64_t)1 << (left_bits + right_bits));
		for (unsigned i = 1; i < histogram.size(); i++) {
			double sum = 0;
			for (unsigned j = 1; j < i; j++) {
				sum += j * p[i - j] * f[j];
			}
			f[i] = -p[i] / (p[0] * log_p0) - sum / (i * p[0]);
		}
		histogram[0] = abs(f[0]);
		for (unsigned i = 1; i < histogram.size(); i++) {
			histogram[i] = abs(f[i] * f[0]);
		}
	}

  protected:
	void update_counts(uint64_t hash_value)
	{
		if (hash_value >> (64 - left_bits) == left_mask) {
			size_t shifted_value = hash_value & (counts.size() - 1);
			++counts[shifted_value];
		}
		++num_total;
	}

  public:
	explicit NtCard(
	    const unsigned kmer_size,
	    const unsigned left_bits = 11,
	    const unsigned right_bits = 27)
	  : kmer_size(kmer_size)
	  , left_bits(left_bits)
	  , right_bits(right_bits)
	  , left_mask((((size_t)1) << (left_bits - 1)) - 1)
	  , counts(1U << right_bits, 0)
	  , num_total(0)
	  , histogram(1 << sizeof(decltype(counts)::value_type) * 8, 0)
	{
	}

	virtual void process(const std::string& seq)
	{
		btllib::NtHash h(seq, 1, kmer_size);
		while (h.roll()) {
			update_counts(h.hashes()[0]);
		}
	}

	std::vector<unsigned> get_histogram(unsigned max_coverage)
	{
		update_estimations();
		std::vector<unsigned> hist;
		hist.reserve(max_coverage + 1);
		hist.push_back(num_total);
		for (size_t i = 0; i <= max_coverage + 1; i++) {
			hist.push_back(histogram[i]);
		}
		return hist;
	}
};

class SeedNtCard : public NtCard
{
  private:
	std::vector<std::string> seed;

  public:
	explicit SeedNtCard(
	    const std::string& seed,
	    const unsigned left_bits = 11,
	    const unsigned right_bits = 27)
	  : NtCard(seed.size(), left_bits, right_bits)
	  , seed(1, seed)
	{
	}

	void process(const std::string& seq) override
	{
		btllib::SeedNtHash h(seq, seed, 1, seed[0].size());
		while (h.roll()) {
			update_counts(h.hashes()[0]);
		}
	}
};

}

#endif // NTCARD_NTCARD_HPP
