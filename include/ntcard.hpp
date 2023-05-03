#ifndef NTCARD_NTCARD_HPP
#define NTCARD_NTCARD_HPP

#include <btllib/nthash.hpp>
#include <string>
#include <vector>

namespace ntcard {

using counter_t = uint16_t;

class NtCard
{
  private:
	unsigned kmer_size;
	unsigned left_bits, right_bits;
	counter_t* t_counter;
	size_t r_buck;
	unsigned left_mask;
	double mean_f0;
	double* mean_f;

	void update_estimations()
	{
		size_t arr_size = 1 << sizeof(counter_t) * 8;
		unsigned* p = new unsigned[arr_size];
		double* mean_p = new double[arr_size];
		for (unsigned i = 0; i < arr_size; i++) {
			p[i] = 0;
			mean_p[i] = 0;
		}
		for (size_t j = 0; j < r_buck; j++) {
			++p[t_counter[j]];
		}
		for (size_t i = 0; i < arr_size; i++) {
			mean_p[i] += p[i];
			mean_p[i] /= 1.0;
		}

		mean_f0 =
		    (ssize_t)((right_bits * log(2) - log(mean_p[0])) * 1.0 * ((size_t)1 << (left_bits + right_bits)));
		for (size_t i = 0; i < arr_size; i++)
			mean_f[i] = 0;
		mean_f[1] = -1.0 * mean_p[1] / (mean_p[0] * (log(mean_p[0]) - right_bits * log(2)));
		for (size_t i = 2; i < arr_size; i++) {
			double sum = 0.0;
			for (size_t j = 1; j < i; j++)
				sum += j * mean_p[i - j] * mean_f[j];
			mean_f[i] = -1.0 * mean_p[i] / (mean_p[0] * (log(mean_p[0]) - right_bits * log(2))) -
			            sum / (i * mean_p[0]);
		}
		for (size_t i = 1; i < arr_size; i++)
			mean_f[i] = abs((ssize_t)(mean_f[i] * mean_f0));
	}

  protected:
	size_t total;

	void update_counts(uint64_t hash_value)
	{
		if (hash_value >> (64 - left_bits) == left_mask) {
			size_t shifted_value = hash_value & (r_buck - 1);
#pragma omp atomic
			++t_counter[shifted_value];
		}
	}

  public:
	explicit NtCard(
	    const unsigned kmer_size,
	    const unsigned left_bits = 11,
	    const unsigned right_bits = 27)
	  : kmer_size(kmer_size)
	  , left_bits(left_bits)
	  , right_bits(right_bits)
	  , r_buck(((size_t)1) << right_bits)
	  , left_mask((((size_t)1) << (left_bits - 1)) - 1)
	  , mean_f0(0.0)
	  , mean_f(new double[1 << sizeof(counter_t) * 8])
	  , total(0)
	{
		t_counter = new counter_t[r_buck]();
	}

	virtual void process(const std::string& seq)
	{
		btllib::NtHash nthash(seq, 1, kmer_size);
		while (nthash.roll()) {
			update_counts(nthash.hashes()[0]);
			++total;
		}
	}

	std::vector<size_t> get_histogram(unsigned max_coverage)
	{
		update_estimations();
		std::vector<size_t> hist;
		hist.reserve(max_coverage + 1);
		hist.push_back(total);
		hist.push_back((size_t)mean_f0);
		for (size_t i = 1; i <= max_coverage; i++) {
			hist.push_back((size_t)mean_f[i]);
		}
		return hist;
	}
};

class SeedNtCard : public NtCard
{
  private:
	std::string seed;

  public:
	explicit SeedNtCard(
	    const std::string& seed,
	    const unsigned left_bits = 7,
	    const unsigned right_bits = 27)
	  : NtCard(seed.size(), left_bits, right_bits)
	  , seed(seed)
	{
	}

	void process(const std::string& seq) override
	{
		btllib::SeedNtHash nthash(seq, { seed }, 1, seed.size());
		while (nthash.roll()) {
			update_counts(nthash.hashes()[0]);
			++total;
		}
	}
};

}

#endif // NTCARD_NTCARD_HPP
