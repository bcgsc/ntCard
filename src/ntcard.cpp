#include "ntcard.hpp"
#include <btllib/nthash.hpp>

void
ntcard::NtCard::process(const std::string& seq)
{
	btllib::NtHash nthash(seq, 1, kmer_size);
	while (nthash.roll()) {
		update_counts(nthash.hashes()[0]);
		++total;
	}
}

void
ntcard::SeedNtCard::process(const std::string& seq)
{
	btllib::SeedNtHash nthash(seq, { seed }, 1, seed.size());
	while (nthash.roll()) {
		update_counts(nthash.hashes()[0]);
		++total;
	}
}

void
ntcard::NtCard::update_counts(const uint64_t hash_value)
{
	if (hash_value >> (64 - left_bits) == left_mask) {
		size_t shifted_value = hash_value & (r_buck - 1);
#pragma omp atomic
		++t_counter[r_buck + shifted_value];
	}
}

void
ntcard::NtCard::update_estimations()
{
	unsigned p[1 << sizeof(counter_t)];
	for (unsigned int& i : p)
		i = 0;
	for (size_t j = 0; j < r_buck; j++)
		++p[t_counter[j]];

	double mean_p[1 << sizeof(counter_t)];
	for (double& i : mean_p)
		i = 0.0;
	for (size_t i = 0; i < (1 << sizeof(counter_t)); i++) {
		mean_p[i] += p[i];
		mean_p[i] /= 1.0;
	}

	mean_f0 =
	    (ssize_t)((right_bits * log(2) - log(mean_p[0])) * 1.0 * ((size_t)1 << (left_bits + right_bits)));
	for (size_t i = 0; i < (1 << sizeof(counter_t)); i++)
		mean_f[i] = 0;
	mean_f[1] = -1.0 * mean_p[1] / (mean_p[0] * (log(mean_p[0]) - right_bits * log(2)));
	for (size_t i = 2; i < (1 << sizeof(counter_t)); i++) {
		double sum = 0.0;
		for (size_t j = 1; j < i; j++)
			sum += j * mean_p[i - j] * mean_f[j];
		mean_f[i] = -1.0 * mean_p[i] / (mean_p[0] * (log(mean_p[0]) - right_bits * log(2))) -
		            sum / (i * mean_p[0]);
	}
	for (size_t i = 1; i < (1 << sizeof(counter_t)); i++)
		mean_f[i] = abs((ssize_t)(mean_f[i] * mean_f0));
}

std::vector<size_t>
ntcard::NtCard::get_histogram(unsigned max_coverage)
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
