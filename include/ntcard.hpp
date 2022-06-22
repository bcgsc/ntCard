#ifndef NTCARD_NTCARD_HPP
#define NTCARD_NTCARD_HPP

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

	void update_estimations();

  protected:
	size_t total;

	void update_counts(uint64_t hash_value);

  public:
	explicit NtCard(
	    const unsigned kmer_size,
	    const unsigned left_bits = 11,
	    const unsigned right_bits = 27)
	  : kmer_size(kmer_size)
	  , left_bits(left_bits)
	  , right_bits(right_bits)
	  , total(0)
	  , r_buck(((size_t)1) << right_bits)
	  , left_mask((((size_t)1) << (left_bits - 1)) - 1)
	  , mean_f0(0.0)
	  , mean_f(new double[1 << sizeof(counter_t) * 8])
	{
		t_counter = new counter_t[r_buck]();
	}

	virtual void process(const std::string& seq);

	std::vector<size_t> get_histogram(unsigned max_coverage);
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

	void process(const std::string& seq) override;
};

}

#endif // NTCARD_NTCARD_HPP
