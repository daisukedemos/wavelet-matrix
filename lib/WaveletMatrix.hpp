#ifndef WAVELET_MATRIX_WAVELET_MATRIX_HPP_
#define WAVELET_MATRIX_WAVELET_MATRIX_HPP_

#include <vector>
#include <stdint.h>
#include <rsdic/RSDic.hpp>

namespace wavelet_matrix{

struct Range{
	Range(uint32_t bpos, uint32_t epos) : bpos(bpos), epos(epos) {}
	uint32_t bpos;
	uint32_t epos;
};

struct ListResult{
	ListResult(uint32_t val, uint32_t count) : val(val), count(count) {}
	uint32_t val;
	uint32_t count;
	int operator < (const ListResult& lr) const{
		if (count != lr.count) return count < lr.count;
		return val < lr.val;
	}
};

class WaveletMatrix{
public:
	WaveletMatrix();
	~WaveletMatrix();
	void Clear();
	uint32_t Length() const;
	uint32_t MaxVal() const;
	uint32_t Lookup(uint32_t pos) const;
	uint32_t Rank(uint32_t val, uint32_t pos) const;
	uint32_t Select(uint32_t val, uint32_t ind) const;
	uint32_t Quantile(uint32_t k, const Range& range) const;
	std::vector<ListResult> ListMode(uint32_t minval, uint32_t maxval, const Range& range, uint32_t num) const;
	std::pair<uint32_t, uint32_t> RankRange(uint32_t val, const Range& range) const;
	std::pair<uint32_t, uint32_t> LookupAndRank(uint32_t pos) const;
	uint32_t GetWorkingSize() const;

	void Save(std::ostream& os) const;
	void Load(std::istream& is);
private:
	uint32_t SelectHelper(uint32_t depth, uint32_t val, uint32_t ind, uint32_t pos) const;
	friend class WaveletMatrixBuilder;
	std::vector<rsdic::RSDic> layers_;
	uint32_t max_val_;
};

}

#endif // WAVELET_MATRIX_WAVELET_MATRIX_HPP_