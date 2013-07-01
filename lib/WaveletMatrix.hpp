#ifndef WAVELET_MATRIX_WAVELET_MATRIX_HPP_
#define WAVELET_MATRIX_WAVELET_MATRIX_HPP_

#include <vector>
#include <stdint.h>
#include <rsdic/RSDic.hpp>

namespace wavelet_matrix{

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
	std::pair<uint32_t, uint32_t> RankRange(uint32_t val, uint32_t bpos, uint32_t epos) const;
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