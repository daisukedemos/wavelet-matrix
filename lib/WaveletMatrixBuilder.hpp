#ifndef WAVLEET_MATRIX_WAVELET_MATRIX_BUILDER_HPP_
#define WAVELET_MATRIX_WAVELET_MATRIX_BUILDER_HPP_

#include <vector>
#include <stdint.h>

namespace wavelet_matrix{

class WaveletMatrix;

class WaveletMatrixBuilder{
public:
	WaveletMatrixBuilder();
	void Add(uint32_t val);
	void Build(WaveletMatrix& wm);
	void BuildWithAVals(WaveletMatrix& wm, std::vector<uint32_t>& avals);
	void Clear();

private:
	std::vector<uint32_t> vals_;
};

}

#endif // #ifndef WAVELET_MATRIX_WAVELET_MATRIX_BUILDER_HPP_