#include <algorithm>
#include <rsdic/RSDicBuilder.hpp>
#include "WaveletMatrix.hpp"
#include "WaveletMatrixBuilder.hpp"

using std::vector;
using rsdic::RSDicBuilder;
using std::copy;

namespace {
uint32_t GetMax(vector<uint32_t>& v){
	uint32_t maxval = 0;
	for (size_t i = 0; i < v.size(); ++i){
		if (v[i] > maxval){
			maxval = v[i];
		}
	}
	return maxval;
}


uint32_t GetLen(uint32_t val){
	uint32_t len = 0;
	while (val > 0){
		val >>= 1;
		++len;
	}
	return len;
}

template <class Vals, class Depth, class NextZeros, class NextOnes, class RSDB>
void Filter(
	const Vals& vals,
	Depth depth,
	NextZeros& next_zeros,
	NextOnes& next_ones,
	RSDB& rsdb){
	for (size_t i = 0; i < vals.size(); ++i){
		uint32_t bit = (vals[i] >> depth) & 1LLU;
		rsdb.PushBack(bit);
		if (!bit){
			next_zeros.push_back(vals[i]);
		} else {
			next_ones.push_back(vals[i]);
		}
	}
}

template <class Vals, class AVals, class Depth, class NextZeros, class NextOnes, class NextZeroAVals, class NextOneAVals, class RSDB>
void FilterWithAVals(
	const Vals& vals,
	const AVals& avals,
	Depth depth,
	NextZeros& next_zeros,
	NextOnes& next_ones,
	NextZeroAVals& next_zero_avals,
	NextOneAVals& next_one_avals,
	RSDB& rsdb){
	for (size_t i = 0; i < vals.size(); ++i){
		uint32_t bit = (vals[i] >> depth) & 1LLU;
		rsdb.PushBack(bit);
		if (!bit){
			next_zeros.push_back(vals[i]);
			next_zero_avals.push_back(avals[i]);
		} else {
			next_ones.push_back(vals[i]);
			next_one_avals.push_back(avals[i]);
		}
	}
}
}

namespace wavelet_matrix{

WaveletMatrixBuilder::WaveletMatrixBuilder(){
}

void WaveletMatrixBuilder::Add(uint32_t val){
	vals_.push_back(val);
}

void WaveletMatrixBuilder::Clear(){
	vals_.clear();
}

void WaveletMatrixBuilder::Build(WaveletMatrix& wm){
	uint32_t max_val = GetMax(vals_);
	uint32_t max_depth = GetLen(max_val);
	wm.max_val_ = max_val;
	wm.layers_.resize(max_depth);

	vector<uint32_t> zeros;
	vector<uint32_t> ones;
	zeros.swap(vals_);
	for (uint32_t depth = 0; depth < max_depth; ++depth){
		vector<uint32_t> next_zeros;
		vector<uint32_t> next_ones;

		RSDicBuilder rsdb;
		Filter(zeros, max_depth - depth - 1, next_zeros, next_ones, rsdb);
		Filter(ones,  max_depth - depth - 1, next_zeros, next_ones, rsdb);
		zeros.swap(next_zeros);
		ones.swap(next_ones);
		rsdb.Build(wm.layers_[depth]);
	}
}

void WaveletMatrixBuilder::BuildWithAVals(WaveletMatrix& wm, vector<uint32_t>& avals){
	uint32_t max_val = GetMax(vals_);
	uint32_t max_depth = GetLen(max_val);
	wm.max_val_ = max_val;
	wm.layers_.resize(max_depth);

	vector<uint32_t> zeros;
	vector<uint32_t> ones;

	vector<uint32_t> zero_avals;
	vector<uint32_t> one_avals;
	zeros.swap(vals_);
	zero_avals = avals;
	for (uint32_t depth = 0; depth < max_depth; ++depth){
		vector<uint32_t> next_zeros;
		vector<uint32_t> next_ones;
		vector<uint32_t> next_zero_avals;
		vector<uint32_t> next_one_avals;

		RSDicBuilder rsdb;
		FilterWithAVals(zeros, zero_avals, max_depth - depth - 1, next_zeros, next_ones, next_zero_avals, next_one_avals, rsdb);
		FilterWithAVals(ones,  one_avals,  max_depth - depth - 1, next_zeros, next_ones, next_zero_avals, next_one_avals, rsdb);
		zeros.swap(next_zeros);
		ones.swap(next_ones);
		zero_avals.swap(next_zero_avals);
		one_avals.swap(next_one_avals);
		rsdb.Build(wm.layers_[depth]);
	}
	copy(zero_avals.begin(), zero_avals.end(), avals.begin());
	copy(one_avals.begin(), one_avals.end(), avals.begin() + zero_avals.size());
}

}