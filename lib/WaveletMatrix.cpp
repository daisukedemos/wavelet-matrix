#include "WaveletMatrix.hpp"

using std::pair;
using std::ostream;
using std::istream;

namespace wavelet_matrix{

WaveletMatrix::WaveletMatrix(): max_val_(0) {
}

WaveletMatrix::~WaveletMatrix(){
}

void WaveletMatrix::Clear() {
	layers_.clear();
	max_val_ = 0;
}

uint32_t WaveletMatrix::Length() const{
	if (layers_.size() == 0) return 0;
	return layers_[0].num();
}

uint32_t WaveletMatrix::MaxVal() const{
	return max_val_;
}

uint32_t WaveletMatrix::Lookup(uint32_t pos) const{
	uint32_t val = 0;
	for (uint32_t depth = 0; depth < layers_.size(); ++depth){
		const rsdic::RSDic& rsdic = layers_[depth];
		if (!rsdic.GetBit(pos)){
			pos = rsdic.Rank(pos, 0);
		} else {
			pos = rsdic.zero_num() + rsdic.Rank(pos, 1);
			val |= (1LLU << (layers_.size() - depth - 1));
		}
	}
	return val;
}

uint32_t WaveletMatrix::Rank(uint32_t val, uint32_t pos) const{
	std::pair<uint32_t, uint32_t> range = RankRange(val, 0, pos);
	return range.second - range.first;
}

std::pair<uint32_t, uint32_t> WaveletMatrix::RankRange(uint32_t val, uint32_t bpos, uint32_t epos) const{
	for (uint32_t depth = 0; depth < layers_.size() && bpos < epos; ++depth){
		uint32_t bit = (val >> (layers_.size() - depth - 1)) & 1LLU;
		const rsdic::RSDic& rsdic = layers_[depth];
		bpos = rsdic.Rank(bpos, bit);
		epos = rsdic.Rank(epos, bit);
		if (bit){
			bpos += rsdic.zero_num();
			epos += rsdic.zero_num();
		}
	}
	return std::make_pair(bpos, epos);
}

std::pair<uint32_t, uint32_t> WaveletMatrix::LookupAndRank(uint32_t pos) const{
	uint32_t val = 0;
	uint32_t bpos = 0;
	uint32_t epos = pos;
	for (uint32_t depth = 0; depth < layers_.size(); ++depth){
		const rsdic::RSDic& rsdic = layers_[depth];
		uint32_t bit = rsdic.GetBit(epos);
		bpos = rsdic.Rank(bpos, bit);
		epos = rsdic.Rank(epos, bit);
		if (bit){
			bpos += rsdic.zero_num();
			epos += rsdic.zero_num();
			val |= (1LLU << (layers_.size() - depth - 1));
		}
	}
	return std::make_pair(val, epos - bpos);
}



uint32_t WaveletMatrix::SelectHelper(uint32_t depth, uint32_t val, uint32_t ind, uint32_t pos) const{
	if (depth == layers_.size()){
		return pos + ind;
	}
	uint32_t bit = (val >> (layers_.size() - depth - 1)) & 1LLU;
	const rsdic::RSDic& rsdic = layers_[depth];
	pos = rsdic.Rank(pos, bit);
	if (!bit){
		ind = SelectHelper(depth+1, val, ind, pos);
		return rsdic.Select(ind, 0);
	} else {
		uint32_t zn = rsdic.zero_num();
		pos += zn;
		ind = SelectHelper(depth+1, val, ind, pos);
		return rsdic.Select(ind - zn, 1);
	}
}

uint32_t WaveletMatrix::Select(uint32_t val, uint32_t ind) const{
	return SelectHelper(0, val, ind, 0);
}

uint32_t WaveletMatrix::GetWorkingSize() const{
	uint32_t ret = 0;
	for (size_t i = 0; i < layers_.size(); ++i){
		ret += layers_[i].GetUsageBytes();
	}
	return ret;
}

void WaveletMatrix::Save(ostream& os) const{
	uint32_t layer_num = layers_.size();
	os.write((const char*)&layer_num, sizeof(layer_num));
	for (uint32_t i = 0; i < layer_num; ++i){
		layers_[i].Save(os);
	}
	os.write(reinterpret_cast<const char*>(&max_val_), sizeof(max_val_));
}

void WaveletMatrix::Load(istream& is){
	uint32_t layer_num = 0;
	is.read(reinterpret_cast<char*>(&layer_num), sizeof(layer_num));
	layers_.resize(layer_num);
	for (uint32_t i = 0; i < layer_num; ++i){
		layers_[i].Load(is);
	}
	is.read(reinterpret_cast<char*>(&max_val_), sizeof(max_val_));

}

}