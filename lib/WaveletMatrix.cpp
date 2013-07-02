#include <queue>
#include "WaveletMatrix.hpp"

using std::pair;
using std::ostream;
using std::istream;
using std::vector;
using std::priority_queue;
using rsdic::RSDic;

#include <iostream>
using namespace std;

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
		val <<= 1;
		if (!rsdic.GetBit(pos)){
			pos = rsdic.Rank(pos, 0);
		} else {
			pos = rsdic.zero_num() + rsdic.Rank(pos, 1);
			val |= 1;
		}
	}
	return val;
}

uint32_t WaveletMatrix::Rank(uint32_t val, uint32_t pos) const{
	std::pair<uint32_t, uint32_t> range = RankRange(val, Range(0, pos));
	return range.second - range.first;
}

std::pair<uint32_t, uint32_t> WaveletMatrix::RankRange(uint32_t val, const Range& range) const{
	uint32_t bpos = range.bpos;
	uint32_t epos = range.epos;
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
		val <<= 1;
		if (bit){
			bpos += rsdic.zero_num();
			epos += rsdic.zero_num();
			val |= 1;
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

uint32_t WaveletMatrix::Quantile(uint32_t k, const Range& range) const{
	uint32_t bpos = range.bpos;
	uint32_t epos = range.epos;
	uint32_t val = 0;
	for (uint32_t depth = 0; depth < layers_.size() && bpos < epos; ++depth){
		const rsdic::RSDic& rsdic = layers_[depth];
		uint32_t nz_bpos = rsdic.Rank(bpos, 0);
		uint32_t nz_epos = rsdic.Rank(epos, 0);
		val <<= 1;
		if (k < nz_epos - nz_bpos){
			bpos = nz_bpos;
			epos = nz_epos;
		} else {
			k -= (nz_epos - nz_bpos);
			bpos = rsdic.zero_num() + bpos - nz_bpos;
			epos = rsdic.zero_num() + epos - nz_epos;
			val |= 1;
		}
	}
	return val;
}

struct WMRange{
	WMRange(uint32_t bpos, uint32_t epos, uint32_t depth, uint32_t prefix_char) :
		bpos(bpos), epos(epos), depth(depth), prefix_char(prefix_char) {}
	uint32_t bpos;
	uint32_t epos;
	uint32_t depth;
	uint32_t prefix_char;
};

class ModeCompararator{
public:
	bool operator() (const WMRange& lhs, const WMRange& rhs) const{
		if (lhs.epos - lhs.bpos != rhs.epos - rhs.bpos){
			return lhs.epos - lhs.bpos < rhs.epos - rhs.bpos;
		} else if (lhs.depth != rhs.depth){
			return lhs.depth < rhs.depth;
		} else if (lhs.prefix_char != rhs.prefix_char){
			return lhs.prefix_char < rhs.prefix_char;
		}
		return lhs.bpos < rhs.bpos;
	}
};

uint32_t PrefixCode(uint32_t x, uint32_t len, uint32_t total_len){
	return x >> (total_len - len);
}

bool CheckPrefix(uint32_t prefix, uint32_t depth, uint32_t total_depth, uint32_t minval, uint32_t maxval){
	if (PrefixCode(minval, depth, total_depth) <= prefix &&
		PrefixCode(maxval-1, depth, total_depth) >= prefix) return true;
	else return false;
}

vector<ListResult> WaveletMatrix::ListMode(uint32_t minval, uint32_t maxval, const Range& range, uint32_t num) const{
	priority_queue<WMRange, vector<WMRange>, ModeCompararator> qs;
	qs.push(WMRange(range.bpos, range.epos, 0, 0));
	vector<ListResult> res;
	while (res.size() < num && !qs.empty()){
		WMRange range = qs.top();
		qs.pop();
		cout << range.epos - range.bpos << " " << range.depth << " " << range.prefix_char << endl;
		if (range.depth == layers_.size()){
			cout << "found " << range.prefix_char << " " << range.epos - range.bpos << endl;
			res.push_back(ListResult(range.prefix_char, range.epos - range.bpos));
			continue;
		}
		const rsdic::RSDic& rsdic = layers_[range.depth];
		uint32_t bpos = range.bpos;
		uint32_t epos = range.epos;
		uint32_t nz_bpos = rsdic.Rank(bpos, 0);
		uint32_t nz_epos = rsdic.Rank(epos, 0);
		uint32_t no_bpos = bpos - nz_bpos + rsdic.zero_num();
		uint32_t no_epos = epos - nz_epos + rsdic.zero_num();
		cout << "count = " << epos - bpos << "\t" << nz_epos - nz_bpos << "\t" << no_epos - no_bpos << endl;
		if (nz_epos - nz_bpos > 0){
			uint32_t next_prefix = range.prefix_char << 1;
			if (CheckPrefix(next_prefix, range.depth+1, layers_.size(), minval, maxval)){
				qs.push(WMRange(nz_bpos, nz_epos, range.depth+1, next_prefix));
			}
		}
		if (no_epos - no_bpos > 0){
			uint32_t next_prefix = (range.prefix_char << 1)+ 1;
			if (CheckPrefix(next_prefix, range.depth+1, layers_.size(), minval, maxval)){
				qs.push(WMRange(no_bpos, no_epos, range.depth+1, next_prefix));
			}
		}
	}
	sort(res.rbegin(), res.rend());
	return res;
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