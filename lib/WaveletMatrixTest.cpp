#include "WaveletMatrix.hpp"
#include "WaveletMatrixBuilder.hpp"
#include <map>
#include <gtest/gtest.h>

using namespace std;
using namespace wavelet_matrix;

TEST(WaveletMatrix, trivial){
	WaveletMatrixBuilder wmb;
	WaveletMatrix wm;
	wmb.Build(wm);
	wmb.Clear();
	wmb.Build(wm);

	EXPECT_EQ(0, wm.Length());
	EXPECT_EQ(0, wm.MaxVal());
	EXPECT_EQ(0, wm.Rank(0, 0));
	EXPECT_EQ(0, wm.Rank(1, 0));
}

vector<uint64_t> SetVector(const string& s){
	vector<uint64_t> ret;
	istringstream is(s);
	for (uint64_t v; is >> v; ){
		ret.push_back(v);
	}
	return ret;
}

TEST(WaveletMatrix, small){
	WaveletMatrixBuilder wmb;
	WaveletMatrix wm;
	wmb.Add(3);
	wmb.Add(2);
	wmb.Add(1);
	wmb.Add(2);
	wmb.Build(wm);

	EXPECT_EQ(4, wm.Length());
	EXPECT_EQ(3, wm.MaxVal());
	EXPECT_EQ(3, wm.Lookup(0));
	EXPECT_EQ(2, wm.Lookup(1));
	EXPECT_EQ(1, wm.Lookup(2));
	EXPECT_EQ(2, wm.Lookup(3));
	EXPECT_EQ(0, wm.Rank(0, 0));
	EXPECT_EQ(0, wm.Rank(1, 0));
	EXPECT_EQ(0, wm.Rank(2, 0));
	EXPECT_EQ(0, wm.Rank(3, 0));
	EXPECT_EQ(0, wm.Rank(0, 1));
	EXPECT_EQ(0, wm.Rank(1, 1));
	EXPECT_EQ(0, wm.Rank(2, 1));
	EXPECT_EQ(1, wm.Rank(3, 1));

	EXPECT_EQ(2, wm.Select(1, 0));
	EXPECT_EQ(1, wm.Select(2, 0));
	EXPECT_EQ(3, wm.Select(2, 1));
	EXPECT_EQ(0, wm.Select(3, 0));
}

TEST(WaveletMatrix, random){
	WaveletMatrixBuilder wmb;
	WaveletMatrix wm;
	vector<uint32_t> orig;
	const int N = 1000;
	const int ALPHA = 100;
	for (uint64_t i = 0; i < N; ++i){
		uint32_t v = rand() % ALPHA;
		wmb.Add(v);
		orig.push_back(v);
	}
	wmb.Build(wm);

	EXPECT_EQ(N, wm.Length());
	vector<uint32_t> cums(ALPHA);
	for (uint64_t i = 0; i < N; ++i){
		uint64_t val = orig[i];
		EXPECT_EQ(val, wm.Lookup(i));
		for (uint j = 0; j < ALPHA; ++j){
			EXPECT_EQ(cums[j], wm.Rank(j, i));
		}
		EXPECT_EQ(i, wm.Select(val, cums[val]));
		++cums[val];
	}

	ostringstream os;
	wm.Save(os);
	istringstream is(os.str());
	WaveletMatrix wm_load;
	wm_load.Load(is);
	fill(cums.begin(), cums.end(), 0);
	EXPECT_EQ(N, wm_load.Length());
	for (uint64_t i = 0; i < N; ++i){
		uint64_t val = orig[i];
		EXPECT_EQ(val, wm_load.Lookup(i));
		for (uint j = 0; j < ALPHA; ++j){
			EXPECT_EQ(cums[j], wm_load.Rank(j, i));
		}
		EXPECT_EQ(i, wm_load.Select(val, cums[val]));
		++cums[val];
	}

	const uint32_t query_N = 1000;
	for (uint32_t i = 0; i < query_N; ++i){
		uint32_t bpos = rand() % (N + 1);
		uint32_t epos = rand() % (N + 1);
		if (bpos == epos) continue;
		if (epos < bpos) swap(bpos, epos);
		uint32_t k = rand() % (epos - bpos);
		uint32_t val = wm.Quantile(k, Range(bpos, epos));
		vector<uint32_t> tmp(epos - bpos);
		copy(orig.begin() + bpos, orig.begin() + epos, tmp.begin());
		sort(tmp.begin(), tmp.end());
		EXPECT_EQ(tmp[k], val);

		uint32_t minval = rand() % ALPHA;
		uint32_t maxval = rand() % ALPHA;
		if (minval == maxval) continue;
		if (maxval < minval) swap(minval, maxval);
		uint32_t num = 1000;
		vector<ListResult> res = wm.ListMode(minval, maxval, Range(bpos, epos), num);

		map<uint32_t, uint32_t> counter;
		for (size_t i = 0; i < tmp.size(); ++i){
			++counter[tmp[i]];
		}
		vector<pair<uint32_t, uint32_t> > ord;
		for (map<uint32_t, uint32_t>::const_iterator it = counter.begin(); it != counter.end(); ++it){
			ord.push_back(make_pair(it->second, it->first));
		}
		sort(ord.rbegin(), ord.rend());
		vector<pair<uint32_t, uint32_t> > filtered_ord;
		for (size_t i = 0; i < ord.size(); ++i){
			if (minval <= ord[i].second  && ord[i].second < maxval){
				filtered_ord.push_back(ord[i]);
				if (filtered_ord.size() == num){
					break;
				}
			}
		}
		ASSERT_EQ(filtered_ord.size(), res.size());
		for (size_t i = 0; i < res.size(); ++i){
			ASSERT_EQ(filtered_ord[i].second, res[i].val);
			ASSERT_EQ(filtered_ord[i].first, res[i].count);
		}
	}
}