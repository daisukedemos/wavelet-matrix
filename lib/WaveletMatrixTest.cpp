#include "WaveletMatrix.hpp"
#include "WaveletMatrixBuilder.hpp"
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

TEST(WaveletMatrix, set){
	WaveletMatrixBuilder wmb;
	WaveletMatrix wm;
	vector<uint64_t> orig;
	const int N = 10;
	const int ALPHA = 4;
	for (uint64_t i = 0; i < N; ++i){
		uint32_t v = rand() % ALPHA;
		wmb.Add(v);
		orig.push_back(v);
	}
	wmb.Build(wm);

	EXPECT_EQ(N, wm.Length());
	vector<uint64_t> cums(ALPHA);
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


}