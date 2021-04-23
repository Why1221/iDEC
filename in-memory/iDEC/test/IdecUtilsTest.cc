#include <gtest/gtest.h>
#include "IdecUtils.h"

#include <bitset>
#include <chrono>
#include <limits>
#include <random>

using ::testing::InitGoogleTest;
using ::testing::Test;

// test counting ones
TEST(IdecUtilsTest, CountingOnes) {
  std::mt19937_64 gen(
      std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_int_distribution<uint64_t> dist(
      0, std::numeric_limits<uint64_t>::max());
  unsigned T = 1e3;
  for (unsigned i = 0; i < T; ++i) {
    uint64_t x = dist(gen);
    std::bitset<64> bits(x);
    EXPECT_EQ(countingOnes(x), (int16_t)bits.count());
  }
}

TEST(IdecUtilsTest, hammingDistance) {
  std::mt19937_64 gen(
      std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_int_distribution<uint64_t> dist(
      0, std::numeric_limits<uint64_t>::max());
  unsigned T = 1e3;
  for (unsigned i = 0; i < T; ++i) {
    uint64_t x = dist(gen);
    uint64_t y = dist(gen);
    std::bitset<64> bits_x(x), bits_y(y);
    auto bits = bits_x ^ bits_y;
    EXPECT_EQ(hammingDistance(x, y), (int16_t)bits.count());
  }
}

int16_t inner_prod(const std::bitset<64>& a, const std::bitset<64>& b) {
  int16_t res = 0;
  for (auto i = (size_t)0; i < a.size(); ++i) {
    res += (a[i] * (b[i] == 0 ? -1 : 1));
  }
  return res;
}

TEST(IdecUtilsTest, calcToWSum) {
  std::mt19937_64 gen(
      std::chrono::system_clock::now().time_since_epoch().count());
  std::uniform_int_distribution<uint64_t> dist(
      0, std::numeric_limits<uint64_t>::max());
  unsigned T = 1e3;
  for (unsigned i = 0; i < T; ++i) {
    uint64_t x = dist(gen);
    uint64_t y = dist(gen);
    std::bitset<64> bits_x(x), bits_y(y);
    EXPECT_EQ(calcToWSum(x, y), inner_prod(bits_x, bits_y));
  }
}

int main(int argc, char** argv) {
  InitGoogleTest(&argc, argv);
  return RUN_ALL_TESTS();
}