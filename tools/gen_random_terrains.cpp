// Generate a fixed batch of random Perlin terrains for Phase 4 property
// tests. Mirrors the random_terrain / random_integer_terrain helpers in
// `Barnes2020-FillSpillMerge/unittests/fsm_tests.cpp`:
//
//   random_terrain(min_size, max_size):
//     perlin(size, seed) where size, seed are PRNG draws
//
//   random_integer_terrain(min_size, max_size):
//     same as random_terrain, then multiply by 100 and truncate to int
//
// We commit the resulting TIFFs to test/test_cases/random/ so the Julia
// tests can iterate over a deterministic, reproducible set of terrains
// without needing to port `perlin` to Julia.
//
// Edge convention mirrors the C++ tests: `dem.setEdges(-1)` so the
// outermost ring is below any interior cell, marking those as ocean
// once the Julia loader does `prepare_label_and_flowdirs(topo, 0.0)`.
//
// Usage:
//   ./gen_random_terrains.exe <out_dir>
//
// Output layout under <out_dir>:
//   small_float_NN.tif    (size 10..30, ~30 files)
//   small_int_NN.tif      (size 10..30, integer-truncated)
//   large_float_NN.tif    (size 100..200, ~5 files)
//   large_int_NN.tif      (size 100..200)
//
// Each TIFF stores Float64; ocean ring cells = -1.0 numerically (NoData
// is set to a value not present in the data; the Julia loader sets
// label = OCEAN by elevation rule, see prepare_label_and_flowdirs).

#include <richdem/terrain_generation.hpp>
#include <richdem/common/Array2D.hpp>

#include <cstdint>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <random>
#include <sstream>
#include <string>

namespace rd = richdem;

// Deterministic seed so re-running this tool produces the same TIFFs.
// Bump if/when we want a fresh batch.
constexpr uint64_t MASTER_SEED = 0xC0FFEEC0FFEE2026ULL;

static const int N_SMALL = 30;
static const int N_LARGE = 5;
static const int SMALL_MIN = 10;
static const int SMALL_MAX = 30;
static const int LARGE_MIN = 100;
static const int LARGE_MAX = 200;

// Truncate to int (mirrors random_integer_terrain).
static void int_truncate(rd::Array2D<double> &dem){
  for(auto i=dem.i0(); i<dem.size(); i++){
    dem(i) *= 100.0;
    dem(i) = static_cast<int>(dem(i));
  }
}

// Mirror C++ `dem.setEdges(-1)` — overwrite the outer ring with -1 so
// the loader can identify the ocean.
static void set_edges_minus_one(rd::Array2D<double> &dem){
  const int W = dem.width();
  const int H = dem.height();
  for(int x=0; x<W; x++){
    dem(x, 0)   = -1.0;
    dem(x, H-1) = -1.0;
  }
  for(int y=0; y<H; y++){
    dem(0,   y) = -1.0;
    dem(W-1, y) = -1.0;
  }
}

int main(int argc, char **argv){
  if(argc != 2){
    std::cerr << "Usage: " << argv[0] << " <out_dir>\n";
    return 1;
  }

  const std::string out_dir = argv[1];
  std::filesystem::create_directories(out_dir);

  std::mt19937_64 gen(MASTER_SEED);
  std::uniform_int_distribution<uint32_t> seed_dist;
  std::uniform_int_distribution<int>      small_size_dist(SMALL_MIN, SMALL_MAX);
  std::uniform_int_distribution<int>      large_size_dist(LARGE_MIN, LARGE_MAX);

  auto emit = [&](const std::string &kind, int idx, rd::Array2D<double> dem){
    set_edges_minus_one(dem);
    // Set NoData to a value not present in the data. The Julia loader
    // converts NoData to NaN, but we want the actual interior elevations
    // and the -1 ocean ring to both survive the load. Use a sentinel
    // below -1.
    dem.setNoData(-9999.0);
    std::ostringstream path;
    path << out_dir << "/" << kind << "_" << std::setw(2) << std::setfill('0') << idx << ".tif";
    dem.saveGDAL(path.str());
    std::cout << "  wrote " << path.str()
              << " (" << dem.width() << "x" << dem.height() << ")\n";
  };

  std::cout << "Generating " << N_SMALL << " small float terrains (size " << SMALL_MIN << ".." << SMALL_MAX << ")\n";
  for(int i=0; i<N_SMALL; i++){
    auto dem = rd::perlin(small_size_dist(gen), seed_dist(gen));
    emit("small_float", i, dem);
  }

  std::cout << "Generating " << N_SMALL << " small integer terrains\n";
  for(int i=0; i<N_SMALL; i++){
    auto dem = rd::perlin(small_size_dist(gen), seed_dist(gen));
    int_truncate(dem);
    emit("small_int", i, dem);
  }

  std::cout << "Generating " << N_LARGE << " large float terrains (size " << LARGE_MIN << ".." << LARGE_MAX << ")\n";
  for(int i=0; i<N_LARGE; i++){
    auto dem = rd::perlin(large_size_dist(gen), seed_dist(gen));
    emit("large_float", i, dem);
  }

  std::cout << "Generating " << N_LARGE << " large integer terrains\n";
  for(int i=0; i<N_LARGE; i++){
    auto dem = rd::perlin(large_size_dist(gen), seed_dist(gen));
    int_truncate(dem);
    emit("large_int", i, dem);
  }

  std::cout << "Done.\n";
  return 0;
}
