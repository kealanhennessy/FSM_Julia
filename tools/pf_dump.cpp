// Priority-Flood oracle dumper for the Julia port's heavy-flooding
// property test.
//
// Runs richdem's Zhou (2016) Priority-Flood depression-filling on a
// terrain GeoTIFF and writes the filled DEM back out as a Float64
// GeoTIFF. The Julia test compares its FSM heavy-flood surface against
// this frozen, trusted-C++ reference (rather than against a Julia
// re-implementation of Priority-Flood), so the cross-algorithm check
// is grounded in the upstream code, consistent with the dephier and
// water-table-depth oracles.
//
// Reference: Zhou, G., Sun, Z., Fu, S., 2016. An efficient variant of
// the Priority-Flood algorithm... Computers & Geosciences 90, 87-96.
//
// The input terrains carry their ocean ring as a numeric -1.0 (the
// declared NoData is an unused -9999), exactly as the Julia loader
// reads them (raw Float64, no NoData->NaN substitution). richdem's
// Zhou2016 does not special-case NoData, so it processes the same grid
// values the Julia side sees — the two stay comparable.
//
// Builds against your own clone of the upstream C++
// Barnes2020-FillSpillMerge (no copy is bundled in this repo); see
// tools/build.sh and the repository root README.md ("Testing" ->
// Option B) for the clone + checkout + patch recipe.
//
// Usage:
//   pf_dump.exe <input.tif> <output.tif>

#include <richdem/depressions/Zhou2016.hpp>
#include <richdem/common/Array2D.hpp>

#include <iostream>
#include <string>

namespace rd = richdem;

int main(int argc, char **argv){
  if(argc != 3){
    std::cerr << "Usage: " << argv[0] << " <input.tif> <output.tif>\n";
    return 1;
  }

  const std::string in_path  = argv[1];
  const std::string out_path = argv[2];

  rd::Array2D<double> dem(in_path);

  // Fill all depressions / remove digital dams in place. This is the
  // same routine the upstream FSM unit test compares against.
  rd::PriorityFlood_Zhou2016<double>(dem);

  // Save as a Float64 GeoTIFF. The Julia side reads this back with
  // ArchGDAL as a Float64 matrix; the round-trip is bit-exact.
  dem.saveGDAL(out_path);

  return 0;
}
