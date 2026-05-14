// Oracle dumper for the Julia port: runs the same ocean-labeling +
// GetDepressionHierarchy steps as fsm.exe (i.e. main.cpp lines 42-80),
// then writes the label grid, flowdirs grid, and full Depression vector
// to a single text file the Julia tests can parse.
//
// Lives here (FSM_Julia/tools/) because it only exists to feed the
// Julia tests. It builds against the vendored C++ snapshot at
// ../vendor/Barnes2020-FillSpillMerge/ — see build.sh next to this file
// for the compile command, and the vendor README for scope.
//
// The vendored dephier.hpp has the deterministic outlet sort comparator
// (tie-break on depa, depb, out_cell) pre-applied, so the dumped
// depression vector matches Julia's get_depression_hierarchy
// bit-for-bit out of the box. See
// FSM_Julia/tools/patches/dephier-deterministic-outlets.patch for the
// underlying diff vs. upstream.
//
// Usage:
//   dh_dump.exe <input.tif> <ocean_level> <output.txt>
//
// Output format:
//
//   WIDTH <W>
//   HEIGHT <H>
//   NDEP <N>
//   LABEL
//   <H rows of W space-separated UInt32>
//   FLOWDIRS
//   <H rows of W space-separated int8>
//   DEPRESSIONS
//   # idx pit_cell out_cell parent odep geolink pit_elev out_elev lchild rchild ocean_parent dep_label cell_count total_elevation dep_vol water_vol ocean_linked
//   <N rows>
//
// `ocean_linked` is `-` if empty, otherwise comma-separated UInt32.
// All flat indices are 0-based (matching C++ in-memory state). The Julia
// reader converts to 1-based on load.

#include <dephier/dephier.hpp>

#include <richdem/common/Array2D.hpp>
#include <richdem/common/constants.hpp>
#include <richdem/misc/misc_methods.hpp>

#include <cstdint>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>

namespace rd = richdem;
namespace dh = richdem::dephier;

int main(int argc, char **argv){
  if(argc != 4){
    std::cerr << "Usage: " << argv[0] << " <input.tif> <ocean_level> <output.txt>\n";
    return 1;
  }

  const std::string topo_path  = argv[1];
  const double      ocean_lvl  = std::stod(argv[2]);
  const std::string out_path   = argv[3];

  rd::Array2D<double> topo(topo_path);
  rd::Array2D<dh::dh_label_t> label   (topo.width(), topo.height(), dh::NO_DEP );
  rd::Array2D<rd::flowdir_t>  flowdirs(topo.width(), topo.height(), rd::NO_FLOW);

  // Mirror main.cpp: bucket-fill ocean from edges, then NoData -> OCEAN.
  rd::BucketFillFromEdges<rd::Topology::D8>(topo, label, ocean_lvl, dh::OCEAN);
  for(unsigned int i=0; i<label.size(); i++){
    if(topo.isNoData(i) || label(i)==dh::OCEAN){
      label(i) = dh::OCEAN;
    }
  }

  auto deps = dh::GetDepressionHierarchy<double, rd::Topology::D8>(topo, label, flowdirs);

  std::ofstream out(out_path);
  if(!out){
    std::cerr << "Could not open " << out_path << " for writing\n";
    return 1;
  }
  // Full Float64 round-trip precision.
  out << std::setprecision(17);

  const int W = topo.width();
  const int H = topo.height();

  out << "WIDTH "  << W << "\n";
  out << "HEIGHT " << H << "\n";
  out << "NDEP "   << deps.size() << "\n";

  out << "LABEL\n";
  for(int y=0; y<H; y++){
    for(int x=0; x<W; x++){
      if(x) out << ' ';
      out << label(x, y);
    }
    out << '\n';
  }

  out << "FLOWDIRS\n";
  for(int y=0; y<H; y++){
    for(int x=0; x<W; x++){
      if(x) out << ' ';
      out << static_cast<int>(flowdirs(x, y));
    }
    out << '\n';
  }

  out << "DEPRESSIONS\n";
  out << "# idx pit_cell out_cell parent odep geolink pit_elev out_elev "
         "lchild rchild ocean_parent dep_label cell_count total_elevation "
         "dep_vol water_vol ocean_linked\n";
  for(std::size_t i=0; i<deps.size(); i++){
    const auto &d = deps[i];
    out << i
        << ' ' << d.pit_cell
        << ' ' << d.out_cell
        << ' ' << d.parent
        << ' ' << d.odep
        << ' ' << d.geolink
        << ' ' << d.pit_elev
        << ' ' << d.out_elev
        << ' ' << d.lchild
        << ' ' << d.rchild
        << ' ' << (d.ocean_parent ? 1 : 0)
        << ' ' << d.dep_label
        << ' ' << d.cell_count
        << ' ' << d.total_elevation
        << ' ' << d.dep_vol
        << ' ' << d.water_vol
        << ' ';
    if(d.ocean_linked.empty()){
      out << '-';
    } else {
      for(std::size_t k=0; k<d.ocean_linked.size(); k++){
        if(k) out << ',';
        out << d.ocean_linked[k];
      }
    }
    out << '\n';
  }

  return 0;
}
