/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#include "input.in"
#include "steenroder/args_parser.hpp"
#include "steenroder/commons.hpp"

#include <steenroder/sparse_matrix.hpp>
#include <steenroder/vector_column.hpp>
#include <steenroder/boundary_matrix.hpp>
#include <steenroder/simplex_matrix.hpp>
#include <steenroder/bars.hpp>
#include <steenroder/reduction.hpp>
#include <steenroder/homology.hpp>
#include <steenroder/steenrod.hpp>
#include <steenroder/sorted_matrix.hpp>
#include <steenroder/sorted_bars.hpp>

using namespace stn;

template<class T>
void read(T& data, const std::string& input_filename, bool use_binary) {
  bool read_successful;

  if(use_binary) {
    read_successful = data.load_binary(input_filename);
  } else {
    read_successful = data.load_ascii(input_filename);
  }

  if(!read_successful) {
    std::cerr << "Error opening file " << input_filename << std::endl;
  }
}

template<class T>
void read_dual(T& data, const std::string& input_filename, bool use_binary) {
  bool read_successful;

  if(use_binary) {
    read_successful = data.load_binary_dual(input_filename);
  } else {
    read_successful = data.load_ascii_dual(input_filename);
  }

  if(!read_successful) {
    std::cerr << "Error opening file " << input_filename << std::endl;
  }
}

template<class T>
void write(T& data, const std::string& name,
           const std::string& output_filename, bool use_binary) {
  if(use_binary) {
    data.save_binary(name, output_filename);
  } else {
    data.save_ascii(name, output_filename);
  }
}

template<typename ColumnType = VectorColumn>
void write_pairs(const ViewFiniteBars<ColumnType>& finite_bars,
                 const ViewInfiniteBars<ColumnType>& infinite_bars,
                 const std::string& output_filename, bool use_binary,
                 const std::string& prefix) {
  std::string filename = output_filename + "_" + prefix +"_pairs.dat";

  if(use_binary) {
    save_pairs_binary(filename, finite_bars, infinite_bars);
  } else {
    save_pairs_ascii(filename, finite_bars, infinite_bars);
  }
}


template<typename ColumnType = VectorColumn>
void write_pairs(const FiniteBars<ColumnType>& finite_bars,
                 const  InfiniteBars<ColumnType>& infinite_bars,
                 const std::string& output_filename, bool use_binary,
                 const std::string& prefix) {
  std::string filename = output_filename + "_" + prefix +"_pairs.dat";

  if(use_binary) {
    save_pairs_binary(filename, finite_bars, infinite_bars);
  } else {
    save_pairs_ascii(filename, finite_bars, infinite_bars);
  }
}

template<typename ColumnType = VectorColumn>
void write_pairs(const Bars<ColumnType>& bars,
                 const std::string& output_filename, bool use_binary,
                 const std::string& prefix) {
  std::string filename = output_filename + "_" + prefix +"_pairs.dat";

  if(use_binary) {
    // save_pairs_binary(filename, bars);
  } else {
    save_pairs_ascii(filename, bars);
  }
}


void compute_steenrod_barcodes(const std::string& input_filename,
                               const std::string& output_filename,
                               const bool use_binary) {
  const dimension_t d = 1;
  const dimension_t k = 1;

  ViewMatrix<VectorColumn> boundary_matrix;
  read(boundary_matrix, input_filename, use_binary);
  write(boundary_matrix, "boundary", output_filename, use_binary);

  SimplexMatrix<VectorColumn> simplex_matrix(boundary_matrix, d, d + k);
  write(simplex_matrix, "simplex", output_filename, use_binary);

  // // Need to delete boundary_matrix to release memory

  // Relative cohomology
  ViewMatrix<VectorColumn> dual_boundary_matrix;
  read_dual(dual_boundary_matrix, input_filename, use_binary);
  write(dual_boundary_matrix, "dual_boundary", output_filename, use_binary);
  index_t n_dimensions = dual_boundary_matrix.get_n_dimensions();
  index_t n_cells = dual_boundary_matrix.get_n_columns();

  ViewInfiniteBars<VectorColumn> dual_infinite_bars_matrix(n_cells, n_dimensions);
  ViewFiniteBars<VectorColumn> dual_finite_bars_matrix(dual_boundary_matrix);

  Homology<TwistReduction<VectorColumn>> dual_homology;
  dual_homology.compute(dual_finite_bars_matrix, dual_infinite_bars_matrix);

  write(dual_finite_bars_matrix, "dual_finite", output_filename, use_binary);
  write(dual_infinite_bars_matrix, "dual_infinite", output_filename, use_binary);

  //dual_finite_bars_matrix.dualize();
  //dual_infinite_bars_matrix.dualize();
  write_pairs(dual_finite_bars_matrix, dual_infinite_bars_matrix,
              output_filename, use_binary, "dual");

  index_t n_finite_bars = dual_finite_bars_matrix.get_n_bars();
  index_t n_infinite_bars = dual_infinite_bars_matrix.get_n_bars();
  Bars<VectorColumn> steenrod_bars_matrix(n_cells);

  Steenrod<StandardReduction<VectorColumn>> steenrod(d, k, n_cells, simplex_matrix);
  steenrod.compute(dual_finite_bars_matrix, dual_infinite_bars_matrix,
                   steenrod_bars_matrix);

  write(steenrod_bars_matrix, "steenrod", output_filename, use_binary);

  //  steenrod_bars_matrix.dualize();
  steenrod_bars_matrix.dualize();
  write_pairs(steenrod_bars_matrix, output_filename, use_binary, "steenrod");

}


int main(int argc, char* argv[]) {
  using namespace stn;
  STN_INSTRUMENT_ON("main", 0)

    ArgsParser args(argc, argv);
  if(args.help){
    args.printUsage();
    return 0;
  }

  bool use_binary = false;
  compute_steenrod_barcodes(args.input_filename,
                            args.output_filename,
                            use_binary);

  return 0;
}
