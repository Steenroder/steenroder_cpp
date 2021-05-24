/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#include "input.in"
#include "steenroder/args_parser.hpp"
#include "steenroder/commons.hpp"
//#include "steenroder/reader.hpp"

#include <steenroder/sparse_matrix.hpp>
#include <steenroder/vector_column.hpp>
#include <steenroder/boundary_matrix.hpp>

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
void write(T& data, const std::string& name,
           const std::string& output_filename, bool use_binary) {
  if(use_binary) {
    data.save_binary(name, output_filename);
  } else {
    data.save_ascii(name, output_filename);
  }
}

void dualize(const std::string& input_filename,
                               const std::string& output_filename,
                               const bool use_binary) {
  BoundaryMatrix<VectorColumn> boundary_matrix;
  read(boundary_matrix, input_filename, use_binary);
  boundary_matrix.dualize();
  write(boundary_matrix, "dualized", output_filename, use_binary);


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
  dualize(args.input_filename, args.output_filename, use_binary);

  return 0;
}
