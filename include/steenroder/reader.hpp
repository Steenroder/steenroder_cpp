/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include <string>
#include <algorithm>

#include "Commons.hpp"
#include "Options.hpp"
#include "SparseMatrix.hpp"


namespace stn {

  template <class T, InputOutputFormat inputOutput>
  class Reader {};

  template <class T>
  class Reader<T, InputOutputFormat::Generic> {
  protected:
    const std::string filename;
    const InputOutputMode mode;

    std::ifstream ifs;

    Reader(const std::string& filename_in,
           const InputOutputMode mode_in)
      : filename(filename_in)
      , mode(mode_in)
    {}

  protected:
    inline bool openFile() {
      if (mode == InputOutputMode::ascii)
        ifs.open(filename, std::ifstream::in);
      else if (mode == InputOutputMode::binary)
        ifs.open(filename, std::ifstream::in | std::ifstream::binary);

      ifs.precision(16);

      if (ifs.fail()) {
        std::cout << "Could not open file " << filename << std::endl;
      }

      return ifs.fail();
    }

    inline bool closeFile() {
      ifs.close();

      if (ifs.fail()) {
        std::cout << "Could not close file " << filename << std::endl;
      }

      return ifs.fail();

    }

  };


  template <class T, InputOutputFormat format, InputOutputMode mode>
  class BoundaryMatrixReader {};

  template <class T>
  class BoundaryMatrixReader<T, InputOutputFormat::PHAT, InputOutputMode::ascii>
    : public Reader<T, InputOutputFormat::Generic> {
  private:
    using Base = Reader<T, InputOutputFormat::Generic>;

    inline bool getLine(std::string& line) {
      if(getline(Base::ifs, line))
        {
          line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
          if(line != "" && line[ 0 ] != '#')
            return true;
          else
            return getLine(line);
        }
      else
        return false;
    }

    inline int getNumberOfLines() {
      int n_lines = 0;

      if(!Base::openFile())
        return 0;

      std::string line;
      while(getLine(Base::ifs, line)) {
        n_lines++;
      }

      if(!Base::closeFile())
        return 0;

      return n_lines;
    }


    inline void readColumn(std::stringstream& line_ss, dimension_t& temp_dimension,
                           column& temp_column) {
        line_ss >> temp_dimension;

        index temp_index;
        temp_column.clear();
        while(line_ss.good()) {
          line_ss >> temp_index;
          temp_column.push_back(temp_index);
        }
        std::sort(temp_column.begin(), temp_column.end());
    }

  public:
    BoundaryMatrixReader(const std::string& filename_in)
      : Base(filename_in, InputOutputMode::ascii) {}

    template<MatrixFormat F>
      bool read(SparseMatrix<T, F>& boundaryMatrix) {
      STN_INSTRUMENT_ON("BoundaryMatrixReader<PHAT, ascii>::read", 3)

      index n_columns = getNumberOfLines();
      if(!n_columns)
        return false;

      boundaryMatrix.set_n_columns(n_columns);

      if(!Base::openFile())
        return false;

      dimension_t temp_dimension;
      column temp_column;
      index idx_column = -1;

      std::string line;
      while(getLine(Base::ifs, line)) {
        idx_column++;
        std::stringstream line_ss(line);

        readColumn(line_ss, temp_dimension, temp_column);
        boundaryMatrix.set_dimension(idx_column, temp_dimension);
        boundaryMatrix.set_column(idx_column, temp_column);
      }

      if(!Base::closeFile())
        return false;

      return true;

    }

  };


}  // namespace stn
