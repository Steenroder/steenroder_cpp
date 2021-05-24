/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include <numeric>
#include <algorithm>

#include "commons.hpp"
#include "sparse_matrix.hpp"
#include "vector_column.hpp"

namespace stn {

  template<typename ColumnType = VectorColumn>
  class ViewMatrix : public SparseMatrix<ColumnType> {
  private:
    using Base = SparseMatrix<ColumnType>;

  protected:
    std::vector<index_t> view;
    dimension_t n_dimensions;
    std::vector<index_t> n_columns_per_dimension;
    std::vector<index_t> start_dimension;

  public:
    ViewMatrix()
      : Base()
      , view()
      , n_dimensions()
      , n_columns_per_dimension()
      , start_dimension()
    {};

    ViewMatrix(const SparseMatrix<ColumnType>& sparse_matrix_in,
               const std::vector<dimension_t>& dimensions_in)
      : Base(sparse_matrix_in)
      , view(sparse_matrix_in.get_n_columns())
      , n_dimensions(*std::max_element(dimensions_in.begin(),
                                       dimensions_in.end()) + 1)
      , n_columns_per_dimension(n_dimensions)
      , start_dimension(n_dimensions)
    {
      create_view(dimensions_in);
    }

    ViewMatrix(const ViewMatrix<ColumnType>& view_matrix_in)
      : Base(view_matrix_in)
      , view(view_matrix_in.get_n_columns())
      , n_dimensions(view_matrix_in.get_n_dimensions())
      , n_columns_per_dimension(view_matrix_in.get_n_columns_per_dimension())
      , start_dimension(view_matrix_in.get_start_dimension())
    {
      for(index_t idx_col = 0; idx_col < view_matrix_in.get_n_columns(); ++idx_col) {
        view[idx_col] = view_matrix_in.get_view(idx_col);
      }
    }

    ViewMatrix(const index_t n_columns_in,
               const dimension_t n_dimensions_in)
      : Base(n_columns_in)
      , view(n_columns_in)
      , n_dimensions(n_dimensions_in)
      , n_columns_per_dimension(n_dimensions_in, n_columns_in)
      , start_dimension(n_dimensions_in, 0)
    {
      std::iota(view.begin(), view.end(), 0);
    }

    // ViewMatrix(const ViewMatrix<ColumnType>& other) {
    //   *this = other;
    // }

    dimension_t get_dimension(index_t idx) const {
      dimension_t dim = 0;
      //idx = view[idx];
      for(index_t n_columns : n_columns_per_dimension) {
        if(idx < n_columns) {
          return dim;
        }
        idx -= n_columns;
        ++dim;
      }
      return -1;
    }

    using Base::get_column;

    index_t get_n_columns() const {
      return Base::get_n_columns();
    }

    void set_n_columns(index_t n_columns) {
      view.resize(n_columns);
      Base::set_n_columns(n_columns);
    }

    dimension_t get_n_dimensions() const {
      return n_dimensions;
    }

    const std::vector<index_t>& get_n_columns_per_dimension() const {
      return n_columns_per_dimension;
    }

    index_t get_n_columns_per_dimension(const dimension_t dim) const {
      return n_columns_per_dimension[dim];
    }

    void set_n_columns_per_dimension(const dimension_t dim,
                                     const index_t n_columns) {
      n_columns_per_dimension[dim] = n_columns;
    }

    void set_n_columns_per_dimension(const std::vector<index_t>& n_columns_per_dimension_in) {
      for(index_t dim = 0; dim < n_columns_per_dimension.size(); ++dim) {
        n_columns_per_dimension[dim] = n_columns_per_dimension_in[dim];
      }
    }


    const std::vector<index_t>& get_start_dimension() const {
      return start_dimension;
    }

    index_t get_start_dimension(const dimension_t dim) const {
      return start_dimension[dim];
    }

    void set_start_dimension(const dimension_t dim,
                             const index_t start) {
      start_dimension[dim] = start;
    }

    void set_start_dimension(const std::vector<index_t>& start_dimension_in) {
      for(index_t dim = 0; dim < start_dimension.size(); ++dim) {
        start_dimension[dim] = start_dimension_in[dim];
      }
    }

    index_t get_view(const index_t idx_view) const {
      return view[idx_view];
    }

    std::vector<index_t>& get_view() {
      return view;
    }

    void set_view(index_t idx_view, index_t idx_col) {
      view[idx_view] = idx_col;
    }

    void set_view(std::vector<index_t>& view_in) {
      view.swap(view_in);
    }


    void create_view(std::vector<dimension_t> dimensions) {
      n_dimensions = *std::max_element(dimensions.begin(),
                                       dimensions.end()) + 1;

      n_columns_per_dimension.resize(n_dimensions);
      std::fill(n_columns_per_dimension.begin(), n_columns_per_dimension.end(), 0);
      for(auto dim : dimensions) {
        ++n_columns_per_dimension[dim];
      }

      start_dimension.resize(n_dimensions);
      start_dimension[0] = 0;
      for(dimension_t dim = 1; dim < n_dimensions; ++dim) {
        start_dimension[dim] = start_dimension[dim-1] + n_columns_per_dimension[dim-1];
      }

      index_t n_columns = get_n_columns();
      std::iota(view.begin(), view.end(), 0);
      std::sort(view.begin(), view.end(),
                [&](index_t idx_a, index_t idx_b) {
                  return dimensions[idx_a] * n_columns + idx_a  <
                    dimensions[idx_b] * n_columns + idx_b;
                });
    }


    // ViewMatrix<ColumnType>&
    // operator=(const ViewMatrix<ColumnType>& other) {
    //   const index_t n_columns = other.get_n_columns();
    //   set_n_columns(n_columns);

    //   VectorColumn temp_col;
    //   for(index_t idx = 0; idx <  n_columns; ++idx) {
    //     set_attribute(idx, other.get_attribute(idx));
    //     other.get_column(idx, temp_col);
    //     Base::set_column(idx, temp_col);
    //   }

    //   return *this;
    // }

    bool load_ascii(std::string filename) {
      // first count number of columns:
      std::string line;
      std::ifstream dummy(filename .c_str());
      if(dummy.fail())
        return false;

      index_t n_columns = 0;
      while(getline(dummy, line)) {
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
        if(line != "" && line[0] != '#')
          n_columns++;

      }
      set_n_columns(n_columns);
      dummy.close();

      std::vector<dimension_t> dimensions(n_columns, -1);

      std::ifstream input_stream(filename.c_str());
      if(input_stream.fail())
        return false;

      VectorColumn temp_col;
      index_t idx = -1;
      while(getline(input_stream, line)) {
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
        if(line != "" && line[0] != '#') {
          idx++;
          std::stringstream ss(line);

          index_t att;
          ss >> att;
          dimensions[idx] = (index_t) att;

          index_t temp_index;
          temp_col.clear();
          while(ss.good()) {
            ss >> temp_index;
            temp_col.push_back((index_t)temp_index);
          }
          std::sort(temp_col.begin(), temp_col.end());
          Base::set_column(idx, temp_col);
        }
      }
      create_view(dimensions);

      input_stream.close();
      return true;
    }

    bool load_ascii_dual(std::string filename) {
      // first count number of columns:
      std::string line;
      std::ifstream dummy(filename .c_str());
      if(dummy.fail())
        return false;

      index_t n_columns = 0;
      while(getline(dummy, line)) {
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
        if(line != "" && line[0] != '#')
          ++n_columns;

      }
      set_n_columns(n_columns);
      dummy.close();

      std::vector<dimension_t> dimensions(n_columns, -1);

      std::ifstream input_stream(filename.c_str());
      if(input_stream.fail())
        return false;

      VectorColumn temp_col;
      index_t idx_col = -1;
      while(getline(input_stream, line)) {
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
        if(line != "" && line[0] != '#') {
          ++idx_col;
          std::stringstream ss(line);

          index_t dim;
          ss >> dim;
          dimensions[n_columns - 1 - idx_col] = dim;

          index_t idx_row;
          temp_col.clear();
          while(ss.good()) {
            ss >> idx_row;
            Base::matrix[n_columns - 1 - idx_row].push_back(n_columns - 1 - idx_col);
          }
        }
      }

      dimension_t max_dimension = *std::max_element(dimensions.begin(),
                                                    dimensions.end());

      for(index_t idx_col = 0; idx_col < n_columns; idx_col++) {
        std::sort(Base::matrix[idx_col].begin(), Base::matrix[idx_col].end());
      }

      create_view(dimensions);

      input_stream.close();
      return true;
    }

    bool load_binary_dual(std::string filename) {
      return true;
    }

    // Format: each line represents a column, first number is attribute, other numbers are the content of the column
    bool save_ascii(const std::string& name, const std::string& output_filename) {
      std::string filename = output_filename + "_" + name + ".dat";
      std::ofstream output_stream(filename.c_str());
      if(output_stream.fail())
        return false;

      VectorColumn temp_col;
      for(dimension_t dim = 0; dim < n_dimensions; ++dim) {
        output_stream << "# dim " << (index_t) dim << std::endl;

        index_t start = start_dimension[dim];
        index_t end = start + n_columns_per_dimension[dim];
        for(index_t view_idx = start; view_idx < end; ++view_idx) {
          output_stream << view[view_idx] << " ";
          Base::get_column(view[view_idx], temp_col);
          output_stream << temp_col;
        }
      }

      output_stream.close();
      return true;
    }

    // Format: n_columns % att1 % N1 % row1 row2 % ...% rowN1 % att2 % N2 % ...
    bool load_binary(std::string filename) {
      std::ifstream input_stream(filename.c_str(), std::ios_base::binary | std::ios_base::in);
      if(input_stream.fail())
        return false;

      index_t n_columns;
      input_stream.read((char*) &n_columns, sizeof(int64_t));
      set_n_columns((index_t) n_columns);

      std::vector<dimension_t> dimensions(n_columns);

      VectorColumn temp_col;
      for(index_t idx = 0; idx < n_columns; idx++) {
        index_t att;
        input_stream.read((char*) &att, sizeof(int64_t));
        dimensions[idx] = (index_t) att;
        index_t n_rows;
        input_stream.read((char*) &n_rows, sizeof(int64_t));
        temp_col.resize((std::size_t) n_rows);
        for(index_t idx = 0; idx < n_rows; ++idx) {
          index_t row;
          input_stream.read((char*)& row, sizeof(int64_t));
          temp_col[idx] = (index_t) row;
        }
        Base::set_column(idx, temp_col);
      }

      create_view(dimensions);

      input_stream.close();
      return true;
    }

    // Format: n_columns % att1 % N1 % row1 row2 % ...% rowN1 % att2 % N2 % ...
    bool save_binary(const std::string& name, const std::string& output_filename) {
      std::string filename = output_filename + "_" + name + ".dat";

      std::ofstream output_stream(filename.c_str(),
                                  std::ios_base::binary | std::ios_base::out);
      if(output_stream.fail())
        return false;

      const index_t n_columns = Base::get_n_columns();
      output_stream.write((char*) &n_columns, sizeof(int64_t));
      VectorColumn temp_col;
      for(index_t idx = 0; idx < n_columns; idx++) {
        index_t att = get_dimension(idx);
        output_stream.write((char*) &att, sizeof(int64_t));
        Base::get_column(idx, temp_col);
        index_t n_rows = temp_col.size();
        output_stream.write((char*) &n_rows, sizeof(int64_t));
        for(index_t row_idx = 0; row_idx < (index_t) temp_col.size(); row_idx++) {
          index_t row = temp_col[row_idx];
          output_stream.write((char*) &row, sizeof(int64_t));
        }
      }

      output_stream.close();
      return true;
    }

  };

} // namespace stn
