/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include "commons.hpp"
#include "sorted_matrix.hpp"
#include "vector_column.hpp"

namespace stn {

  template<typename ColumnType = VectorColumn>
  class SimplexMatrix : public ViewMatrix<ColumnType> {
  private:
    using Base = ViewMatrix<ColumnType>;

    void build_simplex(ColumnType& simplex, const ColumnType& boundary,
                       const dimension_t dim,
                       const ViewMatrix<ColumnType>& boundary_matrix) {
      if(dim == 1) {
        simplex = boundary;
        return;
      }

      VectorColumn temp;
      for(index_t idx_row = 0; idx_row < boundary.size(); ++idx_row) {
        boundary_matrix.get_column(boundary[idx_row], temp);
        simplex |= temp;
      }

      temp.clear();
      build_simplex(temp, simplex, dim - 1, boundary_matrix);
    }

    void init_simplices(const ViewMatrix<ColumnType>& boundary_matrix,
                        const dimension_t dimension_d,
                        const dimension_t dimension_d_k) {
      VectorColumn simplex, boundary;

      index_t start = boundary_matrix.get_start_dimension(dimension_d);
      index_t end = start + boundary_matrix.get_n_columns_per_dimension(dimension_d);
      for(index_t idx_view = start; idx_view < end; ++idx_view) {
        index_t idx_col = boundary_matrix.get_view(idx_view);

        ColumnType simplex;
        boundary_matrix.get_column(idx_col, boundary);
        build_simplex(simplex, boundary, dimension_d, boundary_matrix);
        Base::set_column(idx_col, simplex);
      }

      start = boundary_matrix.get_start_dimension(dimension_d_k);
      end = start + boundary_matrix.get_n_columns_per_dimension(dimension_d_k);
      for(index_t idx_view = start; idx_view < end; ++idx_view) {
        index_t idx_col = boundary_matrix.get_view(idx_view);

        ColumnType simplex;
        boundary_matrix.get_column(idx_col, boundary);
        build_simplex(simplex, boundary, dimension_d_k, boundary_matrix);
        Base::set_column(idx_col, simplex);
      }
    }


  public:
    SimplexMatrix(const ViewMatrix<ColumnType>& boundary_matrix_in,
                  const dimension_t dimension_d_in,
                  const dimension_t dimension_d_k_in)
      : Base(boundary_matrix_in)
    {
      init_simplices(boundary_matrix_in, dimension_d_in, dimension_d_k_in);
    };

    SimplexMatrix(const SimplexMatrix& other) {
      *this = other;
    }

    using Base::get_n_columns;

    index_t is_in(const index_t min_idx, const dimension_t dim,
                  const VectorColumn& candidate) const {
      VectorColumn temp_col;
      index_t start = Base::get_start_dimension(dim);
      index_t end = start + Base::get_n_columns_per_dimension(dim);
      for(index_t idx_view = start; idx_view < end; ++idx_view) {
        index_t idx_col = Base::get_view(idx_view);
        if(idx_col >= min_idx) {
          Base::get_column(idx_col, temp_col);
          if(temp_col == candidate) {
            return idx_col;
          }
        }
      }
      return -1;
    }


    // Format: each line represents a column, first number is attribute, other numbers are the content of the column
    bool save_ascii(const std::string& name, const std::string& output_filename) {
      std::string filename = output_filename + "_" + name + ".dat";
      std::ofstream output_stream(filename.c_str());
      if(output_stream.fail())
        return false;

      VectorColumn temp_col;
      //output_stream << "# dim " << (index_t) dimension << std::endl;

      for(index_t col_idx =0; col_idx < Base::get_n_columns(); ++col_idx) {
        Base::get_column(col_idx, temp_col);
        for(index_t row_idx = 0; row_idx < (index_t) temp_col.size(); ++row_idx)
          output_stream << " " << temp_col[row_idx];
        output_stream << std::endl;
      }

      output_stream.close();
      return true;
    }

    // Format: n_columns % att1 % N1 % row1 row2 % ...% rowN1 % att2 % N2 % ...
    bool save_binary(const std::string& name, const std::string& output_filename) {
      return true;
    }

  };


} // namespace stn
