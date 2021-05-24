/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include "commons.hpp"

namespace stn {

  template<class Matrix>
  class MatrixVector {
  private:
    std::vector<Matrix> matrix_vector;

    void break_down(Matrix& matrix) {
      index_t n_columns = matrix.get_n_columns();
      dimension_t n_matrices = matrix.get_max_dimension();

      for(dimension_t dimension = 0; dimension < n_matrices; ++dimension){
        matrix_vector[dimension].set_n_columns(n_columns);
      }

      dimension_t dimension;
      column temp_col;
      std::vector<index_t> indices(n_matrices, 0);
      for(index_t col = 0; col < n_columns; ++col) {
        dimension = simplex_matrix.get_dimension(idx_cur_col);
        matrix.get_column(col, temp_col);

        matrix_vector[dimension].set_dimension(indices[dimension],
                                               n_columns - 1 - col);
        matrix_vector[dimension].set_column(indices[dimension], temp_col);

        ++indices[dimension];
      }

      for(dimension_t dimension = 0; dimension < n_matrices; ++dimension){
        matrix_vector[dimension].set_n_columns(indices[dimension]);
      }
    }

  public:
    MatrixVector(Matrix& matrix)
    {
      break_down(matrix);
    }

  };

} // namespace stn
