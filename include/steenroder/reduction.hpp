/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include "commons.hpp"
#include "boundary_matrix.hpp"
#include "sorted_matrix.hpp"

namespace stn {

  template<typename ColumnType = VectorColumn>
  class StandardReduction {
  public:
    void operator()(BoundaryMatrix<ColumnType>& boundary_matrix,
                    BoundaryMatrix<ColumnType>& triangular_matrix) {
      const index_t n_columns = boundary_matrix.get_n_columns();
      triangular_matrix.set_n_columns(n_columns);
      std::vector<index_t> pivot_lookup(n_columns, -1);

      for(index_t idx_col = 0; idx_col < n_columns; ++idx_col) {
        index_t pivot = boundary_matrix.get_max_index(idx_col);
        while(pivot != -1 && pivot_lookup[pivot] != -1) {
          boundary_matrix.add(pivot_lookup[pivot], idx_col);
          triangular_matrix.add(pivot_lookup[pivot], idx_col);
          pivot = boundary_matrix.get_max_index(idx_col);
        }
        if(pivot != -1) {
          pivot_lookup[pivot] = idx_col;
        }
      }
    }
  };


  template<typename ColumnType = VectorColumn>
  class TwistReduction {
  public:
    void operator()(ViewMatrix<ColumnType>& boundary_matrix,
                    ViewMatrix<ColumnType>& triangular_matrix ) {
      const index_t n_columns = boundary_matrix.get_n_columns();
      std::vector< index_t > pivot_lookup(n_columns, -1);

      // for(dimension_t dim = boundary_matrix.get_n_dimensions() - 1; dim >= 1 ; --dim) {
      for(dimension_t dim = 0; dim < boundary_matrix.get_n_dimensions() - 1; ++dim) {
        index_t start = boundary_matrix.get_start_dimension(dim);
        index_t end = start + boundary_matrix.get_n_columns_per_dimension(dim);
        for(index_t view_idx = start; view_idx < end; ++view_idx) {
          index_t col_idx = boundary_matrix.get_view(view_idx);
          index_t pivot = boundary_matrix.get_max_index(col_idx);
          while(pivot != -1 && pivot_lookup[pivot] != -1) {
            boundary_matrix.add(pivot_lookup[pivot], col_idx);
            triangular_matrix.add(pivot_lookup[pivot], col_idx);
            pivot = boundary_matrix.get_max_index(col_idx);
          }

          if(pivot != -1) {
            pivot_lookup[pivot] = col_idx;
            boundary_matrix.clear(pivot);
          }
        }
      }
    }
  };

} // namespace stn
