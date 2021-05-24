/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include "commons.hpp"
#include "sorted_matrix.hpp"
#include "sorted_bars.hpp"
#include "reduction.hpp"

namespace stn {

  template<class ReductionAlgorithm>
  class Homology {
  private:
    ReductionAlgorithm reduction;

  public:
    Homology()
      : reduction()
    {}

    template<typename ColumnType = VectorColumn>
    void compute(ViewFiniteBars<ColumnType>& finite_bars,
                 ViewInfiniteBars<ColumnType>& infinite_bars) {
      reduction(finite_bars, infinite_bars);

      const index_t n_columns = finite_bars.get_n_columns();
      std::vector<bool> infinite(n_columns, true);
      for(index_t idx_col = 0; idx_col < n_columns; ++idx_col) {
        if(!finite_bars.is_empty(idx_col)) {
          index_t death = idx_col;
          index_t birth = finite_bars.get_max_index(idx_col);
          infinite[birth] = false;
          infinite[death] = false;
        }
      }

      dimension_t n_dimensions = finite_bars.get_n_dimensions();
      std::vector<index_t> start_infinite_dimension(n_dimensions, 0);
      std::vector<index_t> n_infinite_bars_per_dimension(n_dimensions, 0);
      for(dimension_t dim = 0; dim < n_dimensions; ++dim) {
        if(dim > 0) {
          start_infinite_dimension[dim] = start_infinite_dimension[dim-1]
            + n_infinite_bars_per_dimension[dim-1];
        }

        index_t start = finite_bars.get_start_dimension(dim);
        index_t end = start + finite_bars.get_n_columns_per_dimension(dim);
        for(index_t idx_view = start; idx_view < end; ++idx_view) {
          index_t idx_col = finite_bars.get_view(idx_view);
          if(infinite[idx_col]) {
            infinite_bars.set_birth(idx_col, idx_col);

            infinite_bars.set_view(start_infinite_dimension[dim]
                                   + n_infinite_bars_per_dimension[dim], idx_col);
            ++n_infinite_bars_per_dimension[dim];
          }

          else {
            infinite_bars.clear(idx_col);
          }
        }
      }
      infinite_bars.set_start_dimension(start_infinite_dimension);
      infinite_bars.set_n_columns_per_dimension(n_infinite_bars_per_dimension);

      std::vector<index_t> start_finite_dimension(n_dimensions, 0);
      std::vector<index_t> n_finite_bars_per_dimension(n_dimensions, 0);
      for(dimension_t dim = 0; dim < n_dimensions; ++dim) {
        if(dim > 0) {
          start_finite_dimension[dim] = start_finite_dimension[dim-1]
            + n_finite_bars_per_dimension[dim-1];
        }
        start_finite_dimension[dim] = finite_bars.get_start_dimension(dim);
        index_t end = start_finite_dimension[dim]
          + finite_bars.get_n_columns_per_dimension(dim);
        for(index_t idx_view = start_finite_dimension[dim]; idx_view < end; ++idx_view) {
          index_t idx_col = finite_bars.get_view(idx_view);

          if(!finite_bars.is_empty(idx_col)) {
            index_t death = idx_col;
            index_t birth = finite_bars.get_max_index(idx_col);

            finite_bars.set_birth(idx_col, birth);
            finite_bars.set_death(idx_col, death);

            finite_bars.set_view(start_finite_dimension[dim]
                                 + n_finite_bars_per_dimension[dim], idx_col);
            ++n_finite_bars_per_dimension[dim];
          }

          else {
            finite_bars.clear(idx_col);
          }
        }
      }

      for(dimension_t dim = n_dimensions - 2; dim >= 0 ; --dim) {
        start_finite_dimension[dim+1] = start_finite_dimension[dim];
        n_finite_bars_per_dimension[dim+1] = n_finite_bars_per_dimension[dim];
      }
      start_finite_dimension[0] = 0;
      // This is probably wrong
      n_finite_bars_per_dimension[0] = 0;

      finite_bars.set_start_dimension(start_finite_dimension);
      finite_bars.set_n_columns_per_dimension(n_finite_bars_per_dimension);
    }

  };

} // namespace stn
