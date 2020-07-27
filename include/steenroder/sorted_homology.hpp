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
  class ViewHomology {
  private:
    ReductionAlgorithm reduction;

  public:
    ViewHomology()
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




  // /*  Author: Guillaume Tauzin
  //     License: GPLv3
  // */

  // #pragma once

  // #include "commons.hpp"
  // #include "sorted_matrix.hpp"
  // #include "sorted_bars.hpp"
  // #include "reduction.hpp"

  // namespace stn {

  //   template<class ReductionAlgorithm>
  //   class ViewHomology {
  //   private:
  //     ReductionAlgorithm reduction;
  //     index_t n_finite_bars;
  //     index_t n_infinite_bars;

  //   public:
  //     ViewHomology()
  //       : reduction()
  //       , n_finite_bars(0)
  //       , n_infinite_bars(0)
  //     {}

  //     template<typename ColumnType = VectorColumn>
  //     void compute(ViewFiniteBars<ColumnType>& finite_bars,
  //                  ViewInfiniteBars<ColumnType>& infinite_bars) {
  //       reduction(finite_bars, infinite_bars);


  //       std::map<int, bool> infinite;
  //       for(index_t idx_col = 0; idx_col < finite_bars.get_n_columns(); ++idx_col) {
  //         infinite[idx_col] = true;
  //         if(!finite_bars.is_empty(idx_col)) {
  //           index_t death = idx_col;
  //           index_t birth = finite_bars.get_max_index(idx_col);
  //           infinite[birth] = false;
  //           infinite[death] = false;
  //         }
  //       }

  //       for(std::map<int, bool>::iterator it = infinite.begin(); it != infinite.end(); ++it) {
  //         if(it->second) {
  //           ++n_infinite_bars;
  //         }

  //         if(!finite_bars.is_empty(it->first)) {
  //           ++n_finite_bars;
  //         }
  //       }

  //       // infinite should be an array
  //       // iterate on dim
  //       std::vector<dimension_t> infinite_bars_dimensions(n_infinite_bars);
  //       n_infinite_bars = 0;
  //       for(std::map<int, bool>::iterator it = infinite.begin(); it != infinite.end(); ++it) {
  //         if(it->second) {
  //           infinite_bars.set_birth(it->first, it->first);
  //           infinite_bars_dimensions[n_infinite_bars] =
  //             finite_bars.get_dimension(finite_bars.get_column_index(it->first));

  //           infinite_bars.replace(it->first, n_infinite_bars);
  //           ++n_infinite_bars;
  //         }
  //       }
  //       infinite_bars.set_n_columns(n_infinite_bars);
  //       infinite_bars.create_view(infinite_bars_dimensions);

  //       std::vector<dimension_t> finite_bars_dimensions(n_finite_bars);
  //       n_finite_bars = 0;
  //       // iterate on dim
  //       for(std::map<int, bool>::iterator it = infinite.begin(); it != infinite.end(); ++it) {
  //         if(!finite_bars.is_empty(it->first)) {
  //           index_t death = it->first;
  //           index_t birth = finite_bars.get_max_index(it->first);

  //           finite_bars.set_birth(it->first, birth);
  //           finite_bars.set_death(it->first, death);
  //           finite_bars_dimensions[n_finite_bars] =
  //             finite_bars.get_dimension(finite_bars.get_column_index(birth));

  //           finite_bars.replace(it->first, n_finite_bars);
  //           ++n_finite_bars;
  //         }
  //       }
  //       finite_bars.set_n_columns(n_finite_bars);
  //       finite_bars.create_view(finite_bars_dimensions);
  //     }

  //   };



  // } // namespace stn
