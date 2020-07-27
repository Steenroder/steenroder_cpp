/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include "commons.hpp"
#include "boundary_matrix.hpp"
#include "reduction.hpp"

namespace stn {

  template<class ReductionAlgorithm>
  class Homology {
  private:
    ReductionAlgorithm reduction;
    index_t n_finite_bars;
    index_t n_infinite_bars;

  public:
    Homology()
      : reduction()
      , n_finite_bars(0)
      , n_infinite_bars(0)
    {}

    template<typename ColumnType = VectorColumn>
    void compute(FiniteBars<ColumnType>& finite_bars_matrix,
                 InfiniteBars<ColumnType>& infinite_bars_matrix) {
      reduction(finite_bars_matrix, infinite_bars_matrix);

      std::map<int, bool> infinite;
      for(index_t idx_col = 0; idx_col < finite_bars_matrix.get_n_columns(); ++idx_col) {
        infinite[idx_col] = true;
        if(!finite_bars_matrix.is_empty(idx_col)) {
          index_t death = idx_col;
          index_t birth = finite_bars_matrix.get_max_index(idx_col);
          infinite[birth] = false;
          infinite[death] = false;
        }
      }

      for(std::map<int, bool>::iterator it = infinite.begin(); it != infinite.end(); ++it) {
        if(it->second) {
          infinite_bars_matrix.set_birth(it->first, it->first);
          infinite_bars_matrix.set_dimension(it->first,
                                             finite_bars_matrix.get_dimension(it->first));
          //infinite_bars_matrix.replace(it->first, n_infinite_bars);
          ++n_infinite_bars;
        }
      }
      infinite_bars_matrix.set_n_columns(n_infinite_bars);

      for(std::map<int, bool>::iterator it = infinite.begin(); it != infinite.end(); ++it) {
        // std::cout << it->first << ", "
        //           << (index_t) finite_bars_matrix.get_dimension(it->first) << std::endl;
        if(!finite_bars_matrix.is_empty(it->first)) {
          index_t death = it->first;
          index_t birth = finite_bars_matrix.get_max_index(it->first);

          finite_bars_matrix.set_birth(it->first, birth);
          finite_bars_matrix.set_death(it->first, death);
          finite_bars_matrix.set_dimension(it->first,
                                           finite_bars_matrix.get_dimension(birth));

          //finite_bars_matrix.replace(it->first, n_finite_bars);
          ++n_finite_bars;
        }
      }
      finite_bars_matrix.set_n_columns(n_finite_bars);
    }

  };



} // namespace stn
