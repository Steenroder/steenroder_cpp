/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include "commons.hpp"
#include "sorted_matrix.hpp"
#include "sorted_bars.hpp"
#include "vector_column.hpp"
#include "reduction.hpp"

namespace stn {

  template<typename ReductionAlgorithm, typename ColumnType = VectorColumn>
  class Steenrod {
  private:
    ReductionAlgorithm reduction;
    const dimension_t d;
    const dimension_t k;
    const index_t n_cells;
    const SimplexMatrix<ColumnType>& simplex_matrix;

    bool calculate_index(const index_t idx_vertex, const ColumnType& a_U_b,
                         const ColumnType& bar, const ColumnType& a_bar_U_b_bar) {

      index_t idx_x = std::lower_bound(a_U_b.begin(), a_U_b.end(), bar[idx_vertex])
        - a_U_b.begin();
      index_t idx_x_bar =  std::lower_bound(a_bar_U_b_bar.begin(),
                                            a_bar_U_b_bar.end(), bar[idx_vertex])
        - a_bar_U_b_bar.begin();

      return (idx_x + idx_x_bar ) % 2;
    }

  public:
    Steenrod(const dimension_t d_in, const dimension_t k_in, index_t n_cells_in,
             const SimplexMatrix<ColumnType>& simplex_matrix)
      : reduction()
      , d(d_in)
      , k(k_in)
      , n_cells(n_cells_in)
      , simplex_matrix(simplex_matrix)
    {}

    void compute(ViewFiniteBars<ColumnType>& cohomology_finite_bars,
                 ViewInfiniteBars<ColumnType>& cohomology_infinite_bars,
                 Bars<ColumnType>& steenrod_bars) {
      index_t n_bars = 0;
      index_t start = cohomology_finite_bars.get_start_dimension(d);
      index_t end = start + cohomology_finite_bars.get_n_columns_per_dimension(d);
      for(index_t idx_view = start; idx_view < end; ++idx_view) {
        index_t idx_col = cohomology_finite_bars.get_view(idx_view);
        index_t birth = cohomology_finite_bars.get_birth(idx_col);
        ColumnType cohomology_representative;
        cohomology_finite_bars.get_column(idx_col, cohomology_representative);
        ColumnType steenrod_representative;
        steenrod_square(cohomology_representative, birth, steenrod_representative);
        if(steenrod_representative.size()) {
          steenrod_bars.set_column(n_bars, steenrod_representative);
          steenrod_bars.set_birth(n_bars, birth);
          ++n_bars;
        }
        cohomology_finite_bars.clear(idx_col);
      }

      start = cohomology_infinite_bars.get_start_dimension(d);
      end = start + cohomology_infinite_bars.get_n_columns_per_dimension(d);
      for(index_t idx_view = start; idx_view < end; ++idx_view) {
        index_t idx_col = cohomology_infinite_bars.get_view(idx_view);
        index_t birth = cohomology_infinite_bars.get_birth(idx_col);
        ColumnType cohomology_representative;
        cohomology_infinite_bars.get_column(idx_col, cohomology_representative);
        ColumnType steenrod_representative;
        steenrod_square(cohomology_representative, birth, steenrod_representative);

        if(steenrod_representative.size()) {
          steenrod_bars.set_column(n_bars, steenrod_representative);
          steenrod_bars.set_birth(n_bars, birth);
          ++n_bars;
        }
        cohomology_infinite_bars.clear(idx_col);
      }
      steenrod_bars.set_n_columns(n_bars);
      steenrod_bars.set_n_columns_per_dimension(0, n_bars);

      calculate_deaths(cohomology_finite_bars, steenrod_bars);
    }

    void steenrod_square(const ColumnType& cohomology_representative,
                         const index_t birth,
                         ColumnType& steenrod_representative) {
      ColumnType a, b;
      std::vector<bool> permutations(cohomology_representative.size());
      std::fill(permutations.end() - 2, permutations.end(), true);

      do {
        bool first = true;
        for (int i = 0; i < permutations.size(); ++i) {
          if (permutations[i]) {
            if(first) {
              first = false;
              index_t idx_a = n_cells - 1 - cohomology_representative[i];
              simplex_matrix.get_column(idx_a, a);
            }
            else {
              index_t idx_b = n_cells - 1 - cohomology_representative[i];
              simplex_matrix.get_column(idx_b, b);
              break;
            }
          }
        }

        ColumnType a_U_b = a | b;
        if(k == a_U_b.size() - 1 - d) {
          index_t idx_a_U_b = simplex_matrix.is_in(0, d+k, a_U_b);
          if(idx_a_U_b != -1) {
            ColumnType a_bar = b - a;
            ColumnType b_bar = a - b;
            ColumnType a_bar_U_b_bar = a_bar | b_bar;
            index_t idx_vertex = 0;

            bool pos_a = calculate_index(idx_vertex, a_U_b, a_bar, a_bar_U_b_bar);
            bool pos_b = calculate_index(idx_vertex, a_U_b, b_bar, a_bar_U_b_bar);

            if(pos_a ^ pos_b) {
              bool purity = true;
              for(idx_vertex = 1; idx_vertex < a_bar.size(); ++idx_vertex) {
                if(purity) {
                  bool pos_a_temp = calculate_index(idx_vertex, a_U_b, a_bar, a_bar_U_b_bar);
                  bool pos_b_temp = calculate_index(idx_vertex, a_U_b, b_bar, a_bar_U_b_bar);

                  purity = purity && (pos_a == pos_a_temp) && (pos_b == pos_b_temp);
                }
                else {
                  break;
                }
              }

              if(purity) {
                a_U_b = ColumnType(1, n_cells - idx_a_U_b - 1);
                steenrod_representative |= a_U_b;

              }
            }
          }
        }

      } while (std::next_permutation(permutations.begin(), permutations.end()));

    }

    void calculate_deaths(ViewFiniteBars<ColumnType>& cohomology_finite_bars,
                          Bars<ColumnType>& steenrod_bars) {
      const index_t n_columns_R = cohomology_finite_bars.get_n_columns();
      const index_t n_columns_S = steenrod_bars.get_n_columns();

      std::vector<index_t>& view = cohomology_finite_bars.get_view();
      index_t start = cohomology_finite_bars.get_start_dimension(d + k);
      index_t end = start + cohomology_finite_bars.get_n_columns_per_dimension(d + k);
      view = std::vector<index_t>(view.begin() + start, view.begin() + end);

      view.resize(view.size() + n_columns_S);
      std::iota(view.end() - n_columns_S, view.end(), n_columns_R);

      // sort R based on death
      std::sort(view.begin(), view.end() - n_columns_S,
                [&](index_t idx_a, index_t idx_b) {
                  return cohomology_finite_bars.get_death(idx_a) <=
                    cohomology_finite_bars.get_death(idx_b);
                });

      // sort S based on birth
      std::sort(view.end() - n_columns_S, view.end(),
                [&](index_t idx_a, index_t idx_b) {
                  return steenrod_bars.get_birth(idx_a - n_columns_R) <=
                    steenrod_bars.get_birth(idx_b - n_columns_R);
                });

      const index_t n_columns = cohomology_finite_bars.get_n_columns();
      std::vector<index_t> pivot_lookup(n_columns_R+n_columns_S, -1);

      for(index_t idx_view = 0; idx_view < view.size() - n_columns_S; ++idx_view) {
        index_t idx_col = view[idx_view];
        index_t pivot = cohomology_finite_bars.get_max_index(idx_col);
        if(pivot != -1) {
          pivot_lookup[pivot] = idx_col;
        }
      }

      index_t n_columns_R_birth_S = 0;
      // for each column of S
      for(index_t idx_view_S = view.size() - n_columns_S;
          idx_view_S < view.size(); ++idx_view_S) {
        ColumnType reduced_steenrod_representative;

        steenrod_bars.get_column(view[idx_view_S] - n_columns_R,
                                 reduced_steenrod_representative);

        index_t birth_S = steenrod_bars.get_birth(view[idx_view_S] - n_columns_R);

        bool first_reduction = true;
        while(first_reduction || (view[n_columns_R_birth_S] <= birth_S)) {
          first_reduction = false;

          for(index_t idx_view_S_temp = view.size() - n_columns_S;
              idx_view_S_temp <= idx_view_S; ++idx_view_S_temp) {
            index_t idx_col = view[idx_view_S_temp];

            index_t pivot = steenrod_bars.get_max_index(idx_col - n_columns_R);
            while(pivot != -1 && pivot_lookup[pivot] && pivot_lookup[pivot] != -1) {
              ColumnType temp_col;

              if(pivot_lookup[pivot] < view[n_columns_R_birth_S]) {
                cohomology_finite_bars.get_column(pivot_lookup[pivot], temp_col);
                steenrod_bars.add(temp_col, idx_col - n_columns_R);

              }
              else if(pivot_lookup[pivot] >= n_columns_R && pivot_lookup[pivot] < idx_col) {
                steenrod_bars.get_column(pivot_lookup[pivot] - n_columns_R, temp_col);
                steenrod_bars.add(temp_col, idx_col - n_columns_R);
              }
              else {
                break;
              }
              pivot = steenrod_bars.get_max_index(idx_col - n_columns_R);
            }

            if(pivot != -1 && (pivot_lookup[pivot] == -1 || pivot_lookup[pivot] >= n_columns_R)) {
              pivot_lookup[pivot] = idx_col;
            }
            if(pivot == -1) { // fully reduced
              index_t death = steenrod_bars.get_death(idx_col - n_columns_R);
              if(idx_view_S_temp == idx_view_S) { // last S born dead
                death = steenrod_bars.get_birth(idx_col - n_columns_R);
              }
              else {
                if(death == -1)
                  death = view[n_columns_R_birth_S];
              }
              steenrod_bars.set_death(idx_col - n_columns_R, death);
            }
          }

          ++n_columns_R_birth_S;
        }


      }
    }

  };

} // namespace stn
