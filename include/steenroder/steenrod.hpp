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

      calculate_deaths(cohomology_finite_bars,
                       steenrod_bars);
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

    template<typename T> void print_queue(T q) {
      std::cout << "printing queue" << std::endl;
      while(!q.empty()) {
        std::cout << q.top() << " ";
        q.pop();
      }
      std::cout << '\n';
}

    ColumnType reduce(const ViewFiniteBars<ColumnType>& cohomology_finite_bars,
                      const Bars<ColumnType>& steenrod_bars,
                      const std::vector<index_t>& view, const index_t idx_view_steenrod,
                      const index_t n_columns_R,
                      std::priority_queue<index_t, std::vector<index_t>,
                      std::greater<index_t>>& queue) {
      ColumnType reduced_steenrod_representative;
      steenrod_bars.get_column(view[idx_view_steenrod] - n_columns_R,
                               reduced_steenrod_representative);

      std::cout << "REP: " <<  reduced_steenrod_representative;

      index_t pivot = reduced_steenrod_representative.get_max();

      for(index_t idx_view = idx_view_steenrod - 1; idx_view >= 0; --idx_view) {
        index_t pivot = reduced_steenrod_representative.get_max();
        std::cout << idx_view << std::endl;

        // pas forcement cohomology_finite_bars
        ColumnType temp_col;
        index_t temp_pivot;
        // check if it's not idx_view alone
        if(view[idx_view] < n_columns_R) {
          temp_pivot = cohomology_finite_bars.get_max_index(view[idx_view]);

          if(temp_pivot == pivot) {
            cohomology_finite_bars.get_column(view[idx_view], temp_col);
            reduced_steenrod_representative += temp_col;
            queue.push(view[idx_view]);
            print_queue(queue);
          }
        }
        else {
          temp_pivot = steenrod_bars.get_max_index(view[idx_view] - n_columns_R);

          if(temp_pivot == pivot) {
            steenrod_bars.get_column(view[idx_view], temp_col);
            reduced_steenrod_representative += temp_col;
            queue.push(view[idx_view]);
            print_queue(queue);
          }
        }

        if(!reduced_steenrod_representative.size()) {
          break;
        }
      }

      std::cout << reduced_steenrod_representative;
      return reduced_steenrod_representative;
    }

    void update_steenrod_representatives
    (const ViewFiniteBars<ColumnType>& cohomology_finite_bars,
     Bars<ColumnType>& steenrod_bars,
     const std::vector<index_t>& view, const index_t idx_view_steenrod,
     const index_t n_columns_R, const index_t n_columns_S,
     std::vector<std::priority_queue<index_t, std::vector<index_t>,
     std::greater<index_t>>>& queues,
     const index_t next_birth) {

      for(index_t idx_view = view.size() - n_columns_S; idx_view <= idx_view_steenrod; ++idx_view) {
        index_t idx_col = view[idx_view];

        if(steenrod_bars.get_max_index(idx_col - n_columns_R) != -1) {
          ColumnType steenrod_representative;
          steenrod_bars.get_column(idx_col - n_columns_R,
                                   steenrod_representative);
          std::cout << "REP: " <<  steenrod_representative;

          while(queues[idx_col - n_columns_R].top() <= next_birth
                && !queues[idx_col - n_columns_R].empty()) {
            ColumnType temp_col;
            index_t temp_idx = queues[idx_col - n_columns_R].top();
            cohomology_finite_bars.get_column(temp_idx, temp_col);
            steenrod_representative += temp_col;
            queues[idx_col - n_columns_R].pop();
          }
          steenrod_bars.set_column(idx_col - n_columns_R,
                                   steenrod_representative);
        }
      }
    }

    index_t get_max_index(const index_t idx_col, const index_t n_columns_R,
                          const ViewFiniteBars<ColumnType>& cohomology_finite_bars,
                          const Bars<ColumnType>& steenrod_bars) {
      if(idx_col < n_columns_R) {
        return cohomology_finite_bars.get_max_index(idx_col);
      }
      else {
        return steenrod_bars.get_max_index(idx_col - n_columns_R);
      }
    }

    void sort_by_pivot(const ViewFiniteBars<ColumnType>& cohomology_finite_bars,
                       const Bars<ColumnType>& steenrod_bars,
                       std::vector<index_t>& view, const index_t idx_view_steenrod,
                       const index_t n_columns_R) {
      std::cout << "SORTING" << std::endl;
      for(int i = 0; i < view.size(); ++i)
        std::cout << view[i] <<  " ";
      std::cout << std::endl;
      std::cout << view.size() << std::endl;

      // account for the fact that steenrod reps should be after cohom repsof the same pivot
      std::sort(view.begin(), view.begin() + idx_view_steenrod,
                [&](index_t idx_a, index_t idx_b) {
                  //maybe view[idx] here
                  std::cout << idx_a << ", " << idx_b << std::endl;
                  return get_max_index(idx_a, n_columns_R,
                                       cohomology_finite_bars, steenrod_bars)
                    > get_max_index(idx_b, n_columns_R,
                                    cohomology_finite_bars, steenrod_bars);
                });

      for(int i = 0; i < view.size(); ++i)
        std::cout << view[i] <<  " ";
      std::cout << std::endl;
      std::cout << view.size() << std::endl;

    }

    index_t compute_finite_death(index_t birth, const Bars<ColumnType>& steenrod_bars,
                                 std::priority_queue<index_t, std::vector<index_t>, std::greater<index_t>>>& queue) {
      index_t max = -1;

      while(!queue.empty()) {
            index_t temp = queue.top();
            // get birth from cohom or steenrod
            cohomology_finite_bars.get_birth(temp_idx, temp_col);
            steenrod_representative += temp_col;
            queue.pop();
      if(max <= birth) {
        return birth;
      }
      else {
        return max;
      }
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

      // sort R based on pivot indices i.e. birth
      std::sort(view.begin(), view.end() - n_columns_S,
                [&](index_t idx_a, index_t idx_b) {
                  return cohomology_finite_bars.get_birth(idx_a) >
                    cohomology_finite_bars.get_birth(idx_b);
                });

      // sort S based on pivot indices
      std::sort(view.end() - n_columns_S, view.end(),
                [&](index_t idx_a, index_t idx_b) {
                  return steenrod_bars.get_birth(idx_a - n_columns_R) >
                    steenrod_bars.get_birth(idx_b - n_columns_R);
                });

      for(int i = 0; i < view.size(); ++i)
        std::cout << view[i] <<  " ";
      std::cout << std::endl;
      std::cout << view.size() << std::endl;


      for(index_t idx_view = view.size() - 1 - n_columns_S; idx_view >= 0; --idx_view) {
        std::cout << cohomology_finite_bars.get_max_index(view[idx_view]) <<  " ";
      }
      std::cout << std::endl;
      std::cout << view.size() << std::endl;


      std::vector<std::priority_queue<index_t, std::vector<index_t>,
                                      std::greater<index_t>>> queues(n_columns_S);

      // for each column of S
      for(index_t idx_view = view.size() - n_columns_S; idx_view < view.size(); ++idx_view) {
        index_t idx_col = view[idx_view];

        ColumnType reduced_steenrod_representative =
          reduce(cohomology_finite_bars, steenrod_bars, view, idx_view,
                 n_columns_R, queues[idx_col - n_columns_R]);


        if(reduced_steenrod_representative.size()) {
          std::cout << "IMMORTAL" << std::endl;
          steenrod_bars.set_death(idx_col - n_columns_R, -1);
        }

        else {
          std::cout << "MORTAL" << std::endl;
          index_t birth = steenrod_bars.get_birth(idx_col - n_columns_R);
          index_t death = compute_finite_death(birth, queues[idx_col - n_columns_R]);

          index_t max = queues[idx_col - n_columns_R].empty()
            ? -1 : queues[idx_col - n_columns_R].top();
          std::cout << "birth, max: " << birth << ", " <<  max << std::endl;

          if(death == birth) {
            steenrod_bars.set_death(idx_col - n_columns_R, birth);
            steenrod_bars.clear(idx_col - n_columns_R);
            queues[idx_col - n_columns_R] =
              std::priority_queue<index_t, std::vector<index_t>,
                                  std::greater<index_t>>();
            std::cout << "BORN DEAD" << std::endl;
            // view[view_idx]
          }

          else {
            steenrod_bars.set_death(idx_col - n_columns_R, max);
            // SAVE SOMEWHERE ELSE THE REP
          }

          if(idx_view + 1 < view.size()) {
            index_t next_birth = steenrod_bars.get_birth(view[idx_view + 1] - n_columns_R);
            std::cout << "NEXT BIRTH: " << next_birth << std::endl;

            // update_steenrod_representatives(cohomology_finite_bars, steenrod_bars,
            //                                 view, idx_view, n_columns_R, n_columns_S,
            //                                 queues, next_birth);

            // sort_by_pivot(cohomology_finite_bars, steenrod_bars, view, idx_view,
            //               n_columns_R);
          }
        }
      }

      // std::vector<index_t>view_S(view.end() - n_columns_S, view.end());
      // steenrod_bars.set_view(view_S);
    }

  };

} // namespace stn
