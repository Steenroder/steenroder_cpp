/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include "commons.hpp"
#include "attribute_matrix.hpp"
#include "vector_column.hpp"

namespace stn {

  template<typename ColumnType = VectorColumn>
  class BoundaryMatrix : public AttributeMatrix<ColumnType> {
  private:
    using Base = AttributeMatrix<ColumnType>;

  public:
    using Base::Base;
    using Base::get_column;
    using Base::load_binary;
    using Base::save_ascii;
    using Base::save_binary;
    using Base::set_n_columns;
    using Base::get_n_columns;
    using Base::load_ascii;

    dimension_t get_dimension(index_t idx_col) const {
      return Base::get_attribute(idx_col);
    }

    void set_dimension(index_t idx_col, dimension_t dim) {
      Base::set_attribute(idx_col, dim);
    }

    dimension_t get_max_dimension() const {
      return Base::get_max_attribute();
    }

    void dualize() {
      std::vector<dimension_t> dual_dimensions;
      std::vector<std::vector<index_t>> dual_matrix;

      index_t n_columns = Base::get_n_columns();
      dual_matrix.resize(n_columns);
      dual_dimensions.resize(n_columns);

      std::vector<index_t> dual_sizes(n_columns, 0);

      VectorColumn temp_col;
      for(index_t idx_col = 0; idx_col < n_columns; ++idx_col) {
        Base::get_column(idx_col, temp_col);
        for(index_t idx_col = 0; idx_col < (index_t) temp_col.size(); ++idx_col)
          ++dual_sizes[n_columns - 1 - temp_col[idx_col]];
      }

      #pragma omp parallel for
      for(index_t idx_col = 0; idx_col < n_columns; idx_col++)
        dual_matrix[idx_col].reserve(dual_sizes[idx_col]);

      for(index_t idx_col = 0; idx_col < n_columns; idx_col++) {
        Base::get_column(idx_col, temp_col);
        for(index_t idx_row = 0; idx_row < (index_t) temp_col.size(); idx_row++)
          dual_matrix[n_columns - 1 - temp_col[idx_row]].push_back(n_columns - 1 - idx_col);
      }

      const dimension_t n_dimensions = get_max_dimension() + 1;
      #pragma omp parallel for
      for( index_t idx_col = 0; idx_col < n_columns; ++idx_col) {
        dual_dimensions[n_columns - 1 - idx_col] = n_dimensions - 1 - get_dimension(idx_col);
      }

      #pragma omp parallel for
      for(index_t idx_col = 0; idx_col < n_columns; ++idx_col) {
        std::reverse(dual_matrix[idx_col].begin(), dual_matrix[idx_col].end());
      }

      load_vector_vector(dual_matrix, dual_dimensions);
    }

    template<typename index_type, typename dimemsion_type>
    void load_vector_vector(const std::vector<std::vector<index_type>>& input_matrix,
                            const std::vector<dimemsion_type>& input_dimensions) {
      const index_t n_columns = (index_t) input_matrix.size();
      Base::set_n_columns(n_columns);
      VectorColumn temp_col;
#pragma omp parallel for private(temp_col)
      for(index_t idx_col = 0; idx_col < n_columns; ++idx_col) {
        set_dimension(idx_col, (dimension_t) input_dimensions[idx_col]);

        index_t n_rows = input_matrix[idx_col].size();
        temp_col.resize(n_rows);

        for(index_t idx_row = 0; idx_row < n_rows; ++idx_row) {
          temp_col[idx_row] = (index_t) input_matrix[idx_col][idx_row];
        }
        Base::set_column(idx_col, temp_col);
      }
    }

  };

} // namespace stn
