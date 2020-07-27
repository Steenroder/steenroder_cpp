/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include "commons.hpp"
#include "vector_column.hpp"

namespace stn {

  template<typename ColumnType>
  class SparseMatrix {
  protected:
    std::vector<ColumnType> matrix;
    thread_local_storage<ColumnType> temp_column_buffer;

  public:
    SparseMatrix()
      : matrix()
    {}

    SparseMatrix(const index_t n_columns_in)
      : matrix(n_columns_in)
    {}

    index_t get_n_columns() const {
      return (index_t) matrix.size();
    }

    void set_n_columns(const index_t n_columns) {
      matrix.resize(n_columns);
    }

    void get_column(const index_t idx, ColumnType& col) const {
      col = matrix[idx];
    }
    void set_column(const index_t idx, const ColumnType& col) {
      matrix[idx] = col;
    }

    bool is_empty(const index_t idx) const {
      return matrix[idx].empty();
    }

    index_t get_max_index(const index_t idx) const {
      return matrix[idx].empty() ? -1 : matrix[idx].back();
    }

    void remove_max(const index_t idx) {
      matrix[idx].pop_back();
    }

    void clear(const index_t idx) {
      matrix[idx].clear();
    }

    void swap(const index_t idx_1, const index_t idx_2) {
      std::swap(matrix[idx_1], matrix[idx_2]);
    }

    void erase(const index_t idx) {
      matrix.erase(matrix.begin() + idx);
    }

    void append(const VectorColumn& col) {
      index_t n_columns =  get_n_columns();
      set_n_columns(n_columns + 1);
      set_column(n_columns, col);
    }

    void add(const index_t source, const index_t target) {
      ColumnType& source_col = matrix[source];
      ColumnType& target_col = matrix[target];
      ColumnType& temp_col = temp_column_buffer();

      size_t new_size = source_col.size() + target_col.size();
      if(new_size > temp_col.size()) {
        temp_col.resize(new_size);
      }

      std::vector<index_t>::iterator col_end =
        std::set_symmetric_difference(target_col.begin(), target_col.end(),
                                      source_col.begin(), source_col.end(),
                                      temp_col.begin());

      temp_col.erase(col_end, temp_col.end());
      target_col.swap(temp_col);
    }

    void add(const ColumnType& source_col, const index_t target) {
      ColumnType& target_col = matrix[target];
      ColumnType& temp_col = temp_column_buffer();

      size_t new_size = source_col.size() + target_col.size();
      if(new_size > temp_col.size()) {
        temp_col.resize(new_size);
      }

      std::vector<index_t>::iterator col_end =
        std::set_symmetric_difference(target_col.begin(), target_col.end(),
                                      source_col.begin(), source_col.end(),
                                      temp_col.begin());

      temp_col.erase(col_end, temp_col.end());
      target_col.swap(temp_col);
    }

    index_t get_n_rows(const index_t idx) const {
      VectorColumn col;
      get_column(idx, col);
      return col.size();
    }

    index_t get_max_column_entries() const {
      index_t max_column_entries = -1;
      const index_t n_columns = get_n_columns();
      for(index_t idx = 0; idx < n_columns; idx++)
        max_column_entries = get_n_rows(idx) > max_column_entries ?
          get_n_rows(idx) : max_column_entries;
      return max_column_entries;
    }

    index_t get_max_row_entries() const {
      size_t max_row_entries = 0;
      const index_t n_columns = get_n_columns();
      std::vector< std::vector< index_t > > transposed_matrix(n_columns);
      VectorColumn temp_col;
      for(index_t idx = 0; idx < n_columns; ++idx) {
        get_column(idx, temp_col);
        for(index_t idx = 0; idx < (index_t) temp_col.size(); ++idx)
          transposed_matrix[temp_col[idx]].push_back(idx);
      }

      for(index_t idx = 0; idx < n_columns; ++idx)
        max_row_entries = transposed_matrix[idx].size() > max_row_entries ?
          transposed_matrix[idx].size() : max_row_entries;
      return max_row_entries;
    }

    index_t get_n_entries() const {
      index_t n_nonzero_entries = 0;
      const index_t n_columns = get_n_columns();
      for(index_t idx = 0; idx < n_columns; ++idx)
        n_nonzero_entries += get_n_rows(idx);
      return n_nonzero_entries;
    }

  };

} // namespace stn
