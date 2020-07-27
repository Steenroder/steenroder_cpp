/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include "commons.hpp"
#include "sparse_matrix.hpp"
#include "vector_column.hpp"

namespace stn {

  template<typename ColumnType = VectorColumn, typename AttributeType = dimension_t>
  class AttributeMatrix : public SparseMatrix<ColumnType> {
  private:
    using Base = SparseMatrix<ColumnType>;

    std::vector<AttributeType> attributes;

  public:
    AttributeMatrix()
      : Base()
      , attributes()
    {};


    AttributeMatrix(const index_t n_columns_in)
      : attributes(n_columns_in)
    {
      Base::set_n_columns(n_columns_in);
    };

    AttributeMatrix(const SparseMatrix<ColumnType>& sparse_matrix_in,
                    const std::vector<AttributeType>& attributes_in)
      : Base(sparse_matrix_in)
      , attributes(attributes_in)
    {};

    AttributeMatrix(const AttributeMatrix<ColumnType>& other) {
      *this = other;
    }

    AttributeType get_attribute(index_t idx) const {
      return attributes[idx];
    }

    void set_attribute(index_t idx, AttributeType attribute) {
      attributes[idx] = attribute;
    }

    index_t get_n_columns() const {
      return Base::get_n_columns();
    }

    void set_n_columns(index_t n_columns) {
      attributes.resize(n_columns);
      Base::set_n_columns(n_columns);
    }

    // swaps given columns
    void swap(index_t idx_1, index_t idx_2) {
      Base::swap(idx_1, idx_2);
      std::swap(attributes[idx_1], attributes[idx_2]);
    }

    void erase(index_t idx) {
      Base::erase(idx);
      attributes.erase(attributes.begin() + idx);
    }

    AttributeType get_max_attribute() const {
      AttributeType max_attribute = 0;
      for(index_t idx = 0; idx < get_n_columns(); idx++)
        max_attribute = get_attribute(idx) > max_attribute ?
          get_attribute(idx) : max_attribute;
      return max_attribute;
    }

    AttributeType is_in(index_t start, VectorColumn& candidate){
      VectorColumn temp_col;
      for(index_t idx = start; idx < get_n_columns(); ++idx) {
        Base::get_column(idx, temp_col);
        if(temp_col == candidate) {
          return get_attribute(idx);
        }
      }
      return -1;
    }

    bool operator==(const AttributeMatrix& other) const {
      const index_t n_columns = Base::get_n_columns();

      if(n_columns != other.get_n_columns())
        return false;

      VectorColumn temp_col;
      VectorColumn other_temp_col;
      for(index_t idx = 0; idx < n_columns; ++idx) {
        Base::get_column(idx, temp_col);
        other.get_column(idx, other_temp_col);
        if(temp_col != other_temp_col
           || get_attribute(idx) != other.get_attribute(idx))
          return false;
      }
      return true;
    }

    bool operator!=(const AttributeMatrix& other) const {
      return !(*this == other);
    }

    AttributeMatrix<ColumnType>&
    operator=(const AttributeMatrix<ColumnType>& other) {
      const index_t n_columns = other.get_n_columns();
      set_n_columns(n_columns);

      VectorColumn temp_col;
      for(index_t idx = 0; idx <  n_columns; ++idx) {
        set_attribute(idx, other.get_attribute(idx));
        other.get_column(idx, temp_col);
        Base::set_column(idx, temp_col);
      }

      return *this;
    }

    bool load_ascii(std::string filename) {
      // first count number of columns:
      std::string line;
      std::ifstream dummy(filename .c_str());
      if(dummy.fail())
        return false;

      index_t n_columns = 0;
      while(getline(dummy, line)) {
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
        if(line != "" && line[0] != '#')
          n_columns++;

      }
      set_n_columns(n_columns);
      dummy.close();

      std::ifstream input_stream(filename.c_str());
      if(input_stream.fail())
        return false;

      VectorColumn temp_col;
      index_t idx = -1;
      while(getline(input_stream, line)) {
        line.erase(line.find_last_not_of(" \t\n\r\f\v") + 1);
        if(line != "" && line[0] != '#') {
          idx++;
          std::stringstream ss(line);

          index_t att;
          ss >> att;
          set_attribute(idx, (AttributeType) att);

          index_t temp_index;
          temp_col.clear();
          while(ss.good()) {
            ss >> temp_index;
            temp_col.push_back((index_t)temp_index);
          }
          std::sort(temp_col.begin(), temp_col.end());
          Base::set_column(idx, temp_col);
        }
      }

      input_stream.close();
      return true;
    }

    // Format: each line represents a column, first number is attribute, other numbers are the content of the column
    bool save_ascii(const std::string& name, const std::string& output_filename) {
      std::string filename = output_filename + "_" + name + ".dat";

      std::ofstream output_stream(filename.c_str());
      if(output_stream.fail())
        return false;

      const index_t n_columns = Base::get_n_columns();
      VectorColumn temp_col;
      for(index_t col_idx = 0; col_idx < n_columns; ++col_idx) {
        output_stream << (index_t) get_attribute(col_idx);

        Base::get_column(col_idx, temp_col);
        for(index_t row_idx = 0; row_idx < (index_t) temp_col.size(); ++row_idx)
          output_stream << " " << temp_col[row_idx];
        output_stream << std::endl;
      }

      output_stream.close();
      return true;
    }

    // Format: n_columns % att1 % N1 % row1 row2 % ...% rowN1 % att2 % N2 % ...
    bool load_binary(std::string filename) {
      std::ifstream input_stream(filename.c_str(), std::ios_base::binary | std::ios_base::in);
      if(input_stream.fail())
        return false;

      index_t n_columns;
      input_stream.read((char*) &n_columns, sizeof(int64_t));
      set_n_columns((index_t) n_columns);

      VectorColumn temp_col;
      for(index_t idx = 0; idx < n_columns; idx++) {
        index_t att;
        input_stream.read((char*) &att, sizeof(int64_t));
        set_attribute(idx, (AttributeType) att);
        index_t n_rows;
        input_stream.read((char*) &n_rows, sizeof(int64_t));
        temp_col.resize((std::size_t) n_rows);
        for(index_t idx = 0; idx < n_rows; ++idx) {
          index_t row;
          input_stream.read((char*)& row, sizeof(int64_t));
          temp_col[idx] = (index_t) row;
        }
        Base::set_column(idx, temp_col);
      }

      input_stream.close();
      return true;
    }

    // Format: n_columns % att1 % N1 % row1 row2 % ...% rowN1 % att2 % N2 % ...
    bool save_binary(const std::string& name, const std::string& output_filename) {
      std::string filename = output_filename + "_" + name + ".dat";

      std::ofstream output_stream(filename.c_str(),
                                  std::ios_base::binary | std::ios_base::out);
      if(output_stream.fail())
        return false;

      const index_t n_columns = Base::get_n_columns();
      output_stream.write((char*) &n_columns, sizeof(int64_t));
      VectorColumn temp_col;
      for(index_t idx = 0; idx < n_columns; idx++) {
        index_t att = get_attribute(idx);
        output_stream.write((char*) &att, sizeof(int64_t));
        Base::get_column(idx, temp_col);
        index_t n_rows = temp_col.size();
        output_stream.write((char*) &n_rows, sizeof(int64_t));
        for(index_t row_idx = 0; row_idx < (index_t) temp_col.size(); row_idx++) {
          index_t row = temp_col[row_idx];
          output_stream.write((char*) &row, sizeof(int64_t));
        }
      }

      output_stream.close();
      return true;
    }

  };


  template<typename ColumnType>
  class AttributeMatrixVector {
  protected:
    std::vector<SparseMatrix<ColumnType>> sparse_matrices;
    std::vector<std::vector<index_t>> attribute_vectors;

    void break_down(const AttributeMatrix<ColumnType>& attribute_matrix) {
      const index_t n_columns = attribute_matrix.get_n_columns();
      const dimension_t n_elements = sparse_matrices.size();

      dimension_t elem;
      ColumnType temp_col;
      std::vector<index_t> indices(n_elements, 0);
      for(index_t idx_col = 0; idx_col < n_columns; ++idx_col) {
        elem = attribute_matrix.get_attribute(idx_col);
        attribute_matrix.get_column(idx_col, temp_col);

        attribute_vectors[elem][indices[elem]] = idx_col;
        sparse_matrices[elem].set_column(indices[elem], temp_col);

        ++indices[elem];
      }

      for(dimension_t elem = 0; elem < n_elements; ++elem){
        set_n_columns(elem, indices[elem]);
      }
    }

    void break_down_dualize(const AttributeMatrix<ColumnType>& attribute_matrix) {
      const index_t n_columns = attribute_matrix.get_n_columns();
      const dimension_t n_elements = sparse_matrices.size();

      dimension_t elem;
      ColumnType temp_col;
      std::vector<index_t> indices(n_elements, 0);
      for(index_t idx_col = 0; idx_col < n_columns; ++idx_col) {
        elem = n_elements - 1 - attribute_matrix.get_attribute(idx_col);
        attribute_matrix.get_column(idx_col, temp_col);

        attribute_vectors[elem][indices[elem]] = n_columns - 1 - idx_col;
        sparse_matrices[elem].set_column(indices[elem], temp_col);

        ++indices[elem];
      }

      for(dimension_t elem = 0; elem < n_elements; ++elem){
        set_n_columns(elem, indices[elem]);
      }
    }

  public:
    AttributeMatrixVector(const dimension_t n_elements_in,
                          const index_t n_columns_in)
      : sparse_matrices(n_elements_in)
      , attribute_vectors(n_elements_in)
    {
      for(dimension_t elem = 0; elem < n_elements_in; ++elem){
        set_n_columns(elem, n_columns_in);
      }
    }

    AttributeMatrixVector(const AttributeMatrix<ColumnType>& attribute_matrix_in,
                          const dimension_t n_elements_in,
                          const bool dualize_in)
      : AttributeMatrixVector(n_elements_in,
                              attribute_matrix_in.get_n_columns())

    {
      if(dualize_in) {
        break_down_dualize(attribute_matrix_in);
      }
      else {
        break_down(attribute_matrix_in);
      }
    }

    void set_n_columns(const dimension_t elem, const index_t n_columns) {
      sparse_matrices[elem].set_n_columns(n_columns);
      attribute_vectors[elem].resize(n_columns);
    }

    void set_n_columns(const index_t n_columns) {
      for(dimension_t elem = 0; elem < sparse_matrices.size(); ++elem){
        set_n_columns(elem, n_columns);
      }
    }

    filtration_t get_attribute(const dimension_t elem,
                               const index_t idx_col) {
      return attribute_vectors[elem][idx_col];
    }

    void set_attribute(const dimension_t elem,
                       const index_t idx_col,
                       const index_t att) {
      attribute_vectors[elem][idx_col] = att;
    }

    SparseMatrix<ColumnType>& get_matrix(const dimension_t elem) {
      return sparse_matrices[elem];
    }

    void save_ascii(const std::string& name,
                    const std::string& output_filename) {
      for(dimension_t elem = 0; elem < sparse_matrices.size(); ++elem){
        std::string filename =  name + "_dim_" + std::to_string(elem);
        AttributeMatrix<ColumnType, filtration_t> attribute_matrix(sparse_matrices[elem],
                                                                   attribute_vectors[elem]);
        attribute_matrix.save_ascii(filename, output_filename);
      }
    }

    void save_binary(const std::string& name,
                     const std::string& output_filename) {
      for(dimension_t elem = 0; elem < sparse_matrices.size(); ++elem){
        std::string filename =  name + "_dim_" + std::to_string(elem);

        AttributeMatrix<ColumnType, filtration_t> attribute_matrix(sparse_matrices[elem],
                                                                   attribute_vectors[elem]);
        attribute_matrix.save_binary(filename, output_filename);
      }
    }
  };


} // namespace stn
