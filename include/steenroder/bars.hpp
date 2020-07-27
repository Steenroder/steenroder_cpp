/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include "commons.hpp"
#include "boundary_matrix.hpp"
#include "vector_column.hpp"

namespace stn {

  template<typename ColumnType = VectorColumn>
  class InfiniteBars
    : public BoundaryMatrix<ColumnType> {
  private:
    using Base = BoundaryMatrix<ColumnType>;

  protected:
    const dimension_t n_dimensions;
    const index_t n_cells;
    std::vector<index_t> births;

  public:
    InfiniteBars(const BoundaryMatrix<ColumnType>& boundaryMatrix_in,
                 const dimension_t n_dimensions_in)
      : Base(boundaryMatrix_in)
      , n_dimensions(n_dimensions_in)
      , n_cells(boundaryMatrix_in.get_n_columns())
      , births(boundaryMatrix_in.get_n_columns(), -1)
    {};

    InfiniteBars(const dimension_t n_dimensions_in, const index_t n_cells_in)
      : Base(n_cells_in)
      , n_dimensions(n_dimensions_in)
      , n_cells(n_cells_in)
      , births(n_cells_in, -1)
    {
      for(index_t index = 0; index < n_cells_in; ++index) {
        VectorColumn col = VectorColumn(1, index);
        this->set_column(index, col);
      }
    };

    using Base::set_n_columns;
    using Base::get_n_columns;

    index_t get_birth(index_t idx) const {
      return births[idx];
    }


    void set_birth(index_t idx, index_t birth) {
      births[idx] = birth;
    }

    index_t get_death(index_t idx) const {
      return -1;
    }

    void dualize() {
      for(index_t idx_col = 0; idx_col < get_n_columns(); ++idx_col) {
        dimension_t dim = Base::get_dimension(idx_col);
        index_t birth = get_birth(idx_col);
        Base::set_dimension(idx_col, n_dimensions - 1 - dim);
        set_birth(idx_col, n_cells - 1 - birth);
      }
    }

  };

  template<typename ColumnType = VectorColumn>
  class FiniteBars
    : public InfiniteBars<ColumnType> {
  private:
    using Base = InfiniteBars<ColumnType>;

  protected:
    std::vector<index_t> deaths;

  public:
    FiniteBars(const BoundaryMatrix<ColumnType>& boundaryMatrix_in,
               const dimension_t n_dimensions_in)
      : Base(boundaryMatrix_in, n_dimensions_in)
      , deaths(boundaryMatrix_in.get_n_columns(), -1)
    {};

    using Base::set_n_columns;
    using Base::get_n_columns;
    using Base::set_birth;
    using Base::get_birth;

    index_t get_death(index_t idx) const {
      return deaths[idx];
    }

    void set_death(index_t idx, index_t death) {
      deaths[idx] = death;
    }

    void dualize() {
      for(index_t idx_col = 0; idx_col < get_n_columns(); ++idx_col) {
        dimension_t dim = Base::get_dimension(idx_col);
        index_t birth = get_birth(idx_col);
        index_t death = get_death(idx_col);

        Base::set_dimension(idx_col, Base::n_dimensions - 1 - dim - 1);
        set_birth(idx_col, Base::n_cells - 1 - death);
        set_death(idx_col, Base::n_cells - 1 - birth);
      }
    }

  };


  // Saves the persistence pairs to given file in binary format
  // Format: nr_pairs % newline % dim1 %birth1 % death1 % newline % dim2 % birth2 % death2 % newline ...
  template<typename ColumnType>
  bool save_pairs_ascii(const std::string& filename,
                        const FiniteBars<ColumnType>& finite_bars,
                        const InfiniteBars<ColumnType>& infinite_bars) {
    std::ofstream output_stream(filename.c_str());
    if(output_stream.fail())
      return false;

    index_t n_finite_pairs = finite_bars.get_n_columns();
    index_t n_infinite_pairs = infinite_bars.get_n_columns();
    output_stream << n_finite_pairs + n_infinite_pairs << std::endl;

    for(index_t index = 0; index < n_finite_pairs; ++index) {
      output_stream << (int64_t) finite_bars.get_dimension(index) << " "
                    << finite_bars.get_birth(index) << " "
                    << finite_bars.get_death(index) << std::endl;
    }

    for(index_t index = 0; index < n_infinite_pairs; ++index) {
      output_stream << (int64_t) infinite_bars.get_dimension(index) << " "
                    << infinite_bars.get_birth(index) << " -1"
                    << std::endl;
    }

    output_stream.close();
    return true;
  }

  // Saves the persistence pairs to given file in binary format
  // Format: nr_pairs % birth1 % death1 % birth2 % death2 ...
  template<typename ColumnType>
  bool save_pairs_binary(const std::string& filename,
                         const FiniteBars<ColumnType>& finite_bars,
                         const InfiniteBars<ColumnType>& infinite_bars) {
    std::ofstream output_stream(filename.c_str(), std::ios_base::binary
                                | std::ios_base::out );
    if( output_stream.fail() )
    return false;

    index_t n_finite_pairs = finite_bars.get_n_columns();
    index_t n_infinite_pairs = infinite_bars.get_n_columns();
    index_t n_pairs = n_finite_pairs + n_infinite_pairs;
    output_stream.write((char*) &n_pairs, sizeof(index_t));

    for(index_t index = 0; index < n_finite_pairs; ++index) {
      dimension_t dim = finite_bars.get_dimension(index);
      output_stream.write((char*) &dim, sizeof(dimension_t));
      index_t birth = finite_bars.get_birth(index);
      output_stream.write((char*) &birth, sizeof(index_t));
      index_t death = finite_bars.get_death(index);
      output_stream.write((char*) &death, sizeof(index_t));
    }

    for(index_t index = 0; index < n_infinite_pairs; ++index) {
      dimension_t dim = infinite_bars.get_dimension(index);
      output_stream.write((char*) &dim, sizeof(dimension_t));
      index_t birth = infinite_bars.get_birth(index);
      output_stream.write((char*) &birth, sizeof(index_t));
      index_t death = -1;
      output_stream.write((char*) &death, sizeof(int64_t));
    }

    output_stream.close();
    return true;
  }

} // namespace stn
