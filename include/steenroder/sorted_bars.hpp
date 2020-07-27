/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include "commons.hpp"
#include "sorted_matrix.hpp"
#include "vector_column.hpp"

namespace stn {

  template<typename ColumnType = VectorColumn>
  class ViewInfiniteBars
    : public ViewMatrix<ColumnType> {
  private:
    using Base = ViewMatrix<ColumnType>;

  protected:
    const index_t n_cells;
    std::vector<index_t> births;
    using Base::n_columns_per_dimension;
    using Base::start_dimension;

  public:
    ViewInfiniteBars(const ViewMatrix<ColumnType>& boundaryMatrix_in)
      : Base(boundaryMatrix_in)
      , n_cells(boundaryMatrix_in.get_n_columns())
      , births(boundaryMatrix_in.get_n_columns(), -1)
    {}

    ViewInfiniteBars(const index_t n_cells_in,
                     const dimension_t n_dimensions_in)
      : Base(n_cells_in, n_dimensions_in)
      , n_cells(n_cells_in)
      , births(n_cells_in, -1)
    {
      for(index_t idx_col = 0; idx_col < n_cells_in; ++idx_col) {
        VectorColumn col(1, idx_col);
        Base::set_column(idx_col, col);
      }
    }

    using Base::set_n_columns;
    using Base::get_n_columns;

    index_t get_n_bars() const {
      return std::accumulate(Base::n_columns_per_dimension.begin(),
                             Base::n_columns_per_dimension.end(), 0);
    }

    index_t get_birth(index_t idx) const {
      return births[idx];
    }

    void set_birth(index_t idx, index_t birth) {
      births[idx] = birth;
    }

    index_t get_death(index_t idx) const {
      return -1;
    }

    // Homology needs to affect the view
    void dualize() {
      for(dimension_t dim = 0; dim < Base::n_dimensions; ++dim) {
        index_t start = Base::get_start_dimension(dim);
        index_t end = start + Base::get_n_columns_per_dimension(dim);
        for(index_t idx_view = start; idx_view < end; ++idx_view) {
          index_t idx_col = Base::get_view(idx_view);
          set_birth(idx_col, n_cells - 1 - get_birth(idx_col));
        }
      }
    }

  };

  template<typename ColumnType = VectorColumn>
  class ViewFiniteBars
    : public ViewInfiniteBars<ColumnType> {
  private:
    using Base = ViewInfiniteBars<ColumnType>;

  protected:
    std::vector<index_t> deaths;

  public:
    ViewFiniteBars(const ViewMatrix<ColumnType>& boundaryMatrix_in)
      : Base(boundaryMatrix_in)
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
      for(dimension_t dim = 0; dim < Base::n_dimensions; ++dim) {
        index_t start = Base::get_start_dimension(dim);
        index_t end = start + Base::get_n_columns_per_dimension(dim);
        for(index_t idx_view = start; idx_view < end; ++idx_view) {
          index_t idx_col = Base::get_view(idx_view);

          index_t birth = get_birth(idx_col);
          index_t death = get_death(idx_col);

          set_birth(idx_col, Base::n_cells - 1 - death);
          set_death(idx_col, Base::n_cells - 1 - birth);
        }
      }
    }

  };


  // Saves the persistence pairs to given file in binary format
  // Format: nr_pairs % newline % dim1 %birth1 % death1 % newline % dim2 % birth2 % death2 % newline ...
  template<typename ColumnType>
  bool save_pairs_ascii(const std::string& filename,
                        const ViewFiniteBars<ColumnType>& finite_bars,
                        const ViewInfiniteBars<ColumnType>& infinite_bars) {
    std::ofstream output_stream(filename.c_str());
    if(output_stream.fail())
      return false;

    dimension_t n_dimensions = finite_bars.get_n_dimensions();
    std::vector<index_t> n_finite_bars_per_dimension =
      finite_bars.get_n_columns_per_dimension();

    for(dimension_t dim = 0; dim < n_dimensions; ++dim) {
      output_stream << "# dim " << (index_t) dim << std::endl;

      index_t start_finite = finite_bars.get_start_dimension(dim);
      index_t end_finite = start_finite + finite_bars.get_n_columns_per_dimension(dim);
      index_t start_infinite = infinite_bars.get_start_dimension(dim);
      index_t end_infinite = start_infinite + infinite_bars.get_n_columns_per_dimension(dim);

      output_stream << end_finite - start_finite + end_infinite - start_infinite
                    << std::endl;

      for(index_t idx_view = start_finite; idx_view < end_finite; ++idx_view) {
        index_t idx_col = finite_bars.get_view(idx_view);

        output_stream << finite_bars.get_birth(idx_col) << " "
                      << finite_bars.get_death(idx_col) << std::endl;
      }

      for(index_t idx_view = start_infinite; idx_view < end_infinite; ++idx_view) {
        index_t idx_col = infinite_bars.get_view(idx_view);

        output_stream << infinite_bars.get_birth(idx_col) << " "
                      << -1 << std::endl;
      }
    }

    output_stream.close();
    return true;
  }

  // Saves the persistence pairs to given file in binary format
  // Format: nr_pairs % birth1 % death1 % birth2 % death2 ...
  template<typename ColumnType>
  bool save_pairs_binary(const std::string& filename,
                         const ViewFiniteBars<ColumnType>& finite_bars,
                         const ViewInfiniteBars<ColumnType>& infinite_bars) {
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


  template<typename ColumnType = VectorColumn>
  class Bars : public ViewMatrix<ColumnType> {
  private:
    using Base = ViewMatrix<ColumnType>;

  protected:
    const index_t n_cells;
    std::vector<index_t> births;
    std::vector<index_t> deaths;
    using Base::n_columns_per_dimension;
    using Base::start_dimension;

  public:
    Bars(const ViewMatrix<ColumnType>& boundaryMatrix_in)
      : Base(boundaryMatrix_in)
      , n_cells(boundaryMatrix_in.get_n_columns())
      , births(boundaryMatrix_in.get_n_columns(), -1)
      , deaths(boundaryMatrix_in.get_n_columns(), -1)
    {}

    Bars(const index_t n_cells_in)
      : Base(n_cells_in, 1)
      , n_cells(n_cells_in)
      , births(n_cells_in, -1)
      , deaths(n_cells_in, -1)
    {
      for(index_t idx_col = 0; idx_col < n_cells_in; ++idx_col) {
        VectorColumn col(1, idx_col);
        Base::set_column(idx_col, col);
      }
    }

    using Base::set_n_columns;
    using Base::get_n_columns;

    index_t get_n_bars() const {
      return std::accumulate(Base::n_columns_per_dimension.begin(),
                             Base::n_columns_per_dimension.end(), 0);
    }

    index_t get_birth(index_t idx) const {
      return births[idx];
    }

    void set_birth(index_t idx, index_t birth) {
      births[idx] = birth;
    }

    index_t get_death(index_t idx) const {
      return deaths[idx];
    }

    void set_death(index_t idx, index_t death) {
      deaths[idx] = death;
    }

    // Homology needs to affect the view
    void dualize() {
      for(dimension_t dim = 0; dim < Base::n_dimensions; ++dim) {
        index_t start = Base::get_start_dimension(dim);
        index_t end = start + Base::get_n_columns_per_dimension(dim);
        for(index_t idx_view = start; idx_view < end; ++idx_view) {
          index_t idx_col = Base::get_view(idx_view);
          std::cout << idx_view << ", " <<  idx_col << std::endl;
          if(get_death(idx_col) == -1) {
            set_birth(idx_col, n_cells - 1 - get_birth(idx_col));
          }
          else {
            index_t birth = get_birth(idx_col);
            index_t death = get_death(idx_col);

            set_birth(idx_col, n_cells - 1 - death);
            set_death(idx_col, n_cells - 1 - birth);
          }
        }
      }
    }
  };

  // Saves the persistence pairs to given file in binary format
  // Format: nr_pairs % newline % dim1 %birth1 % death1 % newline % dim2 % birth2 % death2 % newline ...
  template<typename ColumnType>
  bool save_pairs_ascii(const std::string& filename,
                        const Bars<ColumnType>& bars) {
    std::ofstream output_stream(filename.c_str());
    if(output_stream.fail())
      return false;

    dimension_t n_dimensions = bars.get_n_dimensions();

    for(dimension_t dim = 0; dim < n_dimensions; ++dim) {
      output_stream << "# dim " << (index_t) dim << std::endl;
      index_t start = bars.get_start_dimension(dim);
      index_t end = start + bars.get_n_columns_per_dimension(dim);

      output_stream << end - start << std::endl;

      for(index_t idx_view = start; idx_view < end; ++idx_view) {
        index_t idx_col = bars.get_view(idx_view);

        output_stream << bars.get_birth(idx_col) << " "
                      << bars.get_death(idx_col) << std::endl;
      }
    }

    output_stream.close();
    return true;
  }

} // namespace stn
