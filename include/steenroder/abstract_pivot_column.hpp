/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include "sparse_matrix.hpp"

namespace stn {

    // Note: We could even make the rep generic in the underlying Const representation
    //       But I cannot imagine that anything else than vector<vector<index_t>> would
    //       make sense
    template<typename PivotColumn>
    class AbstractPivotColumn : public SparseMatrix {

    protected:
        using Base = SparseMatrix;

        // For parallization purposes, it could be more than one full column
        mutable thread_local_storage< PivotColumn > pivot_cols;
        mutable thread_local_storage< index_t > idx_of_pivot_cols;

        PivotColumn& get_pivot_col() const {
            return pivot_cols();
        }

        bool is_pivot_col( index_t idx ) const {
            return idx_of_pivot_cols() == idx;
        }

        void release_pivot_col() {
            index_t idx = idx_of_pivot_cols();
            if( idx != -1 ) {
                this->matrix[ idx ].clear();
                pivot_cols().get_col_and_clear( this->matrix[ idx ] );
            }
            idx_of_pivot_cols() = -1;
        }

        void make_pivot_col( index_t idx ) {
            release_pivot_col();
            idx_of_pivot_cols() = idx;
            get_pivot_col().add_col( matrix[ idx ] );
        }

    public:
        void _set_num_cols( index_t nr_of_cols ) {
            #pragma omp parallel for
            for( int tid = 0; tid < omp_get_num_threads(); tid++ ) {
                pivot_cols[ tid ].init( nr_of_cols );
                idx_of_pivot_cols[ tid ] = -1;
            }
            Base::_set_num_cols( nr_of_cols );
        }

        void _add( index_t source, index_t target ) {
            if( !is_pivot_col( target ) )
                make_pivot_col( target );
            get_pivot_col().add_col( matrix[source] );
        }

        void _add( column& source_col, index_t target ) {
            if( !is_pivot_col( target ) )
                make_pivot_col( target );
            get_pivot_col().add_col( source_col );
        }


        void _sync() {
            #pragma omp parallel for
            for( int tid = 0; tid < omp_get_num_threads(); tid++ )
                release_pivot_col();
        }

        void _get_col( index_t idx, column& col  ) const {
          is_pivot_col( idx ) ? get_pivot_col().get_col( col ) :
            Base::_get_col( idx, col );
        }

        bool _is_empty( index_t idx ) const {
          return is_pivot_col( idx ) ? get_pivot_col().is_empty() :
            Base::_is_empty( idx );
        }

        index_t _get_max_index( index_t idx ) const {
          return is_pivot_col( idx ) ? get_pivot_col().get_max_index()
            : Base::_get_max_index( idx );
        }

        void _clear( index_t idx ) {
          is_pivot_col( idx ) ? get_pivot_col().clear() :
            Base::_clear( idx );
        }

        void _set_col( index_t idx, const column& col  ) {
          is_pivot_col( idx ) ? get_pivot_col().set_col( col ) :
            Base::_set_col( idx, col );
        }

        void _remove_max( index_t idx ) {
          is_pivot_col( idx ) ? get_pivot_col().remove_max() :
            Base::_remove_max( idx );
        }

        void finalize( index_t idx ) { Base::_finalize( idx ); }
    };

} // namespace stn
