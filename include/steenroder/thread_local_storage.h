/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include "commons.hpp"

// should ideally be equal to the cache line size of the CPU
#define PHAT_TLS_SPACING_FACTOR 64

namespace stn {

  // ThreadLocalStorage with some spacing to avoid "false sharing" (see wikipedia)
  template< typename T >
  class thread_local_storage
  {
  public:

    thread_local_storage() : per_thread_storage( omp_get_max_threads() * PHAT_TLS_SPACING_FACTOR ) {};

    T& operator()() {
      return per_thread_storage[ omp_get_thread_num() * PHAT_TLS_SPACING_FACTOR ];
    }

    const T& operator()() const {
      return per_thread_storage[ omp_get_thread_num() * PHAT_TLS_SPACING_FACTOR ];
    }

    T& operator[]( int tid ) {
      return per_thread_storage[ tid * PHAT_TLS_SPACING_FACTOR ];
    }

    const T& operator[]( int tid ) const {
      return per_thread_storage[ tid * PHAT_TLS_SPACING_FACTOR ];
    }

  protected:
    std::vector< T > per_thread_storage;
  };

} // namespace stn
