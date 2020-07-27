#pragma once

// STL includes
#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <set>
#include <list>
#include <map>
#include <algorithm>
#include <queue>
#include <cassert>
#include <sstream>
#include <algorithm>
#include <iomanip>
#include <cmath>
#include <cstdlib>
#include <iterator>
#include <stdint.h>


#ifdef USE_SCOREP
  #include <scorep/SCOREP_User.h>

  #define STN_INSTRUMENT_ON(name, colorID) \
    { SCOREP_USER_REGION(name, SCOREP_USER_REGION_TYPE_FUNCTION) }

#else  // USE_SCOREP
  #define STN_INSTRUMENT_ON(name, colorID)
#endif  // USE_SCOREP

#define STN_INSTRUMENT_OFF(name, colorID)


// basic types. index_t can be changed to int32_t to save memory on small instances
namespace stn {
    typedef int64_t index_t;
    typedef int64_t filtration_t;
    typedef int8_t dimension_t;
} // namespace stn;

// OpenMP (proxy) functions
#if defined _OPENMP
    #include <omp.h>
#else
    #define omp_get_thread_num() 0
    #define omp_get_max_threads() 1
    #define omp_get_num_threads() 1
	void omp_set_num_threads( int ) {};
    #include <time.h>
    #define omp_get_wtime() (float)clock() / (float)CLOCKS_PER_SEC
#endif

#include <steenroder/thread_local_storage.h>
