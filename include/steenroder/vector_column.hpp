/*  Author: Guillaume Tauzin
    License: GPLv3
*/

#pragma once

#include <iostream>
#include <vector>

#include "commons.hpp"

namespace stn {

  class VectorColumn :
    public std::vector<index_t> {
  private:
    using Base = std::vector<index_t>;
    using Base::Base;

  public:
    index_t get_max() {
      if(size()) {
        return Base::operator[](size() - 1);
      }
      else {
        return -1;
      }
    }

    VectorColumn& operator+=(const VectorColumn& right) {
      VectorColumn temp;

      size_t new_size = right.size() + size();
      if(new_size > temp.size()) {
        temp.resize(new_size);
      }

      std::vector<index_t>::iterator col_end =
        std::set_symmetric_difference(begin(), end(),
                                      right.begin(), right.end(),
                                      temp.begin());

      temp.erase(col_end, temp.end());
      swap(temp);
      return *this;
    }

    VectorColumn& operator-=(const VectorColumn& right) {
      VectorColumn temp;

      size_t new_size = right.size() + size();
      if(new_size > temp.size()) {
        temp.resize(new_size);
      }

      std::vector<index_t>::iterator col_end =
        std::set_difference(begin(), end(),
                            right.begin(), right.end(),
                            temp.begin());

      temp.erase(col_end, temp.end());
      swap(temp);
      return *this;
    }

    VectorColumn& operator |=(const VectorColumn& right) {
      VectorColumn temp;

      size_t new_size = right.size() + size();
      if(new_size > temp.size()) {
        temp.resize(new_size);
      }

      std::vector<index_t>::iterator col_end =
        std::set_union(begin(), end(),
                       right.begin(), right.end(),
                       temp.begin());

      temp.erase(col_end, temp.end());
      swap(temp);
      return *this;
    }

  };

  std::ostream& operator<<(std::ostream& os, const VectorColumn& col) {
    os << "[";
    for (int i = 0; i < col.size(); ++i) {
      os << col[i];
      if (i != col.size() - 1)
        os << ", ";
    }
    os << "]\n";
    return os;
  }


  VectorColumn operator+(const VectorColumn& right, const VectorColumn& left) {
    VectorColumn res = left;
    res += right;
    return res;
  }

  VectorColumn operator-(const VectorColumn& right, const VectorColumn& left) {
    VectorColumn res = left;
    res -= right;
    return res;
  }

  VectorColumn operator|(const VectorColumn& right, const VectorColumn& left) {
    VectorColumn res = left;
    res |= right;
    return res;
  }

} // namespace stn
