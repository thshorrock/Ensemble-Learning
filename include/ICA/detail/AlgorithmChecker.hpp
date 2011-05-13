#pragma once

namespace std{
  template<class InputIterator, class T>
  inline
  T
  accumulate_checked(InputIterator first, InputIterator last, T init )
  {
    return std::accumulate(first,last, init);
  }


  template<class InputIterator>
  inline
  int
  accumulate_checked(InputIterator first, InputIterator last, int init );

}
