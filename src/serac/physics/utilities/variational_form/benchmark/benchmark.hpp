#pragma once

#include <sstream>

namespace serac {
  namespace profiling {

    template<typename ...T>
    std::string concat(T... args) {
      std::stringstream ss;
      (ss << ... << args);
      return ss.str();
    }
  
  }
}
