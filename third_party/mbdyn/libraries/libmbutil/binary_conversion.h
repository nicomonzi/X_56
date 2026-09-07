#ifndef ___BINARY_CONVERSION_H_INCLUDED___
#define ___BINARY_CONVERSION_H_INCLUDED___

#include <array>
#include <iostream>
#include <type_traits>
#include <vector>

#include "except.h"

namespace BinaryConversion {
     template <typename T, typename std::enable_if_t<std::is_fundamental<T>::value, bool> = true>
     std::ostream& WriteBinary(std::ostream& os, const T& value) {
          return os.write(reinterpret_cast<const char*>(&value), sizeof(value));
     }
     
     template <typename T, typename std::enable_if_t<std::is_fundamental<T>::value, bool> = true>
     std::istream& ReadBinary(std::istream& is, T& value) {
          return is.read(reinterpret_cast<char*>(&value), sizeof(value));
     }

     std::ostream& WriteBinary(std::ostream& os, const std::string& value);

     std::istream& ReadBinary(std::istream& is, std::string& value);
          
     template <typename T>
     std::istream& ReadBinary(std::istream& is, std::vector<T>& v) {
          size_t n;

          ReadBinary(is, n);

          v.reserve(n);
          v.clear();
          
          T a;

          while (n--) {
               ReadBinary(is, a);
               v.emplace_back(a);
          }

          return is;
     }

     template <typename T>
     std::ostream& WriteBinary(std::ostream& os, const std::vector<T>& v) {
          const size_t n = v.size();
          
          WriteBinary(os, n);

          for (const auto& a: v) {
               WriteBinary(os, a);
          }

          return os;
     }


     template <typename T, size_t N>
     std::istream& ReadBinary(std::istream& is, std::array<T, N>& v) {
          size_t n;

          ReadBinary(is, n);

          if (n != N) {
               throw ErrFile(MBDYN_EXCEPT_ARGS);
          }

          for (size_t i = 0; i < n; ++i) {
               ReadBinary(is, v[i]);
          }

          return is;
     }

     template <typename T, size_t N>
     std::ostream& WriteBinary(std::ostream& os, const std::array<T, N>& v) {
          const size_t n = v.size();
          
          WriteBinary(os, n);

          for (const auto& a: v) {
               WriteBinary(os, a);
          }

          return os;
     }
}
#endif
