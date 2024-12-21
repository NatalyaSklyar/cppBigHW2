#ifndef BIGHW2__FIXED_H_
#define BIGHW2__FIXED_H_

#include <cassert>
#include <cmath>
#include <cstdint>
#include <iostream>

template <typename T, typename T2, size_t N, size_t K, typename R>
struct FixedT {
  static constexpr size_t n = N;
  static constexpr size_t k = K;

  constexpr FixedT(int v) : v(v << K) {}
  constexpr FixedT(long long v) : v(v << K) {}
  constexpr FixedT(float f) : v(f * (1 << K)) {}
  constexpr FixedT(double f) : v(f * (1 << K)) {}
  constexpr FixedT() : v(0) {}
  template <template <size_t, size_t> class FixedT2, size_t N2, size_t K2>
  explicit constexpr FixedT(const FixedT2<N2, K2> &x) {
    if constexpr (K2 < K) {
      v = (T)x.v << (K - K2);
    } else {
      v = (T)x.v >> (K2 - K);
    }
  }

  static constexpr R from_raw(T x) {
    R ret;
    ret.v = x;
    return ret;
  }

  T v;

  auto operator<=>(const FixedT &) const = default;
  bool operator==(const FixedT &) const = default;

  operator double() const { return v / (double)(1 << K); }
  operator float() const { return v / (float)(1 << K); }

  friend R operator+(R a, R b) { return R::from_raw(a.v + b.v); }

  friend R operator-(R a, R b) { return R::from_raw(a.v - b.v); }

  friend R operator*(R a, R b) { return R::from_raw(((T2)a.v * b.v) >> K); }

  friend R operator/(R a, R b) { return R::from_raw(((T2)a.v << K) / b.v); }

  friend R &operator+=(R &a, R b) { return a = a + b; }

  friend R &operator-=(R &a, R b) { return a = a - b; }

  friend R &operator*=(R &a, R b) { return a = a * b; }

  friend R &operator/=(R &a, R b) { return a = a / b; }

  friend R operator-(R x) { return R::from_raw(-x.v); }

  friend R abs(R x) {
    if (x.v < 0) {
      x.v = -x.v;
    }
    return x;
  }

  friend std::ostream &operator<<(std::ostream &out, R x) {
    return out << x.v / (double)(1 << K);
  }
};

template <size_t N, size_t K> struct Fixed;

template <size_t K>
struct Fixed<64, K> : public FixedT<int64_t, __int128, 64, K, Fixed<64, K>> {
  using FixedT<int64_t, __int128, 64, K, Fixed<64, K>>::FixedT;
};

template <size_t K>
struct Fixed<32, K> : public FixedT<int32_t, int64_t, 32, K, Fixed<32, K>> {
  using FixedT<int32_t, int64_t, 32, K, Fixed<32, K>>::FixedT;
};

template <size_t K>
struct Fixed<16, K> : public FixedT<int16_t, int32_t, 16, K, Fixed<16, K>> {
  using FixedT<int16_t, int32_t, 16, K, Fixed<16, K>>::FixedT;
};

template <size_t K>
struct Fixed<8, K> : public FixedT<int8_t, int16_t, 8, K, Fixed<8, K>> {
  using FixedT<int8_t, int16_t, 8, K, Fixed<8, K>>::FixedT;
};

#endif // BIGHW2__FIXED_H_
