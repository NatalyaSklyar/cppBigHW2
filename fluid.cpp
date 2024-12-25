#include "fixed.h"
#include "simulation.h"
#include <cassert>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>

using namespace std;

#define FLOAT float
#define DOUBLE double
#define FIXED(N, K) Fixed<N, K>

#ifndef TYPES
#define TYPES FIXED(32, 16), FIXED(16, 8), DOUBLE
#endif

#define S(N, M) Size<N, M>

#ifndef SIZES
#define SIZES S(36, 84), S(10, 10)
#endif

constexpr size_t T = 1'000'000;
constexpr string_view FIELD_FILE = "field.txt";

string get_field() {
  ifstream file(FIELD_FILE);
  string field, line;
  while (getline(file, line)) {
    field += line + '\n';
  }
  return field;
}

template <typename T> struct is_fixed : std::false_type {};

template <size_t N, size_t K> struct is_fixed<Fixed<N, K>> : std::true_type {};

template <typename T> struct is_fast_fixed : std::false_type {};

template <typename T> static bool type_matches(std::string type) {
  if constexpr (std::is_same_v<T, double>) {
    return type == "DOUBLE";
  } else if constexpr (std::is_same_v<T, float>) {
    return type == "float";
  } else if constexpr (is_fixed<T>::value) {
    if (!type.starts_with("FIXED(")) {
      return false;
    }
    int N, K;
    sscanf(type.c_str(), "FIXED(%d,%d)", &N, &K);
    return T::n == N && T::k == K;
  }
  return false;
}

template <typename PreasureT, typename VelocityT, typename VelocityFlowT,
          size_t N, size_t M>
void run_simulation(std::string field) {

  FluidSimulation<PreasureT, VelocityT, VelocityFlowT, N, M> simulation(
      (typename FluidSimulation<PreasureT, VelocityT, VelocityFlowT, N,
                                M>::State(field)));

  for (size_t i = 0; i < T; ++i) {
    simulation.tick(i);
  }
}

template <typename PreasureT, typename VelocityT, typename VelocityFlowT,
          typename T, typename... Ts>
void size_to_template(const std::string &p_type, const std::string &v_type,
                      const std::string &v_flow_type,
                      const std::string &field) {
  size_t N = 0, M = 0;
  for (size_t i = 0; i < field.size(); ++i) {
    if (field[i] == '\n') {
      M = i;
      break;
    }
  }
  N = field.size() / (M + 1);

  if (T::n == N && T::m == M) {
    run_simulation<PreasureT, VelocityT, VelocityFlowT, T::n, T::m>(field);
  } else if constexpr (sizeof...(Ts) > 0) {
    size_to_template<PreasureT, VelocityT, VelocityFlowT, Ts...>(
        p_type, v_type, v_flow_type, field);
  } else {
    throw std::runtime_error("No matching template found size");
  }
}

template <typename PreasureT, typename VelocityT, typename T, typename... Ts>
void v_flow_type_to_template(const std::string &p_type,
                             const std::string &v_type,
                             const std::string &v_flow_type,
                             const std::string &field) {
  if (type_matches<T>(v_flow_type)) {
    size_to_template<PreasureT, VelocityT, T, SIZES>(p_type, v_type,
                                                     v_flow_type, field);
  } else if constexpr (sizeof...(Ts) > 0) {
    v_flow_type_to_template<PreasureT, VelocityT, Ts...>(p_type, v_type,
                                                         v_flow_type, field);
  } else {
    throw std::runtime_error("No matching template found v-flow");
  }
}

template <typename PreasureT, typename T, typename... Ts>
void v_type_to_template(const std::string &p_type, const std::string &v_type,
                        const std::string &v_flow_type,
                        const std::string &field) {
  if (type_matches<T>(v_type)) {
    v_flow_type_to_template<PreasureT, T, TYPES>(p_type, v_type, v_flow_type,
                                                 field);
  } else if constexpr (sizeof...(Ts) > 0) {
    v_type_to_template<PreasureT, Ts...>(p_type, v_type, v_flow_type, field);
  } else {
    throw std::runtime_error("No matching template found velocity");
  }
}

template <typename T, typename... Ts>
void p_type_to_template(const std::string &p_type, const std::string &v_type,
                        const std::string &v_flow_type,
                        const std::string &field) {
  if (type_matches<T>(p_type)) {
    v_type_to_template<T, TYPES>(p_type, v_type, v_flow_type, field);
  } else if constexpr (sizeof...(Ts) > 0) {
    p_type_to_template<Ts...>(p_type, v_type, v_flow_type, field);
  } else {
    throw std::runtime_error("No matching template found p-type");
  }
}

void run_with_params(const std::string &p_type, const std::string &v_type,
                     const std::string &v_flow_type, const std::string &field) {
  p_type_to_template<TYPES>(p_type, v_type, v_flow_type, field);
}

int main(int argc, char *argv[]) {
  std::string p_type = "", v_type = "", v_flow_type = "";

  for (int i = 1; i < argc; ++i) {
    if (std::string(argv[i]) == "--p-type") {
      p_type = argv[++i];
    } else if (std::string(argv[i]) == "--v-type") {
      v_type = argv[++i];
    } else if (std::string(argv[i]) == "--v-flow-type") {
      v_flow_type = argv[++i];
    }
  }
  if (p_type.empty() || v_type.empty() || v_flow_type.empty()) {
    return 1;
  }
  run_with_params(p_type, v_type, v_flow_type, get_field());

  return 0;
}
