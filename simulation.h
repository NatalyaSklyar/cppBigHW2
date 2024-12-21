#ifndef BIGHW2__SIMULATION_H_
#define BIGHW2__SIMULATION_H_

#include <algorithm>
#include <cstring>
#include <array>
#include <cassert>
#include <cmath>
#include <iostream>
#include <random>
#include <tuple>
#include <utility>

using namespace std;

constexpr std::array<pair<int, int>, 4> deltas{
    {{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};

template <size_t N, size_t M> struct Size {
  static constexpr size_t n = N;
  static constexpr size_t m = M;
};

mt19937 rnd(1337);

template <typename PreasureT, typename VelocityT, typename VelocityFlowT,
          size_t N, size_t M>
class FluidSimulation {
private:
  template <typename T> struct VectorField {
    array<T, deltas.size()> v[N][M];
    T &add(int x, int y, int dx, int dy, T dv) {
      return get(x, y, dx, dy) += dv;
    }

    T &get(int x, int y, int dx, int dy) {
      size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
      assert(i < deltas.size());
      return v[x][y][i];
    }
  };

private:
  struct ParticleParams {
    ParticleParams(FluidSimulation &sim) : sim(sim) {}
    char type;
    PreasureT cur_p;
    array<VelocityT, deltas.size()> v;
    FluidSimulation &sim;

    void swap_with(int x, int y) {
      swap(sim.field[x][y], type);
      swap(sim.p[x][y], cur_p);
      swap(sim.velocity.v[x][y], v);
    }
  };

public:
  struct State {
    PreasureT rho[256];
    PreasureT p[N][M]{}, old_p[N][M];
    VelocityT g;
    char field[N][M + 1];
    VectorField<VelocityT> velocity{};
    VectorField<VelocityFlowT> velocity_flow{};
    int last_use[N][M]{};
    int UT = 0;
    int dirs[N][M]{};

    State(std::string field) {
      rho[' '] = 0.01;
      rho['.'] = 1000;
      g = 0.1;

      for (size_t i = 0; i < N; ++i) {
        for (size_t j = 0; j < M; ++j) {
          this->field[i][j] = field[i * (M + 1) + j];
        }
      }

      for (size_t x = 0; x < N; ++x) {
        for (size_t y = 0; y < M; ++y) {
          if (this->field[x][y] == '#')
            continue;
          for (auto [dx, dy] : deltas) {
            dirs[x][y] += (this->field[x + dx][y + dy] != '#');
          }
        }
      }
    }
  };

public:
  FluidSimulation(const State &state) {
    memcpy(rho, state.rho, sizeof(rho));
    memcpy(p, state.p, sizeof(p));
    memcpy(old_p, state.old_p, sizeof(old_p));
    memcpy(field, state.field, sizeof(field));
    memcpy(dirs, state.dirs, sizeof(dirs));
    memcpy(last_use, state.last_use, sizeof(last_use));
    memcpy(velocity.v, state.velocity.v, sizeof(velocity.v));
    memcpy(velocity_flow.v, state.velocity_flow.v, sizeof(velocity_flow.v));

    g = state.g;
    UT = state.UT;
  }

  State export_state() const {
    State state = {};
    memcpy(state.rho, rho, sizeof(rho));
    memcpy(state.p, p, sizeof(p));
    memcpy(state.old_p, old_p, sizeof(old_p));
    memcpy(state.field, field, sizeof(field));
    memcpy(state.dirs, dirs, sizeof(dirs));
    memcpy(state.last_use, last_use, sizeof(last_use));
    memcpy(state.velocity.v, velocity.v, sizeof(velocity.v));
    memcpy(state.velocity_flow.v, velocity_flow.v, sizeof(velocity_flow.v));

    state.g = g;
    state.UT = UT;

    return state;
  }

private:
  char field[N][M + 1];

  PreasureT rho[256];
  PreasureT p[N][M]{}, old_p[N][M];
  VelocityT g;
  VectorField<VelocityT> velocity{};
  VectorField<VelocityFlowT> velocity_flow{};
  int last_use[N][M]{};
  int UT = 0;
  int dirs[N][M]{};

  tuple<VelocityT, bool, pair<int, int>> propagate_flow(int x, int y,
                                                        VelocityT lim) {
    last_use[x][y] = UT - 1;
    VelocityT ret = 0;
    for (auto [dx, dy] : deltas) {
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
        auto cap = velocity.get(x, y, dx, dy);
        auto flow = velocity_flow.get(x, y, dx, dy);
        if (flow == VelocityFlowT(cap)) {
          continue;
        }
        auto vp = min(lim, cap - VelocityT(flow));
        if (last_use[nx][ny] == UT - 1) {
          velocity_flow.add(x, y, dx, dy, VelocityFlowT(vp));
          last_use[x][y] = UT;
          return {vp, 1, {nx, ny}};
        }
        auto [t, prop, end] = propagate_flow(nx, ny, vp);
        ret += t;
        if (prop) {
          velocity_flow.add(x, y, dx, dy, VelocityFlowT(t));
          last_use[x][y] = UT;
          return {t, prop && end != pair(x, y), end};
        }
      }
    }
    last_use[x][y] = UT;
    return {ret, 0, {0, 0}};
  }

  VelocityT random01() {
    return VelocityT((double)((rnd() & ((1 << 16) - 1))) / (1 << 16));
  }

  void propagate_stop(int x, int y, bool force = false) {
    if (!force) {
      bool stop = true;
      for (auto [dx, dy] : deltas) {
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 &&
            velocity.get(x, y, dx, dy) > VelocityT(0)) {
          stop = false;
          break;
        }
      }
      if (!stop) {
        return;
      }
    }
    last_use[x][y] = UT;
    for (auto [dx, dy] : deltas) {
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] == '#' || last_use[nx][ny] == UT ||
          velocity.get(x, y, dx, dy) > VelocityT(0)) {
        continue;
      }
      propagate_stop(nx, ny);
    }
  }

  VelocityT move_prob(int x, int y) {
    VelocityT sum = 0;
    for (size_t i = 0; i < deltas.size(); ++i) {
      auto [dx, dy] = deltas[i];
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
        continue;
      }
      auto v = velocity.get(x, y, dx, dy);
      if (v < VelocityT(0)) {
        continue;
      }
      sum += v;
    }
    return sum;
  }

  bool propagate_move(int x, int y, bool is_first) {
    last_use[x][y] = UT - is_first;
    bool ret = false;
    int nx = -1, ny = -1;
    do {
      std::array<VelocityT, deltas.size()> tres;
      VelocityT sum = 0;
      for (size_t i = 0; i < deltas.size(); ++i) {
        auto [dx, dy] = deltas[i];
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
          tres[i] = sum;
          continue;
        }
        auto v = velocity.get(x, y, dx, dy);
        if (v < VelocityT(0)) {
          tres[i] = sum;
          continue;
        }
        sum += v;
        tres[i] = sum;
      }

      if (sum == VelocityT(0)) {
        break;
      }

      VelocityT p = random01() * sum;
      size_t d = std::ranges::upper_bound(tres, p) - tres.begin();

      auto [dx, dy] = deltas[d];
      nx = x + dx;
      ny = y + dy;
      assert(velocity.get(x, y, dx, dy) > VelocityT(0) &&
             field[nx][ny] != '#' && last_use[nx][ny] < UT);

      ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
    } while (!ret);
    last_use[x][y] = UT;
    for (size_t i = 0; i < deltas.size(); ++i) {
      auto [dx, dy] = deltas[i];
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 &&
          velocity.get(x, y, dx, dy) < VelocityT(0)) {
        propagate_stop(nx, ny);
      }
    }
    if (ret) {
      if (!is_first) {
        ParticleParams pp(*this);
        pp.swap_with(x, y);
        pp.swap_with(nx, ny);
        pp.swap_with(x, y);
      }
    }
    return ret;
  }

public:
  void tick(size_t i) {
    PreasureT total_delta_p = 0;
    // Apply external forces
    for (size_t x = 0; x < N; ++x) {
      for (size_t y = 0; y < M; ++y) {
        if (field[x][y] == '#')
          continue;
        if (field[x + 1][y] != '#')
          velocity.add(x, y, 1, 0, g);
      }
    }

    // Apply forces from p
    memcpy(old_p, p, sizeof(p));
    for (size_t x = 0; x < N; ++x) {
      for (size_t y = 0; y < M; ++y) {
        if (field[x][y] == '#')
          continue;
        for (auto [dx, dy] : deltas) {
          int nx = x + dx, ny = y + dy;
          if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
            auto delta_p = old_p[x][y] - old_p[nx][ny];
            auto force = delta_p;
            auto &contr = velocity.get(nx, ny, -dx, -dy);
            if (contr * VelocityT(rho[(int)field[nx][ny]]) >=
                VelocityT(force)) {
              contr -= VelocityT(force / rho[(int)field[nx][ny]]);
              continue;
            }
            force -= PreasureT(contr) * rho[(int)field[nx][ny]];
            contr = 0;
            velocity.add(x, y, dx, dy,
                         VelocityT(force / rho[(int)field[x][y]]));
            p[x][y] -= force / PreasureT(dirs[x][y]);
            total_delta_p -= force / PreasureT(dirs[x][y]);
          }
        }
      }
    }

    // Make flow from velocities
    velocity_flow = {};
    bool prop = false;
    do {
      UT += 2;
      prop = 0;
      for (size_t x = 0; x < N; ++x) {
        for (size_t y = 0; y < M; ++y) {
          if (field[x][y] != '#' && last_use[x][y] != UT) {
            auto [t, local_prop, _] = propagate_flow(x, y, 1);
            if (t > VelocityT(0)) {
              prop = 1;
            }
          }
        }
      }
    } while (prop);

    // Recalculate p with kinetic energy
    for (size_t x = 0; x < N; ++x) {
      for (size_t y = 0; y < M; ++y) {
        if (field[x][y] == '#')
          continue;
        for (auto [dx, dy] : deltas) {
          auto old_v = velocity.get(x, y, dx, dy);
          auto new_v = velocity_flow.get(x, y, dx, dy);
          if (old_v > VelocityT(0)) {
            assert(VelocityT(new_v) <= old_v);
            velocity.get(x, y, dx, dy) = VelocityT(new_v);
            auto force = (VelocityFlowT(old_v) - new_v) *
                         VelocityFlowT(rho[(int)field[x][y]]);
            if (field[x][y] == '.')
              force *= VelocityFlowT(0.8);
            if (field[x + dx][y + dy] == '#') {
              p[x][y] += PreasureT(force) / PreasureT(dirs[x][y]);
              total_delta_p += PreasureT(force) / PreasureT(dirs[x][y]);
            } else {
              p[x + dx][y + dy] +=
                  PreasureT(force) / PreasureT(dirs[x + dx][y + dy]);
              total_delta_p +=
                  PreasureT(force) / PreasureT(dirs[x + dx][y + dy]);
            }
          }
        }
      }
    }

    UT += 2;
    prop = false;
    for (size_t x = 0; x < N; ++x) {
      for (size_t y = 0; y < M; ++y) {
        if (field[x][y] != '#' && last_use[x][y] != UT) {
          if (random01() < move_prob(x, y)) {
            prop = true;
            propagate_move(x, y, true);
          } else {
            propagate_stop(x, y, true);
          }
        }
      }
    }

    if (prop) {
      cout << "Tick " << i << ":\n";
      for (size_t x = 0; x < N; ++x) {
        cout << field[x] << "\n";
      }
    }
  }
};

#endif // BIGHW2__SIMULATION_H_
