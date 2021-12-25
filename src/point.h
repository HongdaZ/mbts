#ifndef __POINT__
#define __POINT__

#include <algorithm> // std::min std::max
#include <cassert> // assert macro
#include <cstddef> // std::size_t
#include <limits> // std::numeric_limits
#include <string> // std::string std::to_string

template<typename T>
class point {
public:
  
  using self = point<T>;
  using coordinate = T;
  
  T x;
  T y;
  
  explicit point()
    : x(0), y(0) {
  }
  
  point(T x_coordinate, T y_coordinate)
    : x(x_coordinate), y(y_coordinate) {
  }
  
  std::string to_string() const {
    return "(" + std::to_string(x) + ", " + std::to_string(y) + ")";
  }
  
  bool operator==(point<T> const& other) const {
    return x == other.x and y == other.y;
  }
  
  bool operator!=(point<T> const& other) const {
    return not (*this == other);
  }
  
  bool operator<(point<T> const& other) const {
    return (x < other.x) or (x == other.x and y < other.y);
  }
  
  bool operator>(point<T> const& other) const {
    return other < *this;
  }
  
  friend std::ostream& operator<<(std::ostream& os, point<T> const& p) {
    return os << p.to_string();
  }
};

namespace double_precision {

template<typename P>
bool left_turn(P const& p, P const& q, P const& r) {
  
#ifdef MEASURE_TURNS
  
  ++turns;
  
#endif
  
  using T = typename P::coordinate;
  using Z = long long;
  Z mid_x = Z(q.x);
  Z mid_y = Z(q.y);
  Z dx_1 = mid_x - Z(p.x);
  Z dx_2 = Z(r.x) - mid_x;
  Z dy_1 = mid_y - Z(p.y);
  Z dy_2 = Z(r.y) - mid_y;
  
  /* We have 31 bits + 1 for the sign in int, and
   we may need (31 (+1 from minus)) * 2 (from mult) = 64 bits
   in total, but we only have 63 bits + 1 for the sign.
   If we overflow the same direction, the order is preserved. */
  bool check1 = Z(T(dx_1)) != dx_1 and Z(T(dy_2)) != dy_2;
  bool check2 = Z(T(dx_2)) != dx_2 and Z(T(dy_1)) != dy_1;
  if (check1 or check2) {
    // overflow happening
    if (check1) {
      if ((dx_1 > 0) == (dy_2 > 0)) {
        // first part overflows above
        if (check2 and ((dx_2 > 0) == (dy_1 > 0))) {
          // both overflow above
          return dx_1 * dy_2 > dx_2 * dy_1;
        }
        else {
          return true; // only first term overflows above
        }
      }
      else {
        // first term overflows below
        if (check2 and ((dx_2 > 0) != (dy_1 > 0))) {
          // second term also overflows below
          return dx_1 * dy_2 > dx_2 * dy_1;
        }
        else {
          // only first term overflows below
          return false;
        }
      }
    }
    else {
      // second term overflows; result is inverse of the direction
      return (dx_2 > 0) != (dy_1 > 0);
    }
    return dx_1 * dy_2 > dx_2 * dy_1;
  }
  else {
    // standard setting
    return dx_1 * dy_2 > dx_2 * dy_1;
  }
}

template<typename P>
bool no_turn(P const& p, P const& q, P const& r) {
  
#ifdef MEASURE_TURNS
  
  ++turns;
  
#endif
  
  using T = typename P::coordinate;
  using Z = long long;
  Z mid_x = Z(q.x);
  Z mid_y = Z(q.y);
  Z dx_1 = mid_x - Z(p.x);
  Z dx_2 = Z(r.x) - mid_x;
  Z dy_1 = mid_y - Z(p.y);
  Z dy_2 = Z(r.y) - mid_y;
  
  bool check1 = Z(T(dx_1)) != dx_1 and Z(T(dy_2)) != dy_2;
  bool check2 = Z(T(dx_2)) != dx_2 and Z(T(dy_1)) != dy_1;
  if (check1 or check2) {
    // if they do not both overflow, they must be different
    if (check1 and check2) {
      // if overflow direction not the same, they are different
      if (((dx_1 > 0) == (dy_2 > 0)) == ((dx_2 > 0) == (dy_1 > 0))) {
        // the non-overflow part must also be the same
        return dx_1 * dy_2 == dx_2 * dy_1;
      }
    }
    return false;
  }
  else {
    return dx_1 * dy_2 == dx_2 * dy_1;
  }
}

// may do controlled overflow on very large values

template<typename P>
unsigned long long signed_area(P const& p, P const& q, P const& r) {
  assert(not right_turn(p, q, r));
  
#ifdef MEASURE_TURNS
  
  ++turns;
  
#endif
  
  using Z = long long;
  Z mid_x = Z(q.x);
  Z mid_y = Z(q.y);
  Z dx_1 = mid_x - Z(p.x);
  Z dx_2 = Z(r.x) - mid_x;
  Z dy_1 = mid_y - Z(p.y);
  Z dy_2 = Z(r.y) - mid_y;
  return (unsigned long long) (dx_1 * dy_2 - dx_2 * dy_1);
}
}

using namespace double_precision;

template<typename P>
bool right_turn(P const& p, P const& q, P const& r) {
  return left_turn(r, q, p);
}

// is p on the line segment {q, r}?

template<typename P>
bool on_line_segment(P const& p, P const& q, P const& r) {
  P left = std::min(q, r);
  P right = std::max(q, r);
  if (p < left) {
    return false;
  }
  if (p > right) {
    return false;
  }
  return no_turn(p, q, r);
}

#endif