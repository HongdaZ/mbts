#ifndef __VALIDATION__
#define __VALIDATION__

#include <algorithm> // std::copy std::sort std::equal std::rotate
#include <cassert> // assert macro
#include <cstdlib> // std::size_t
#include <iostream> // std streams
#include <iterator> // std::iterator_traits
#include "point.hpp" // left_turn right_turn
#include <vector> // std::vector

namespace validation {

  template<typename I, typename J>
  bool same_multiset(I p, I q, J r, J s) {
    using P = typename std::iterator_traits<I>::value_type;
    using S = std::vector<P>;
    std::size_t n = q - p;
    std::size_t m = s - r;
    if (m != n) {
      return false;
    }
    S backup;
    backup.resize(n);
    std::copy(p, q, backup.begin());
    std::sort(backup.begin(), backup.end(),
      [](P const& a, P const& b) -> bool {
        return (a.x < b.x) or (a.x == b.x and a.y < b.y);
      });
    S other;
    other.resize(n);
    std::copy(r, s, other.begin());
    std::sort(other.begin(), other.end(),
      [](P const& a, P const& b) -> bool {
        return (a.x < b.x) or (a.x == b.x and a.y < b.y);
      });
    return std::equal(backup.begin(), backup.end(), other.begin(),
      [](P const& a, P const& b) -> bool {
        return (a.x == b.x) and (a.y == b.y);
      });
  }

  // [p, r) circular chain of vertices

  template<typename I>
  bool left_turns_only(I p, I r) {
    std::size_t h = r - p;
    if (h < 3) {
      return true;
    }
    for (I q = p; q != r; ++q) {
      I next = q + 1;
      I next_after = next + 1;
      if (q == r - 2) {
        next_after = p;
      }
      else if (q == r - 1) {
        next = p;
        next_after = p + 1;
      }
      if (not left_turn(*q, *next, *next_after)) {
        std::cerr << "Left-turn: " << *q << " " << *next << " "
                  << *next_after << " failed\n";
        return false;
      }
    }
    return true;
  }

  // [p, r) circular chain of vertices

  template<typename I>
  bool right_turns_only(I p, I r) {
    std::size_t h = r - p;
    if (h < 3) {
      return true;
    }
    for (I q = p; q != r; ++q) {
      I next = q + 1;
      I next_after = next + 1;
      if (q == r - 2) {
        next_after = p;
      }
      else if (q == r - 1) {
        next = p;
        next_after = p + 1;
      }
      if (not right_turn(*q, *next, *next_after)) {
        std::cerr << "Right-turn: " << *q << " " << *next << " "
                  << *next_after << " failed\n";
        return false;
      }
    }
    return true;
  }

  // [p, r) half-circular chain of vertices

  template<typename I>
  bool monotone(I p, I r) {
    assert(p != r);
    using P = typename std::iterator_traits<I>::value_type;
    std::size_t h = r - p;
    if (h == 0 or h == 1) {
      return true;
    }
    I q = p;
    while (q != r - 1) {
      P a = *q;
      P b = *(q + 1);
      if (not (a.x < b.x or (a.x == b.x and a.y < b.y))) {
        std::cerr << "Convex: " << a << " " << b
                  << ": not monotone\n";
        return false;
      }
      ++q;
    }
    return true;
  }

  // [p, r) circular chain of vertices

  template<typename I>
  bool convex_polygon(I p, I r) {
    using P = typename std::iterator_traits<I>::value_type;
    using S = std::vector<P>;
    using J = typename S::iterator;
    std::size_t h = r - p;
    if (h == 0 or h == 1) {
      return true;
    }
    if (h == 2) {
      if (*p != *(p + 1)) {
        return true;
      }
      else {
        std::cerr << "Convex: h = 2: two equal points\n";
        return false;
      }
    }
    // h >= 3
    bool clockwise = right_turn(*p, *(p + 1), *(p + 2));
    if (clockwise and (not right_turns_only(p, r))) {
      std::cerr << "Convex: h = " << h
                << ": not a spiral (cw)\n";
      return false;
    }
    if (not clockwise and (not left_turns_only(p, r))) {
      std::cerr << "Convex: h = " << h
                << ": not a spiral (ccw)\n";
      return false;
    }

    using P = typename std::iterator_traits<I>::value_type;
    using S = std::vector<P>;
    using J = typename S::iterator;
    S hull;
    hull.resize(h);
    std::copy(p, r, hull.begin());
    J west = std::min_element(hull.begin(), hull.end(),
      [](P const& a, P const& b) -> bool {
        return (a.x < b.x) or (a.x == b.x and a.y < b.y);
      });
    (void) std::rotate(hull.begin(), west, hull.end());
    hull.push_back(*west);
    west = hull.begin();
    J east = std::max_element(hull.begin(), hull.end() - 1,
      [](P const& a, P const& b) -> bool {
        return (a.x < b.x) or (a.x == b.x and a.y < b.y);
      });
    assert(west != east);
    if (clockwise) {
      if (not monotone(west, east + 1)) {
        std::cerr << "Convex: h = " << h
                  << ": upper hull not monotone (cw)\n";
        return false;
      }
      std::reverse(east, hull.end());
      if (not monotone(east, hull.end())) {
        std::cerr << "Convex: h = " << h
                  << ": lower hull not monotone (cw)\n";
        return false;
      }
      return true;
    }
    // counterclockwise
    if (not monotone(west, east + 1)) {
      std::cerr << "Convex: h = " << h
                << ": lower hull not monotone (ccw)\n";
      return false;
    }
    std::reverse(east, hull.end());
    if (not monotone(east, hull.end())) {
      std::cerr << "Convex: h = " << h
                << ": upper hull not monotone (ccw)\n";
      return false;
    }
    return true;
  }

  // [p, r) upper hull from left to right

  template<typename I, typename P>
  bool above(I p, I r, P const& u) {
    std::size_t h = r - p;
    assert(h > 1);
    I last = r - 1;
    assert(u.x >= (*p).x and u.x <= (*last).x);
    if ((*p).x == (*(p + 1)).x) {
      ++p;
    }
    if (h != 1 and (*last).x == (*(last - 1)).x) {
      --r;
    }
    I q = std::lower_bound(p, r, u,
      [](P const& a, P const& b) -> bool {
        return (a.x < b.x);
      });
    if (q == p) {
      return u.y > (*p).y;
    }
    return left_turn(*(q - 1), *q, u);
  }

  // [p, r) lower hull from left to right

  template<typename I, typename P>
  bool below(I p, I r, P const& u) {
    std::size_t h = r - p;
    assert(h > 1);
    I last = r - 1;
    assert(u.x >= (*p).x and u.x <= (*last).x);
    if ((*p).x == (*(p + 1)).x) {
      ++p;
    }
    if (h != 1 and (*last).x == (*(last - 1)).x) {
      --r;
    }
    I q = std::lower_bound(p, r, u,
      [](P const& a, P const& b) -> bool {
        return (a.x < b.x);
      });
    if (q == p) {
      return u.y < (*p).y;
    }
    return right_turn(*(q - 1), *q, u);
  }

  // [r, s) multiset of points; [p, q) convex polygon

  template<typename I, typename J>
  bool all_inside(I r, I s, J p, J q) {
    std::size_t h = q - p;
    if (h == 0) {
      if (r != s) {
        std::cerr << "Inside: h = 0: no points can be inside\n";
        return false;
      }
      return true;
    }
    if (h == 1) {
      for (J i = r; i != s; ++i) {
        if (*i != *p) {
          std::cerr << "Inside: h = 1: all points not equal\n";
          return false;
        }
      }
      return true;
    }
    if (h == 2) {
      for (J i = r; i != s; ++i) {
        if (not on_line_segment(*i, *p, *(p + 1))) {
          std::cerr << "Inside: h = 2: all points not collinear\n";
          return false;
        }
      }
      return true;
    }
    using P = typename std::iterator_traits<I>::value_type;
    using S = std::vector<P>;
    using K = typename S::iterator;
    S hull;
    hull.resize(h);
    std::copy(p, q, hull.begin());
    K west = std::min_element(hull.begin(), hull.end(),
      [](P const& a, P const& b) -> bool {
        return (a.x < b.x) or (a.x == b.x and a.y < b.y);
      });
    (void) std::rotate(hull.begin(), west, hull.end());
    hull.push_back(*west);
    west = hull.begin();
    K east = std::max_element(hull.begin(), hull.end() - 1,
      [](P const& a, P const& b) -> bool {
        return (a.x < b.x) or (a.x == b.x and a.y < b.y);
      });
    assert(west != east);
    for (I i = r; i != s; ++i) {      
      if ((*i).x < (*west).x) {
        std::cerr << "Inside: " << *i
                  << " on left of the output\n";
        return false;
      }
      if ((*i).x > (*east).x) {
        std::cerr << "Inside: " << *i
                  << " on right of the output\n";
        return false;
      }
    }
    bool clockwise = right_turn(*west, *(west + 1), *(west + 2));
    if (clockwise) {
      for (I i = r; i != s; ++i) {      
        if (above(west, east + 1, *i)) {
          std::cerr << "Inside: " << *i
                    << " above the upper hull (cw)\n";
          return false;
        }
      }
      std::reverse(east, hull.end());
      for (I i = r; i != s; ++i) {      
        if (below(east, hull.end(), *i)) {
          std::cerr << "Inside: " << *i
                    << " below the lower hull (cw)\n";
          return false;
        }
      }
    }
    else { // counterclockwise
      for (I i = r; i != s; ++i) {      
        if (below(west, east + 1, *i)) {
          std::cerr << "Inside: " << *i
                    << " below the lower hull (ccw)\n";
          return false;
        }
      }
      std::reverse(east, hull.end());
      for (I i = r; i != s; ++i) {      
        if (above(east, hull.end(), *i)) {
          std::cerr << "Inside: " << *i
                    << " above the upper hull (ccw)\n";
          return false;
        }
      }
    }
    return true;
  }
}

#endif