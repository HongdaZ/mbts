#ifndef __QUICKHULL__
#define __QUICKHULL__

#include "point.hpp" // signed_area left_turn

#include <algorithm> // std::minmax_element std::partitioning ...
#include <cassert> // assert macro
#include <cstdlib> // std::size_t
#include <iterator> // std::iterator_traits
#include <utility> // std::pair std::make_pair ...
#include "validation.hpp" // same_multiset convex_polygon all_inside
#include <vector> // std::vector

namespace quickhull {

  // move *st <- *rd and *nd <- *th by swaps in parallel

  template<typename I>
  void parallel_iter_swap(I st, I nd, I rd, I th) {
    assert(st != nd and rd != th);
    std::iter_swap(st, rd);
    if (th == st) {
      std::iter_swap(nd, rd);
    }
    else {
      std::iter_swap(nd, th);
    }
  }

  template<typename I>
  std::pair<I, I> find_poles(I first, I past) {
    using P = typename std::iterator_traits<I>::value_type;
    auto pair = std::minmax_element(first, past,
      [](P const& a, P const& b) -> bool {
        return (a.x < b.x) or (a.x == b.x and a.y < b.y);
      });
    return pair;
  }

  template<typename I>
  I find_furthest(I first, I past, I antipole) {
    assert(first != past);
    I pole = first;
    I answer = pole;
    using U = decltype(signed_area(*pole, *pole, *pole));
    U best = 0;
    if ((*antipole).x == (*pole).x) { // vertical
      for (I i = first + 1; i != past; ++i) {
        U _Phi = signed_area(*pole, *antipole, *i);
        if (_Phi > best or (_Phi == best and (*i).y < (*answer).y)) {
          answer = i;
          best = _Phi;
        }
      }
    }
    else {
      for (I i = first + 1; i != past; ++i) {
        U _Phi = signed_area(*pole, *antipole, *i);
        if (_Phi > best or (_Phi == best and (*i).x < (*answer).x)) {
          answer = i;
          best = _Phi;
        }
      }
    }
    return answer;
  }

  template<typename I>
  I partition_left_right(I first, I past, I antipole) {
    assert(first != past);
    using P = typename std::iterator_traits<I>::value_type;
    I pole = first;
    I middle = std::partition(pole + 1, past,
      [&](P const& q) -> bool {
        return not left_turn(*pole, q, *antipole);
      });
    return middle;
  }

  template <typename I>
  void swap_blocks(I source, I past_the_end, I target) {
    if (source == target or source == past_the_end) {
      return;
    }
    using P = typename std::iterator_traits<I>::value_type;
    I hole = target;
    P p = *target;
    I const last = past_the_end - 1;
    while (true) {
      *hole = *source;
      ++hole;
      if (source == last) {
        break;
      }
      *source = *hole;
      ++source;
    }
    *source = p;
  }

  template <typename I>
  void move_away(I here, I rest, I past) {
    if (here == rest or rest == past) {
      return;
    }
    if (rest - here < past - rest) {
      swap_blocks(here, rest, past - (rest - here));
    }
    else {
      swap_blocks(rest, past, here);
    }
  }

  template<typename I>
  I recurse(I pole, I past, I antipole) {
    std::size_t n = std::distance(pole, past);
    if (n == 1) {
      return past;
    }
    if (n == 2) {
      if (no_turn(*(pole + 1), *pole, *antipole)) {
        return pole + 1;
      }
      else {
        return past;
      }
    }
    I pivot = find_furthest(pole, past, antipole);
    if (no_turn(*pivot, *pole, *antipole)) {
      return pole + 1;
    }
    I last = past - 1;
    std::iter_swap(pivot, last); // pivot at the end
    I mid = partition_left_right(pole, last, last);
    I eliminated = recurse(pole, mid, last);
    std::iter_swap(mid, last);
    std::iter_swap(eliminated, mid); // pivot at its final place
    pivot = eliminated;
    std::size_t m = past - mid;
    ++mid;
    ++eliminated;
    move_away(eliminated, mid, past);
    I interior = partition_left_right(pivot, pivot + m, antipole); 
    eliminated = recurse(pivot, interior, antipole);
    return eliminated;
  }

  template<typename I>
  I solve(I first, I past) {
    std::size_t n = past - first;
    if (n < 2 or (n == 2 and *first != *(first + 1))) {
      return past;
    }
    std::pair<I, I> pair = find_poles(first, past);
    I west = first;
    I east = past - 1;
    parallel_iter_swap(west, east, std::get<0>(pair), std::get<1>(pair));
    if (*west == *east) {
      return first + 1;
    }
    I middle = partition_left_right(west, east, east);
    std::size_t m = past - middle;
    I eliminated = recurse(first, middle, east);
    std::iter_swap(middle, east);
    std::iter_swap(eliminated, middle); // east at its final place
    east = eliminated;
    ++middle;
    ++eliminated;
    move_away(eliminated, middle, past);
    eliminated = recurse(east, east + m, west); // downunder
    return eliminated;
  }

  template<typename I>
  bool check(I first, I past) {
    using P = typename std::iterator_traits<I>::value_type;
    using S = std::vector<P>;
    using J = typename S::iterator;
    S data;
    std::size_t n = past - first;
    data.resize(n);
    std::copy(first, past, data.begin());
    J rest = solve(data.begin(), data.end());
    bool ok = validation::same_multiset(data.begin(), data.end(), first, past) and validation::convex_polygon(data.begin(), rest) and validation::all_inside(rest, data.end(), data.begin(), rest);
    return ok;
  }
}

#endif
