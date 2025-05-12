#pragma once

#include <cassert>
#include <cstddef>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <vector>

namespace sliding_window_convex_hull {

constexpr long double EPS = 1e-25;

// This is a custom deque implementation that uses a circular buffer
// with a size that is a power of two. It provides O(1) amortized
// complexity for push and pop operations at both ends.
template <typename T, typename Alloc = std::allocator<T>>
struct CustomizeDeque {
  size_t LowBit(size_t x) { return x & -x; }

  CustomizeDeque(size_t n) {
    size_t buffer_size = n + 1;
    while (buffer_size != LowBit(buffer_size)) {
      buffer_size += LowBit(buffer_size);
    }
    buffer_size_mask_ = buffer_size - 1;
    assert(buffer_size > 0 && (buffer_size & buffer_size_mask_) == 0);

    printf("buffer_size = %ld\n", buffer_size);

    front_ = 0;
    back_ = 0;
    buffer_size_ = buffer_size;
    data_ = (Alloc().allocate(buffer_size));
    assert(data_ != nullptr);
  }

  void clear() {
    front_ = 0;
    back_ = 0;
  }

  ~CustomizeDeque() {
    if (data_ != nullptr) {
      Alloc().deallocate(data_, buffer_size_);
      data_ = nullptr;
    }
  }

  void push_back(const T &value) {
    data_[back_] = value;
    back_ = (back_ + 1) & buffer_size_mask_;
    assert(front_ != back_);
  }

  void push_front(const T &value) {
    front_ = (front_ - 1 + buffer_size_) & buffer_size_mask_;
    data_[front_] = value;
    assert(front_ != back_);
  }
  void pop_front() { front_ = (front_ + 1) & buffer_size_mask_; }
  void pop_back() { back_ = (back_ - 1 + buffer_size_) & buffer_size_mask_; }

  size_t size() const {
    return (back_ + buffer_size_ - front_) & buffer_size_mask_;
  }

  T &front() { return data_[front_]; }
  const T &front() const { return data_[front_]; }

  T &back() {
    size_t prev_back = (back_ - 1) & buffer_size_mask_;
    return data_[prev_back];
  }
  const T &back() const {
    size_t prev_back = (back_ - 1) & buffer_size_mask_;
    return data_[prev_back];
  }

  T &operator[](size_t index) {
    assert(index < size());
    return data_[(front_ + index) & buffer_size_mask_];
  }
  const T &operator[](size_t index) const {
    assert(index < size());
    return data_[(front_ + index) & buffer_size_mask_];
  }

  size_t front_, back_;
  size_t buffer_size_, buffer_size_mask_;
  T *data_;
};

struct Slope {
  long double dx, dy;
  Slope operator*(long double k) const { return Slope{dx * k, dy * k}; }
  bool operator<(const Slope &p) const { return dy * p.dx < dx * p.dy; }
  bool operator>(const Slope &p) const { return dy * p.dx > dx * p.dy; }
  bool operator<=(const Slope &p) const { return !(*this > p); }
  bool operator>=(const Slope &p) const { return !(*this < p); }
  bool operator==(const Slope &p) const { return dy * p.dx == dx * p.dy; }

  explicit operator long double() const { return dy / dx; }
  long double decode() { return (long double)(*this); }

  static long double CrossProduct(const Slope &a, const Slope &b) {
    return a.dx * b.dy - a.dy * b.dx;
  }
};

struct Point {
  long double x;
  long double y;
  Point(long double x = 0, long double y = 0) : x(x), y(y) {}
  Slope operator-(const Point &p) const { return Slope{x - p.x, y - p.y}; }
  static long double CrossProduct(const Point &O, const Point &A,
                                  const Point &B) {
    Slope OA = A - O;
    Slope OB = B - O;
    return Slope::CrossProduct(OA, OB);
  }
};

template <template <typename> class Alloc = std::allocator>
struct PopConvexHull {
  const Point *p_;
  size_t n_;
  bool is_up_;
  size_t l_, r_;
  CustomizeDeque<size_t, Alloc<size_t>> hull_pos_;
  size_t *prev_hull_pos_;
  size_t *next_hull_pos_;
  std::vector<size_t> tmp_hull_pos_;

  PopConvexHull(const Point *p, size_t n, bool up, size_t *prev_hull_pos,
                size_t *next_hull_pos)
      : p_(p),
        n_(n),
        is_up_(up),
        l_(0),
        r_(0),
        hull_pos_(n),
        prev_hull_pos_(prev_hull_pos),
        next_hull_pos_(next_hull_pos) {}

  size_t operator[](size_t i) const { return hull_pos_[i]; }

  size_t size() const { return hull_pos_.size(); }

  void Clear() {
    hull_pos_.clear();
    l_ = r_;
  }

  void Initialize(size_t l, size_t r) {
    hull_pos_.clear();

    l_ = l;
    r_ = r;
    for (size_t i = r; i-- > l;) {
      Point p = p_[i];
      size_t hull_size = hull_pos_.size();
      while (hull_size >= 2) {
        size_t pos1 = hull_pos_[0];
        size_t pos2 = hull_pos_[1];
        assert(i < pos1 && pos1 < pos2);

        double cross_product = Point::CrossProduct(p, p_[pos1], p_[pos2]);
        bool should_pop = is_up_ ? cross_product > EPS : cross_product < -EPS;
        if (should_pop) {
          hull_pos_.pop_front();
        } else {
          hull_size = hull_pos_.size();
          break;
        }
        hull_size = hull_pos_.size();
      }
      next_hull_pos_[i] = hull_pos_.size() ? hull_pos_.front() : i;
      prev_hull_pos_[next_hull_pos_[i]] = i;
      hull_pos_.push_front(i);
    }
  }

  void PopHullFront() {
    assert(l_ == hull_pos_.front());
    hull_pos_.pop_front();
    if (hull_pos_.size() == 0) {
      l_ = r_;
    } else {
      l_ = hull_pos_.front();
      prev_hull_pos_[l_] = l_;
    }
  }

  void PopFront() {
    tmp_hull_pos_.clear();

    assert(l_ == hull_pos_.front());
    hull_pos_.pop_front();
    if (hull_pos_.size() == 0) {
      return;
    }

    size_t next_l = l_ + 1;
    while (next_l != hull_pos_.front()) {
      tmp_hull_pos_.push_back(next_l);
      next_l = next_hull_pos_[next_l];
    }
    for (size_t i = tmp_hull_pos_.size(); i-- > 0;) {
      prev_hull_pos_[hull_pos_[0]] = tmp_hull_pos_[i];
      hull_pos_.push_front(tmp_hull_pos_[i]);
    }
    prev_hull_pos_[hull_pos_[0]] = hull_pos_[0];

    l_++;

    assert(l_ == hull_pos_.front());
  }

  void print() {
    printf("%s PopHull:\n", is_up_ ? "Lower" : "Upper");
    size_t s = size();
    for (size_t i = 0; i < s; i++) {
      printf("%ld ", operator[](i));
    }
    printf("\n");
  }

  size_t PrevHullPos(size_t i) const {
    assert(i >= l_ && i < r_);
    return prev_hull_pos_[i];
  }

  size_t NextHullPos(size_t i) const {
    assert(i >= l_ && i < r_);
    return next_hull_pos_[i];
  }

  void verify_hull() {
    size_t hull_size = size();
    for (size_t i = 0; i < hull_size; i++) {
      size_t pos = operator[](i);
      size_t prev = (i == 0) ? pos : operator[](i - 1);
      size_t next = (i == hull_size - 1) ? pos : operator[](i + 1);
      assert(pos >= l_ && pos < r_);
      assert(PrevHullPos(pos) == prev);
      assert(NextHullPos(pos) == next);
    }

    for (size_t i = 0; i + 2 < hull_size; i++) {
      size_t pos1 = operator[](i);
      size_t pos2 = operator[](i + 1);
      size_t pos3 = operator[](i + 2);
      assert(pos1 < pos2 && pos2 < pos3);
      assert(pos1 >= l_ && pos2 >= l_ && pos3 >= l_);
      double cross_product = Point::CrossProduct(p_[pos1], p_[pos2], p_[pos3]);
      if (is_up_) {
        assert(cross_product < EPS);
      } else {
        assert(cross_product > -EPS);
      }
    }

    for (size_t i = 0; i + 1 < hull_size; i++) {
      size_t pos1 = operator[](i);
      size_t pos2 = operator[](i + 1);
      assert(pos1 < pos2);
      assert(pos1 >= l_ && pos2 >= l_);
      for (size_t j = pos1 + 1; j < pos2; j++) {
        double cross_product = Point::CrossProduct(p_[pos1], p_[pos2], p_[j]);
        if (is_up_) {
          assert(cross_product < 0);
        } else {
          assert(cross_product > 0);
        }
      }
    }
  }
};

template <template <typename> class Alloc = std::allocator>
struct PushConvexHull {
  const Point *p_;
  size_t n_;
  bool is_up_;
  size_t l_, r_;
  CustomizeDeque<size_t, Alloc<size_t>> hull_pos_;
  size_t *prev_hull_pos_;
  size_t *next_hull_pos_;

  PushConvexHull(const Point *p, size_t n, bool up, size_t *prev_hull_pos,
                 size_t *next_hull_pos)
      : p_(p),
        n_(n),
        is_up_(up),
        l_(0),
        r_(0),
        hull_pos_(n),
        prev_hull_pos_(prev_hull_pos),
        next_hull_pos_(next_hull_pos) {}

  size_t operator[](size_t i) const {
    assert(i < size());
    return hull_pos_[i];
  }

  size_t size() const { return hull_pos_.size(); }

  void Clear() {
    hull_pos_.clear();
    l_ = r_;
  }

  size_t PopHullFront() {
    hull_pos_.pop_front();
    if (hull_pos_.size() == 0) {
      // run out of elements
      l_ = r_;
      return l_;
    } else {
      l_ = hull_pos_.front();
      prev_hull_pos_[l_] = l_;
      return l_;
    }
  }

  void PushBack() {
    assert(r_ < n_);
    size_t hull_size = hull_pos_.size();
    while (hull_size >= 2) {
      size_t pos1 = hull_pos_[hull_size - 2];
      size_t pos2 = hull_pos_[hull_size - 1];
      double cross_product = Point::CrossProduct(p_[pos1], p_[pos2], p_[r_]);
      bool should_pop = is_up_ ? cross_product > EPS : cross_product < -EPS;
      if (should_pop) {
        hull_pos_.pop_back();
      } else {
        hull_size = hull_pos_.size();
        break;
      }
      hull_size = hull_pos_.size();
    }
    if (hull_size == 0) {
      prev_hull_pos_[r_] = r_;
    } else {
      next_hull_pos_[hull_pos_.back()] = r_;
      prev_hull_pos_[r_] = hull_pos_.back();
    }
    next_hull_pos_[r_] = r_;
    hull_pos_.push_back(r_);
    r_++;
  }

  void print() {
    printf("%s PushHull:\n", is_up_ ? "Lower" : "Upper");
    size_t s = size();
    for (size_t i = 0; i < s; i++) {
      printf("%ld ", operator[](i));
    }
    printf("\n");
  }

  size_t PrevHullPos(size_t i) const {
    assert(i >= l_ && i < r_);
    return prev_hull_pos_[i];
  }

  size_t NextHullPos(size_t i) const {
    assert(i >= l_ && i < r_);
    return next_hull_pos_[i];
  }

  void verify_hull() {
    size_t hull_size = size();
    for (size_t i = 0; i < hull_size; i++) {
      size_t pos = operator[](i);
      size_t prev = (i == 0) ? pos : operator[](i - 1);
      size_t next = (i == hull_size - 1) ? pos : operator[](i + 1);
      assert(pos >= l_ && pos < r_);
      assert(PrevHullPos(pos) == prev);
      assert(NextHullPos(pos) == next);
    }

    for (size_t i = 0; i + 2 < hull_size; i++) {
      size_t pos1 = operator[](i);
      size_t pos2 = operator[](i + 1);
      size_t pos3 = operator[](i + 2);
      assert(pos1 < pos2 && pos2 < pos3);
      assert(pos1 >= l_ && pos2 >= l_ && pos3 >= l_);
      double cross_product = Point::CrossProduct(p_[pos1], p_[pos2], p_[pos3]);
      if (is_up_) {
        assert(cross_product < EPS);
      } else {
        assert(cross_product > -EPS);
      }
    }

    for (size_t i = 0; i + 1 < hull_size; i++) {
      size_t pos1 = operator[](i);
      size_t pos2 = operator[](i + 1);
      assert(pos1 < pos2);
      assert(pos1 >= l_ && pos2 >= l_);
      for (size_t j = pos1 + 1; j < pos2; j++) {
        double cross_product = Point::CrossProduct(p_[pos1], p_[pos2], p_[j]);
        if (is_up_) {
          assert(cross_product < 0);
        } else {
          assert(cross_product > 0);
        }
      }
    }
  }
};

template <template <typename> class Alloc = std::allocator>
struct SlidingWindowConvexHull {
  const Point *p_;
  size_t n_;
  bool is_up_;
  size_t l_, r_, mid_;
  size_t tan_l_, tan_r_;
  size_t tan_l_pos_, tan_r_pos_;

  size_t *prev_hull_pos_;
  size_t *next_hull_pos_;
  PopConvexHull<Alloc> pop_hull_;
  PushConvexHull<Alloc> push_hull_;

  SlidingWindowConvexHull(const Point *p, size_t n, bool up)
      : p_(p),
        n_(n),
        is_up_(up),
        l_(0),
        r_(0),
        mid_(0),
        tan_l_(0),
        tan_r_(0),
        tan_l_pos_(0),
        tan_r_pos_(0),
        prev_hull_pos_(Alloc<size_t>().allocate(n)),
        next_hull_pos_(Alloc<size_t>().allocate(n)),
        pop_hull_(p, n, up, prev_hull_pos_, next_hull_pos_),
        push_hull_(p, n, up, prev_hull_pos_, next_hull_pos_) {}

  ~SlidingWindowConvexHull() {
    if (prev_hull_pos_ != nullptr) {
      Alloc<size_t>().deallocate(prev_hull_pos_, n_);
      prev_hull_pos_ = nullptr;
    }
    if (next_hull_pos_ != nullptr) {
      Alloc<size_t>().deallocate(next_hull_pos_, n_);
      next_hull_pos_ = nullptr;
    }
  }

  size_t size() const {
    size_t size = push_hull_.size() - tan_r_pos_;
    if (pop_hull_.size()) {
      size += tan_l_pos_ + 1;
    }
    return size;
  }

  size_t operator[](size_t i) const {
    assert(i < size());
    if (pop_hull_.size()) {
      return (i <= tan_l_pos_) ? pop_hull_[i]
                               : push_hull_[i - tan_l_pos_ - 1 + tan_r_pos_];
    } else {
      return push_hull_[i];
    }
  }

  Point p(size_t i) const { return p_[operator[](i)]; }

  size_t GetNextHullPos(size_t i) const {
    assert(i >= l_ && i < r_);
    return next_hull_pos_[i];
  }

  size_t GetHullPos(size_t i) const {
    assert(i < size());
    return (i < tan_l_pos_) ? pop_hull_[i]
                            : push_hull_[i - tan_l_pos_ + tan_r_pos_];
  }

  bool move_tan_l_left() {
    if (tan_l_pos_ == 0) {
      return false;
    }
    assert(pop_hull_.size() > 0 && push_hull_.size() > 0);

    size_t pos1 = pop_hull_[tan_l_pos_ - 1];
    size_t pos2 = pop_hull_[tan_l_pos_];
    double cross_product = Point::CrossProduct(p_[pos1], p_[pos2], p_[tan_r_]);
    // try to avoid moving by EPS
    bool should_move = is_up_ ? cross_product > EPS : cross_product < -EPS;
    if (should_move) {
      tan_l_pos_--;
      tan_l_ = pop_hull_[tan_l_pos_];
      return true;
    }
    return false;
  }

  bool move_tan_r_left() {
    if (tan_r_pos_ == 0) {
      return false;
    }
    assert(pop_hull_.size() > 0 && push_hull_.size() > 0);

    size_t pos1 = push_hull_[tan_r_pos_ - 1];
    size_t pos2 = push_hull_[tan_r_pos_];
    double cross_product = Point::CrossProduct(p_[tan_l_], p_[pos1], p_[pos2]);
    // try to move by EPS
    bool should_move = is_up_ ? cross_product < EPS : cross_product > -EPS;
    if (should_move) {
      tan_r_pos_--;
      tan_r_ = push_hull_[tan_r_pos_];
      return true;
    }
    return false;
  }

  void PushBack() {
    assert(r_ < n_);
    if (pop_hull_.size() == 0) {
      push_hull_.PushBack();
      tan_r_pos_ = 0;
      tan_r_ = push_hull_[0];
    } else {
      double cross_product = Point::CrossProduct(p_[tan_l_], p_[tan_r_], p_[r_]);
      bool should_reset = is_up_ ? cross_product > EPS : cross_product < -EPS;
      push_hull_.PushBack();
      if (should_reset) {
        tan_r_ = r_;
        tan_r_pos_ = push_hull_.size() - 1;
        while (move_tan_l_left());
      }
    }
    r_++;
  }

  void PopHullFront() {
    assert(size() > 0);

    if (size() == 1) {
      push_hull_.r_ = pop_hull_.r_ = l_ = mid_ = r_;
      push_hull_.Clear();
      pop_hull_.Clear();
      tan_l_pos_ = tan_r_pos_ = 0;
      return;
    } else {
      size_t pos1 = operator[](1);
      if (pos1 >= mid_) {
        pop_hull_.r_ = l_ = mid_ = pos1;
        pop_hull_.Clear();
        while (push_hull_[0] != pos1) {
          push_hull_.PopHullFront();
        }
        tan_l_pos_ = tan_r_pos_ = 0;
        tan_r_ = push_hull_[0];
      } else {
        pop_hull_.PopHullFront();
        if (tan_l_pos_ == 0) {
          tan_l_ = pop_hull_[0];
        } else {
          tan_l_pos_ --;
          assert(tan_l_ == pop_hull_[tan_l_pos_]);
        }
  
        while (true) {
          if (move_tan_l_left()) {
            continue;
          }
          if (move_tan_r_left()) {
            continue;
          }
          break;
        }
      }
    }

    l_ = (size() > 0) ? operator[](0) : r_;
    mid_ = push_hull_.l_;
    return;
  }

  void PopFront() {
    if (l_ < operator[](0)) {
      // do nothing
      goto final;
    }
    assert(l_ == operator[](0));

    assert(l_ != r_);
    assert(size() > 0);

    if (pop_hull_.size() == 0) {
      mid_ = push_hull_.PopHullFront();
      pop_hull_.Initialize(l_, mid_);
      tan_l_pos_ = tan_r_pos_ = 0;
      tan_l_ = pop_hull_[0];
      if (push_hull_.size() > 0) {
        tan_r_ = push_hull_[0];
        prev_hull_pos_[mid_] = mid_;
      }
    }

    if (pop_hull_.size() == 1) {
      pop_hull_.PopFront();
      tan_l_pos_ = tan_r_pos_ = 0;
      // no valid tan_l_
      if (push_hull_.size() > 0) {
        tan_r_ = push_hull_[0];
      }
      goto final;
    } else {
      if (tan_l_pos_ == 0) {
        tan_l_pos_ = 1;
        tan_l_ = pop_hull_[1];
      }

      size_t l_size = pop_hull_.size();
      pop_hull_.PopFront();
      tan_l_pos_ = tan_l_pos_ + pop_hull_.size() - l_size;
      assert(tan_l_ == pop_hull_[tan_l_pos_]);

      while (true) {
        if (move_tan_l_left()) {
          continue;
        }
        if (move_tan_r_left()) {
          continue;
        }
        break;
      }
    }

  final:
    l_++;
    return;
  }

  size_t PrevHullPos(size_t i) const {
    assert(i >= l_ && i < r_);
    if (push_hull_.size() > 0 && pop_hull_.size() > 0 && i == tan_r_) {
      return tan_l_;
    } else {
      return prev_hull_pos_[i];
    }
  }

  size_t NextHullPos(size_t i) const {
    assert(i >= l_ && i < r_);
    if (push_hull_.size() > 0 && pop_hull_.size() > 0 && i == tan_l_) {
      return tan_r_;
    } else {
      return next_hull_pos_[i];
    }
  }

  void print() {
    printf("%s ConvexHull:\n", is_up_ ? "Lower" : "Upper");
    size_t s = size();
    for (size_t i = 0; i < s; i++) {
      printf("%ld ", operator[](i));
    }
    printf("\n");
    pop_hull_.print();
    push_hull_.print();
    printf("l_ = %ld, r_ = %ld, mid_ = %ld\n", l_, r_, mid_);
    printf("tan_l_ = %ld, tan_r_ = %ld\n", tan_l_, tan_r_);
    printf("tan_l_pos_ = %ld, tan_r_pos_ = %ld\n", tan_l_pos_, tan_r_pos_);
  }

  void verify_hull() {
    size_t hull_size = size();
    for (size_t i = 0; i < hull_size; i++) {
      size_t pos = operator[](i);
      size_t prev = (i == 0) ? pos : operator[](i - 1);
      size_t next = (i == hull_size - 1) ? pos : operator[](i + 1);
      assert(pos >= l_ && pos < r_);
      assert(PrevHullPos(pos) == prev);
      assert(NextHullPos(pos) == next);
    }

    for (size_t i = 0; i + 2 < hull_size; i++) {
      size_t pos1 = operator[](i);
      size_t pos2 = operator[](i + 1);
      size_t pos3 = operator[](i + 2);
      assert(pos1 < pos2 && pos2 < pos3);
      assert(pos1 >= l_ && pos2 >= l_ && pos3 >= l_);
      double cross_product = Point::CrossProduct(p_[pos1], p_[pos2], p_[pos3]);
      if (is_up_) {
        assert(cross_product < EPS);
      } else {
        assert(cross_product > -EPS);
      }
    }

    for (size_t i = 0; i + 1 < hull_size; i++) {
      size_t pos1 = operator[](i);
      size_t pos2 = operator[](i + 1);
      assert(pos1 < pos2);
      assert(pos1 >= l_ && pos2 >= l_);
      for (size_t j = pos1 + 1; j < pos2; j++) {
        double cross_product = Point::CrossProduct(p_[pos1], p_[pos2], p_[j]);
        if (is_up_) {
          assert(cross_product < 0);
        } else {
          assert(cross_product > 0);
        }
      }
    }
  }

  void verify() {
    size_t hull_size = size();
    if (hull_size < 2) {
      return;
    }
    assert(hull_size <= n_);
    assert(l_ < r_);
    assert(l_ >= 0 && r_ <= n_);

    if (l_ < mid_) {
      assert(pop_hull_.size() > 0);
      assert(tan_l_pos_ < pop_hull_.size());
      assert(tan_l_ == pop_hull_[tan_l_pos_]);
      assert(tan_l_ >= l_ && tan_l_ < mid_);

    } else {
      assert(l_ == mid_);
      assert(pop_hull_.size() == 0);
    }

    if (r_ > mid_) {
      assert(push_hull_.size() > 0);
      assert(tan_r_pos_ < push_hull_.size());
      assert(tan_r_ == push_hull_[tan_r_pos_]);
      assert(tan_r_ >= mid_ && tan_r_ < r_);
    } else {
      assert(r_ == mid_);
      assert(push_hull_.size() == 0);
    }

    verify_hull();
    pop_hull_.verify_hull();
    push_hull_.verify_hull();
  }
};

}  // namespace sliding_window_convex_hull