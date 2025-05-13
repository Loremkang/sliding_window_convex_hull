#include <iostream>
#include <vector>
#include <algorithm>
#include "sliding_window_convex_hull.h"

using namespace std;

uint64_t hash64(uint64_t x) {
    x = (~x) + (x << 21); // x = (x << 21) - x - 1;
    x = x ^ (x >> 24);
    x = (x + (x << 3)) + (x << 8); // x * 0x000000A9
    x = x ^ (x >> 14);
    x = (x + (x << 2)) + (x << 4); // x * 0x000000A9
    x = x ^ (x >> 28);
    x = x + (x << 31);
    return x;
}

int main(int argc, char *argv[]) {

    size_t size = 1000000;
    uint64_t seed = 0;
    uint64_t limit = 100000000;

    if (argc >= 2) {
        size = atoi(argv[1]);
    }

    if (argc >= 3) {
        seed = atoi(argv[2]);
    }

    if (argc >= 4) {
        limit = atoi(argv[3]);
    }

    std::vector<uint64_t> data(size);
    for (size_t i = 0; i < size; i++) {
        data[i] = hash64(hash64(i) + hash64(seed)) % limit;
    }

    std::sort(data.begin(), data.end());

    sliding_window_convex_hull::Point *points = new sliding_window_convex_hull::Point[size];
    for (size_t i = 0; i < size; i++) {
        points[i] = sliding_window_convex_hull::Point(data[i], i);
    }

    sliding_window_convex_hull::SlidingWindowConvexHull<> swch(points, size, true);

    enum action {
        PUSH,
        POP,
        POPHULL,
        STOP
    };

    auto to_str_action = [&](action act) {
        switch (act) {
            case PUSH:
                return "PUSH";
            case POP:
                return "POP";
            case POPHULL:
                return "POPHULL";
            case STOP:
                return "STOP";
        }
        return "UNKNOWN";
    };

    auto act = [&](size_t rnd) {
        if (swch.l_ == size) {
            return STOP;
        }
        if (swch.l_ == swch.r_) {
            return PUSH;
        }
        if (swch.r_ == size) {
            return POP;
        }
        uint64_t action = hash64(rnd) % 100;
        if (action == 0) {
            return POPHULL;
        } else if (action < 20) {
            return POP;
        } else {
            return PUSH;
        }
    };

    for (size_t i = 0; true; i ++) {
        auto action = act(i);
        if (action == STOP) {
            break;
        }
        size_t pos = (action == PUSH) ? swch.r_ : swch.l_;
        size_t l = swch.l_;
        size_t r = swch.r_;
        printf("i: %ld Action: %s Pos: %ld Size: %ld Dist: (%ld, %ld) \n", i, to_str_action(action), pos, swch.size(), l, r);
        fflush(stdout);
        if (action == PUSH) {
            swch.PushBack();
            assert(swch.l_ == l && swch.r_ == r + 1);
            if (swch.size() > 0) {
            assert(swch[0] == swch.l_ && swch[swch.size() - 1] == swch.r_ - 1);
            }
        } else if (action == POP){
            swch.PopFront();
            assert(swch.l_ == l + 1 && swch.r_ == r);
            if (swch.size() > 0) {
                assert(swch[0] == swch.l_ && swch[swch.size() - 1] == swch.r_ - 1);
            }
        } else if (action == POPHULL) {
            size_t size = swch.size();
            if (size > 1) {
                size_t nxt = swch[1];
                swch.PopHullFront();
                assert(swch.l_ == nxt && swch.r_ == r);
                assert(swch[0] == swch.l_ && swch[swch.size() - 1] == swch.r_ - 1);
            } else {
                swch.PopHullFront();
                assert(swch.size() == 0);
                assert(swch.l_ == swch.r_);
            }
            assert(swch.size() + 1 == size);
        }
        // swch.print();
        swch.verify();
    }
    return 0;    
}