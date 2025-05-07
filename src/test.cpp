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
        STOP
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
        uint64_t action = hash64(rnd) % 3;
        return action == 0 ? POP : PUSH;
    };

    for (size_t i = 0; true; i ++) {
        auto action = act(i);
        if (action == STOP) {
            break;
        }
        size_t pos = (action == PUSH) ? swch.r_ : swch.l_;
        size_t l = swch.l_;
        size_t r = swch.r_;
        printf("i: %ld Action: %s Pos: %ld Size: %ld Dist: (%ld, %ld) \n", i, (action == PUSH ? "PUSH" : "POP"), pos, swch.size(), l, r);
        fflush(stdout);
        if (action == PUSH) {
            swch.PushBack();
        } else {
            swch.PopFront();
        }
        // swch.print();
        swch.verify();
    }
    return 0;    
}