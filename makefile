all: test

test: test.cpp sliding_window_convex_hull.h
	g++ -o test test.cpp -g -Wall -Wextra -Wpedantic -std=c++17 -O0

.PHONY: clean

clean:
	rm -f test