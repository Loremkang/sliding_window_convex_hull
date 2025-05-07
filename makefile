all: ./build/test

./build/test: ./src/test.cpp ./include/sliding_window_convex_hull.h
	@echo "Current directory: $(CURDIR)/"
	@echo "Makefile directory: $(dir $(abspath $(lastword $(MAKEFILE_LIST))))"
	@if [ "$(CURDIR)/" != "$(dir $(abspath $(lastword $(MAKEFILE_LIST))))" ]; then \
		echo "Error: You are not in the Makefile directory."; \
		exit 1; \
	fi
	mkdir -p build
	g++ -o ./build/test ./src/test.cpp -g -Wall -Wextra -Wpedantic -std=c++17 -O0 -Iinclude

.PHONY: clean

clean:
	@echo "Current directory: $(CURDIR)/"
	@echo "Makefile directory: $(dir $(abspath $(lastword $(MAKEFILE_LIST))))"
	@if [ "$(CURDIR)/" != "$(dir $(abspath $(lastword $(MAKEFILE_LIST))))" ]; then \
		echo "Error: You are not in the Makefile directory."; \
		exit 1; \
	fi
	rm -rf ./build