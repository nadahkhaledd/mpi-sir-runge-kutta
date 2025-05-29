# Compiler and flags
MPICC = mpic++
CFLAGS = -Wall -O2

# Source files
SRCS_CORE = $(wildcard src/core/*.cpp)
SRCS_TESTS = $(wildcard src/tests/*.cpp)
SRCS_MAIN = main.cpp

# Test sources
TEST_SOURCES = src/tests/TestRunner.cpp src/tests/TestSuite.cpp src/tests/test_main.cpp src/tests/TestConfig.cpp

# Object files
OBJS_CORE = $(patsubst src/core/%.cpp,output/core/%.o,$(SRCS_CORE))
OBJS_TESTS = $(patsubst src/tests/%.cpp,output/tests/%.o,$(SRCS_TESTS))
OBJS_MAIN = $(patsubst %.cpp,output/%.o,$(SRCS_MAIN))

# Executables
TARGET = sir_simulation
TEST_TARGET = sir_test_suite

# Default target
all: $(TARGET)

# Test target
test: $(OBJS_CORE) $(OBJS_TESTS)
	$(MPICC) $(CFLAGS) -o $(TEST_TARGET) $^

# Main executable
$(TARGET): $(OBJS_CORE) $(OBJS_MAIN)
	$(MPICC) $(CFLAGS) -o $@ $^

# Test executable
$(TEST_TARGET): $(OBJS_CORE) $(OBJS_TESTS)
	$(MPICC) $(CFLAGS) -o $@ $^

# Compile core source files
output/core/%.o: src/core/%.cpp
	@mkdir -p output/core
	$(MPICC) $(CFLAGS) -c $< -o $@

# Compile test source files
output/tests/%.o: src/tests/%.cpp
	@mkdir -p output/tests
	$(MPICC) $(CFLAGS) -c $< -o $@

# Compile main
output/%.o: %.cpp
	@mkdir -p output
	$(MPICC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf output $(TARGET) $(TEST_TARGET)

analyze: test
	@mkdir -p data/analysis
	python scripts/analyze_tests.py
	python scripts/visualize_comparison.py

plot: analyze
	@echo "Generating plots in data/analysis/"
	python scripts/PlottingSIRModelResults.py

.PHONY: all test clean analyze plot