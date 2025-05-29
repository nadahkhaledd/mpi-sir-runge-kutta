# Compiler and flags
MPICC = mpic++
CFLAGS = -Wall -O2

# Directories
MAIN_DIR = src/main
TEST_DIR = src/test
OUTPUT_DIR = output

# Source files
SRCS_MAIN = $(wildcard $(MAIN_DIR)/*.cpp)
SRCS_TEST = $(wildcard $(TEST_DIR)/*.cpp)
MAIN_CPP = main.cpp

# Object files
OBJS_MAIN = $(patsubst $(MAIN_DIR)/%.cpp,$(OUTPUT_DIR)/main/%.o,$(SRCS_MAIN))
OBJS_TEST = $(patsubst $(TEST_DIR)/%.cpp,$(OUTPUT_DIR)/test/%.o,$(SRCS_TEST))
OBJ_MAIN_CPP = $(OUTPUT_DIR)/main_exe.o

# Executables
TARGET = sir_simulation
TEST_TARGET = sir_test_suite

# Default target
all: $(TARGET)

# Test target
test: $(TEST_TARGET)

# Main executable
$(TARGET): $(OBJS_MAIN) $(OBJ_MAIN_CPP)
	$(MPICC) $(CFLAGS) -o $@ $^

# Test executable
$(TEST_TARGET): $(OBJS_MAIN) $(OBJS_TEST)
	$(MPICC) $(CFLAGS) -o $@ $^

# Compile main.cpp
$(OBJ_MAIN_CPP): $(MAIN_CPP)
	@mkdir -p $(OUTPUT_DIR)
	$(MPICC) $(CFLAGS) -c $< -o $@

# Compile main source files
$(OUTPUT_DIR)/main/%.o: $(MAIN_DIR)/%.cpp
	@mkdir -p $(OUTPUT_DIR)/main
	$(MPICC) $(CFLAGS) -c $< -o $@

# Compile test source files
$(OUTPUT_DIR)/test/%.o: $(TEST_DIR)/%.cpp
	@mkdir -p $(OUTPUT_DIR)/test
	$(MPICC) $(CFLAGS) -c $< -o $@

clean:
	rm -rf $(OUTPUT_DIR) $(TARGET) $(TEST_TARGET)

.PHONY: all test clean