# Your existing Makefile content...
# CC = g++
MPICC = mpic++
# CFLAGS = -Wall -g -O0 # For debugging
CFLAGS = -Wall -O2     # For release/optimization
# LDFLAGS = -lm # Example for math library

# Define source files
SRCS_MPI_HANDLER = src/MPIHandler.cpp
SRCS_GRID_SIM = src/GridSimulation.cpp
SRCS_CSV_PARSER = src/CSVParser.cpp
SRCS_SIR_CELL = src/SIRCell.cpp
SRCS_SIR_MODEL = src/SIRModel.cpp
SRCS_TIMING_UTILS = src/TimingUtils.cpp # <<< ADDED
SRCS_MAIN = main.cpp

# Define object files based on source files
OBJS_MPI_HANDLER = $(SRCS_MPI_HANDLER:.cpp=.o)
OBJS_MPI_HANDLER := $(patsubst src/%,output/%,$(OBJS_MPI_HANDLER))

OBJS_GRID_SIM = $(SRCS_GRID_SIM:.cpp=.o)
OBJS_GRID_SIM := $(patsubst src/%,output/%,$(OBJS_GRID_SIM))

OBJS_CSV_PARSER = $(SRCS_CSV_PARSER:.cpp=.o)
OBJS_CSV_PARSER := $(patsubst src/%,output/%,$(OBJS_CSV_PARSER))

OBJS_SIR_CELL = $(SRCS_SIR_CELL:.cpp=.o)
OBJS_SIR_CELL := $(patsubst src/%,output/%,$(OBJS_SIR_CELL))

OBJS_SIR_MODEL = $(SRCS_SIR_MODEL:.cpp=.o)
OBJS_SIR_MODEL := $(patsubst src/%,output/%,$(OBJS_SIR_MODEL))

OBJS_TIMING_UTILS = $(SRCS_TIMING_UTILS:.cpp=.o) # <<< ADDED
OBJS_TIMING_UTILS := $(patsubst src/%,output/%,$(OBJS_TIMING_UTILS)) # <<< ADDED

OBJS_MAIN = $(SRCS_MAIN:.cpp=.o)
OBJS_MAIN := $(patsubst src/%,output/%,$(OBJS_MAIN))


# Output directory
OUTPUT_DIR = output/

# Target executable
TARGET = sir_simulation

# All objects
ALL_OBJS = $(OBJS_CSV_PARSER) $(OBJS_GRID_SIM) $(OBJS_MPI_HANDLER) $(OBJS_SIR_CELL) $(OBJS_SIR_MODEL) $(OBJS_TIMING_UTILS) $(OBJS_MAIN) # <<< ADDED OBJS_TIMING_UTILS

# Default target
all: $(TARGET)

$(TARGET): $(ALL_OBJS)
	$(MPICC) $(CFLAGS) -o $@ $^ $(LDFLAGS)

# Rule to compile .cpp files in src/ to .o files in output/
$(OUTPUT_DIR)%.o: src/%.cpp | $(OUTPUT_DIR)
	$(MPICC) $(CFLAGS) -c $< -o $@

# Create output directory if it doesn't exist
$(OUTPUT_DIR):
	mkdir -p $(OUTPUT_DIR)

# Clean rule
clean:
	rm -rf $(OUTPUT_DIR) $(TARGET)

.PHONY: all cleans