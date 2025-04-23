CXX = mpic++    
CXXFLAGS = -Wall -O2    

# Dynamically find all .cpp files in src/ and include main.cpp explicitly
SRCS = $(wildcard src/*.cpp) main.cpp        
OBJS = $(patsubst src/%.cpp,output/%.o,$(SRCS))  
OBJS := $(patsubst main.cpp,output/main.o,$(OBJS))  # Handle main.cpp separately
EXEC = sir_simulation   

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)   

output/%.o: src/%.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

output/main.o: main.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: all clean