CXX = mpic++    
CXXFLAGS = -Wall -O2    

SRCS = src/SIRCell.cpp src/CSVParser.cpp src/SIRModel.cpp src/GridSimulation.cpp src/MPIHandler.cpp main.cpp        
OBJS = $(patsubst src/%.cpp,output/%.o,$(SRCS))  
EXEC = sir_simulation   

all: $(EXEC)

$(EXEC): $(OBJS)
	$(CXX) $(CXXFLAGS) -o $@ $(OBJS)   

output/%.o: src/%.cpp
	mkdir -p $(dir $@)
	$(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	rm -f $(OBJS) $(EXEC)

.PHONY: all clean