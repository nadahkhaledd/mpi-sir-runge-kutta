CXX = mpic++    
CXXFLAGS = -Wall -O2    

SRCS = SIRCell.cpp CSVParser.cpp SIRModel.cpp GridSimulation.cpp MPIHandler.cpp main.cpp        
OBJS = $(SRCS:.cpp=.o)  
EXEC = sir_simulation   

all: $(EXEC)

$(EXEC): $(OBJS)
	    $(CXX) $(CXXFLAGS) -o $@ $(OBJS)   

%.o: %.cpp
	    $(CXX) $(CXXFLAGS) -c $< -o $@

clean:
	    rm -f $(OBJS) $(EXEC)

.PHONY: all clean
