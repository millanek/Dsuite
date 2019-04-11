
CXXFLAGS=-std=c++11
CXX=g++
BIN := Build
LDFLAGS=-lz

all: $(BIN)/Dsuite

$(BIN)/evo: $(BIN)/Dsuite.o $(BIN)/Dsuite_utils.o $(BIN)/D.o $(BIN)/gzstream.o $(BIN)/Dmin.o $(BIN)/Dmin_combine.o
	$(CXX) $(CXXFLAGS) $(LDFLAGS) $^ -o $@

$(BIN)/%.o: %.cpp
	$(CXX) -c $(CXXFLAGS) $< -o $@

$(BIN):
	mkdir -p $@

# Dependencies
$(BIN)/Dsuite: $(BIN)/Dsuite.o $(BIN)/Dsuite_utils.o $(BIN)/D.o $(BIN)/gzstream.o $(BIN)/Dmin.o $(BIN)/Dmin_combine.o | $(BIN)
